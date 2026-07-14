Domain definition file (``domain.yaml``)
========================================

The domain definition file tells TOPO how to scale the **sidechain–sidechain
(SS) contact energies** of the structure-based non-bonded potential on a
per-domain basis. It is passed through ``domain_def`` in ``md.ini`` and consumed
by :func:`topo.utils.nonbonded.read_yaml_config` and
:func:`~topo.utils.nonbonded.get_scaling_ss_matrix`.

It is **optional**: if you omit ``domain_def``, every SS contact is scaled by
1.0 (i.e. no domain scaling at all).

.. tip::

   **Where to get domain boundaries.** If your simulation uses a native
   structure from the AlphaFold Database, you can obtain ready-made domain
   assignments (residue ranges) from `TED — The Encyclopedia of Domains
   <https://ted.cathdb.info/>`_, which provides consensus domain definitions
   for the AlphaFold predicted structures. Look up your model by its UniProt
   accession and use the reported domain residue ranges to populate the
   ``intra_domains`` blocks below.


YAML syntax in 60 seconds
-------------------------

If you have not used YAML before, these are the only rules you need for a
``domain.yaml``. YAML is a plain-text format for nested key/value data.

* **Key/value pairs** are written ``key: value`` — note the **space after the
  colon** (``nscale:1.1`` is wrong; ``nscale: 1.1`` is right).
* **Nesting is by indentation**, using **spaces only — never tabs**. Use a
  consistent number of spaces per level (2 is used throughout this page). A
  block indented further "belongs to" the key above it::

      intra_domains:      # this key ...
        A:                #   ... has child A (indented 2 spaces)
          nscale: 1.1     #     ... which has child nscale (indented 4 spaces)

* **Lists** are written either one item per line, each starting with a dash
  and a space (``-``)::

      residues:
        - 1-117
        - 166-214

  or inline in square brackets on one line — ``residues: [1-117, 166-214]``.
  Both forms are identical; this page uses the inline ``[...]`` form because
  residue lists are short.
* **Inline mappings** collapse a whole block onto one line with braces:
  ``A: { residues: [1-164], nscale: 1.114 }`` is exactly the same as the
  three-line indented block for ``A``. Use whichever reads better.
* **Comments** start with ``#`` and run to the end of the line. They are
  ignored by the parser.
* **Numbers vs. text.** ``1.1556`` and ``214`` are read as numbers.
  A residue *range* like ``1-117`` contains a dash, so YAML reads it as the
  text string ``"1-117"`` — which is exactly what TOPO wants (see
  ``residues`` below). You normally do **not** need quotes; add them only if
  you want to be explicit (``"1-117"`` and ``1-117`` behave identically here).

.. tip::

   Indentation mistakes (a stray tab, or a child not indented under its
   parent) are the most common cause of a YAML file that "looks right" but
   fails to parse. If you get a parse error, check indentation first. Any
   online "YAML lint" validator, or ``python -c "import yaml,sys;
   print(yaml.safe_load(open('domain.yaml')))"``, will confirm the file is
   syntactically valid before you feed it to TOPO.


Quick reference
---------------

.. code-block:: yaml

    n_residues: 214                      # REQUIRED: total residue count (1..n_residues)
    intra_domains:                       # REQUIRED: at least one domain
      A:
        residues: [1-117, 166-214]       # list of ranges ("a-b", inclusive) and/or ints
        class: alpha-beta                # OPTIONAL: structural class (optimizer only)
        nscale: 1.1556                 # SS-contact scale factor WITHIN domain A
      B:
        residues: [118-165]
        class: alpha
        nscale: 1.6871
    inter_domains:                       # OPTIONAL: scale factor BETWEEN domains
      A-B: 1.8611                         # key is "<domain1>-<domain2>"


Field reference
---------------

Every key TOPO reads from ``domain.yaml``, its type, and its allowed values.
Any key **not** in this table is silently ignored by the reader.

.. list-table::
   :header-rows: 1
   :widths: 22 12 12 54

   * - Key
     - Required?
     - Type
     - Meaning / allowed values
   * - ``n_residues``
     - **yes**
     - integer
     - Total residue count of the chain. Defines the full range ``1..n_residues``
       used to detect unassigned residues. Must match the structure being
       simulated.
   * - ``intra_domains``
     - **yes**
     - mapping
     - One entry per domain. Each key is a **domain name** (see naming rules
       below); each value is a mapping with ``residues`` and ``nscale`` (and
       optionally ``class``).
   * - ``intra_domains.<name>.residues``
     - **yes**
     - list
     - Residues belonging to this domain. Each list item is either a **bare
       integer** (``5``) or a **range string** ``"start-end"`` (``1-117``,
       inclusive of both ends). Mix freely: ``[1, 2, "5-10", 166-214]``. A
       domain may list several ranges to cover non-contiguous segments.
   * - ``intra_domains.<name>.nscale``
     - **yes**
     - float
     - Multiplicative scale factor for SS contact energy **within** this domain.
       ``1.0`` = unscaled. Any positive float is allowed; ``0.0`` removes the
       domain's internal native contacts.
   * - ``intra_domains.<name>.class``
     - no
     - string
     - Structural class, used **only** by the nscale optimizer, ignored by the
       runner. One of ``alpha``, ``beta``, ``alpha-beta`` (or ``a`` / ``b`` /
       ``c``).
   * - ``intra_domains.<name>.strength``
     - no
     - float
     - **Deprecated** alias for ``nscale``. Still accepted (prints a one-time
       notice); use ``nscale`` instead.
   * - ``inter_domains``
     - no
     - mapping
     - Scale factors **between** domains. Each key is ``"<name1>-<name2>"``;
       each value is a float. Symmetric — ``A-B`` also covers ``B-A``. Omit for
       single-domain systems. Any pair you do not list defaults to ``1.0``.

.. important::

   **Domain-name rules.** A domain name may be any text (``A``, ``B``,
   ``NTD``, ``core`` …) with two hard constraints:

   * It must **not contain a hyphen** (``-``). The ``inter_domains`` keys are
     split on ``-`` to recover the two domain names, so a name like ``N-term``
     would break that parsing. Use ``Nterm`` or ``N_term`` instead.
   * The name **``X`` is reserved.** TOPO auto-creates a domain called ``X``
     for any residues you leave unassigned (see the note under *How it is
     interpreted*). If you define your own ``X`` **and** also leave some
     residues unassigned, your ``X`` will be overwritten. Do not name a domain
     ``X``.
   * Names are **case-sensitive and matched exactly** between blocks. The name
     used in an ``inter_domains`` key must be spelled identically to the
     domain defined in ``intra_domains`` — if you define ``NTD`` and ``CTD``,
     the interface key must be ``NTD-CTD`` (not ``ntd-ctd``, ``Ntd-Ctd``, or a
     different label). A mismatched name is **not** an error: the pair you
     wrote is stored under a name that never appears in ``intra_domains``, so
     it is never applied, and the real interface silently stays at the default
     scale of ``1.0``.


How it is interpreted
---------------------

* **Residue numbering is 1-based** and matches the residue order of the input
  PDB. Ranges are written as the **string** ``"start-end"`` and are
  **inclusive** of both ends (``1-117`` → residues 1, 2, …, 117). Single
  residues may be given as bare integers (e.g. ``[5, 9, 20-25]``).
* A single domain may span **multiple, non-contiguous segments** — just list
  several ranges (see the multi-segment scenario below). The contact matrix is
  ordered by sorted residue number internally, so the segments do not need to be
  contiguous or in order.
* ``nscale`` is the multiplicative scale factor applied to the SS contact
  **energy** for native contacts whose two residues are in that domain. It is a
  *scaling factor*, not an absolute energy — at ``nscale = 1.0`` the contacts keep
  their unscaled sidechain–sidechain well depths.

  .. note::

     ``nscale`` was previously called ``strength``. The old key is still accepted
     as a **deprecated alias** (the reader prints a one-time notice and keeps
     working), but new ``domain.yaml`` files should use ``nscale``.
* ``inter_domains`` gives the scale factor for native contacts whose two
  residues are in **different** domains. The key ``"A-B"`` is split on ``-`` and
  stored symmetrically, so ``A-B`` also covers ``B-A``.
* Only the **SS** part of the contact energy is scaled. The hydrogen-bond and
  backbone–sidechain contributions are not affected by these factors.
* ``class`` is an **optional** per-domain field (``alpha``, ``beta`` or
  ``alpha-beta``; ``a``/``b``/``c`` accepted). It is used **only** by the nscale
  optimizer (:doc:`../tutorials/05_opt_nscal`) to pick which *n*\ :sub:`scale`
  ladder a domain climbs, and is **ignored** by the runner — the YAML reader only
  reads ``residues`` and ``nscale``. Include it on domains you intend to
  optimize; omit it otherwise.

.. important::

   **If an inter-domain pair is not listed in ``inter_domains``, its scale
   factor defaults to 1.0** (the identity — no scaling). The scaling matrix only
   *modulates* contacts that already exist in the native contact map, so an
   unspecified interface leaves those native contacts at full nscale rather
   than removing them. To intentionally **decouple** two domains (remove their
   inter-domain native contacts), set the pair explicitly to ``0.0``.

.. note::

   **Unassigned residues** (any residue in ``1..n_residues`` not listed in a
   domain) are automatically collected into a domain named ``'X'`` with intra
   nscale ``1.0`` and inter nscale ``1.0`` to every other domain. This is a
   convenience fallback; for reproducible runs it is best to assign every
   residue explicitly so no residue silently lands in ``X``.


Scenarios
---------

1. Single domain (whole protein, no special scaling)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Either omit ``domain_def`` entirely (everything scaled by 1.0), or define one
domain covering the whole chain with the desired uniform factor:

.. code-block:: yaml

    n_residues: 106
    intra_domains:
      A:
        residues: [1-106]
        nscale: 2.5044

No ``inter_domains`` block is needed because there is only one domain.


2. Two contiguous domains with a coupled interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    n_residues: 283
    intra_domains:
      A:
        residues: [1-164]
        nscale: 1.114
      B:
        residues: [165-283]
        nscale: 1.114
    inter_domains:
      A-B: 2.124

Native A–A and B–B contacts are scaled by 1.114; A–B interface contacts by
2.124.


3. Multi-segment (discontiguous) domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A domain that is split in sequence (e.g. an N-terminal lobe that wraps around a
second domain). Adenylate kinase (1AKE) is a classic example: domain A is
residues 1–117 **and** 166–214, with domain B (118–165) inserted between them.

.. code-block:: yaml

    n_residues: 214
    intra_domains:
      A:
        residues: [1-117, 166-214]     # both segments belong to ONE domain A
        nscale: 1.1556
      B:
        residues: [118-165]
        nscale: 1.6871
    inter_domains:
      A-B: 1.8611

All contacts within A (including segment-1↔segment-2 contacts) use 1.1556;
within B use 1.6871; across the A/B interface use 1.8611.


4. Two domains, intentionally decoupled
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make two domains independent (remove their inter-domain native contacts),
set the pair to ``0.0`` **explicitly**:

.. code-block:: yaml

    n_residues: 283
    intra_domains:
      A: { residues: [1-164],   nscale: 1.0 }
      B: { residues: [165-283], nscale: 1.0 }
    inter_domains:
      A-B: 0.0                  # explicitly remove A-B interface contacts

(The compact ``{ ... }`` inline-mapping form shown here is equivalent to the
indented block form.) Note: simply *omitting* ``A-B`` would leave the interface
at the default scale of 1.0 (contacts kept), not remove it.


5. Three or more domains with selective coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

List each inter-domain pair whose nscale you want to change from the default
of 1.0. Any pair you omit keeps its native contacts at scale 1.0; to remove an
interface, set it to ``0.0``.

.. code-block:: yaml

    n_residues: 400
    intra_domains:
      A: { residues: [1-120],   nscale: 1.20 }
      B: { residues: [121-260], nscale: 1.35 }
      C: { residues: [261-400], nscale: 1.10 }
    inter_domains:
      A-B: 1.50
      B-C: 1.45
      A-C: 0.0     # explicitly decouple A and C (omitting it would keep them at 1.0)


6. Partial assignment (auto ``X`` domain)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you only assign part of the chain, the remaining residues are gathered into
domain ``X`` (intra 1.0, inter 1.0 to all). Useful for a flexible linker or tail
you want left at the default scale:

.. code-block:: yaml

    n_residues: 214
    intra_domains:
      A: { residues: [1-117, 166-214], nscale: 1.1556 }
      B: { residues: [118-160],        nscale: 1.6871 }
      # residues 161-165 are unassigned -> domain X (nscale 1.0,
      # X-A = X-B = 1.0)


Common pitfalls
---------------

* **Empty / missing ``nscale``** — every domain must have a numeric
  ``nscale``. A blank value parses as ``None`` and raises an error when the
  scaling matrix is built.
* **Omitting an ``inter_domains`` pair** leaves that interface at the default
  scale of 1.0 (native contacts kept, unscaled). If you instead want to *remove*
  an interface, you must set the pair to ``0.0`` explicitly.
* **``n_residues`` too small** — must equal the true number of residues; if it
  is smaller than the highest residue you list, the contact matrix and your
  domain assignment will be inconsistent.
* **Quoting ranges** — ``1-117`` is parsed as a YAML string, which is what the
  parser expects. Quotes are optional (``1-117`` and ``"1-117"`` behave
  identically). Mixed ints and range strings — ``[1, 2, "5-10"]`` — are also
  valid. Only single ``start-end`` ranges are supported; do **not** write
  strides or open-ended ranges (``"1-117-2"`` or ``"117-"`` will error).
* **Hyphen in a domain name** — domain names must not contain ``-`` (it clashes
  with the ``inter_domains`` ``"A-B"`` key syntax). See the domain-name rules in
  *Field reference*.
* **Naming a domain ``X``** — ``X`` is reserved for auto-collected unassigned
  residues; a user-defined ``X`` can be silently overwritten. Pick another name.
* **Name mismatch between blocks** — an ``inter_domains`` key must spell the
  domain names exactly as in ``intra_domains`` (case-sensitive). A typo such as
  ``ntd-ctd`` for domains ``NTD``/``CTD`` does not error; it is simply never
  applied, and that interface silently stays at the default scale of ``1.0``.
* **Tabs / bad indentation** — YAML requires spaces, not tabs, and consistent
  indentation per nesting level. This is the most common parse failure; see
  *YAML syntax in 60 seconds*.
* **Chain / residue order** — assignments follow the residue order of the input
  PDB (1-based). Make sure your numbering matches the structure actually being
  simulated.

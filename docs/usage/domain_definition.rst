Domain definition file (``domain.yaml``)
========================================

The domain definition file tells TOPO how to scale the **sidechain–sidechain
(SS) contact energies** of the structure-based non-bonded potential on a
per-domain basis. It is passed through ``domain_def`` in ``md.ini`` and consumed
by :func:`topo.utils.nonbonded.read_yaml_config` and
:func:`~topo.utils.nonbonded.get_scaling_ss_matrix`.

It is **optional**: if you omit ``domain_def``, every SS contact is scaled by
1.0 (i.e. no domain scaling at all).


Quick reference
---------------

.. code-block:: yaml

    n_residues: 214                      # REQUIRED: total residue count (1..n_residues)
    intra_domains:                       # REQUIRED: at least one domain
      A:
        residues: [1-117, 166-214]       # list of ranges ("a-b", inclusive) and/or ints
        class: alpha-beta                # OPTIONAL: structural class (optimizer only)
        strength: 1.1556                 # SS-contact scale factor WITHIN domain A
      B:
        residues: [118-165]
        class: alpha
        strength: 1.6871
    inter_domains:                       # OPTIONAL: scale factor BETWEEN domains
      A-B: 1.8611                         # key is "<domain1>-<domain2>"


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
* ``strength`` is the multiplicative scale factor applied to the SS contact
  **energy** for native contacts whose two residues are in that domain.
* ``inter_domains`` gives the scale factor for native contacts whose two
  residues are in **different** domains. The key ``"A-B"`` is split on ``-`` and
  stored symmetrically, so ``A-B`` also covers ``B-A``.
* Only the **SS** part of the contact energy is scaled. The hydrogen-bond and
  backbone–sidechain contributions are not affected by these factors.
* ``class`` is an **optional** per-domain field (``alpha``, ``beta`` or
  ``alpha-beta``; ``a``/``b``/``c`` accepted). It is used **only** by the strength
  optimizer (:doc:`../tutorials/05_opt_nscal`) to pick which *n*\ :sub:`scale`
  ladder a domain climbs, and is **ignored** by the runner — the YAML reader only
  reads ``residues`` and ``strength``. Include it on domains you intend to
  optimize; omit it otherwise.

.. important::

   **If an inter-domain pair is not listed in ``inter_domains``, its scale
   factor defaults to 1.0** (the identity — no scaling). The scaling matrix only
   *modulates* contacts that already exist in the native contact map, so an
   unspecified interface leaves those native contacts at full strength rather
   than removing them. To intentionally **decouple** two domains (remove their
   inter-domain native contacts), set the pair explicitly to ``0.0``.

.. note::

   **Unassigned residues** (any residue in ``1..n_residues`` not listed in a
   domain) are automatically collected into a domain named ``'X'`` with intra
   strength ``1.0`` and inter strength ``1.0`` to every other domain. This is a
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
        strength: 2.5044

No ``inter_domains`` block is needed because there is only one domain.


2. Two contiguous domains with a coupled interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    n_residues: 283
    intra_domains:
      A:
        residues: [1-164]
        strength: 1.114
      B:
        residues: [165-283]
        strength: 1.114
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
        strength: 1.1556
      B:
        residues: [118-165]
        strength: 1.6871
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
      A: { residues: [1-164],   strength: 1.0 }
      B: { residues: [165-283], strength: 1.0 }
    inter_domains:
      A-B: 0.0                  # explicitly remove A-B interface contacts

(The compact ``{ ... }`` inline-mapping form shown here is equivalent to the
indented block form.) Note: simply *omitting* ``A-B`` would leave the interface
at the default scale of 1.0 (contacts kept), not remove it.


5. Three or more domains with selective coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

List each inter-domain pair whose strength you want to change from the default
of 1.0. Any pair you omit keeps its native contacts at scale 1.0; to remove an
interface, set it to ``0.0``.

.. code-block:: yaml

    n_residues: 400
    intra_domains:
      A: { residues: [1-120],   strength: 1.20 }
      B: { residues: [121-260], strength: 1.35 }
      C: { residues: [261-400], strength: 1.10 }
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
      A: { residues: [1-117, 166-214], strength: 1.1556 }
      B: { residues: [118-160],        strength: 1.6871 }
      # residues 161-165 are unassigned -> domain X (strength 1.0,
      # X-A = X-B = 1.0)


Common pitfalls
---------------

* **Empty / missing ``strength``** — every domain must have a numeric
  ``strength``. A blank value parses as ``None`` and raises an error when the
  scaling matrix is built.
* **Omitting an ``inter_domains`` pair** leaves that interface at the default
  scale of 1.0 (native contacts kept, unscaled). If you instead want to *remove*
  an interface, you must set the pair to ``0.0`` explicitly.
* **``n_residues`` too small** — must equal the true number of residues; if it
  is smaller than the highest residue you list, the contact matrix and your
  domain assignment will be inconsistent.
* **Quoting ranges** — ``1-117`` is parsed as a YAML string, which is what the
  parser expects. Do *not* write it as a number; ``[1, 2, "5-10"]`` (mixed ints
  and range strings) is also valid.
* **Chain / residue order** — assignments follow the residue order of the input
  PDB (1-based). Make sure your numbering matches the structure actually being
  simulated.

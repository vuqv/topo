Native-contact analysis (the *Q* score)
=======================================

The standard observable for a structure-based model is **Q, the fraction of
native contacts formed** — a number between 0 (fully unfolded) and 1 (fully
folded) that tells you, frame by frame, how native-like a conformation is.
TOPO ships a ready-made *Q* calculator,
:mod:`topo.analysis.native_contacts`, that reports *Q* for the **whole protein**
and for **every domain and domain–domain interface** declared in your
``domain.yaml``.

It is the same scorer the nscale optimizer uses internally
(:doc:`../tutorials/A5_opt_nscal`), exposed as a standalone command and a small
Python API.


What *Q* measures
-----------------

A **native contact** is a pair of residues that are close together in the
**folded reference structure**. *Q* is the fraction of those pairs that are still
close in a given simulation frame:

.. math::

   Q(\text{frame}) = \frac{1}{N_\mathrm{native}}
   \sum_{(i,j)\,\in\,\text{native}}
   \Theta\!\big(\,\tau\, d^\mathrm{native}_{ij} - d_{ij}(\text{frame})\,\big)

where :math:`\Theta` is the step function (1 if the contact is *formed*, else 0),
:math:`d^\mathrm{native}_{ij}` is the reference Cα–Cα distance of the contact,
:math:`d_{ij}` is the instantaneous Cα–Cα distance, and :math:`\tau` is a
tolerance factor allowing some stretch.

The method has two stages:

1. **Define native contacts once, from the all-atom reference.** Residues *i* and
   *j* form a native contact when

   * they are more than ``local_separation`` residues apart in sequence
     (:math:`|i-j| >` ``local_separation``, default **3**), **and**
   * at least one **heavy-atom** pair (one atom per residue) lies within
     ``cutoff`` (default **4.5 Å**).

   The reference Cα–Cα distance :math:`d^\mathrm{native}_{ij}` is stored for each.

2. **Score each CG frame.** A native contact is counted as *formed* when

   .. math::

      d_{ij}(\text{frame}) \le \tau\, d^\mathrm{native}_{ij},
      \qquad \tau = \texttt{tolerance}\ (\text{default } 1.2)

   i.e. the Cα–Cα distance is within 20 % of its native value. *Q* is the
   fraction of that group's native contacts that are formed.

Contacts are then **grouped by domain**: an *intra-domain* group holds contacts
with both residues in the same domain; an *interface* group holds contacts with
one residue in each of two domains. A group with no native contacts (e.g. two
domains that never touch) yields ``Q = NaN``.

.. note::

   **This native-contact set is defined more simply than the force field's.**
   The *Q* scorer calls a contact "any heavy-atom pair within 4.5 Å"
   (with :math:`|i-j| > 3`). The model's *energy* contacts (:ref:`theory
   <theory-contacts>`) are built from hydrogen bonds + backbone–sidechain +
   sidechain–sidechain terms with :math:`|i-j| > 2`. The two sets overlap heavily
   but are **not identical** — *Q* is a geometric folding metric, not a readout of
   the exact pairs that carry energy.

   The sequence-separation cutoffs also differ **by design**: :math:`|i-j| > 2`
   (i.e. :math:`\ge 3`) is part of the **model parameterization** — the separation
   at which native pairs are given energy (``LOCAL_SEPARATION = 2``) — whereas the
   *Q* default :math:`|i-j| > 3` (i.e. :math:`\ge 4`) is our **analysis convention**
   for the folding metric. Consequently *Q*'s default does **not** score the
   :math:`|i-j| = 3` pairs the model does carry energetically; pass
   ``--local-separation 2`` if you want *Q* to use the model's separation instead.
   Both are reasonable definitions of "native"; just don't expect the *Q* contact
   count to equal the ``native contacts:`` number printed when the model is built.


Command-line use
----------------

.. code-block:: bash

   python -m topo.analysis.native_contacts \
       -d domain.yaml \          # domain definition (defines the Q groups)
       -r reference.pdb \        # all-atom reference (defines native contacts)
       -p traj/traj.psf \        # CG topology
       -f traj/traj.dcd \        # CG trajectory to score
       -o Q.csv                  # output CSV

Arguments:

.. list-table::
   :header-rows: 1
   :widths: 22 12 12 54

   * - Flag
     - Type
     - Default
     - Meaning
   * - ``-d``, ``--domain``
     - str
     - *required*
     - ``domain.yaml`` whose domains/interfaces define the *Q* groups.
   * - ``-r``, ``--reference``
     - str
     - *required*
     - **All-atom** reference structure that defines the native contacts (typically the same PDB you simulated).
   * - ``-p``, ``--psf``
     - str
     - *required*
     - CG (Cα) topology, e.g. ``traj/traj.psf``.
   * - ``-f``, ``--dcd``
     - str
     - *required*
     - CG trajectory to score, e.g. ``traj/traj.dcd``.
   * - ``-o``, ``--output``
     - str
     - *required*
     - Output CSV path (parent dirs created if needed).
   * - ``--cutoff``
     - float [Å]
     - ``4.5``
     - Heavy-atom distance defining a native contact.
   * - ``--local-separation``
     - int
     - ``3``
     - Minimum sequence separation (``|i-j| >`` this).
   * - ``--tolerance``
     - float
     - ``1.2``
     - Stretch factor :math:`\tau`; a contact is formed when ``d <= tolerance * d_native``.
   * - ``--start-frame``
     - int
     - ``0``
     - First frame to score (inclusive).
   * - ``--end-frame``
     - int
     - all
     - Last frame to score (exclusive).

.. important::

   The ``--reference`` must be the **all-atom** structure (it needs sidechain
   heavy atoms to define contacts), while ``--psf``/``--dcd`` are the **CG**
   topology and trajectory (one Cα bead per residue). The CG model must have
   exactly one bead per reference residue, or the tool errors out.


Output CSV
----------

One row per scored frame, one column per *Q* group:

.. code-block:: text

   Frame,Q_protein,Q_A,Q_B,Q_A-B
   0,0.981,0.992,0.965,0.900
   1,0.978,0.990,0.961,0.880
   ...

* ``Frame`` — 0-based frame index (respects ``--start-frame``).
* ``Q_protein`` — whole-protein *Q* (all native contacts).
* ``Q_<domain>`` — one column per domain (intra-domain *Q*).
* ``Q_<a>-<b>`` — one column per domain pair (interface *Q*; ``NaN`` if the pair
  shares no native contacts).

The tool also prints the native-contact count of each group and the mean *Q* per
group to the console.


Python API
----------

The same routines are importable for custom analysis (this is what the optimizer
uses):

.. code-block:: python

   import MDAnalysis as mda
   from topo.analysis import (
       load_domains, reference_residue_geometry,
       build_native_contacts, fraction_native_contacts,
   )
   import numpy as np

   # 1. Reference geometry + heavy atoms
   u_ref = mda.Universe("reference.pdb")
   ca_positions, resindex_to_pos = reference_residue_geometry(u_ref)
   n_res = ca_positions.shape[0]
   heavy = u_ref.select_atoms("protein and not name H*")
   heavy_res = np.fromiter((resindex_to_pos[ri] for ri in heavy.resindices),
                           dtype=int, count=heavy.n_atoms)

   # 2. Native contacts for the whole protein (group1=group2=None)
   pairs, dnat = build_native_contacts(
       None, None, heavy.positions, heavy_res, ca_positions,
       cutoff=4.5, local_separation=3)

   # 3. Score a CG trajectory frame by frame
   u_cg = mda.Universe("traj/traj.psf", "traj/traj.dcd")
   for ts in u_cg.trajectory:
       q = fraction_native_contacts(u_cg.atoms.positions, pairs, dnat, tolerance=1.2)
       print(ts.frame, q)

* :func:`~topo.analysis.native_contacts.build_native_contacts` takes
  ``(group1, group2, ...)``: ``(None, None)`` = whole protein, ``(set, None)`` =
  intra-domain, ``(set, set)`` = interface.
* :func:`~topo.analysis.native_contacts.load_domains` reads the domain residue
  sets (as 0-based indices) from ``domain.yaml``.


Typical workflow
----------------

1. **Run a simulation** (:doc:`../tutorials/A1_single_domain`) producing
   ``traj/traj.psf`` and ``traj/traj.dcd``.
2. **Score it** against the structure you built from::

      python -m topo.analysis.native_contacts -d domain.yaml \
          -r P0CX28_clean.pdb -p traj/traj.psf -f traj/traj.dcd -o Q.csv

3. **Plot / interpret.** A stable folded run sits near ``Q ≈ 1``; raising
   ``ref_t`` toward the melting temperature drives *Q* down as native contacts
   break. In an annealing run (:doc:`../tutorials/A6_anneal`), score the
   ``_quench`` trajectory to confirm *Q* → 0 (the protein really unfolded) and
   the production trajectory to watch *Q* climb back as it refolds.

For multi-copy runs, split the combined trajectory first
(:func:`topo.split_chains`) and score each per-copy DCD independently — exactly
what the optimizer does to gather statistics over many trajectories.

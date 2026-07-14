.. TOPO citation / how-to-cite page.

How to cite TOPO
================

If TOPO contributed to work you are publishing, please cite it. Citing the
software lets others find the exact model and version you used, and credits the
people who maintain it. If your work is published with the help of TOPO, please
also give it a star on `GitHub <https://github.com/vuqv/topo>`_ in addition to
citing it — it helps others find the project.

.. note::

   No paper or preprint describing TOPO exists yet, so the software entry below
   is the primary reference. Once a paper is available, cite that as the primary
   reference and keep the software entry as a secondary, version-specific
   citation. Keep :download:`CITATION.cff <../CITATION.cff>` at the repository
   root in sync with this page.


Primary citation
----------------

Cite the software itself:

   Vu, Q. (2026). *TOPO: a topology-based coarse-grained model for folded
   proteins* (Version 2026.1) [Computer software]. Zenodo.
   https://doi.org/10.5281/zenodo.21360706

The DOI above is the **concept DOI**, which always resolves to the latest
release. To pin the exact version you ran, use the version-specific DOI for
2026.1 instead: https://doi.org/10.5281/zenodo.21360707

BibTeX
~~~~~~

A software entry you can drop into your ``.bib`` file:

.. code-block:: bibtex

   @software{topo,
     author  = {Vu, Quyen},
     title   = {{TOPO}: a topology-based coarse-grained model for folded proteins},
     year    = {2026},
     version = {2026.1},
     doi     = {10.5281/zenodo.21360706},
     url      = {https://github.com/vuqv/topo},
     note     = {Built on OpenMM}
   }

When the accompanying article is published, add an ``@article`` entry and cite
both:

.. code-block:: bibtex

   @article{topo_paper,
     author  = {[AUTHORS]},
     title   = {[PAPER TITLE]},
     journal = {[JOURNAL]},
     year    = {[YEAR]},
     volume  = {[VOL]},
     pages   = {[PAGES]},
     doi     = {[DOI]}
   }


Cite the version you used
-------------------------

Different TOPO versions can produce different numbers, so record the exact
version and, if possible, the commit:

* **Version** — the release tag, or ``import topo; print(topo.__version__)`` /
  the value in ``pyproject.toml``.
* **Commit** — ``git rev-parse --short HEAD`` inside the source tree.

If you archive a release on `Zenodo <https://zenodo.org/>`_ (or a similar
service), it mints a citable, version-specific DOI. Put the *concept* DOI (all
versions) in the template above and the *version* DOI in your paper's methods.


The models TOPO implements
--------------------------

TOPO was designed specifically to reproduce the O'Brien-lab coarse-grained and
co-translational-synthesis models. **After the software itself, these are the most
important references to cite** — they define the physics TOPO runs. Cite the
coarse-grained model whenever you use TOPO at all, and add the CSP reference when
you run synthesis (``topo-csp`` / ``topo-cylinder``).

Coarse-grained (structure-based / Gō-like) model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The one-bead-per-residue Cα potential — bonds, angles, sequence-dependent
dihedrals, and native-contact wells:

* O'Brien, E. P., Christodoulou, J., Vendruscolo, M. & Dobson, C. M. Trigger
  factor slows co-translational folding through kinetic trapping while
  sterically protecting the nascent chain from aberrant cytosolic interactions.
  *J. Am. Chem. Soc.* **134**\ (26):10920–10932 (2012).
  https://doi.org/10.1021/ja302305u
* Jiang, Y. *et al.* How synonymous mutations alter enzyme structure and
  function over long timescales. *Nat. Chem.* **15**:308–318 (2023).
  https://doi.org/10.1038/s41557-022-01091-z

Component potentials of that model, depending on which terms your run uses:

* **Native-contact wells (LJ 12-10-6).** Karanicolas, J. & Brooks III, C. L. The
  origins of asymmetry in the folding transition states of protein L and protein
  G. *Protein Sci.* **11**\ (10):2351–2361 (2002).
  https://doi.org/10.1110/ps.0205402
* **Pairwise contact energies (BT reference-state potential).** Betancourt, M. R.
  & Thirumalai, D. Pair potentials for protein folding: choice of reference
  states and sensitivity of predicted native states to variations in the
  interaction schemes. *Protein Sci.* **8**\ (2):361–369 (1999).
  https://doi.org/10.1110/ps.8.2.361
* **Double-well angle potential.** Best, R. B., Chen, Y.-G. & Hummer, G. Slow
  protein conformational dynamics from multiple experimental structures: the
  helix/sheet transition of arc repressor. *Structure* **13**\ (12):1755–1763
  (2005). https://doi.org/10.1016/j.str.2005.08.009

Co-translational synthesis (CSP) model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Used by the ``topo-csp`` / ``topo-cylinder`` runners. The per-codon, three-stage
continuous-synthesis protocol reproduces the O'Brien-lab elongation scheme,
described in the same **Jiang et al. (2023)** paper cited above — cite it here too
when you run synthesis.


Other methods and datasets
--------------------------

Depending on which parts of the software you use, please also cite the following.

**MD engine (always).** All dynamics run on OpenMM:

* Eastman, P. *et al.* OpenMM 7: Rapid development of high performance algorithms
  for molecular dynamics. *PLoS Comput. Biol.* **13**\ (7):e1005659 (2017).
  https://doi.org/10.1371/journal.pcbi.1005659

**Codon dwell-time datasets.** If you run synthesis with a species dwell-time
table (see :doc:`usage/codon_dwell_times`), cite the source dataset:

* *E. coli* — Fluitt, A., Pienaar, E. & Viljoen, H. *Comput. Biol. Chem.*
  **31**:335–346 (2007). https://doi.org/10.1016/j.compbiolchem.2007.07.003
* Yeast — Gardin, J. *et al.* *eLife* **3**:e03735 (2014).
  https://doi.org/10.7554/eLife.03735
* Human — Gobet, C. *et al.* *PNAS* **117**\ (17):9630–9641 (2020).
  https://doi.org/10.1073/pnas.1918145117
* *N. crassa* — Yang, Q. *et al.* *Nucleic Acids Res.* **47**\ (17):9243–9258
  (2019). https://doi.org/10.1093/nar/gkz710


Machine-readable metadata
-------------------------

The repository ships a `Citation File Format
<https://citation-file-format.github.io/>`_ file,
:download:`CITATION.cff <../CITATION.cff>`, at its root. GitHub reads it to show
a **“Cite this repository”** button that exports formatted APA and BibTeX
automatically. Keep ``CITATION.cff`` and this page in step whenever the citation
details change.


Questions
---------

For anything not covered here — collaboration, a preprint DOI, or how to cite a
specific analysis — open an issue on the
`GitHub repository <https://github.com/vuqv/topo/>`_.

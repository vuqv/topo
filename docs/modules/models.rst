Models
======

The models module provides a single entry point for building the TOPO coarse-grained system: ``buildCoarseGrainModel``. It initializes a system with the necessary force field parameters and optional structure-based non-bonded interactions.

Coarse-grained alpha-carbon (CA) model
++++++++++++++++++++++++++++++++++++++++

The TOPO model represents the protein as beads centered at the alpha carbons of each residue. It uses:

- **Bonded**: harmonic bonds, Gaussian-sum angles, and periodic torsion potentials to maintain chain connectivity and local geometry.
- **Non-bonded**: structure-based (native) contacts built from the input structure, STRIDE hydrogen-bond output, and optional domain definition; electrostatics via Yukawa potential.

To create a TOPO CA model, call:

.. code-block:: python

   topo.models.buildCoarseGrainModel(structure_file, minimize=False, model='topo',
                                     domain_def=None, stride_output_file=None, box_dimension=None)

- **structure_file**: path to the PDB/CIF structure.
- **minimize**: if ``True``, run energy minimization after building the system.
- **model**: model name; currently only ``'topo'`` is supported.
- **domain_def**: optional path to domain definition file (e.g. ``domain.yaml``) for contact scaling between domains.
- **stride_output_file**: optional path to STRIDE output file for hydrogen-bond-based contact energies. If ``None``, the non-bonded builder may attempt to run STRIDE on the structure.
- **box_dimension**: optional float or list of three floats [x, y, z] in nm to enable periodic boundary conditions (cubic or rectangular box).

The Hamiltonian has the form:

.. math::

   H = \sum_{bonds} V_{bond} + \sum_{angle} V_{angle} + \sum_{torsion} V_{torsion}
       + \sum_{i,j} \Phi_{ij}^{nb} + \sum_{i,j} \Phi_{ij}^{el}

where :math:`\Phi_{ij}^{nb}` is the structure-based non-bonded contact potential and :math:`\Phi_{ij}^{el}` is the Yukawa electrostatics.

Bonded potential
++++++++++++++++

**Harmonic bond**

.. math::

   V_{bond} = \frac{k_b}{2}(r - r_0)^2

Default values: :math:`k_b = 20920\,\mathrm{kJ/(mol\,nm^2)}`, :math:`r_0 = 0.381\,\mathrm{nm}`.

Angle potential (Gaussian-sum)
+++++++++++++++++++++++++++++

.. math::

   U_{angle}(\theta) = -\frac{1}{\gamma}
   \ln \left[ e^{-\gamma[ k_{\alpha} (\theta-\theta_{\alpha})^2 + \epsilon_{\alpha} ]}
        + e^{-\gamma k_{\beta} (\theta-\theta_{\beta})^2} \right]

Parameters: :math:`\gamma = 0.1\,\mathrm{mol/kcal}`,
:math:`\epsilon_{\alpha} = 4.3\,\mathrm{kcal/mol}`,
:math:`\theta_{\alpha} = 1.6\,\mathrm{rad}`,
:math:`\theta_{\beta} = 2.27\,\mathrm{rad}`.

Torsion potential (periodic)
++++++++++++++++++++++++++++

.. math::

   U_{torsion}(\theta) = -\ln\left[ U_{torsion,\alpha}(\theta, \epsilon_d)
       + U_{torsion,\beta}(\theta, \epsilon_d) \right]

with

.. math::

   U_{torsion,\alpha}(\theta, \epsilon_d)
   = e^{-k_{\alpha,1}(\theta-\theta_{\alpha,1})^2 - \epsilon_d}
   + e^{-k_{\alpha,2}(\theta-\theta_{\alpha,2})^4 + e_0}
   + e^{-k_{\alpha,2}(\theta-\theta_{\alpha,2}+2\pi)^4 + e_0}

.. math::

   U_{torsion,\beta}(\theta, \epsilon_d)
   = e^{-k_{\beta,1}(\theta-\theta_{\beta,1})^2 + e_1 + \epsilon_d}
   + e^{-k_{\beta,1}(\theta-\theta_{\beta,1}-2\pi)^2 + e_1 + \epsilon_d}
   + e^{-k_{\beta,2}(\theta-\theta_{\beta,2})^4 + e_2}
   + e^{-k_{\beta,2}(\theta-\theta_{\beta,2}-2\pi)^4 + e_2}

Parameters (from Karanicolas–Brooks style dihedral set): e.g.
:math:`k_{\alpha,1} = 11.4\,\mathrm{kcal/(mol\,rad^2)}`,
:math:`k_{\alpha,2} = 0.15\,\mathrm{kcal/(mol\,rad^4)}`,
:math:`\theta_{\alpha,1} = 0.9\,\mathrm{rad}`,
:math:`\theta_{\alpha,2} = 1.02\,\mathrm{rad}`,
:math:`e_0 = 0.27\,\mathrm{kcal/mol}`,
and corresponding :math:`k_{\beta}`, :math:`\theta_{\beta}`, :math:`e_1`, :math:`e_2` (see ``topo.parameters``).

Non-bonded potentials
+++++++++++++++++++++

**Structure-based contacts**

Native contacts are built by ``topo.utils.build_nonbonded_interaction.build_nonbonded_interaction`` from:

- Hydrogen bonds (STRIDE): 0.75 kcal/mol per H-bond; pairs with 2+ H-bonds are capped at 1.5 kcal/mol.
- Backbone–sidechain and sidechain–sidechain contacts from distances (e.g. 4.5 Å cutoff) and optional domain scaling (YAML).

The result is a distance matrix and an energy matrix (in nm and kJ/mol) passed to the system as a custom non-bonded force. Non-native pairs use a repulsive well (sigma from structure).

**Yukawa electrostatics**

.. math::

   \Phi_{ij}^{el}(r) = \frac{q_i q_j}{4\pi\epsilon_0 D r} e^{-\kappa r}

:math:`q_i, q_j` are residue charges; :math:`D` is dielectric constant; :math:`\kappa` is inverse Debye length. Cut-off is typically 3.5 nm when using PBC.

Bonded exclusion rule for non-bonded: **1–2** (nearest-neighbor bonds excluded from non-bonded).

.. autoclass:: topo.core.models

   .. automethod:: __init__
   .. automethod:: buildCoarseGrainModel

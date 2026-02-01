Simulation control options
==========================

An example of how the simulation config file (e.g. ``md.ini``) looks:

.. code-block::

        [OPTIONS]
        md_steps = 500_000   ; number of steps (underscores allowed)
        dt = 0.01 ; time step in ps
        nstxout = 1000 ; number of steps to write checkpoint = nstxout
        nstlog = 1000 ; number of steps to print log
        nstcomm = 100 ; frequency for center of mass motion removal
        ; TOPO model (only option currently)
        model = topo

        ; temperature coupling
        tcoupl = yes
        ref_t = 310          ; Kelvin
        tau_t = 0.01         ; ps^-1

        ; pressure coupling
        pcoupl = no
        ref_p = 1
        frequency_p = 25

        ; periodic boundary condition
        pbc = yes
        ; box_dimension = x (cubic) or [x, y, z] (rectangular), in nm
        box_dimension = 30   ; or [30, 30, 60]

        ; input
        protein_code = 2ww4
        pdb_file = 2ww4.pdb
        ; optional: for structure-based non-bonded (TOPO)
        domain_def = domain.yaml
        stride_output_file = stride.dat
        ; output
        checkpoint = 2ww4.chk
        ; GPU/CPU
        device = GPU
        ppn = 4
        ; restart
        restart = no
        minimize = no

General information
+++++++++++++++++++
Simulation parameters are read from an `.ini` file (e.g. ``md.ini``) using Python's :mod:`configparser`. The section title ``[OPTIONS]`` is required.

* Comments: inline or new line, start with ``;`` or ``#``
* Keyword and value separated by ``=`` or ``:``

Run control
+++++++++++

::

    md_steps:   (int)
                Total number of integration steps. Underscores are allowed (e.g. 500_000).
    ------------------------------------------------------------------------------------
    dt:         (float)
                (0.01) [ps] Time step for integration
    ------------------------------------------------------------------------------------
    nstxout:    (int)
                Steps between writing coordinates and checkpoint to output files
    ------------------------------------------------------------------------------------
    nstlog:     (int)
                Steps between writing energies to the log file
    ------------------------------------------------------------------------------------
    nstcomm:    (int)
                (100) Frequency for center-of-mass motion removal

Model parameter
+++++++++++++++
Currently only the **topo** model is supported: topology-based coarse-grained CA model with structure-based non-bonded contacts and Yukawa electrostatics.

::

    model:      (string)
                topo: TOPO model (default). Uses domain_def and stride_output_file when provided for contact-based non-bonded interactions.


Temperature coupling
+++++++++++++++++++++

::

    tcoupl:     (bool)
                yes (default) : The only available option for now, we don't care about NVE ensemble.
    ------------------------------------------------------------------------------------
    ref_t:      (double)
                (300) [K] : Reference temperature in unit of Kelvin
    ------------------------------------------------------------------------------------
    tau_t:      (double)
                [ps^-1] : The friction coefficient which couples the system to the heat bath (in inverse picoseconds)

Pressure coupling
+++++++++++++++++++

::

    pcoupl      (bool)
                yes : Using pressure coupling

                no (default) : Run on NVT ensemble only
    ------------------------------------------------------------------------------------
    ref_p       (double)
                 (1) [bar] The default pressure acting on the system.
    ------------------------------------------------------------------------------------
    frequency_p (int)
                (25) [steps] the frequency at which Monte Carlo pressure changes should be attempted

Periodic boundary condition:
+++++++++++++++++++++++++++++
if pcoupl is yes then pbc must be yes.

::

    pbc         (bool)
                yes : Using periodic boundary condition.
                        If this option is chosen, then it will affect to non-bonded forces in the system,
                        and the coordinate writen in PDB and DCD file as well. No worries since I have handled these.

                no (default) : Without periodic boundary condition.
    ------------------------------------------------------------------------------------
    box_dimension   (float or list of float)
                [nm] An example of box dimension:
                If you want a cubic box of 30x30x30 nm^3, put: 30 or [30, 30, 30]
                If you want a rectangular box? Put:  [30, 30, 60]

File input/output
+++++++++++++++++

::

    protein_code    (string)
                    Output prefix, e.g. {protein_code}.dcd, {protein_code}.log
    ------------------------------------------------------------------------------------
    pdb_file        (string)
                    [.pdb, .cif] Input structure for topology and initial coordinates
    ------------------------------------------------------------------------------------
    domain_def      (string, optional)
                    Path to domain definition YAML (e.g. domain.yaml) for TOPO structure-based non-bonded scaling. Omit for single-domain or if not using domain scaling.
    ------------------------------------------------------------------------------------
    stride_output_file  (string, optional)
                    Path to STRIDE output file for hydrogen-bond-based contact energies in TOPO. If omitted, the builder may attempt to run STRIDE on the structure.
    ------------------------------------------------------------------------------------
    checkpoint      (string)
                    [.chk] Checkpoint file name; required for restart and for saving state.

Simulation platform
+++++++++++++++++++++
Simulation can be run on CPU with number of threads is control by `ppn` or using GPU.
If `device=CPU` then ppn need to be specify, otherwise simulation will run on 1 core

::

    device          (string)
                    GPU : Use gpu to run simulation

                    CPU (default) : use cpu to run simulation, if you specify cpu, you should modify ppn option, it control
                            how many cores will be used to run simulation, if not, default is 1.
    ------------------------------------------------------------------------------------
    ppn             (int)
                    (1) [threads] Number of threads used to run simulation on CPU. When using GPU,
                                performance is boosted a lot so ppn in that case is set to 1.

Restart simulation
++++++++++++++++++++

::

    restart         (bool)
                    yes : restart simulation from checkpoint file. This can be True, 1 or whatever are not (FALSE)
                            in python condition. If this option is selected, minimize will be force to False.

                    no (default) : Run simulation from beginning, if this option is selected, you can choose if you want to minimize your
                        system before running simulation.
    ------------------------------------------------------------------------------------

    minimize        (bool)
                    yes (default) : perform energy minimization before run molecular dynamics.

                    no : Not running energy minimization. This is default option when restart option is set to yes.
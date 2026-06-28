# Tasks
- [ ] Auto-split per-chain trajectories — **multi-copy runs only** (`n_copies > 1`, non-interacting copies):
  - [ ] At end of run, optionally emit per-chain DCDs `traj_<k>.dcd` (+ appendix) so the manual `split_chains` step isn't needed. Make it a config flag (e.g. `split_chains = yes|no`), not unconditional.
  - [ ] Keep the combined `traj.dcd` as the canonical/source-of-truth output; per-chain files are *derived* (don't replace the combined one — it's needed for some analyses and is the restart/append target).
  - [ ] Per-chain centering applies to independent copies only (current `split_chains(center=True)`).
  - [ ] Do **NOT** auto-split or auto-center genuine *interacting* multichain systems/complexes (see "interacting chains" below) — the inter-chain arrangement is the physics; splitting + independent centering would destroy it.
  - [ ] Decide naming convention: code/optimizer currently use 0-based `traj_0..N-1`; note proposed 1-based `_1..N`. Pick one and document the mapping (e.g. `traj_1.dcd` = copy index 0).
  - [ ] Note: if split runs at finalize, a restart re-splits the whole grown DCD (idempotent but redundant) — acceptable, or split incrementally.
- [ ] interacting chains

## Optional / later
- [ ] OpenMM XML serialization support (discussed; revisit later):
  - [ ] **System XML export** (`XmlSerializer.serialize(system)`): optional flag (e.g. `write_system_xml = yes`) to dump the fully-built System (all forces + tabulated contact R/eps matrices + exclusions) to `<outname>_system.xml`. Value: reproducibility/provenance (pins every force-field number), decouples build from run, and lets reruns load the model with no STRIDE installed.
  - [ ] **State XML** as a *portable* restart alternative to the binary `.chk` (OpenMM checkpoints are not portable across GPU/OpenMM build/platform). Highest practical payoff for long resubmitted runs.
  - [ ] *(optional)* run-from-XML import path in the runner: `system_xml = path` -> skip build, load forces from XML (still build the cheap CA topology from PDB for the Topology + positions).
  - [ ] Caveats to handle: Topology not stored in System XML (pair with .psf/.pdb); n×n tabulated-function size for large complexes (few MB now, tens of MB for ribosome-nascent-chain); reconstruct topoReporter force-group *names* from group order after deserialize.
  - Skip ForceField-template XML (amber14.xml-style): residue templates don't fit a structure-based contact model.

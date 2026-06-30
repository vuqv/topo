topo package
============

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   topo.core
   topo.parameters
   topo.reporter
   topo.utils
   topo.analysis

Submodules
----------

topo.mdrun module
-----------------

The canonical simulation runner (``topo-mdrun`` / ``python -m topo.mdrun``).

.. automodule:: topo.mdrun.mdrun
   :members:
   :undoc-members:
   :show-inheritance:

topo.optimize module
--------------------

The strength (n_scale) optimizer (``topo-optimize`` / ``python -m topo.optimize``).

.. automodule:: topo.optimize.optimize
   :members:
   :undoc-members:
   :show-inheritance:

topo.csp module
-----------------------

The co-translational synthesis package. Its user-facing runner is the
**codon-resolved Continuous Synthesis Protocol** (``topo-csp`` /
``python -m topo.csp``) — see the :doc:`usage/continuous_synthesis` reference and
Tutorials 12/13. The protocol is built on a shared, low-level MD engine
(:mod:`topo.csp.core`) and rigid-ribosome scenery (:mod:`topo.csp.ribosome`), which
the Tutorial-9 cylinder runner also reuses.

Continuous Synthesis Protocol (kinetics + 3-stage loop)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: topo.csp.protocol
   :members:
   :show-inheritance:

.. automodule:: topo.csp.kinetics
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: ribosome_traffic_times

Core MD engine & ribosome scenery (shared machinery)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: topo.csp.core
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: topo.csp.ribosome
   :members:
   :undoc-members:
   :show-inheritance:

Ribosome preparation & movie tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: topo.csp.cg_ribosome
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: topo.csp.truncate_ribosome
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: topo.csp.movie
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: topo
   :members:
   :undoc-members:
   :show-inheritance:

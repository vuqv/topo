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

topo.translation module
-----------------------

The protein synthesis (nascent-chain elongation) runner
(``topo-elongate`` / ``python -m topo.translation``) and the rigid-ribosome
scenery it uses. See the :doc:`usage/protein_synthesis` reference and
Tutorial 7.

.. automodule:: topo.translation.elongate
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: topo.translation.ribosome
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: topo
   :members:
   :undoc-members:
   :show-inheritance:

Pade_PyCpp
==========

Padé approximation in Python and C++ (with Python link)

Source Files
++++++++++++

``src/expm.py``
---------------

This file contains an implementation of the Padé approximation method.
The standard scipy version (``scipy.linalg.expm``) will be replaced with this version in the future.

Testing
+++++++

``src/test/demo_expm.py``
-------------------------
This file executes the Python and C++ versions of the Padé approximation and does some timings.
It produced the output plots ``src/test/results/expmtimings_*.pdf``.

``src/test/demo_cores.py``
--------------------------
This file does some timings for the multithreaded C++ implementation.
It stores results in the files ``src/test/results/<N>threads_gcc.npz`` where the values for ``<N>`` is determined by the ``OMP_NUM_THREADS`` environment variable.

If multiple runs are needed (better results), this file can just be executed multiple times.
It then reads the existing timing information and stores the minimum of that and the new timings.



Intersystem crossing (isc) rate of NV- center in Diamond
===========

Warning!! This tutorial is not up-to-date as of 11/2/20, contact tjsmart@ucsc.edu for help

1. Obtain relaxed geometry of the triplet excited state (`relax-cdftdn1/`).

2. Obtain relaxed geometry of the singlet ground state (`relax-s0/`).

3. Generate a series of linear interpolated structure (e.g. by running `mix_structure_linear.py`) and run SCF calculation for both ground state and excited state (see calculations in `lin1/` and `lin2/`).

4. Calculate isc using `calc_auto.py` (see `isc/`).

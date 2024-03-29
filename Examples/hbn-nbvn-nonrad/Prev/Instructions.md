Nonradiative lifetime of NBVN (in-plane) 6x6 h-BN
===========

Warning!! This tutorial is not up-to-date as of 11/2/20, contact tjsmart@ucsc.edu for help

Follow along on kairay! :)

1. Obtain relaxed geometry of the ground state (a normal QE relax) and that of the excited state (QE relax with constrained occupations wherein a state which is occupied in the ground state has occupation 0, and a state which is unoccupied in the ground state set to 1, to simulate an excited state).

```bash
cd relax-gs
sbatch job
cd ../relax-cdftup1
sbatch job
cd ..
```

2. Generate a series of linear interpolated structure (e.g. by running `mix_structure_linear.py`) and run SCF calculation for both ground state and excited state.

```bash
cd lin-gs
./run.sh
./submit.sh
cd ../lin-cdftup1
./run.sh
./submit.sh
cd ..
```

3. Compute the nonradiative lifetime accordingly using `calc_auto.py`.

```bash
cd ../nonrad
./run.sh
```


4. (Bonus) calculate intersystem crossing rate with a fictitous setup. This is not physical in this case but just demonstrates this calculation (see NV diamond example for realistic example).

```bash
cd ../isc
./run.sh
```

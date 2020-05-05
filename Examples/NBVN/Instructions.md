Nonradiative lifetime of NBVN (in-plane) 6x6 h-BN
===========

1. Obtain relaxed geometry of the ground state electrong state (a normal QE relax) and that of the excited state (QE relax with constrained occupations wherein a state which is occupied in the ground state has occupation 0, and a state which is unoccupied in the ground state set to 1, to simulate an excited state).

```bash
cd relax-gs
sbatch job
cd ../relax-cdftup1
sbatch job
cd ..
```

2. Generate a series of linear interpolated structure (by running `mix_structure_linear.py`) and run SCF calculation for both ground state and excited state.

```bash
cd lin-gs
./run.sh
./submit.sh
cd ../lin-cdftup1
./run.sh
./submit.sh
cd ..
```

3. Compute the nonradiative lifetime accordingly using calc_auto.py.

```bash
cd ../calcauto
./run.sh
```

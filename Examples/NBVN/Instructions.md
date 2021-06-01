Nonradiative lifetime of NBVN (in-plane) 6x6 h-BN
===========

Follow along on kairay! :)

1. Obtain relaxed geometry of the ground state (a normal QE relax) and that of the excited state (QE relax with constrained occupations wherein a state which is occupied in the ground state has occupation 0, and a state which is unoccupied in the ground state set to 1, to simulate an excited state).

```bash
cd relax-gs
sbatch relax.job            # you will need to create your own job script
cd ../relax-cdftup1
sbatch relax.job
cd ..
```

2. Generate a series of linear interpolated structure (e.g. by running `./run_lin.sh`) and run SCF calculation for both ground state (`lin1`) and excited state (`lin2`).

```bash
./run_lin.sh
```

3. Compute the nonradiative lifetime accordingly using `calc_auto.py`.

```bash
cd nonrad
./run.sh
```

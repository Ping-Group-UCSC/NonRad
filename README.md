Non-Radiative Code
===================================

Description
------------------------------------
This code calculates the non-raditive recombination capture coefficient and lifetime as in the articles:

[F. Wu, T. J. Smart, J. Xu, and Y. Ping, *Physical Review B (Rapid Communication)* **100**, 081407 (2019).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.081407 "Carrier recombination mechanism at defects in wide band gap two-dimensional materials from first principles")

[A. Alkauskas, Q. Yan, and C. G. Van de Walle, *Physical Review B* **90**, 075202 (2014).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.075202 "First-principles theory of nonradiative carrier capture via multiphonon emission")

[L. Shi, K. Xu, and L.-W. Wang, *Physical Review B* **91**, 205315 (2015).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.205315 "Comparative study of ab initio nonradiative recombination rate calculations under different formalisms")

Prerequisites:
------------------------------------
* [Python](https://www.python.org/downloads)
* [pyyaml](https://pypi.org/project/PyYAML/), [numpy](https://pypi.org/project/numpy/), [scipy](https://pypi.org/project/scipy/), [mpmath](https://pypi.org/project/mpmath/), [matplotlib](https://pypi.org/project/matplotlib/)
 > pip install pyyaml numpy scipy mpmath matplotlib

Examples:
------------------------------------
Try out the example calculations under the directory `./Examples/`


Author(s)
------------------------------------
This code was written by Feng Wu  
Revised and released by Tyler J. Smart  
Python3 implementation by Tyler J. Smart  
SOC extension implemented by Tyler J. Smart  


Setup for running the code on kairay (with Python2)
------------------------------------
1. Load numpy

```bash
module load numpy/1.11.1
```

2. Install pip

```bash
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py --user
```

Note! This installs pip to your local directory `~/.local/bin`. If this is not in your environmental `$PATH` then you will need to add the line `export PATH=$PATH:~/.local/bin` to your `~/.bashrc`, followed by `source ~/.bashrc`.

Check!! If the above is done correctly running `pip --version` should print something similiar to the below:
 > kairay:[~]$ pip --version  
 > pip 20.1 from /home/tjsmart/.local/lib/python2.7/site-packages/pip (python 2.7)

3. Install needed packages (may take a second)

```bash
pip install pyyaml scipy mpmath matplotlib
```

4. Done! :) You should be good to go.

```bash
python2 /path/to/nonrad/code/calc_auto.py
```
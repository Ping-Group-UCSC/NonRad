#!/usr/bin/env bash
for d in ratio-* ; do cd $d ; sbatch ../job ; cd .. ; done

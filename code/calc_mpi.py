#!/usr/bin/env python


from mpi4py import MPI
# main serial routine still written in calc_auto
from calc_auto import serial_routine, test


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

is_root = (rank == 0)

if ( is_root ):
    print('I am root, my rank is ',rank)
    test()
else:
    print('I am child, my rank is ',rank)


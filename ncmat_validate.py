#!/usr/bin/env python3
from swanfile2 import SwnData
import xarray as xr
import matplotlib.pyplot as plt
import ipdb
import os
import numpy as np

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(\
    description='''Validate swan matlab files converted into netcdf.'''\
    , usage='%(prog)s [options]')

    parser.add_argument('-o',  metavar='out_file', \
        help='output file: (optional, default: ./in_file[PTvar].png)')

    parser.add_argument('-s', metavar='swan in_file', required = True, \
        help='''swan mat file (required).''')

    parser.add_argument('-n', metavar='netcdf in_file', required = True, \
        help='''converted swan netcdf file (required).''')

    parser.add_argument('-w', metavar='warmup', action = 'store_const',
        const = True, default = False,
        help='''The swan file has warm up data starting prior to
        its first supposed output''')

    in_args = parser.parse_args()

    if in_args.o:
        fout = in_args.o

    else:
        fout = in_args.s.replace('.mat','.png')

    sdata = SwnData()

    sdataset = sdata.readswan(in_args.s, in_args.w)
    ndataset = xr.open_dataset(in_args.n)

    for key in sdataset.keys():
        tmp = ndataset[key] - sdataset[key]
        tmp.max(axis = 0).plot()
        plt.title('MAX = {}'.format(tmp.max().data))
        plt.savefig(fout.replace('.png','diff-{}.png'.format(key)))
        plt.clf()

        ndataset[key].isel(time=0).plot()
        plt.savefig(fout.replace('.png','T0-{}.png'.format(key)))
        plt.clf()
        
        ndataset[key].isel(time=-1).plot()
        plt.savefig(fout.replace('.png','T1-{}.png'.format(key)))
        plt.clf()

        t1 = ndataset[key].isel(longitude=10, latitude=10)
        t2 = sdataset[key].isel(longitude=10, latitude=10)
        p1 = np.abs((t1-t2).values).sum()
        ((t1-t2).values/t1).plot()
        plt.savefig(fout.replace('.png','diffpct-{}.png'.format(key)))
        plt.clf()

        match = (
                (ndataset.time.isel(time=0) == \
                sdataset.time.isel(time=0)) &
                (ndataset.time.isel(time=-1) == \
                sdataset.time.isel(time=-1))
                ).values

        print('{}, absolute errors p1: {:1.5f} time match: {} '.format(os.path.basename(in_args.s),p1, match))

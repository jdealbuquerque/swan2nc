#!/usr/bin/env python3
from swanfile2 import SwnData
import xarray as xr
import matplotlib.pyplot as plt

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
        plt.savefig(fout.replace('.png','{}.png'.format(key)))
        plt.close()

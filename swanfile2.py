#!/usr/bin/env python3
import scipy.io
import datetime as dt
import xarray as xr
from collections import OrderedDict as od
import os
import re
import numpy as np
import ipdb

class SwnData:

    # swan vars to be found in file
    __svars_lst = ['Depth',
        'Dir',
        'Dspr',
        'FSpr',
        'Hsig',
        'Hswell',
        'Lwavp',
        'PkDir',
        'Qp',
        'RTm01',
        'RTm_10',
        'RTpeak',
        'Tm02',
        'TPsmoo',
        'Windv_x',
        'Windv_y',
        'Wlen',
        'DrPT',
        'DsPT'
        'HsPT',
        'StPT',
        'TpPT',
        'WfPT',
        'WlPT']


    # ignore variables
    __ivars_lst = ['__header__', '__version__', '__globals__']

    # coordinates variables
    __cvars_lst = ['Xp', 'Yp']

    # netcdf attributes
    _TimeAtt = {
            'long_name': 'Time',
            'standard_name' : 'time'}

    _TimeEnc = {'units': "hours since 1980-01-01 00:00 UTC",
            'dtype' : 'f8'}

    _LonAtt = {'units' : "degrees_east",
            'long_name' : "Longitude",
            'standard_name' : "longitude"}

    _LonEnc = {'dtype' : 'f4'}

    _latAtt = {'units' : "degrees_north",
            'long_name' : "Latitude",
            'standard_name' : "latitude"}

    _LatEnc = {'dtype' : 'f4'}

    _DepthAtt = {'units' : "m", # max: 3826.213, min: 0.062288705, 2^15 = 32768: too small
            'long_name' : "Depth",
            'standard_name' : "sea_floor_depth_below_sea_level"}

    _DepthEnc = {'dtype' : 'f4',
            '_FillValue' : -99999.,
            'scale_factor': 1.0}

    _DirAtt = {'units': 'degrees_true',
            'long_name' : "Mean Wave Direction",
            'standard_name' : "sea_surface_wave_from_direction"}

    _DirEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _DsprAtt = {'units' : 'degrees',
            'long_name' : "Directional Spreading",
            'standard_name' : "sea_surface_wave_directional_spreading"}

    _DsprEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.01}

    _FsprAtt = {'units' : 'dimensionless', #TODO  Double Check 
            'long_name' : "Frequency Spectral Width (Kappa)",
            'standard_name' : "sea_surface_wave_normalized_width_of_frequency_spectrum"}

    _FsprEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.0001}

    _HsigAtt = {'units' : 'm',
            'long_name' : "Significant Wave Height",
            'standard_name' : "sea_surface_wave_significant_height"}

    _HsigEnc = {'dtype' : 'i2', # Corrected
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _HswellAtt = {'units' : 'm',
            'long_name' : "Wave Height of Swell Part",
            'standard_name' : "sea_surface_swell_wave_significant_height"}

    _HswellEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _LwavpAtt = {'units' : 'm',
            'long_name' : "Peak Wave Length",
            'standard_name' : "sea_surface_wave_length_peak"}

    _LwavpEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _PkDirAtt = {'units' : 'degrees_true', 
            'long_name' : "Peak Wave Spectrum Direction",
            'standard_name' : "sea_surfave_wave_peak_direction"}

    _PkDirEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _QpAtt = {'units':'dimensionless',
            'long_name' : "Peakedness of the Wave Spectrum",  
            'standard_name' : "sea_surface_wave_peakedness"}

    _QpEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _RTm01Att = {'units':'s',
            'long_name' : "Mean Relative Wave Period",
            'standard_name' : "sea_surface_wave_mean_relative_period"}

    _RTm01Enc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _RTm_10Att = {'units':'s',
            'long_name' : "Mean Relative Wave Period",
            'standard_name' : "sea_surface_wave_mean_relative_period_from_variance"}

    _RTm_10Enc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _RTpeakAtt = {'units':'s',
            'long_name' : "Peak Period of the Variance Density Spectrum",
            'standard_name' : "sea_surface_wave_relative_peak_period"}

    _RTpeakEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _Tm02Att = {'units':'s',
            'long_name' : "Zero Crossing Period of the Variance Density Spectrum",
            'standard_name' : "sea_surface_wave_zero_crossing_period"}

    _Tm02Enc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _TPsmooAtt = {'units':'s',
            'long_name' : "Smoothed Peak Period of the Variance Density Spectrum",
            'standard_name' : "smoothed_sea_surface_wave_relative_peak_period"}

    _TPsmooEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _Windv_xAtt = {'units' : 'm s-1',
            'long_name' : "Eastward wind 10 m above ground",
            'standard_name' : "eastward_wind"}

    _Windv_xEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _Windv_yAtt = {'units' : 'm s-1',
            'long_name' : "Northward wind 10 m above ground",
            'standard_name' : "Northward_wind"}

    _Windv_yEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _WlenAtt = {'units' : 'm',
            'long_name' : "mean Wave Length of the Variance Density Spectrum",
            'standard_name' : "sea_surface_wave_mean_length"}

    _WlenEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _DrPTAtt = {'units': 'degrees_true',
            'long_name' : "Wave Direction of Partition {}",
            'standard_name' : "partition_{}_wave_from_direction"}

    _DrPTEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _DsPTAtt = {'units': 'degrees',
            'long_name' : "Directional Spreading of Partition {}",
            'standard_name' : "partition_{}_wave_directional_spreading"}

    _DsPTEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.01}

    _HsPTAtt = {'units' : 'm',
            'long_name' : "Significant Wave Height of Partition {}",
            'standard_name' : "partition_{}_wave_significant_height"}

    _HsPTEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _StPTAtt = {'units' : 'dimensionless', # TODO: double check
            'long_name' : "Wave steepness of Partition {}",
            'standard_name' : "partition_{}_wave_steepness"}

    _StPTEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _TpPTAtt = {'units' : 's',
            'long_name' : "Peak Period of Partition {}",
            'standard_name' : "partition_{}_wave_peak_period"}

    _TpPTEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.001}

    _WfPTAtt = {'units' : 'dimensionless', # vmax = 1 #TODO: dcheck units
            'long_name' : "Wind Fraction of Partition {}",
            'standard_name' : "partition_{}_wind_fraction"}

    _WfPTEnc = {'dtype' : 'i2',
            '_FillValue' : -9999,
            'scale_factor': 0.0001}

    _WlPTAtt = {'units' : 'm',
            'long_name' : "Mean Wave Length of Partition {}",
            'standard_name' : "partition_{}_mean"}

    _WlPTEnc = {'dtype' : 'f4',
            '_FillValue' : -9999.,
            'scale_factor': 1.0}

    _AttDict = {'time' : _TimeAtt,
            'longitude' : _LonAtt,
            'latitude'  : _latAtt,
            'Depth'   : _DepthAtt,
            'Dir'     : _DirAtt,
            'Dspr'    : _DsprAtt,
            'FSpr'    : _FsprAtt,
            'Hsig'    : _HsigAtt,
            'Hswell'  : _HswellAtt,
            'Lwavp'   : _LwavpAtt,
            'PkDir'   : _PkDirAtt,
            'Qp'      : _QpAtt ,
            'RTm01'   : _RTm01Att,
            'RTm_10'  : _RTm_10Att,
            'RTpeak'  : _RTpeakAtt,
            'Tm02'    : _Tm02Att,
            'TPsmoo'  : _TPsmooAtt,
            'Windv_x' : _Windv_xAtt,
            'Windv_y' : _Windv_yAtt,
            'Wlen'    : _WlenAtt,
            'DrPT'    : _DrPTAtt,
            'DsPT'    : _DsPTAtt,
            'HsPT'    : _HsPTAtt,
            'StPT'    : _StPTAtt,
            'TpPT'    : _TpPTAtt,
            'WfPT'    : _WfPTAtt,
            'WlPT'    : _WlPTAtt}

    _EncDict = {'time' : _TimeEnc,
            'longitude' : _LonEnc,
            'latitude'  : _LatEnc,
            'Depth'   : _DepthEnc,
            'Dir'     : _DirEnc,
            'Dspr'    : _DsprEnc,
            'FSpr'    : _FsprEnc,
            'Hsig'    : _HsigEnc,
            'Hswell'  : _HswellEnc,
            'Lwavp'   : _LwavpEnc,
            'PkDir'   : _PkDirEnc,
            'Qp'      : _QpEnc ,
            'RTm01'   : _RTm01Enc,
            'RTm_10'  : _RTm_10Enc,
            'RTpeak'  : _RTpeakEnc,
            'Tm02'    : _Tm02Enc,
            'TPsmoo'  : _TPsmooEnc,
            'Windv_x' : _Windv_xEnc,
            'Windv_y' : _Windv_yEnc,
            'Wlen'    : _WlenEnc,
            'DrPT'    : _DrPTEnc,
            'DsPT'    : _DsPTEnc,
            'HsPT'    : _HsPTEnc,
            'StPT'    : _StPTEnc,
            'TpPT'    : _TpPTEnc,
            'WfPT'    : _WfPTEnc,
            'WlPT'    : _WlPTEnc}


    def __init__(self):
        return


    def iocheck(self,fname):
        io = os.path.isfile(fname)
        if io:
            pass
        else:
            raise IOError('File {0} not found.'.format(fname))


    def _set_vars(self, mat):

        tvars = 0
        fvars = []

        for var in mat.keys():
            if var in self.__ivars_lst or var in self.__cvars_lst:
                continue
        
            if 'Time' in var: tvars+=1

            else: fvars.append(var[:-16])

            if tvars > 1:
                break

        self.vars = fvars

        return

    def _set_longitude(self, mat):

        self.longitude = mat['Xp'][0]
        return

    def _set_latitude(self, mat):

        self.latitude = mat['Yp'].T[0]
        return

    def _get_tmstp(self, mat):

        keys = mat.keys()
        keys = list(keys)
        keys.sort()
        x=re.compile('Time_')
        new_list = list(filter(x.match,keys))

        t_stamps = []

        for i in new_list:
            t_stamp = re.sub('Time_','',i)
            t_stamps.append(t_stamp)

        return t_stamps

    def _set_time(self, mat, warmup):

        tstamps = self._get_tmstp(mat)
        idx = 0

        if warmup:
            for idx, tstamp in enumerate(tstamps):
                if tstamp[4:] == '0101_000000' or tstamp[4:] == '0701_000000':
                    break

                print('Ignoring warm-up time {}'.format(tstamp))

        self.time = [dt.datetime.strptime(ts, '%Y%m%d_%H%M%S') for ts in tstamps[idx:]]
        self.start_time = self.time[0]
        self.end_time = self.time[-1]

        return

    def _gen_data_dict(self):

        ddict = od()

        for var in self.vars:
            ddict[var] = np.zeros((len(self.time),
                                  len(self.latitude),
                                  len(self.longitude)))

        return ddict

    def _set_coord_dict(self):

        self.coord_dict = od({'time' : self.time,
                            'latitude' : self.latitude,
                            'longitude' : self.longitude})

        return
        
    def _load_swn_data(self, swan_file, warmup = False):

        self.iocheck(swan_file)
        mat = scipy.io.loadmat(swan_file)

        self._set_longitude(mat)
        self._set_latitude(mat)
        self._set_vars(mat)
        self._set_time(mat, warmup)
        self._set_coord_dict()

        d_dict = self._gen_data_dict()

        idx = 0
        last_time = self.start_time

        for var in mat.keys():
            vname = var[:-16]
            vtime = var[-15:]

            if vname in self.vars \
            and self.start_time <= \
            dt.datetime.strptime(vtime, '%Y%m%d_%H%M%S'):

                cur_time = dt.datetime.strptime(vtime, '%Y%m%d_%H%M%S')

                if cur_time != last_time:
                    idx += 1
                    last_time = cur_time

                d_dict[vname][idx] = mat[var]

        self._set_attr_dict(d_dict.keys())
        self._set_enc_dict(d_dict.keys())
        d_dict = self._set_var_dict(d_dict)

        return d_dict

    def _set_var_dict(self, d_dict):

        var_dict = od()

        for key in d_dict.keys():

            tmp_var = xr.DataArray(
                    data = d_dict[key],
                    dims = ['time','latitude','longitude'],
                    coords = self.coord_dict,
                    attrs = self.attr_dict[key]
                    )

            tmp_var.encoding = self.enc_dict[key]

            tmp_var.longitude.attrs = self.attr_dict['longitude']
            tmp_var.latitude.attrs = self.attr_dict['latitude']
            tmp_var.time.attrs = self.attr_dict['time']
            tmp_var.longitude.encoding = self.enc_dict['longitude']
            tmp_var.latitude.encoding = self.enc_dict['latitude']
            tmp_var.time.encoding = self.enc_dict['time']

            var_dict[key] = tmp_var


        return var_dict

    def _set_enc_dict(self, dict_keys):
        '''Sets the encoding metadata for the variables and dimensions in self.enc_dict'''

        enc_dict = od()

        for key in list(dict_keys) \
                   + \
                   list(self.coord_dict.keys()):

            if 'PT' in key:
                vkey = key[:-2]
                tmp = self._EncDict[vkey].copy()

            else:
                tmp = self._EncDict[key].copy()

            enc_dict[key] = tmp

        self.enc_dict = enc_dict

    def _set_attr_dict(self, dict_keys):

        attr_dict = od()

        for key in list(dict_keys) \
                + \
                list(self.coord_dict.keys()):

            if 'PT' in key:
                pnum = key[-2:]
                vkey = key[:-2]
                tmp = self._AttDict[vkey].copy()
                tmp['long_name'] = tmp['long_name'].format(pnum)
                tmp['standard_name'] = tmp['standard_name'].format(pnum)
            else:
                tmp = self._AttDict[key].copy()

            attr_dict[key] = tmp

        self.attr_dict = attr_dict


    def _2Xarray(self, d_dict):

        ds = xr.Dataset(d_dict)

        return ds

    def readswan(self, f_in, warmup = False):    

        swandata = self._load_swn_data(f_in, warmup = warmup)
        swandset = self._2Xarray(swandata)

        return swandset


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(\
    description='''Converts swan matlab files into netcdf.'''\
    , usage='%(prog)s [options]')

    parser.add_argument('-o',  metavar='out_file', \
        help='output file: (optional, default: ./in_file.nc)')

    parser.add_argument('-i', metavar='in_file', required = True, \
        help='''swan mat file (required).''')

    parser.add_argument('-w', metavar='warmup', action = 'store_const',
        const = True, default = False, 
        help='''The file has warm up data starting prior to
        its first supposed output''')

    in_args = parser.parse_args()

    if in_args.o:
        fout = in_args.o

    else:
        fout = in_args.i.replace('.mat','.nc')

    sdata = SwnData()

    sdataset = sdata.readswan(in_args.i, in_args.w)
    sdataset.to_netcdf(fout, format = 'NETCDF4_CLASSIC')


#!/usr/bin/env python3
import scipy.io

class SwnData:

    __svars_lst = ['Depth',
        'Dir',
        'Dspr',
        'Fspr',
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

    _TimeAtt = {'units': "hours since 1990-01-01 00:00 UTC", # arbitrary
            'long_name': 'Time',
            'standard_name' : 'time',
            'type' : 'f8'}

    _LonAtt = {'units' : "degrees_east",
            'long_name' : "Longitude",
            'standard_name' : "longitude",
            'type' : 'f4'}

    _latAtt = {'units' : "degrees_north",
            'long_name' : "Latitude",
            'standard_name' : "latitude",
            'type' : 'f4'}

    _DepthAtt = {'units' : "m" # max: 3826.213, min: 0.062288705, 2^15 = 32768: too small
            'long_name' : "Depth",
            'standard_name' : "sea_floor_depth_below_sea_level",
            'type' : 'f4',
            'scale_factor': 1.0}

    _DirAtt = {'units': 'degrees_true',
            'long_name' : "Mean Wave Direction",
            'standard_name' : "sea_surface_wave_from_direction",
            'type' : 'f4',
            'scale_factor': 1.0}

    _DsprAtt = {'units' : 'degrees',
            'long_name' : "Directional Spreading",
            'standard_name' : "sea_surface_wave_directional_spreading",
            'type' : 'i2',
            'scale_factor': 100.0}

    _FsprAtt = {'units' : 'dimensionless', #TODO  Double Check 
            'long_name' : "Frequency Spectral Width (Kappa)",
            'standard_name' : "sea_surface_wave_normalized_width_of_frequency_spectrum",
            'type' : 'i2',
            'scale_factor': 10000.0}

    _HsigAtt = {'units' : 'm',
            'long_name' : "Significant Wave Height",
            'standard_name' : "sea_surface_wave_significant_height",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _HswellAtt = {'units' : 'm',
            'long_name' : "Wave Height of Swell Part",
            'standard_name' : "sea_surface_swell_wave_significant_height",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _LwavpAtt = {'units' : 'm',
            'long_name' : "Peak Wave Length",
            'standard_name' : "sea_surface_wave_length_peak",
            'type' : 'f4',
            'scale_factor': 1.0}

    _PkDirAtt = {'units' : 'degrees_true', 
            'long_name' : "Peak Wave Spectrum Direction",
            'standard_name' : "sea_surfave_wave_peak_direction",
            'type' : 'f4',
            'scale_factor': 1.0}

    _QpAtt = {'units':'dimensionless',
            'long_name' : "Peakedness of the Wave Spectrum",  
            'standard_name' : "sea_surface_wave_peakedness",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _RTm01Att = {'units':'s',
            'long_name' : "Mean Relative Wave Period",
            'standard_name' : "sea_surface_wave_mean_relative_period",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _RTm_10Att = {'units':'s',
            'long_name' : "Mean Relative Wave Period",
            'standard_name' : "sea_surface_wave_mean_relative_period_from_variance",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _RTpeakAtt = {'units':'s',
            'long_name' : "Peak Period of the Variance Density Spectrum",
            'standard_name' : "sea_surface_wave_relative_peak_period",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _Tm02Att = {'units':'s',
            'long_name' : "Zero Crossing Period of the Variance Density Spectrum",
            'standard_name' : "sea_surface_wave_zero_crossing_period",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _TPsmooAtt = {'units':'s',
            'long_name' : "Smoothed Peak Period of the Variance Density Spectrum",
            'standard_name' : "smoothed_sea_surface_wave_relative_peak_period",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _Windv_xAtt = {'units' : 'm s-1'
            'long_name' : "Eastward wind 10 m above ground",
            'standard_name' : "eastward_wind",
            'type' : 'f4',
            'scale_factor': 1.0}

    _Windv_yAtt = {'units' : 'm s-1'
            'long_name' : "Northward wind 10 m above ground",
            'standard_name' : "Northward_wind",
            'type' : 'f4',
            'scale_factor': 1.0}

    _WlenAtt = {'units' : 'm'
            'long_name' : "mean Wave Length of the Variance Density Spectrum",
            'standard_name' : "sea_surface_wave_mean_length",
            'type' : 'f4',
            'scale_factor': 1.0}

    _DrPTAtt = {'units': 'degrees_true',
            'long_name' : "Wave Direction of Partition {}",
            'standard_name' : "partition_{}_wave_from_direction",
            'type' : 'f4',
            'scale_factor': 1.0}

    _DsPTAtt = {'units': 'degrees',
            'long_name' : "Directional Spreading of Partition {}",
            'standard_name' : "partition_{}_wave_directional_spreading",
            'type' : 'i2',
            'scale_factor': 100.0}

    _HsPTAtt = {'units' : 'm',
            'long_name' : "Significant Wave Height of Partition {}",
            'standard_name' : "partition_{}_wave_significant_height",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _StPTAtt = {'units' : 'dimensionless', # TODO: double check
            'long_name' : "Wave steepness of Partition {}",
            'standard_name' : "partition_{}_wave_steepness",
            'type' : 'f4',
            'scale_factor': 1.0}

    _TpPTAtt = {'units' : 's',
            'long_name' : "Peak Period of Partition {}",
            'standard_name' : "partition_{}_wave_peak_period",
            'type' : 'i2',
            'scale_factor': 1000.0}

    _WfPTAtt = {'units' : 'dimensionless', # vmax = 1 #TODO: dcheck units
            'long_name' : "Wind Fraction of Partition {}",
            'standard_name' : "partition_{}_wind_fraction",
            'type' : 'i2',
            'scale_factor': 10000.0}

    _WlPTAtt = {'units' : 'm',
            'long_name' : "Mean Wave Length of Partition {}",
            'standard_name' : "partition_{}_mean",
            'type' : 'f4',
            'scale_factor': 1.0}

    def __init__(self,):


    def iocheck(self,fname):
        io = os.path.isfile(fname)
        if io:
            pass
        else:
            raise IOError('File {0} not found.'.format(fname))


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


    def _load_swn_data(self, swan_file, s_vars):

        self.iocheck(swan_file)
        mat = scipy.io.loadmat(swan_file)

        self._get_swnvr(mat)
        ts_list = self._get_tmstp(mat)

        t_beg, t_end = self._get_time_of_file(swan_file)

        if not s_vars.__contains__('dates'):
            s_vars['dates'] = np.asarray([])

            s_vars['lon'] = mat['Xp'][0]
            s_vars['lat'] = mat['Yp'].T[0]

# TODO: Consider the approach of creating a 3d array with np.ones
#       and just add each grid to it instead of appending. It should
#       be faster
        for t_stmp in ts_list:
            time = dt.datetime.strptime(t_stmp ,'%Y%m%d_%H%M%S' )

            if time >= t_beg:
                if not s_vars['dates'].__contains__(time):
                    s_vars['dates'] = np.append(s_vars['dates'],time)

                for s_var in self.swn_vars:
                    if not s_vars.__contains__(s_var):
                        # adding another dimension to stack up the data
                        if self.DEBUG:
                            tmp = np.asarray(mat[s_var+'_'+t_stmp][:20,:20],\
                                    dtype = np.float32)
                        else:
                            tmp = np.asarray(mat[s_var+'_'+t_stmp],\
                                    dtype = np.float32)

                        tmp.shape = (1, tmp.shape[0], tmp.shape[1])

                        s_vars[s_var] = tmp
                    else:
                        # adding another dimension to stack up the data
                       tmp = np.asarray(mat[s_var+'_'+t_stmp],\
                                dtype = np.float32)


                        tmp.shape = (1, tmp.shape[0], tmp.shape[1])

                        s_vars[s_var] = np.append(s_vars[s_var], tmp, axis = 0)

        return s_vars



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 16:59:12 2020

@author: romerowo
"""

from pylab import *
import pandas

from Star import Star
from Zodi import Zodi
from Observatory_Position_and_Time import Observatory_Position_and_Time
from Exozodi import Exozodi
from Detector import Detector
from Completeness import Completeness_MC_Calculator

class Orbit_Simulator():
    def __init__(self, HIP):
        print('Initializing Orbit_Simulator')
        
        '''
        Retrieve ExoCat entry for the input star Hipparcos ID (HIP)
        '''
        exocat = pandas.read_excel('ExoCat1.xls', header=0)
        i_entry = np.where(exocat['HIP'] == HIP)[0][0]
        
        exocat_star =  exocat.loc[i_entry]

        self.star_name = exocat_star['COMMON']
        print( self.star_name)

        '''
        Initizalize the target models. 
        '''
        self.star_model  = Star(L_bol = exocat_star['Lbol'],
                                d_pc = exocat_star['d(pc)'],
                                Temperature=exocat_star['Teff'],
                                Mass = exocat_star['M*(Msun)'],
                                MV = exocat_star['Mv'])
        self.zodi_model    = Zodi()
        self.exozodi_model = Exozodi()   
        
        '''
        Initizalize the observatory position and time calculator. 
        '''
        self.obs_PT = Observatory_Position_and_Time(exocat_star['RA(ICRS)'],
                                                    exocat_star['DE(ICRS)'])
                
        # Note: may want to consider not calling the function above during 
        # initialization
 
    #########################################################################
    #########################################################################
    #########################################################################

    def sample_orbits(self, # Set the range of orbital inclination angles. (0 deg is face-on, 90 deg is edge-on)
                      num_samples = 1000,       # number of instances to simulate
                      SNR_cut     = 7.,
                      inc_deg_min = None,       # for full range, = None
                      inc_deg_max = None,       # for full range, = None
                      nzod = 4.5,               # set the number of zodi strength for background.
                      T_int = 1.*24.*60.*60.,   # 1 day
                      wl_1_i = 615.e-9,         # lower bound for detection band
                      wl_2_i = 800.e-9,         # upper bound for detection band 
                      geometric_albedo=0.2,     # planet albedo
                      e2e_efficiency   = 0.035, # end-to-end efficiency for imaging
                      contrast = 4.e-11,        # end-to-end efficiency for imaging
                      reference_inner_orbital_radius_AU = 0.95,
                      reference_outer_orbital_radius_AU = 1.67,
                      min_planet_radius_km    = 0.8*6371., # nottom of range of planet radii
                      max_planet_radius_km    = 1.4*6371., # top of range of planet radii
                      planet_radius_km = 6371., #planet radius for Earth
                      IWA = 100.                # Inner working angle 
                      ):
        
        '''
        Initialize completeness calculator (this is used to extract orbits)
        '''
        CMC_calc      = Completeness_MC_Calculator()

        
        '''
        Define the imaging detector model
        '''
        
        imaging_detector_model  = Detector(Diameter=2.4, 
                                           throughput=e2e_efficiency, 
                                           quantum_efficiency=1., 
                                           contrast=contrast, 
                                           band=[wl_1_i, wl_2_i], 
                                           PSF_diameter_mas = 65.)
        '''
        Run the completeness calculator to sample the orbits.
        '''
        C = CMC_calc.Completeness(self.star_model,
                                  imaging_detector_model,
                                  self.zodi_model,
                                  self.exozodi_model,
                                  self.obs_PT,
                                  n_zodi=nzod,
                                  reference_inner_orbital_radius_AU = reference_inner_orbital_radius_AU,
                                  reference_outer_orbital_radius_AU = reference_outer_orbital_radius_AU,
                                  sampling_mode = 'HabEx', # 'SAG13', 'HabEx'
                                  min_planet_radius_km    = min_planet_radius_km,
                                  max_planet_radius_km    = max_planet_radius_km,
                                  planet_radius_km = planet_radius_km,
                                  geometric_albedo = geometric_albedo,
                                  T_int            = T_int, # seconds
                                  t_obs_days       = self.t_obs_array,
                                  IWA_lim_mas      = IWA, # we cycle through this in the loop below
                                  SNR_cut          = SNR_cut,# we cycle through this in the loop below
                                  num_samples      = num_samples,
                                  scale_radius     = True, # scale the orbit radius by star luminosity
                                  diagnostic_plots = False)
        
        '''
        This is the number of times the planet was detected for each instance
        '''
        hits = np.sum(CMC_calc.cut, axis=0)
        
        
        '''
        Compute the positions of the planet as seen in the plane of the sky.
        '''
        # Check that the point moves by a PSF
        projected_x_mas   = np.arcsin(CMC_calc.r_pos_AU[:,0,:] * CMC_calc.AU_in_meters / (self.star_model.d_pc * CMC_calc.parsec_in_meters)) / CMC_calc.mas_in_radians
        projected_y_mas   = np.arcsin(CMC_calc.r_pos_AU[:,1,:] * CMC_calc.AU_in_meters / (self.star_model.d_pc * CMC_calc.parsec_in_meters)) / CMC_calc.mas_in_radians
        
        HZ_AU = np.array([CMC_calc.r_inner_AU,CMC_calc.r_outer_AU])
        HZ_mas = np.arcsin(HZ_AU * CMC_calc.AU_in_meters / (self.star_model.d_pc * CMC_calc.parsec_in_meters)) / CMC_calc.mas_in_radians
        orb_radius_mas = np.arcsin(CMC_calc.r_orb_AU * CMC_calc.AU_in_meters / (self.star_model.d_pc * CMC_calc.parsec_in_meters)) / CMC_calc.mas_in_radians

        '''
        Save the parts that are used for orbit fitting
        '''
        self.SNR_cut      = SNR_cut 
        self.hits         = hits
        self.a_SM         = orb_radius_mas
        self.t_0          = CMC_calc.K_t_0_seconds/(24.*60.**2)
        self.Period       = CMC_calc.K_Period_seconds/(24.*60.**2)
        self.eccentricity = CMC_calc.K_eccentricity
        self.LAN          = CMC_calc.K_LAN
        self.inclination  = CMC_calc.K_inclination
        self.arg_peri     =  CMC_calc.K_arg_peri
        self.projected_x_mas = projected_x_mas 
        self.projected_y_mas = projected_y_mas 
        self.orbital_period_days = CMC_calc.K_Period_seconds/(24.*60.**2), # in days
        self.orb_radius_mas = orb_radius_mas
        self.HZ_mas = HZ_mas # habitable zone bounds
        self.SNR = CMC_calc.SNR



    #########################################################################
    #########################################################################
    #########################################################################
    
    '''
    Get the target observation times. This function finds that start of 
    each observing window if there are two observing windows per year. 
    If there is a single observing windows per year, the function will 
    identify two equally-spaced windows to target.
    '''

    def get_observation_time_array(self, 
                                   init_date = '2029/1/1', # format 'YYYY/MM/DD'
                                   duration_yrs=2):
        '''
        Initialize parameters
        '''
        self.init_date = init_date
        self.duration_days = np.int(duration_yrs*365.25)

        '''
        Establish target timeline
        The first step is to get the sun angles with respect to the 
        boresight of the telescope as a function of time.
        '''

        self.ang = self.obs_PT.get_sun_angles(initial_date = init_date, 
                                    num_days = self.duration_days)

        self.days = np.arange(0,len(self.ang))

        '''
        Set a cut for the solar keep-out angle of 54 degrees of RST and
        83 degrees of the Starshade
        '''

        self.obs_cut = np.logical_and(self.ang>54., self.ang<83.)

        '''
        Find the indices correspoding to the start and stop times of each 
        observing window
        '''
        se_indices = self.get_window_start_stop()
        print( 'Time Windows (%d): '%len(se_indices), se_indices )
        
        window_durations = []
        for k in range(len(se_indices)): 
            window_durations.append(se_indices[k][1]-se_indices[k][0])
        print( 'Window Durations (days): ', window_durations )
        #print( se_indices[0], se_indices[1])
        #print( se_indices[0][0], se_indices[0][1], se_indices[1][0], se_indices[1][1])
        
        '''
        Make an array of days where the observation is desired to be taken.
        There are two cases:
            1) when there are 2 observation windows per year
            2) when there is  1 observation window per year 
        '''
        t_obs_array = []
        if len(se_indices) == duration_yrs*2:
            for k in range(duration_yrs*2):
                t_obs_array.append( (se_indices[k][0]) )
            
        if len(se_indices) == duration_yrs*1:
            for k in range(duration_yrs*1):
                t_obs_array.append( se_indices[k][0] ) 
                t_obs_array.append( (se_indices[k][0] + se_indices[k][1])//2 ) 

        '''
        This is the case where there are the mission observing period
        starts in the middle of a window.
        '''
        if len(se_indices) == duration_yrs*2+1:
            for k in range(1, duration_yrs*2):
                t_obs_array.append( (se_indices[k][0]) )
            first_window_duration = se_indices[0][1]-se_indices[0][0]
            last_window_duration  = se_indices[-1][1]-se_indices[-1][0]
            if first_window_duration < last_window_duration:
                t_obs_array.append( se_indices[-1][0] ) 
            if first_window_duration >= last_window_duration:
                t_obs_array.append( se_indices[0][0] ) 

        if len(se_indices) == duration_yrs*1 + 1:
            for k in range(1, duration_yrs*1):
                t_obs_array.append( se_indices[k][0] ) 
                t_obs_array.append( (se_indices[k][0] + se_indices[k][1])//2 ) 
            first_window_duration = se_indices[0][1]-se_indices[0][0]
            last_window_duration  = se_indices[-1][1]-se_indices[-1][0]
            flag = False
            if first_window_duration < int(1.5*last_window_duration):
                t_obs_array.append( se_indices[-1][0] ) 
                t_obs_array.append( (se_indices[-1][0] + se_indices[-1][1])//2 ) 
                Flag = True
            if first_window_duration > int(1.5*last_window_duration):
                t_obs_array.append( se_indices[0][0] ) 
                t_obs_array.append( (se_indices[0][0] + se_indices[0][1])//2 ) 
                Flag = True
            if Flag==False:
                t_obs_array.append( se_indices[0][0] ) 
                t_obs_array.append( se_indices[-1][0] ) 

            t_obs_array = sort(t_obs_array)
                
        self.t_obs_array = np.array(t_obs_array)
        print( 'Observation Times: ', self.t_obs_array )


    #########################################################################
    #########################################################################
    #########################################################################

    def get_window_start_stop(self):
        windows = np.array(self.obs_cut).astype('float')
        diff_win = np.diff(windows)
        starts = np.where(diff_win>1.-1.e-10)[0] + 1 # this is because diff(a)[n] = a[n+1] - a[n]
        ends   = np.where(diff_win<-1.+1.e-10)[0] + 1 # this is to make the point inclusive a[n:m] ranges from a[m] to a[n-1]
        if self.obs_cut[0] == 1: 
            starts = np.concatenate([[0], starts])
            ends   = np.concatenate([ends, [self.days[-1]]])
        #print ('starts', starts)
        #print ('ends', ends)
        start_end_array=[]
        for k in range(0,len(starts)):
            start_end_array.append([starts[k], ends[k]])
        return start_end_array


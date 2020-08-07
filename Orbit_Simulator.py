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
    def __init__(self,
                 HIP,
                 duration_yrs = 2, # must be an integer number!
                 obs_density = 1, # number of visits per window
                 init_date = '2029/1/1' # format 'YYYY/MM/DD'
                ):
        
        '''
        Initialize parameters
        '''
        self.init_date = init_date
        self.duration_days = np.int(duration_yrs*365.25)

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
        self.obs_PT = Observatory_Position_and_Time(exocat['RA(ICRS)'][i_entry],
                                                    exocat['DE(ICRS)'][i_entry])

        self.CMC_calc      = Completeness_MC_Calculator()
        #detector_model = Detector(Diameter=2.4, throughput=0.2, quantum_efficiency=0.9, contrast=1.e-10, band=[wl_1, wl_2], PSF_diameter_mas = 65.)

        print('Initializing Orbit_Simulator')

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
        Get the start and end of each observing window.
        The start and stop dates are encoded in self.se_indices
        '''
        se_indices = self.get_window_start_stop()
        print( 'se_indices', se_indices, len(se_indices))
        print( se_indices[0], se_indices[1])
        print( se_indices[0][0], se_indices[0][1], se_indices[1][0], se_indices[1][1])
        
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
        print( 't_obs_array', self.t_obs_array )

    def get_window_start_stop(self):
        # returns the indices of where the observation windows start and end.
        windows = np.array(self.obs_cut).astype('float')
        diff_win = np.diff(windows)
        starts = np.where(diff_win>1.-1.e-10)[0] + 1 # this is because diff(a)[n] = a[n+1] - a[n]
        ends   = np.where(diff_win<-1.+1.e-10)[0] + 1 # this is to make the point inclusive a[n:m] ranges from a[m] to a[n-1]
        if self.obs_cut[0] == 1: 
            starts = np.concatenate([[0], starts])
            ends   = np.concatenate([ends, [self.days[-1]]])
        print ('starts', starts)
        print ('ends', ends)
        #print( starts, ends )
        #print( self.days[starts], self.days[ends] )
        #print( len(starts), len(ends) )
        # format into pairs
        start_end_array=[]
        for k in range(0,len(starts)):
            start_end_array.append([starts[k], ends[k]])
        #print( start_end_array )
        # will probably need some if statements in the case where the observing window starts or ends away from 2 yr window.
        #figure()
        #plot(days,obs_cut,'.')
        #plot(days,windows,'.')
        
        #plot(days[1:],np.diff(windows), 'o')
        #ylim(0.,1.1)
        return start_end_array





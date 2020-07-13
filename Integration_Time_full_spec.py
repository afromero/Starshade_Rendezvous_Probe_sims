#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 12:19:11 2018

@author: romerowo
"""

# Command line interface.
import matplotlib

matplotlib.use('Agg')

import argparse
parser=argparse.ArgumentParser(description='Integration Time Calculator')
parser.add_argument("-hip",    "--HIP",              help="Input Star HIP",          type=int)
parser.add_argument("-diam",   "--diameter",         help="Telescope Diameter",      default = [2.37], type=float)
parser.add_argument("-isc",    "--imaging_SNR_cut",  help="Imaging  SNR cut",        default = [5.,7.],  type=float)
parser.add_argument("-ssc",    "--spectral_SNR_cut", help="Spectral SNR cut",        default = [20.,20.], type=float)
parser.add_argument("-spcR",   "--spc_R",            help="Spectral Resolution",     default = [50., 50.], type=float)
parser.add_argument('-e2e_eff',"--e2e_efficiency", nargs='+', help='List of end to end efficiencies', default = [0.024, 0.035], type=float)
parser.add_argument('-contr',  "--inst_contrast",  nargs='+', help='List of instrument contrast', default = [1.e-10, 4.e-11], type=float)
parser.add_argument('-iwa',    "--inner_working_angle",  nargs='+', help='List of instrument contrast', default = [103.], type=float)
parser.add_argument("-tp",     "--Table_Path",       help="path for look-up tables", default = './', type=str)
parser.add_argument("-slw",     "--Slit_Width",      help="set slit width (in mas) for spectroscopy", default = 120., type=float)
parser.add_argument("-sg",     "--Solar_Glint",      help="switch to activate solar glint", default = 0., type=float)
parser.add_argument("-outdir", "--outdir",           help="output directry",         default='.', type=str)
parser.add_argument("-outtag", "--outtag",           help="output tag",              default='imgspc', type=str)
args=parser.parse_args()


# This script calculates the integration time required to achieve a desired SNR.
# This is done as a function of 
# 1. zodiacal dust brightness, 
# 2. planet radius (at fixed albedo),
# 3. zodiacal dust brightness, 
# 4. telescope throughput, 
# 5. contrast
# 6. Integration Time

# Create matrix of completeness values by parameters

from pylab import *
from scipy.special import erf
import time

'''
Load EXOCAT and Computational classes
'''
import pandas
if args.Table_Path[-1]=='/': args.Table_Path=args.Table_Path[:-1]
exocat = pandas.read_excel('%s/ExoCat1.xls'%args.Table_Path, header=0)
#gp = pandas.read_excel("%s/EKGPS.xls"%args.Table_Path)
print('args.Table_Path',args.Table_Path)
'''
Initialize Star model.
'''
from Star import Star
i_entry = np.where(exocat['HIP'] == args.HIP)[0][0]
exocat_star =  exocat.iloc[i_entry]

print('exocat_star[\'COMMON\']: ', exocat_star['COMMON'])

star_model = Star(L_bol = exocat_star['Lbol'], d_pc = exocat_star['d(pc)'], Temperature=exocat_star['Teff'], Mass = exocat_star['M*(Msun)'], MV = exocat_star['Mv'])

'''
Initialize solar system zodi model.
'''
from Zodi import Zodi
zodi_model = Zodi(args.Table_Path)

''' 
Initialize observatory model.
'''
from Observatory_Position_and_Time import Observatory_Position_and_Time
obs_PT = Observatory_Position_and_Time(exocat_star['RA(ICRS)'], exocat_star['DE(ICRS)'])

'''
Initialize exozodi model.
'''
from Exozodi import Exozodi
exozodi_model = Exozodi()

'''
Import detector model (initialize in calculation loop).
'''
from Detector import Detector

'''
Initialize completeness (Monte Carlo) calculator
'''
from Completeness import Completeness_MC_Calculator
CMC_calc = Completeness_MC_Calculator(LUT_path = args.Table_Path)

'''
Initialize Solar Glint model
'''
from Solar_Glint import Solar_Glint


print('Finished Importing Libraries')
######################################################
'''
Set ranges of input parameters
'''

# Diameter and psf are looped under the same index
diameter_array = np.array(args.diameter)
psf_mas_array  =  65. * 2.37 / diameter_array# scaling from WFIRST PSF FWHM 

# End-to-end efficiency of imaging and spectroscopy are looped under the same index
e2e_IMG_throughput_array = np.array(args.e2e_efficiency)
# IFS was descoped. This is now a slit prism spectrometer. 
# Transmission through the prism expected to be ~97% (TBC)
e2e_IFS_throughput_array =  0.97*np.array(args.e2e_efficiency)

# Each of the following has its own loop index.

contrast_array       = np.array(args.inst_contrast)

IWA_array            = np.array(args.inner_working_angle)

nzodi_array          = np.arange(0.,10.1, 0.5)

#Planet radii are Monte Carlo sampled. This is set for indexing purposes. 
planet_radius_array  = [6371.] # km

# A range of integration times are used to estimate the best one for imaging. 
T_int_array          = 10**np.arange(-0.4, 2.1, 0.2)*24.*60.*60. # in seconds

# ALL THESE GO TOGETHER WITH THE SAME INDEX
spc_R_array          = np.array(args.spc_R)
img_SNR_cuts         = np.array(args.imaging_SNR_cut)
spc_SNR_cuts         = np.array(args.spectral_SNR_cut)

'''
Initialize Completeness Matrix
'''
# loop order
# e2e_throughput, Contrast, IWA (defines detector)
# nzodi, planet radius (assumptions about nature)
# SNR_cuts, integration time
N_results = 7
Completeness_matrix = -1*np.ones((len(img_SNR_cuts), # coupled to spcR and spc_SNR_cuts
                                  len(diameter_array), # coupled to psf_mas
                                  len(e2e_IFS_throughput_array),
                                  len(contrast_array),
                                  len(IWA_array),
                                  len(nzodi_array),
                                  len(planet_radius_array),
                                  len(T_int_array),
                                  N_results # single-visit, at least 1, 2, 3, 4 successes
                                  ))
print('Completeness_matrix.shape',Completeness_matrix.shape)
print('')
print('=============================================')
print('')
print('spc_R_array', spc_R_array)
print('img_SNR_cuts', img_SNR_cuts) 
print('spc_SNR_cuts', spc_SNR_cuts)
print('')
print('diameter_array', diameter_array)
print('psf_mas_array', psf_mas_array)
print('')
print('e2e_IFS_throughput_array', e2e_IFS_throughput_array)
print('e2e_IMG_throughput_array', e2e_IMG_throughput_array)
print('')
print('contrast_array', contrast_array)
print('')
print('IWA_array', IWA_array)
print('')
print('nzodi_array', nzodi_array)
print('')
print('planet_radius_array', planet_radius_array)
print('')
print('T_int_array',  T_int_array)
print('')
print('=============================================')
print('')

'''
Get observing windows
'''
init_date = '2029/1/1'
duration_days = 2*366 # two years and change. The code expects this to be an integer.
obs_PT.get_observing_windows(init_date, duration_days, sun_ang_min=54., sun_ang_max=83.)
print('observing time windows and durations', obs_PT.starts, obs_PT.durations)
obs_time_array = obs_PT.starts
obs_time_durations = obs_PT.durations

'''
In the line below, usable windows have a minimum duration of 16 days. 
This to deal with windows that land within the start or end time of the 
observation campaign or with sources that are too close to the ecliptic 
poles and do not have observation windows that are long enough to make 
meaningful measurements.
'''
id_good = np.where(obs_time_durations>16)[0]
obs_time_array     = obs_time_array[id_good] 
obs_time_durations = obs_time_durations[id_good]
print('* observing time array', obs_time_array, obs_time_durations)
slit_factor = args.Slit_Width/65. # slit width divided by PSF
jitter_factor = (erf(args.Slit_Width/np.sqrt(2.)/(65.*0.55)))**1.3 # parameterization of monte carlo results with 13 mas RMS Gaussian jitter

'''
For sources with ecliptic latitudes above 54 deg, there is only one window
per year, although it is typically very long.
The loop below finds these cases splits each yearly window into two parts.
'''
while(len(obs_time_array)<4):
    idx = np.argmax(obs_time_durations)
    # print obs_time_array[idx], obs_time_durations[idx]
    obs_time_array = np.concatenate([obs_time_array[0:idx], 
                                      [obs_time_array[idx]], 
                                      [obs_time_array[idx]+obs_time_durations[idx]//2], 
                                      obs_time_array[idx+1:] ])
    obs_time_durations = np.concatenate([obs_time_durations[0:idx], 
                                         [obs_time_durations[idx]-obs_time_durations[idx]//2], 
                                         [obs_time_durations[idx]//2], 
                                         obs_time_durations[idx+1:]])
#print '*observing time array', obs_time_array, obs_time_durations

'''
Initialize Imaging and Reference Spectroscopy Bands
THESE WILL GET UPDATED WHEN WE LOOP THROUGH SPC_R
'''
img_wl_1 = 615.e-9 # imaging band lower frequency
img_wl_2 = 800.e-9 # imaging band upper frequency
spc_wl_1 = 725.e-9*(1. - 1./2./50.) # spectral sub-band lower frequency, R=50
spc_wl_2 = 725.e-9*(1. + 1./2./50.) # spectral sub-band upper frequency, R=50

'''
Loop Through Parameters
'''
count = 0

# note: imaging e2e througput is  6.9% (based on CGI requirements), spectral  4.5% (based on CGI requirements)
imaging_detector_model = Detector(Diameter=2.37, 
                                  throughput=0.045*(0.44/0.29),               
                                  quantum_efficiency=1.,  # throughput value includes this already.
                                  contrast=1.e-10,        # current requirement 
                                  band=[img_wl_1, img_wl_2], 
                                  PSF_diameter_mas = 65.)

spectrum_detector_model = Detector(Diameter=2.37,
                                   throughput=0.045,
                                   quantum_efficiency=1.,
                                   contrast=1.e-10,
                                   band=[spc_wl_1, spc_wl_2],
                                   PSF_diameter_mas = 65.)
 
SG = Solar_Glint(lower_wavelength = img_wl_1, upper_wavelength = img_wl_2, 
                 fnm1 = "%s/glint_620_800_53deg.txt"%args.Table_Path,
                 fnm2 = "%s/glint_620_800_63deg.txt"%args.Table_Path,
                 fnm3 = "%s/glint_620_800_73deg.txt"%args.Table_Path,
                 fnm4 = "%s/glint_620_800_83deg.txt"%args.Table_Path)

print('obs_time_array', obs_time_array)

'''
We run the Completeness calculator once to obtain the the planet orbit 
simulation part (should eventually make this part its own class so we don't 
have to do this). 
The instrument response is then corrected with the input values in the loops 
that follow. 
'''

print('Running Planet Monte Carlo')
C = CMC_calc.Completeness(
        star_model,
        imaging_detector_model,
        zodi_model,
        exozodi_model,
        obs_PT,
        n_zodi=4.5, # current best estimate
        reference_inner_orbital_radius_AU = 0.95,
        reference_outer_orbital_radius_AU = 1.67,
        sampling_mode = 'HabEx', # 'SAG13', 'HabEx'
        min_planet_radius_km    = 0.8*6371.,
        max_planet_radius_km    = 1.4*6371.,
        planet_radius_km = 6371., # Earth
        geometric_albedo = 0.20,
        #geometric_albedo_min = 0.025, # sets geom albedo to fixed value of 0.2.
        #geometric_albedo_max = 0.43,
        T_int            = 24.*60.**2, # seconds
        t_obs_days       = obs_time_array,
        IWA_lim_mas      = 100, # we cycle through this in the loop below
        SNR_cut          = 5.,# we cycle through this in the loop below
        num_samples      = 100000,
        scale_radius     = True, # scale the orbit radius by star luminosity
        diagnostic_plots = False)
npoints = float(CMC_calc.cut.shape[1])


# get the sun angles at the day of observation
# these do not change in the loop below. 
solar_lat, solar_lon = obs_PT.get_star_sun_lat_lon(initial_date = '2029/1/1', num_days=366*2)
ang = obs_PT.get_sun_angles(initial_date = '2029/1/1', num_days=365)
img_zodi_flux = np.zeros((len(obs_time_array), int(npoints)))
spc_zodi_flux = np.zeros((len(obs_time_array), int(npoints)))
#print len(obs_time_array), obs_time_array
start_time = time.time() # start timer

''' LOOP THROUGH START DAYS OF EACH OBSERVING WINDOW '''
print('Starting Detection Loop')
for k in range(0,len(obs_time_array)):
    print('k', k, 'of', len(obs_time_array))

    for i_spcR in range(0,len(spc_R_array)): # coupled with SNR_threshold
        print('\ti_spcR', i_spcR, 'of', len(spc_R_array))
        # SNR thresholds are now coupled to the spectral resolution R. 
        i_snr = i_spcR

        ''' GET THE REFERENCE SPECTRAL BAND AND FLUXES IN THAT BAND'''
        spc_wl_1 = 725.e-9*(1. - 1./2./spc_R_array[i_spcR]) # spectral sub-band lower frequency, R=50
        spc_wl_2 = 725.e-9*(1. + 1./2./spc_R_array[i_spcR]) # spectral sub-band upper frequency, R=50
        
        img_star_flux = star_model.get_Flux(img_wl_1, img_wl_2)        
        spc_star_flux = star_model.get_Flux(spc_wl_1, spc_wl_2)        

        img_planet_flux = CMC_calc.F_p_F_star * img_star_flux
        spc_planet_flux = CMC_calc.F_p_F_star * spc_star_flux

        phi_rad = np.arctan2(CMC_calc.r_pos_AU[k,1,:], CMC_calc.r_pos_AU[k,0,:])

        solar_glint_flux_img = []

        for ii in range(len(phi_rad)):
            solar_glint_flux_img.append(SG.get_flux(CMC_calc.projected_angle_mas[k,ii], phi_rad[ii], ang[obs_time_array[k]%365]))
        solar_glint_flux_img = np.array(solar_glint_flux_img)            
        solar_glint_flux_spc = solar_glint_flux_img*(spc_wl_2-spc_wl_1)/(img_wl_2-img_wl_1)

        for i_diam in range(0,len(diameter_array)): # coupled with PSF
            print('\t\ti_diam', i_diam, 'of', len(diameter_array))
            
            for i_thr in range(0,len(e2e_IFS_throughput_array)):
                print('\t\t\ti_thr', i_thr, 'of', len(e2e_IFS_throughput_array))
                
                for i_con in range(0,len(contrast_array)):

                    ''' GET THE IMAGING AND SPECTRAL DETECTOR MODELS '''
                    # assume detector QE = 1 since it is included in the end-to-end efficiency
                    imaging_detector_model = Detector(Diameter=diameter_array[i_diam], 
                                                      throughput=e2e_IMG_throughput_array[i_thr], # imaging has more throughput than IFS by this ratio 
                                                      quantum_efficiency=1., # set to 1 because this is already accounted for in the e2e efficiency
                                                      contrast=contrast_array[i_con], 
                                                      band=[img_wl_1, img_wl_2], 
                                                      PSF_diameter_mas = psf_mas_array[i_diam])
                    
                    spectrum_detector_model = Detector(Diameter=diameter_array[i_diam], 
                                                      throughput=e2e_IFS_throughput_array[i_thr], 
                                                      quantum_efficiency=1., # set to 1 because this is already accounted for in the e2e efficiency
                                                      contrast=contrast_array[i_con], 
                                                      band=[spc_wl_1, spc_wl_2], 
                                                      PSF_diameter_mas = psf_mas_array[i_diam])


                    ''' GET THE SOLAR SYSTEM ZODI FLUXES AT THE TIMES OF OBSERVATION'''
                    img_zodi_flux[k,:] = zodi_model.get_zodi_flux(imaging_detector_model,  solar_lat[obs_time_array[k]], solar_lon[obs_time_array[k]])
                    spc_zodi_flux[k,:] = zodi_model.get_zodi_flux(spectrum_detector_model, solar_lat[obs_time_array[k]], solar_lon[obs_time_array[k]])

                    ''' GET THE EXOZODI FLUXES'''
                    # Reference Exozodi levels
                    # these get rescaled later. 
                    # Done here because there is a relatively expensive calculation performed over a Planck distribution
                    n_zodi = 1.
                    img_ref_exozodi_flux = exozodi_model.scaled_exozodi(star_model, imaging_detector_model, n_zodi)
                    spc_ref_exozodi_flux = exozodi_model.scaled_exozodi(star_model, spectrum_detector_model, n_zodi)

                    for i_exz in range(0,len(nzodi_array)):

                        # rescale exozodi values
                        img_exozodi_flux = nzodi_array[i_exz] * img_ref_exozodi_flux
                        spc_exozodi_flux = nzodi_array[i_exz] * spc_ref_exozodi_flux
            
                        for i_rpl in range(0, len(planet_radius_array)): 

                            # rescale by planet radius here 
                            img_planet_flux = (planet_radius_array[i_rpl]/6371.)**2 * CMC_calc.F_p_F_star * img_star_flux
                            spc_planet_flux = (planet_radius_array[i_rpl]/6371.)**2 * CMC_calc.F_p_F_star * spc_star_flux

                            for i_tint in range(0,len(T_int_array)):
                                img_SNR = imaging_detector_model.SNR(img_planet_flux, 
                                                                     img_star_flux, 
                                                                     img_exozodi_flux, 
                                                                     img_zodi_flux,
                                                                     T_int_array[i_tint],
                                                                     solar_glint_flux = args.Solar_Glint*solar_glint_flux_img)
                                spc_SNR = spectrum_detector_model.SNR(jitter_factor * spc_planet_flux, 
                                                                      slit_factor * spc_star_flux, 
                                                                      slit_factor * spc_exozodi_flux, 
                                                                      slit_factor * spc_zodi_flux, 
                                                                      25.*24.*60.**2, # assume full 25 days
                                                                      solar_glint_flux = slit_factor*args.Solar_Glint*solar_glint_flux_img)
                                                                      #T_int_array[i_tint])
                                #img_SNR_SG = imaging_detector_model.SNR(img_planet_flux, 
                                #                                     img_star_flux, 
                                #                                     img_exozodi_flux, 
                                #                                     img_zodi_flux,
                                #                                     T_int_array[i_tint],
                                #                                     solar_glint_flux = 1.*solar_glint_flux_img)
                                
                                #spc_SNR_SG = spectrum_detector_model.SNR(spc_planet_flux, 
                                #                                      spc_star_flux, 
                                #                                      spc_exozodi_flux, 
                                #                                      spc_zodi_flux, 
                                #                                      25.*24.*60.**2, 
                                #                                      solar_glint_flux = 1.*solar_glint_flux_spc)
                                
            
                                for i_iwa in range(0,len(IWA_array)):
                                    
                                    # SNR is now coupled to spc_R, i_snr redefined above
                                    #for i_snr in range(0,len(img_SNR_cuts)):

                                    img_cut = np.logical_and(CMC_calc.projected_angle_mas>IWA_array[i_iwa], img_SNR>=img_SNR_cuts[i_snr])
                                    spc_cut = np.logical_and(CMC_calc.projected_angle_mas>IWA_array[i_iwa], spc_SNR>=spc_SNR_cuts[i_snr])

                                    # compute the various search completeness estimates
                                    img_C_single_visit = float(np.sum(img_cut[0,:]))/npoints
                                    img_C1_all_visits = float(np.sum(np.sum(img_cut, axis=0)>=1))/npoints
                                    img_C2_all_visits = float(np.sum(np.sum(img_cut, axis=0)>=2))/npoints
                                    img_C3_all_visits = float(np.sum(np.sum(img_cut, axis=0)>=3))/npoints
                                    img_C4_all_visits = float(np.sum(np.sum(img_cut, axis=0)>=4))/npoints
                                    spc_C1_all_visits = float(np.sum(np.sum(spc_cut, axis=0)>=1))/npoints
                                    
                                    num_dets_3 = np.sum(img_cut, axis=0)>=3
                                    num_spec_1 = np.sum(spc_cut, axis=0)>=1
                                    
                                    # combined completeness of at least 3 imaging detections and at least 1 spectral observation.
                                    C_al3d_1sp = np.logical_and(num_dets_3, num_spec_1)                             
                                    img_3obs_spc_1obs = float(np.sum(C_al3d_1sp))/npoints



                                    # put the estimates in the results array
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 0] = img_C_single_visit
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 1] = img_C1_all_visits
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 2] = img_C2_all_visits
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 3] = img_C3_all_visits
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 4] = img_C4_all_visits
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 5] = spc_C1_all_visits
                                    Completeness_matrix[i_snr, i_diam, i_thr, i_con, i_iwa, i_exz, i_rpl, i_tint, 6] = img_3obs_spc_1obs
                            
                                    count += 1

        #if(args.mode==1): tag = 'spectrum'
        np.savez('%s/star_%d_%s.npz'%(args.outdir, args.HIP, args.outtag),
                 solar_glint_level = args.Solar_Glint,
                 e2e_IFS_throughput_array = e2e_IFS_throughput_array, 
                 e2e_IMG_throughput_array = e2e_IMG_throughput_array, 
                 diameter_array           = diameter_array,
                 psf_mas_array            = psf_mas_array,
                 contrast_array           = contrast_array,
                 IWA_array                = IWA_array,
                 nzodi_array              = nzodi_array,
                 planet_radius_array      = planet_radius_array,
                 spc_R_array              = spc_R_array,
                 img_SNR_cuts             = img_SNR_cuts,
                 spc_SNR_cuts             = spc_SNR_cuts,
                 T_int_array              = T_int_array,
                 Completeness_matrix      = Completeness_matrix
                 )
            
                                

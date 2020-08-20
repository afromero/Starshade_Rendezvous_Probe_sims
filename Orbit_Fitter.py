# SET UP emcee
import emcee
from Kepler import Kepler
import numpy as np
import pickle
import sys

pi = np.pi

class Orbit_likelihood(object):
    # initialize the likelihood function estimator
    def __init__(self, KPC, time, data, data_uncertainty,
                 a_SM_min=0., 
                 a_SM_max=1000., 
                 t_0_min=0., 
                 t_0_max=1000., 
                 star_mass_Msun = 1.,
                 star_mass_Msun_fractional_uncertainty=0.1,
                 #Period_min = 0., 
                 #Period_max = 1000., 
                 eccentricity_min = 0., 
                 eccentricity_max = 1., 
                 LAN_min = 0., 
                 LAN_max = 2.*pi, 
                 inclination_min = 0., 
                 inclination_max = 2.*pi,
                 arg_peri_min = 0.,
                 arg_peri_max = 2.*pi,
                 distance_to_star_pc = 1.,
                 distance_to_star_pc_fractional_uncertainty=0.01
                ):
        print('Initializing Orbit_likelihood')

        self.KPC              = KPC
        self.time             = time
        self.data             = data
        self.data_uncertainty = data_uncertainty
        
        self.a_SM_min         = a_SM_min 
        self.a_SM_max         = a_SM_max 
        self.t_0_min          = t_0_min 
        self.t_0_max          = t_0_max 
        #self.Period_min       = Period_min 
        #self.Period_max       = Period_max 
        self.star_mass_Msun    = star_mass_Msun
        self.star_mass_Msun_uncertainty = self.star_mass_Msun * star_mass_Msun_fractional_uncertainty
        self.eccentricity_min = eccentricity_min 
        self.eccentricity_max = eccentricity_max 
        self.LAN_min          = LAN_min 
        self.LAN_max          = LAN_max 
        self.inclination_min  = inclination_min 
        self.inclination_max  = inclination_max
        self.arg_peri_min     = arg_peri_min
        self.arg_peri_max     = arg_peri_max
        self.distance_to_star_pc = distance_to_star_pc
        self.distance_to_star_pc_uncertainty = self.distance_to_star_pc*distance_to_star_pc_fractional_uncertainty

        # Conversion factors
        self.mAU_in_meters  = 1.49598e+8 # 1 AU is 1.49597e+11 meters. When convertain mas * pc = mAU = 10^-3 AU.
        self.G_Newton       = 6.67408e-11 # m3 kg-1 s-2
        self.Msun_in_kg     = 1.989e30 #kg
        self.day_in_seconds = 24.*60.**2 # 1 day in seconds.
    def Period(self, a_in_mas, d_pc, Mass_in_Msun):
        a_in_meters = a_in_mas * d_pc * self.mAU_in_meters
        Mass_in_kg = Mass_in_Msun * self.Msun_in_kg
        Period_in_seconds = 2.*pi * np.sqrt( a_in_meters**3 / ( self.G_Newton * Mass_in_kg) )
        Period_in_days = Period_in_seconds / self.day_in_seconds
        return Period_in_days
    
    def logprior(self, _parms, _Period_days):
        a_SM, t_0, Mass_Msun, eccentricity, LAN, inclination, arg_peri, d_pc = _parms
        pass_cond = a_SM>self.a_SM_min
        
        pass_cond = np.logical_and(a_SM<self.a_SM_max, pass_cond)

        pass_cond = np.logical_and(t_0>self.t_0_min, pass_cond)

        pass_cond = np.logical_and(t_0<self.t_0_max, pass_cond)
        
        #pass_cond = np.logical_and(Period>self.Period_min, pass_cond)
        #pass_cond = np.logical_and(Period<self.Period_max, pass_cond)
        #print('Period', Period, pass_cond)

        pass_cond = np.logical_and(eccentricity>self.eccentricity_min, pass_cond)
        pass_cond = np.logical_and(eccentricity<self.eccentricity_max, pass_cond)
        
        pass_cond = np.logical_and(LAN>self.LAN_min, pass_cond)
        pass_cond = np.logical_and(LAN<self.LAN_max, pass_cond)

        pass_cond = np.logical_and(inclination>self.inclination_min, pass_cond)
        pass_cond = np.logical_and(inclination<self.inclination_max, pass_cond)

        pass_cond = np.logical_and(arg_peri>self.arg_peri_min, pass_cond)
        pass_cond = np.logical_and(arg_peri<self.arg_peri_max, pass_cond)

        pass_cond = np.logical_and(t_0<_Period_days, pass_cond)

        if(pass_cond==True):
            distance_LL = -0.5*(d_pc      - self.distance_to_star_pc)**2 / self.distance_to_star_pc_uncertainty**2
            mass_LL     = -0.5*(Mass_Msun - self.star_mass_Msun)**2      / self.star_mass_Msun_uncertainty**2
            return distance_LL + mass_LL
        return -np.inf

    def loglhood(self,_Kepler_parms):
        data_x, data_y = self.data
        model_x, model_y = self.KPC.get_Orbit_points(self.time, *_Kepler_parms)
        model_x = model_x[:,0]
        model_y = model_y[:,0]
        unc_x, unc_y = self.data_uncertainty
        # currently not using time uncertainty (integration time)
        arr_x = (data_x - model_x)/unc_x
        arr_y = (data_y - model_y)/unc_y
        return -0.5*np.sum(arr_x**2 + arr_y**2)

    
    def __call__(self, _parms):
        a_SM, t_0, Mass_Msun, eccentricity, LAN, inclination, arg_peri, d_pc = _parms
        Period_days = self.Period(a_SM, d_pc, Mass_Msun)
        lprior = self.logprior(_parms, Period_days)
        # don't waste time calculating likelihoods for parameters that are not allowed.
        if( np.isinf(lprior)):
            return lprior

        # return the prior plus the calculated log likelihood
        if(np.isnan(lprior)):
            print('lprior is nan')
            return -np.inf
        
        Kepler_parms = [a_SM, t_0, Period_days, eccentricity, LAN, inclination, arg_peri]
        likelihood = self.loglhood(Kepler_parms)

        if(np.isnan(likelihood)):
            print('loglhood is nan')
            return -np.inf

        return  lprior + likelihood

#########################################################################################
#########################################################################################
#########################################################################################

#def dither_parms(_parms):
#    pos = [parms + 1e-2*np.random.uniform(-1.,1.,len(_parms)) for i in range(nwalkers)]
#    # Dither the initial positions
#    parms       = parms*(1. + 1.e-8*np.random.randn(len(parms)))
#    parms[3:7] += np.random.uniform(0.,1.e-8,4)
#    return parms

def Orbit_emcee(KPC, # an instance of the UHE_fluence calculator class
                time, data, data_uncertainty,
                initial_parm_vals = np.array([0.,0.,0.,0.,0.,0.,0.,0.]),
                ndim = 8,
                nwalkers = 100, # ndim * 10
                niterations = 1000, # number of steps in mcmc sampler
                interval = 1000, # the interval over which we save the data.
                out_tag = '',
                out_dir = './',
                a_SM_min=0., 
                a_SM_max=1000., 
                t_0_min=0., 
                t_0_max=1000., 
                star_mass_Msun = 1.,
                star_mass_Msun_fractional_uncertainty=0.1,
                #Period_min = 0., 
                #Period_max = 1000., 
                eccentricity_min = 0., 
                eccentricity_max = 1., 
                LAN_min = 0., 
                LAN_max = 2.*pi, 
                inclination_min = 0., 
                inclination_max = 2.*pi,
                arg_peri_min = 0.,
                arg_peri_max = 2.*pi,
                distance_to_star_pc = 1.,
                distance_to_star_pc_fractional_uncertainty=0.01
                ):
    # initialize the likelihood function
    logposterior = Orbit_likelihood(KPC, time, data, data_uncertainty,
                                    a_SM_min   = a_SM_min,
                                    a_SM_max   = a_SM_max,
                                    t_0_min    = t_0_min,
                                    t_0_max    = t_0_max,
                                    star_mass_Msun = star_mass_Msun,
                                    star_mass_Msun_fractional_uncertainty=star_mass_Msun_fractional_uncertainty,
                                    #Period_min = Period_min,
                                    #Period_max = Period_max,
                                    eccentricity_min = eccentricity_min,
                                    eccentricity_max = eccentricity_max,
                                    LAN_min          = LAN_min,
                                    LAN_max          = LAN_max,
                                    inclination_min  = inclination_min,
                                    inclination_max  = inclination_max,
                                    arg_peri_min     = arg_peri_min,
                                    arg_peri_max     = arg_peri_max,
                                    distance_to_star_pc = distance_to_star_pc,
                                    distance_to_star_pc_fractional_uncertainty=distance_to_star_pc_fractional_uncertainty
                                   )

    #print('initial_parm_vals', initial_parm_vals)
    #print('Initial Likelihood', logposterior(initial_parm_vals))

    # Set up the emcee Ensemble sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logposterior, threads=1)
    
    # Dither the initial positions
    pos = [initial_parm_vals + 1.e-3*np.random.uniform(-1.,1.,ndim) for i in range(nwalkers)]

    sys.stdout.flush()
    # Set up an interval of iterations over which to output the chains.
    N_loops = int(niterations / interval)
    print(interval, N_loops, N_loops*interval, niterations)

    for N in range(0,N_loops):
        print('N=', N, 'of', N_loops-1)
        sys.stdout.flush()
        # SAVING THE 
        #pickle.dump(sampler, open("uhe_emcee_sampler.p", "wb"))
        pos, prob, state = sampler.run_mcmc(pos, interval)
        '''
        out_chain  = sampler.chain[:,N*interval:(N+1)*interval,:]
        out_lnprob = sampler.lnprobability[:,N*interval:(N+1)*interval]

        #pickle.dump(sampler.chain, open("uhe_emcee_chain_%d.p"%N, "wb"))
        #pickle.dump(sampler.lnprobability, open("uhe_emcee_lnprobability_%d.p"%N, "wb"))
        print('\tsampler.chain.shape', sampler.chain.shape, sampler.lnprobability.shape)

        print('\tout_chain.shape', out_chain.shape, out_lnprob.shape)
        pickle.dump(out_chain, open("%s/uhe_emcee_chain_%s_%d.p"%(out_dir,out_tag,N), "wb"))
        pickle.dump(out_lnprob, open("%s/uhe_emcee_lnprob_%s_%d.p"%(out_dir,out_tag,N), "wb"))
        print('\t%1.2e'%(prob[0]))
        print('\t%1.2e\t%1.1f\t%1.1f\t%1.2e\t%1.2e\t%1.2e\t%1.2e'%(pos[0][0], pos[0][1], pos[0][2], pos[0][3], pos[0][4], pos[0][5], pos[0][6]))
        '''
        # Print out the mean acceptance fraction. In general, acceptance_fraction
        print("\tMean acceptance fraction: %1.2e"%np.mean(sampler.acceptance_fraction))
        # Estimate the integrated autocorrelation time for the time series in each
        # parameter.
        #try:
        #    print("\tAutocorrelation time:", sampler.get_autocorr_time())
        #except:
        #    print("Autocorrelation time estimate failed")
        sys.stdout.flush()
    return sampler

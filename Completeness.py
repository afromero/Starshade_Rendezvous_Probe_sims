from pylab import *
from scipy.interpolate import interp1d
from Kepler import Kepler

class Completeness_MC_Calculator:
    '''
    Estimate the Brown Completeness Via MC 
    Return a function with the following inputs:
        Stellar Properties: Star Mass, Bolometric Luminosity, Distance from Sun
        Planet Properties: Inner and Outer radial regions, planet radius, geometric albedo
    Ouptputs:
        Completeness
        Others to come
        
    Note: Characterization requires star will remain detectable or in view for the integration time required.
    '''
    def __init__(self, LUT_path ='./'):
        '''
        Constants
        '''
        self.AU_in_km = 1.496e+8
        self.mas_in_radians = 4.84814e-9
        self.Grav_constant = 6.67408e-11 #m^3 kg^-1 s^-2
        self.mass_of_Sun_in_kg = 1.989e30 #kg
        self.AU_in_meters     = 1.496e+11
        self.parsec_in_meters = 3.086e+16
        self.KPC = Kepler(force_new_LUT=False, LUT_path=LUT_path)
        self.year_in_seconds = 365.25*24.*60.*60.

    def Completeness(self,
                     star_model,
                     detector_model,
                     zodi_model,
                     exozodi_model,
                     obs_PT,
                     n_zodi = 1.,
                     reference_inner_orbital_radius_AU = 0.75, 
                     reference_outer_orbital_radius_AU = 1.77, 
                     planet_radius_km = 6371., 
                     geometric_albedo = 0.2,
                     geometric_albedo_min = None,
                     geometric_albedo_max = None,
                     T_int            = 3600., # seconds
                     t_obs_days       = np.array([0.]),
                     IWA_lim_mas      = 100.,
                     mag_lim          = 30.,
                     SNR_cut          = 5.,
                     num_samples      = 10000,
                     inclination_deg_min = None,
                     inclination_deg_max = None,
                     eccentricity_min = None,
                     eccentricity_max = None,
                     scale_radius = True, 
                     sampling_mode = 'log-uniform', # 'SAG13', 'HabEx'
                     min_planet_radius_km    = None,
                     max_planet_radius_km    = None,
                     diagnostic_plots = False):
        # Flow of this function
        # 1. Randomly sample planet orbital parameters. 
        # 2. For each instance get the planet position for each day in the observing schedule.
        # 3. For each sample, calculate the illumination phase angle and planet flux. zodi flux, and exozodi flux.
        # 4. Use the detector model to calculate SNR as a function of integration time.
        self.t_obs = np.array(t_obs_days).astype('float') * 24.*60.*60.
        #print 't_obs_days', t_obs_days
        # rescale inner and outer orbital radii to star luminosity
        self.r_inner_AU = reference_inner_orbital_radius_AU 
        self.r_outer_AU = reference_outer_orbital_radius_AU 
        if(scale_radius):
            self.r_inner_AU *= np.sqrt(star_model.L_bol)
            self.r_outer_AU *= np.sqrt(star_model.L_bol)
            
        # randomly sample keplerian parameters
        #a_SM, t_0, Period, eccentricity, LAN, inclination, arg_peri
        if sampling_mode == 'log-uniform':
            self.K_a_SM_AU        = self.get_random_orbit_radius(self.r_inner_AU, self.r_outer_AU, num_samples)
            self.K_Period_seconds = self.get_orbital_period(self.K_a_SM_AU, star_model.Mass) # Using Kepler's law
            self.planet_radius_km      = planet_radius_km * np.ones(num_samples) 
        if sampling_mode == 'SAG13':
            # SAG13 values. Gamma not used for sampling.
            Gamma = 0.38
            alpha = -0.19
            beta  = +0.26
            Period_min_yr = self.get_orbital_period(self.r_inner_AU, star_model.Mass)/self.year_in_seconds # Using Kepler's law
            Period_max_yr = self.get_orbital_period(self.r_outer_AU, star_model.Mass)/self.year_in_seconds # Using Kepler's law
            Period_yr = np.random.uniform(
                    Period_min_yr**(beta), 
                    Period_max_yr**(beta), 
                    num_samples)**(1./beta)
            self.K_Period_seconds = Period_yr * self.year_in_seconds
            self.K_a_SM_AU = self.get_orbital_semimajor_axis_AU(Period_yr, star_model.Mass)
            self.planet_radius_km = np.random.uniform(
                    (min_planet_radius_km/6371.)**(alpha),
                    (max_planet_radius_km/6371.)**(alpha),
                    num_samples)**(1./alpha) * 6371.
                    
        if sampling_mode == 'HabEx':
            # SAG13 values. Gamma not used for sampling.
            Gamma = 0.38
            alpha = -0.19
            beta  = +0.26
            Period_min_yr = self.get_orbital_period(self.r_inner_AU, star_model.Mass)/self.year_in_seconds # Using Kepler's law
            Period_max_yr = self.get_orbital_period(self.r_outer_AU, star_model.Mass)/self.year_in_seconds # Using Kepler's law
            Period_yr = np.random.uniform(
                    Period_min_yr**(beta), 
                    Period_max_yr**(beta), 
                    num_samples)**(1./beta)
            self.K_Period_seconds = Period_yr * self.year_in_seconds
            self.K_a_SM_AU = self.get_orbital_semimajor_axis_AU(Period_yr, star_model.Mass)
            self.planet_radius_km = np.random.uniform(
                    (min_planet_radius_km/6371./np.sqrt(self.K_a_SM_AU/np.sqrt(star_model.L_bol)))**(alpha),
                    (max_planet_radius_km/6371.)**(alpha),
                    num_samples)**(1./alpha) * 6371.


        self.K_t_0_seconds    = np.random.uniform(0., self.K_Period_seconds, num_samples)
        if eccentricity_min==None:
            self.K_eccentricity   = np.zeros(num_samples)
        if eccentricity_min!=None:
            self.K_eccentricity   = np.random.uniform(eccentricity_min, eccentricity_max, num_samples)
        self.K_LAN            = np.random.uniform(0.,2.*pi, num_samples)  
        if(inclination_deg_min==None):
            self.K_inclination    = np.arccos(np.random.uniform(-1.,1., num_samples))
        if(inclination_deg_min!=None):
            self.K_inclination    = np.arccos(np.random.uniform(np.cos(inclination_deg_max*pi/180.),
                                                 np.cos(inclination_deg_min*pi/180.), 
                                                 num_samples))
            
        self.K_arg_peri       = np.random.uniform(0.,2.*pi, num_samples)
        # get keplerian position vectors
        self.K_x = np.zeros((len(self.t_obs), num_samples))
        self.K_y = np.zeros((len(self.t_obs), num_samples))
        self.K_z = np.zeros((len(self.t_obs), num_samples))
        self.r_pos_AU = np.zeros((len(self.t_obs), 3, num_samples))
        for ii in range(0, num_samples):
            vals = self.KPC.get_Orbit_vectors(self.t_obs/(24.*60.**2),                # in days
                                              self.K_a_SM_AU[ii],
                                              self.K_t_0_seconds[ii]/(24.*60.**2),    # in days
                                              self.K_Period_seconds[ii]/(24.*60.**2), # in days
                                              self.K_eccentricity[ii],
                                              self.K_LAN[ii],
                                              self.K_inclination[ii],
                                              self.K_arg_peri[ii])
            #self.K_x[:,ii], self.K_y[:,ii], self.K_z[:,ii] = vals[:,0]
            self.K_x[:,ii] = vals[0,:,0]
            self.K_y[:,ii] = vals[1,:,0]
            self.K_z[:,ii] = vals[2,:,0]
            '''
            if ii==100: 
                print 'vals.shape',vals.shape
                print 'vals',vals
                print 'vals[:,0]',vals[:,0]
                print 'vals[0]',vals[0]
                print 'vals[1]',vals[1]
                print 'vals[2]',vals[2]
                print 'self.K_x[:,ii]', self.K_x[:,ii]
                print 'self.K_y[:,ii]', self.K_y[:,ii]
                print 'self.K_z[:,ii]', self.K_z[:,ii]
            '''
            self.r_pos_AU[:,0,ii] = self.K_x[:,ii]
            self.r_pos_AU[:,1,ii] = self.K_y[:,ii]
            self.r_pos_AU[:,2,ii] = self.K_z[:,ii]
        #print 'self.t_obs/(24.*60.**2)', self.t_obs/(24.*60.**2)
        #print 'self.r_pos_AU.shape', self.r_pos_AU.shape
        #self.r_orb_AU           = self.K_a_SM_AU.copy()
        self.r_orb_AU           = np.sqrt(  self.r_pos_AU[:,0,:]**2 
                                          + self.r_pos_AU[:,1,:]**2 
                                          + self.r_pos_AU[:,2,:]**2)
        '''
        for ii in range(0,10):
            print 'XxX', self.r_orb_AU[ii], np.sqrt(self.r_pos_AU[0,0,ii]**2+self.r_pos_AU[0,1,ii]**2+self.r_pos_AU[0,2,ii]**2)
            print 'XxX', self.r_orb_AU[ii], np.sqrt(self.r_pos_AU[1,0,ii]**2+self.r_pos_AU[1,1,ii]**2+self.r_pos_AU[1,2,ii]**2)
            print 'XxX', self.r_orb_AU[ii], np.sqrt(self.r_pos_AU[2,0,ii]**2+self.r_pos_AU[2,1,ii]**2+self.r_pos_AU[2,2,ii]**2)
            print 'XxX', self.r_orb_AU[ii], np.sqrt(self.r_pos_AU[3,0,ii]**2+self.r_pos_AU[3,1,ii]**2+self.r_pos_AU[3,2,ii]**2)
        #self.orb_period_seconds = self.K_Period_seconds.copy()
        #self.inclination        = self.K_inclination.copy()
        '''
        '''
        self.r_orb_AU   = self.get_random_orbit_radius(self.r_inner_AU, self.r_outer_AU, num_samples)
        
        
        # Get orbital periods
        self.orb_period_seconds = self.get_orbital_period(self.r_orb_AU, star_model.Mass) # Using Kepler's law

        # Sample orbital planes
        self.rot_x_vec_orb, self.rot_y_vec_orb, self.rot_z_vec_orb = self.get_random_orbital_planes(num_samples, 
                                                                                     inclination_deg_min = inclination_deg_min, 
                                                                                     inclination_deg_max = inclination_deg_max)
        self.orb_phase = np.random.uniform(0.,2.*pi, num_samples)

        # Get positions for sr_pos_AUampled exoplanets, orbits, and observations
        self.r_pos_AU = self.get_planet_positions(self.t_obs, self.orb_period_seconds, self.orb_phase, self.r_orb_AU, self.rot_x_vec_orb, self.rot_y_vec_orb)
        '''
        #print 'self.r_pos_AU.shape', self.r_pos_AU.shape
        #print 'r_pos_AU.shape', r_pos_AU.shape

        self.geometric_albedo = geometric_albedo*np.ones(num_samples)
        if geometric_albedo_min!=None:
            self.geometric_albedo = np.random.uniform(geometric_albedo_min, geometric_albedo_max, num_samples)
        #print r_pos_AU.shape
        # Get illumination phase angles
        self.beta = self.get_beta(self.r_pos_AU)
        #print 'self.beta.shape', self.beta.shape

        #dmag = self.get_dmag(geometric_albedo, self.beta,  planet_radius_km, self.r_orb_AU)
        self.projected_distance_AU = self.r_orb_AU * np.sin(self.beta)
        #self.projected_angle_mas   = np.arcsin(self.projected_distance_AU * self.AU_in_meters / (star_distance_pc * self.parsec_in_meters)) / self.mas_in_radians
        self.projected_angle_mas   = np.arcsin(self.projected_distance_AU * self.AU_in_meters / (star_model.d_pc * self.parsec_in_meters)) / self.mas_in_radians
        #self.mag = dmag + star_model.MV        
        self.F_p_F_star = self.get_planet_flux_ratio(self.geometric_albedo, self.beta,  planet_radius_km, self.r_orb_AU)

        self.star_flux = star_model.get_Flux(detector_model.band[0], detector_model.band[1])        
        self.planet_flux = self.F_p_F_star*self.star_flux
        self.exozodi_flux = exozodi_model.scaled_exozodi(star_model, detector_model, n_zodi)
        
        # get the worst case zodi NEED TO UPDATE THIS TO USE THE ZODI AT OBSERVATION.
        solar_lat, solar_lon = obs_PT.get_star_sun_lat_lon(initial_date = '2029/1/1', num_days=365)
        ang = obs_PT.get_sun_angles(initial_date = '2029/1/1', num_days=365)
        min_idx = np.argmin(np.abs(ang-54.))
        self.zodi_flux = zodi_model.get_zodi_flux(detector_model, solar_lat[min_idx], solar_lon[min_idx])
        #print 'min_idx, solar_lat[min_idx], solar_lon[min_idx], ang[min_idx], zodi_flux',min_idx, solar_lat[min_idx], solar_lon[min_idx], ang[min_idx], zodi_flux

        SNR = detector_model.SNR(self.planet_flux, self.star_flux, self.exozodi_flux, self.zodi_flux, T_int)
        SNR_shape = SNR.shape
        SNR_ravel_cut = np.ravel(SNR)
        SNR_ravel_cut[np.ravel(self.projected_angle_mas)<IWA_lim_mas] = 1.e-9
        SNR = SNR_ravel_cut.reshape(SNR_shape)
        #print 'CMC, SNR.shape, self.projected_angle_mas.shape',SNR.shape, self.projected_angle_mas.shape
        #for k in range(0,len(SNR)):
        #    (SNR[k])[self.projected_angle_mas[k]>IWA_lim_mas] = 1.e-9
        
        self.SNR = SNR.copy()
        #print SNR
        #cut = np.logical_and(self.projected_angle_mas[0,:]>IWA_lim_mas, self.mag[0,:]<mag_lim)
        self.cut = np.logical_and(self.projected_angle_mas>IWA_lim_mas, SNR>=SNR_cut)
        completeness = float(np.sum(self.cut[0,:]))/(len(self.cut[0,:]))
        
        if (diagnostic_plots):
            figure(1)
            semilogy(self.beta[0,:]*180./pi, SNR[0,:], 'k.')
            grid(True)
            ylim(0.1, 10**np.ceil(np.log10(np.max(SNR[0,:]))))
            ylabel('SNR')
            xlabel(r'Illumination Phase Angle $\beta$, deg')
            
            figure(2)
            ax = subplot(111)
            norm = float(np.sum( self.projected_angle_mas[0,:]>IWA_lim_mas )) / float(len(self.projected_angle_mas[0,:]))
            u = np.linspace(0.,norm, len(SNR[0,:]))
            semilogx(sorted(SNR[0,:]), u[::-1], '.')
            plot([SNR_cut, SNR_cut], [0.,1.], 'r--')
            x0, x1 = ax.get_xlim()
            yticks(np.arange(0.,1.05,0.1))
            plot([x0,x1], [completeness, completeness], 'r--')
            grid(True)
            
            print ('Fluxes (W/m^2)')
            print ('max planet_flux %1.1e'%np.max(self.planet_flux[0,:][self.cut]))
            print ('exozodi_flux    %1.1e'%self.exozodi_flux )
            print ('zodi_flux       %1.1e'%self.zodi_flux )
            print ('star_flux       %1.1e'%(self.star_flux*detector_model.contrast))

            print ('\nElectrons')
            print ('max planet      %d'%detector_model.get_photons(np.max(self.planet_flux[0,:][self.cut]), T_int))
            print ('exozodi         %d'%detector_model.get_photons(self.exozodi_flux, T_int) )
            print ('zodi            %d'%detector_model.get_photons(self.zodi_flux, T_int))
            print ('star leakage    %d'%detector_model.get_photons(self.star_flux*detector_model.contrast, T_int))
            print ('dark_counts     %d'%int(np.ceil(detector_model.dark_current_counts(T_int))))
            print ('CIC_counts      %d'%int(np.ceil(detector_model.clock_induced_charge_counts(T_int))))

            
        return completeness
        
    def get_random_orbit_radius(self, r_IHZ, r_OHZ, num_samples):
        u_IHZ = np.log(r_IHZ)
        u_OHZ = np.log(r_OHZ)
        u = np.random.uniform(u_IHZ, u_OHZ, num_samples)
        r = np.exp(u)
        return r

    def get_orbital_period(self, orbital_radius_AU, star_mass_Msun):
        # returns period in seconds
        star_mass_kg = star_mass_Msun * self.mass_of_Sun_in_kg
        orbital_radius_m = orbital_radius_AU * self.AU_in_meters
        period = 2.*pi * np.sqrt(orbital_radius_m**3 / (self.Grav_constant * star_mass_kg))
        return period

    def get_orbital_semimajor_axis_AU(self, orbital_Period_years, star_mass_Msun):
        # returns the semimajor axis in AU
        star_mass_kg = star_mass_Msun * self.mass_of_Sun_in_kg
        orbital_Period_s = orbital_Period_years * 365.25 * 24. * 60. * 60.
        orbital_radius_m = ( self.Grav_constant * star_mass_kg * (orbital_Period_s / (2.*pi) )**2 )**(1./3.)
        #period = 2.*pi * np.sqrt(orbital_radius_m**3 / (self.Grav_constant * star_mass_kg))
        return orbital_radius_m / self.AU_in_meters


    def get_random_orbital_planes(self,num_samples, inclination_deg_min=None, inclination_deg_max=None):
        ph_orb = np.random.uniform(0.,2.*pi, num_samples)
        if(inclination_deg_min==None):
            th_orb = np.arccos(np.random.uniform(-1.,1., num_samples))
        if(inclination_deg_min!=None):
            th_orb = np.arccos(np.random.uniform(np.cos(inclination_deg_max*pi/180.),
                                                 np.cos(inclination_deg_min*pi/180.), 
                                                 num_samples))
            #th_orb = inclination_deg*pi/180.*np.ones(num_samples)
        # rotated orbit coordinates
        r_hat_orb  = np.array([np.sin(th_orb) * np.cos(ph_orb), np.sin(th_orb)* np.sin(ph_orb),  np.cos(th_orb)])
        th_hat_orb = np.array([np.cos(th_orb) * np.cos(ph_orb), np.cos(th_orb)* np.sin(ph_orb), -np.sin(th_orb)])
        ph_hat_orb = np.array([-np.sin(ph_orb)                , np.cos(ph_orb)                , np.zeros(len(th_orb))])

        self.inclination = th_orb
        self.orbit_rot = ph_orb
        # check orthogonality
        #print r_hat_orb.shape, th_hat_orb.shape, ph_hat_orb.shape
        #print r_hat_orb[0,:]*th_hat_orb[0,:]  + r_hat_orb[1,:]*th_hat_orb[1,:]  + r_hat_orb[2,:]*th_hat_orb[2,:]
        #print r_hat_orb[0,:]*ph_hat_orb[0,:]  + r_hat_orb[1,:]*ph_hat_orb[1,:]  + r_hat_orb[2,:]*ph_hat_orb[2,:]
        #print th_hat_orb[0,:]*ph_hat_orb[0,:] + th_hat_orb[1,:]*ph_hat_orb[1,:] + th_hat_orb[2,:]*ph_hat_orb[2,:]
        #print 'th_orb, r_hat_orb[2,:]', th_orb*180./pi, r_hat_orb[2,:], np.arccos(r_hat_orb[2,:])*180./pi
        return th_hat_orb, ph_hat_orb, r_hat_orb # orbit "x", "y", and "z" directions


    def get_planet_positions(self, t_obs, orbital_period, orbit_phase, orbit_radius, orbit_x, orbit_y):

        planet_pos = np.zeros((len(t_obs), orbit_x.shape[0], orbit_x.shape[1]))

        for i_t_obs in range(0,len(t_obs)):
            planet_pos[i_t_obs,:,:] = np.cos(2.*pi*t_obs[i_t_obs]/orbital_period + orbit_phase) * orbit_x  \
                                    + np.sin(2.*pi*t_obs[i_t_obs]/orbital_period + orbit_phase) * orbit_y 
            planet_pos[i_t_obs,:,:] *= orbit_radius
        return planet_pos

    def get_beta(self, vec_planet_pos):
        # assume observation is along the z-axis
        cos_beta = vec_planet_pos[:,2,:]/np.sqrt(vec_planet_pos[:,0,:]**2 +vec_planet_pos[:,1,:]**2  + vec_planet_pos[:,2,:]**2 )
        return arccos(cos_beta)
    
    def Lambertian_Phase(self, beta):
        # range of beta is 0., pi
        return (np.sin(beta) + (np.pi-beta)*np.cos(beta)) / np.pi

    def get_planet_flux_ratio(self, geom_albedo, beta,  planet_radius_km, star_planet_radius_AU):
        Phi_L = self.Lambertian_Phase(beta)
        refl = np.zeros(beta.shape)
        star_planet_radius_km = star_planet_radius_AU * self.AU_in_km
        for k in range(0,len(refl)):
            refl[k] = geom_albedo * Phi_L[k,:] * ( planet_radius_km / star_planet_radius_km[k,:])**2 
        return refl

    def get_dmag(self, geom_albedo, beta,  planet_radius_km, star_planet_radius_AU):
        Phi_L = self.Lambertian_Phase(beta)
        dmag = np.zeros(beta.shape)
        star_planet_radius_km = star_planet_radius_AU * self.AU_in_km
        for k in range(0,len(dmag)):
            dmag[k] = -2.5 * np.log10( geom_albedo * Phi_L[k,:] * ( planet_radius_km / star_planet_radius_km)**2 )
        return dmag



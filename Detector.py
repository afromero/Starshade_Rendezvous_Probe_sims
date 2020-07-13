import numpy as np
class Detector:
    '''
    NOTES: 
    '''
    def __init__(self, Diameter=2.4, throughput=0.2, quantum_efficiency=0.9, contrast=1.e-10, band=[615.e-9, 800.e-9], PSF_diameter_mas=65.):
        '''
        Initialize detector parameters
        '''
        self.band = band
        self.A_tau_eta = np.pi*(Diameter/2.)**2 * throughput * quantum_efficiency
        lam_c = np.mean(self.band) # meters
        hc = 1.9864e-25 # Watts * seconds / meter
        self.photon_energy = hc/lam_c
        #print 'photon_energy, Watt*s', photon_energy
        self.contrast = contrast
        self.PSF_diameter_mas = PSF_diameter_mas

        self.mas_in_radians = 4.84813681e-9
        self.steradian_in_arcsec2 = 4.25e10 #arcsec2
        self.d_Omega_PSF_sr = 2.*np.pi * ( 1. - np.cos(self.PSF_diameter_mas*self.mas_in_radians/2.))
        self.d_Omega_PSF_arcsec2 = self.d_Omega_PSF_sr * self.steradian_in_arcsec2
        
        # detector noise parameters
        self.t_frame = 1000. # seconds
        self.i_dark = 7.e-4  # electrons/second/pix
        self.m_pix = 4.      # number of pixels encompassing the PSF (from Nemati, 2014)
        self.q_CIC = 0.02    # electrons / pixel / frame
        self.read_noise = 0. # from https://wfirst.ipac.caltech.edu/sims/Param_db.html, asssume it is negligible
        
    def get_photons(self, flux, T_int):
        return flux * T_int * self.A_tau_eta / self.photon_energy
    
    def dark_current_counts(self,T_int):
        return self.m_pix * self.i_dark * T_int

    def clock_induced_charge_counts(self,T_int):
        return self.m_pix * self.q_CIC  * np.ceil(T_int / self.t_frame)
    
    def get_instrument_noise_counts(self, T_int):
        counts   = 0
        counts += self.dark_current_counts(T_int)
        counts += self.clock_induced_charge_counts(T_int)
        return counts
        
    def SNR(self, planet_flux, star_flux, exozodi_flux, zodi_flux, T_int, solar_glint_flux = []):
        star_leakage        = star_flux * self.contrast
        planet_photons      = self.get_photons(planet_flux, T_int)
        star_leakge_photons = self.get_photons(star_leakage, T_int)
        exozodi_photons     = self.get_photons(exozodi_flux, T_int)
        zodi_photons        = self.get_photons(zodi_flux, T_int)
        
        instrument_counts   = self.get_instrument_noise_counts(T_int)
        Signal_photons      = planet_photons
        Noise_Sq_photons    = planet_photons + star_leakge_photons + zodi_photons + exozodi_photons + instrument_counts

        #if solar_glint_flux != None:
        if len(solar_glint_flux) > 0:
            solar_glint_photons = self.get_photons(solar_glint_flux, T_int)
            Noise_Sq_photons   += solar_glint_photons
        
        return Signal_photons / np.sqrt(Noise_Sq_photons)
        

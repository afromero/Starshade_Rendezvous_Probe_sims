from pylab import *
import matplotlib.colors as colors
rcParams['font.size']=15


class Exozodi:
    def __init__(self):
        self.zodi_unit = 22. # mag arcsec^-2 (Absolute V-band)
        self.V_band_zero_point = 3.917e-8 # Watts m^-2 um^-1
        self.zodi_absolute_brightness = 10**(-0.4*self.zodi_unit) * self.V_band_zero_point # W m^-2 um^-1 arcsec^-2 at 10 pc
        self.MV_sun = 4.83 # absolute magnitude of Sun
        self.V_band_wavelength = 550.e-9
        self.Teff_Sun = 5777.
        
    def scaled_exozodi(self, star_model, detector_model, n_zodi):
        
        # set to zodi
        sc_zod = self.zodi_absolute_brightness
        
        # DOES NOT scale to distance, it's a specific intensity, not a flux !
        
        # scale to star intensity
        sc_zod *= 10.**(-0.4 * (star_model.MV - self.MV_sun) )
        #print '10.**(-0.4 * (star_model.MV - self.MV_sun) )', 10.**(-0.4 * (star_model.MV - self.MV_sun) )
        # scale to habitable zone (1/r^2 Earth-equivalent distance) 
        # NOTE: This really should use the radius of the planet scaled to 1 AU.
        sc_zod *= 1./star_model.L_bol
        #print '1./star_model.L_bol', 1./star_model.L_bol
        
        # scale to the wavelength of observation and temperature(value is for V-band)
        # for the equation below, you already accounted for the temperature difference of the Sun and the star.
        #sc_zod *= star_model.planck_val(np.mean(detector_model.band), star_model.Temperature) / star_model.planck_val(self.V_band_wavelength, self.Teff_Sun)
        # only shift within the spectrum of the star.
        sc_zod *= star_model.planck_val(np.mean(detector_model.band), star_model.Temperature) / star_model.planck_val(self.V_band_wavelength, star_model.Temperature)

        #print 'star_model.planck_val(np.mean(detector_model.band), star_model.Temperature) / star_model.planck_val(self.V_band_wavelength, star_model.Temperature)', star_model.planck_val(np.mean(detector_model.band), star_model.Temperature) / star_model.planck_val(self.V_band_wavelength, star_model.Temperature)
        
        band_0 = detector_model.band[0]
        band_1 = detector_model.band[1]
        #sc_zod *= star_model.get_Flux(band_0, band_1) / star_model.planck_val(np.mean(detector_model.band), star_model.Temperature)
        
        # the long way
        #wl_array = np.linspace(band_0, band_1, 10000)
        #intg = sum(star_model.planck_val(wl_array, star_model.Temperature))*(wl_array[1]-wl_array[0])/star_model.planck_val(np.mean(detector_model.band), star_model.Temperature) * 1.e6
  
        sc_zod *= (band_1 - band_0) * 1.e6 # This is now in units of W/m^2/arcsec^2
    
        #print '(band_1 - band_0) * 1.e6', (band_1 - band_0) * 1.e6
        sc_zod *= detector_model.d_Omega_PSF_arcsec2 # this is in units of W/m^2
        
        return sc_zod * n_zodi
        
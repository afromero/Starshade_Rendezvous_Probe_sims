from pylab import *
from scipy.interpolate import interp1d

stefan_boltzmann_constant = 5.670367e-8 # Watts m^-2 Kelvin^-4
boltzmann_constant        = 1.38064852e-23  # Watts / Hz / K-1 
planck_constant           = 6.62607004e-34 #m^2 kg / s
speed_of_light            = 299792458. # m / s
L_Sun_Watts               = 3.846e26 # Watts
pc_in_meters              = 3.086e+16 # meters

def planck_law(wavelength_meters, temperature_Kelvin):
    val  =  2.*planck_constant*speed_of_light**2 / wavelength_meters**5
    exp_arg = planck_constant*speed_of_light / (wavelength_meters*boltzmann_constant*temperature_Kelvin)
    val  /= np.exp( exp_arg ) - 1.
    return val # in Watts / sr / m^2 / m

def Stefan_Boltzmann(temperature_kelvin):
    return stefan_boltzmann_constant * temperature_kelvin**4 # Watts / meter^2


class Star:
    '''
    Class for stellar modeling
    '''
    
    def __init__(self, L_bol=None, d_pc=None, Temperature=None, Mass=None, MV=None, exocat_init=None):
        '''
        INITIALIZE VALUES
        '''
        if(L_bol!=None):
            self.L_bol       = L_bol
            self.d_pc        = d_pc
            self.Temperature = Temperature

            self.L_bol_Watts = self.L_bol * L_Sun_Watts
            self.d_meters    = self.d_pc*pc_in_meters
            self.Flux_bol    = self.L_bol_Watts / 4. / pi / self.d_meters**2 # Watts / m^2
            self.Mass        = Mass
            self.MV          = MV

    def get_band_fractional_irradiance(self, wavelength_meters_min, wavelength_meters_max):
        wl_array = np.linspace(wavelength_meters_min, wavelength_meters_max, 10000)
        intg = pi*sum(planck_law(wl_array, self.Temperature))*(wl_array[1]-wl_array[0])
        return intg/Stefan_Boltzmann(self.Temperature)
    
    def get_Flux(self, wavelength_meters_min, wavelength_meters_max):
        fractional_irradiance = self.get_band_fractional_irradiance(wavelength_meters_min, wavelength_meters_max)
        #print self.Flux_bol, fractional_irradiance, fractional_irradiance*self.Flux_bol
        return fractional_irradiance*self.Flux_bol
    
    def planck_val(self, wavelength_meters, temperature_Kelvin):
        return planck_law(wavelength_meters, temperature_Kelvin)

'''
# check against values and plots in wikipedia entry for planck's law
temp = 10000.

figure()
wl = 10**np.arange(-7., -5.5, 0.01)
#print planck_law(wl, 5000.)
plot(wl*1.e6, planck_law(wl, temp))

figure()
# now check for some relevant values for Starshade Rendezvous
# "optical" band
wl = np.linspace(1.e-8, 1.e-4, 1000)
#print planck_law(wl, 5000.)
semilogx(wl*1.e6, planck_law(wl, temp))
figure()
semilogx(wl*1.e6, cumsum(planck_law(wl, temp)))

print pi*sum(planck_law(wl, temp))*(wl[1]-wl[0]), Stefan_Boltzmann(temp), pi*sum(planck_law(wl, temp))*(wl[1]-wl[0])/Stefan_Boltzmann(temp)

figure()
wl = np.linspace(0.1e-6, 1.e-6, 1000)
plot(wl*1.e6,pi*cumsum(planck_law(wl, temp))*(wl[1]-wl[0])/Stefan_Boltzmann(temp))

# in-band black body radiation
def f_bb(wl_lower_meters, wl_upper_meters, temp):
    wl = np.linspace(wl_lower_meters, wl_upper_meters, 1000)
    # could use a numpy integrate to set precision, but this is good enough for now.
    return pi*sum(planck_law(wl, temp))*(wl[1]-wl[0])/Stefan_Boltzmann(temp)

print '%1.3f'%f_bb(0.6e-6, 0.8e-6, 5790.)
print '%1.3f'%f_bb(0.6e-6, 0.8e-6, 5260.)
print '%1.3f'%f_bb(0.6e-6, 0.8e-6, 6530.)
print '%1.3f'%f_bb(0.6e-6, 0.8e-6, 5084.)
print '%1.3f'%f_bb(0.6e-6, 0.8e-6, 6299.)
'''

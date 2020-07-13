from pylab import *
from Star import Star
from scipy import interpolate

#import matplotlib.colors as colors
#rcParams['font.size']=15

class Solar_Glint:
    def __init__(self, lower_wavelength = 620.e-9, upper_wavelength = 800.e-9, fnm1 = None, fnm2 = None, fnm3 = None, fnm4 = None):
        
        if fnm1 == None: fnm1 = "/Users/romerowo/Desktop/StarShade/ExoCat/glint_620_800_53deg.txt"
        if fnm2 == None: fnm2 = "/Users/romerowo/Desktop/StarShade/ExoCat/glint_620_800_63deg.txt"
        if fnm3 == None: fnm3 = "/Users/romerowo/Desktop/StarShade/ExoCat/glint_620_800_73deg.txt"        
        if fnm4 == None: fnm4 = "/Users/romerowo/Desktop/StarShade/ExoCat/glint_620_800_83deg.txt"

        self.ang_array = np.array([53., 63., 73., 83.])
        
        self.solar_mag = -26.74 # apparent magnitude of the Sun
        self.wl_1 = lower_wavelength
        self.wl_2 = upper_wavelength
        
        # use star model for the Sun
        Sun = Star(L_bol = 1., d_pc = 4.84813681e-6, Temperature=5778.)
        flux = Sun.get_Flux(self.wl_1, self.wl_2) # total flux in W/m^2
        #print 'flux, star1.Flux_bol', flux, star1.Flux_bol # in Watts / m^2

        self.data1 = np.loadtxt(fname=fnm1)
        self.data2 = np.loadtxt(fname=fnm2)
        self.data3 = np.loadtxt(fname=fnm3)
        self.data4 = np.loadtxt(fname=fnm4)
        
        self.flux1 = flux * 10**(-(self.data1-self.solar_mag)/2.5) # total flux in W/m^2
        self.flux2 = flux * 10**(-(self.data2-self.solar_mag)/2.5) # total flux in W/m^2
        self.flux3 = flux * 10**(-(self.data3-self.solar_mag)/2.5) # total flux in W/m^2
        self.flux4 = flux * 10**(-(self.data4-self.solar_mag)/2.5) # total flux in W/m^2
        
        self.pixel_mas = 3. 
        self.image_size = len(self.flux1)
        self.extent_mas = (self.image_size-1)/2 * self.pixel_mas
        

    def get_pixel(self, theta_mas, phi_rad):
        y_val = int(np.round(theta_mas*np.cos(phi_rad) / self.pixel_mas)) + (self.image_size+1)/2
        x_val = int(np.round(theta_mas*np.sin(phi_rad) / self.pixel_mas)) + (self.image_size+1)/2
        return int(x_val), int(y_val)
            #def plot_image(fnum=0):
        
    def get_flux(self, theta_mas, phi_rad, solar_angle):
        x, y = self.get_pixel(theta_mas, phi_rad)
        if x<0 or x>self.image_size-1 or y<0 or y>self.image_size-1: 
            return 0
        flux_array = np.array([self.flux1[x,y], self.flux2[x,y], self.flux3[x,y], self.flux4[x,y]])
        fun = interpolate.interp1d(self.ang_array, flux_array, fill_value = "extrapolate")
        return fun(solar_angle)
        
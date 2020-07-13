from pylab import *
import matplotlib.colors as colors
rcParams['font.size']=15

from scipy.interpolate import interp2d
from scipy.interpolate import interp1d

class Zodi:
    def __init__(self, Zodi_table_path='./'):
        
        self.lam = np.array([0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.0, 1.2, 2.2, 3.5, 4.8, 12, 25, 60, 100, 140]) # um
        self.Blam = np.array([2.5e-8, 5.3e-7, 2.2e-6, 2.6e-6, 2.0e-6, 1.3e-6, 1.2e-6, 8.1e-7, 1.7e-7, 5.2e-8, 1.2e-7, 7.5e-7, 3.2e-7, 1.8e-8, 3.2e-9, 6.9e-10]) # W/m2/sr/um

        # table 17 in Leinert et al. (1998)
        # Zodiacal Light brightness function of solar LON (rows) and LAT (columns)
        # values given in W m^-2 sr^-1 um^-1 for a wavelength of 500 nm
        #self.Izod = np.loadtxt('/Users/romerowo/EXOSIMS/EXOSIMS/ZodiacalLight/Leinert98_table17.txt')*1e-8 # W/m2/sr/um
        if Zodi_table_path[-1]=='/': Zodi_table_path = Zodi_table_path[:-1] 
        self.Izod = np.loadtxt('%s/Leinert98_table17.txt'%Zodi_table_path)*1e-8 # W/m2/sr/um
        # create data point coordinates
        self.lon_pts = np.array([0., 5, 10, 15, 20, 25, 30, 35, 40, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]) # deg
        self.lat_pts = np.array([0., 5, 10, 15, 20, 25, 30, 45, 60, 75, 90]) # deg
        self.y_pts, self.x_pts = np.meshgrid(self.lat_pts, self.lon_pts)
        self.points = np.array(zip(np.concatenate(self.x_pts), np.concatenate(self.y_pts)))
        
        #zodi_500nm = griddata(self.points, self.Izod, (self.x_pts, self.y_pts), method='linear')
        self.zodi_500nm = interp2d(self.lat_pts, self.lon_pts, np.log10(self.Izod), kind='linear')
        # function is (lat, lon) *NOT* (lon,lat)
        self.zodi_spec_shape = interp1d(self.lam, np.log10(self.Blam))
        
        
    def get_zodi_500nm(self, solar_lat = 0., solar_lon = 0.):
        return 10**self.zodi_500nm(np.abs(solar_lat),solar_lon)
    
    def get_zodi_spectrum(self, band_0, band_1, num_points=1000):
        lam_array = np.linspace(band_0*1.e-6, band_1*1.e6, num_points)
        d_lam = lam_array[1] - lam_array[0]
        return 10**self.zodi_spec_shape(lam_array)

    def get_zodi_integrated_spectrum(self, band_0, band_1, num_points=1000):
        lam_array = np.linspace(band_0*1.e6, band_1*1.e6, num_points) # band changed from meters to microns
        d_lam = lam_array[1] - lam_array[0]
        return np.sum(10**self.zodi_spec_shape(lam_array)) * d_lam 

    def get_zodi_flux(self, detector_model, solar_lat, solar_lon, num_points=1000):
        
        zodi_int = self.get_zodi_integrated_spectrum(detector_model.band[0], detector_model.band[1], num_points=num_points)
        zodi_norm = self.get_zodi_500nm(solar_lat, solar_lon) / self.Izod[12,0]
        #arcsecond_in_radians = 4.84814e-6
        #d_Omega_sr = 2.*pi*( 1. - np.cos(psf_diameter_arcsec/2. * arcsecond_in_radians) )

        val = zodi_int * zodi_norm * detector_model.d_Omega_PSF_sr # in W / m^-2
        return val[0]
    
    def plot_input_data(self):
        figure()
        pcolormesh(self.x_pts, self.y_pts, self.Izod,  norm=colors.LogNorm())
        xlabel('Lon, deg')
        ylabel('Lat, deg')
        colorbar(label='W m$^{-2}$ sr$^{-1}$ $\mu$m$^{-1}$')
        
        figure()
        loglog(self.lam, self.Blam)
        plot([np.min(self.lam), np.max(self.lam)], [self.Izod[12,0], self.Izod[12,0]])
        xlabel(r'$\lambda$, $\mu$m')
        ylabel(r'Zodi Brightness, W m$^{-2}$ sr$^{-1}$ $\mu$m$^{-1}$')
        grid(True)
        title('Brightness Spectrum at lon, lat = (%1.0f$^\circ$, %1.0f$^\circ$)'%(self.x_pts[12,0], self.y_pts[12,0]), fontsize=16)

        figure()
        loglog(self.lam, self.Blam*self.Izod[0,5]/self.Izod[12,0])
        plot([np.min(self.lam), np.max(self.lam)], [self.Izod[0,5], self.Izod[0,5]])
        xlabel(r'$\lambda$, $\mu$m')
        ylabel(r'Zodi Brightness, W m$^{-2}$ sr$^{-1}$ $\mu$m$^{-1}$')
        grid(True)
        title('Brightness Spectrum at lon, lat = (%1.0f$^\circ$, %1.0f$^\circ$)'%(self.x_pts[0,5], self.y_pts[0,5]), fontsize=16)
        print ('Izod[12,0]',self.Izod[12,0], self.x_pts[12,0], self.y_pts[12,0])
        print ('Izod[0,5]', self.Izod[0,5],  self.x_pts[0,5],  self.y_pts[0,5])
        
        figure()
        print ('self.get_zodi_500nm(solar_lat = 0., solar_lon =  0.)', self.get_zodi_500nm(solar_lat = 0.,   solar_lon =  0.))
        print ('self.get_zodi_500nm(solar_lat =  0., solar_lon = 90.)', self.get_zodi_500nm(solar_lat =  0., solar_lon = 90.))
        print ('self.get_zodi_500nm(solar_lat =  0., solar_lon = 90.)', self.get_zodi_500nm(solar_lat =  0., solar_lon = 180.))
        for k in range(0,10):
            lat = np.random.uniform(0.,90.)
            lon = np.random.uniform(0.,180.)
            print ('self.get_zodi_500nm(solar_lat =  %1.1f, solar_lon =  %1.1f)'%(lat,lon), self.get_zodi_500nm(solar_lat = lat, solar_lon =  lon))

        
    def exozodi_extensions(self):
        figure()
        pcolormesh(self.x_pts, self.y_pts, self.Izod,  norm=colors.LogNorm(), cmap='jet_r')
        xlabel('Solar Lon, deg')
        ylabel('Solar Lat, deg')
        colorbar(label='W m$^{-2}$ sr$^{-1}$ $\mu$m$^{-1}$')

        V_band_zero_point = 3.917e-8 # Watts m^-2 um^-1
        parsec_in_AU = 206265.
        sr_in_arcsec2 = 4.25e10
        abs_mag = -2.5*np.log10(self.Izod / sr_in_arcsec2 / V_band_zero_point * (550./500)**2)
        figure()
        pcolormesh(self.x_pts, self.y_pts, abs_mag, vmin = np.max(abs_mag), vmax=np.min(abs_mag))
        xlabel('Lon, deg')
        ylabel('Lat, deg')
        cb1 = colorbar(label='$M_V$, mag arcsec$^{-2}$')
        cb_ticks = np.arange( np.ceil(np.min(abs_mag)), np.floor(np.max(abs_mag))+0.1, 1.)
        cb1.set_ticks(cb_ticks)
        cb1.set_ticklabels(cb_ticks[::-1])
        figure()
        ax = subplot(111)
        plot(self.y_pts[0,:], abs_mag[0,:], 'o-')
        xlabel('Solar Latitude, deg')
        ylabel('$M_V$, mag arcsec$^{-2}$')
        yticks(arange(13,25,1.))
        y0, y1 = ax.get_ylim()
        ylim(y1,y0)
        grid(True)
        
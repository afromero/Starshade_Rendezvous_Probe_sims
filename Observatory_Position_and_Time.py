from pylab import *
import ephem

#from scipy.interpolate import interp1d

#stefan_boltzmann_constant = 5.670367e-8 # Watts m^-2 Kelvin^-4
#boltzmann_constant        = 1.38064852e-23  # Watts / Hz / K-1 
#planck_constant           = 6.62607004e-34 #m^2 kg / s
#speed_of_light            = 299792458. # m / s
#L_Sun_Watts               = 3.846e26 # Watts
#pc_in_meters              = 3.086e+16 # meters




'''
def planck_law(wavelength_meters, temperature_Kelvin):
    val  =  2.*planck_constant*speed_of_light**2 / wavelength_meters**5
    exp_arg = planck_constant*speed_of_light / (wavelength_meters*boltzmann_constant*temperature_Kelvin)
    val  /= np.exp( exp_arg ) - 1.
    return val # in Watts / sr / m^2 / m

def Stefan_Boltzmann(temperature_kelvin):
    return stefan_boltzmann_constant * temperature_kelvin**4 # Watts / meter^2
'''

class Observatory_Position_and_Time:
    '''
    Class for stellar modeling
    '''
    
    def __init__(self, RA=None, Dec=None):
        '''
        INITIALIZE VALUES
        '''
        target_coord = ephem.Equatorial(RA*pi/180., Dec*pi/180., epoch=ephem.J2000) # input is in radians
        #print 'target_coord.ra, target_coord.dec', target_coord.ra, target_coord.dec
        self.star = ephem.FixedBody()
        self.star._ra    = target_coord.ra
        self.star._dec   = target_coord.dec
        self.star._epoch = ephem.J2000
        self.star.compute()
        self.ecliptic = ephem.Ecliptic(target_coord)
        self.ecliptic_lat_deg, self.ecliptic_lon_deg  = float(self.ecliptic.lat)*180./pi, float(self.ecliptic.lon)*180./pi
        #print 'star.ra, star.dec', self.star.ra, self.star.dec, self.star._epoch

    def get_star_sun_lat_lon(self, initial_date = '2025/1/1', num_days=1100):
        d_lat_array = []
        d_lon_array = []
        for k in range(0,num_days,1):
            #print target_coord.ra, target_coord.dec
            #star = ephem.FixedBody(target_coord)
            #print 'xx', star.ra, star.dec
            day_0 = ephem.Date(initial_date) 
            Julian_day = float(day_0)+float(k)
            #print '*', ephem.Date(Julian_day)
            sun = ephem.Sun(Julian_day)
            self.star.compute(Julian_day)
            #print '\t',ephem.degrees(ephem.separation(sun, star)), 180./pi*float(ephem.degrees(ephem.separation(sun, star)))
            # IN DEGREES
            # get ecliptic for star and sun in degrees
            sun_ecliptic  =  ephem.Ecliptic(sun)
            star_ecliptic =  ephem.Ecliptic(self.star)
            sun_elat, sun_elon   = float(sun_ecliptic.lat)*180./pi,  float(sun_ecliptic.lon)*180./pi
            star_elat, star_elon = float(star_ecliptic.lat)*180./pi, float(star_ecliptic.lon)*180./pi

            d_lat = np.abs(star_elat - sun_elat)
            d_lon = np.abs(star_elon - sun_elon)
            while(d_lon>360.): d_lon-=360. 
            while(d_lon<0.): d_lon+=360. 
            if(d_lon>180.): d_lon = 360. - d_lon
                
            d_lat_array.append(d_lat)
            d_lon_array.append(d_lon)
            
            #print '%1.2f %1.2f %1.2f %1.2f :: %1.2f %1.2f'%(sun_elat, sun_elon, star_elat, star_elon, d_lat, d_lon)
            #lat.append(180./pi*float(ephem.degrees(ephem.separation(sun, self.star))))
        return np.array(d_lat_array), np.array(d_lon_array)

    def get_sun_angles(self, initial_date = '2025/1/1', num_days=1100):
        ang = []
        for k in range(0,num_days,1):
            #print target_coord.ra, target_coord.dec
            #star = ephem.FixedBody(target_coord)
            #print 'xx', star.ra, star.dec
            day_0 = ephem.Date(initial_date) 
            Julian_day = float(day_0)+float(k)
            #print '*', ephem.Date(Julian_day)
            sun = ephem.Sun(Julian_day)
            self.star.compute(Julian_day)
            #print '\t',ephem.degrees(ephem.separation(sun, star)), 180./pi*float(ephem.degrees(ephem.separation(sun, star)))
            # IN DEGREES
            ang.append(180./pi*float(ephem.degrees(ephem.separation(sun, self.star))))

        return np.array(ang)
    
    def get_observing_windows(self, init_date, duration_days, sun_ang_min=54., sun_ang_max=83.):
        # returns the indices of where the observation windows start and end.
        ang = self.get_sun_angles(initial_date = init_date, num_days=duration_days) # pad to get the wrap-around
        #print 'ang', ang
        self.days = np.arange(0,len(ang))
        self.obs_cut = np.logical_and(ang>sun_ang_min, ang<sun_ang_max) # these are acceptable

        windows = np.array(self.obs_cut).astype('float')
        diff_win = np.diff(windows)
        self.starts = np.where(diff_win>1.-1.e-10)[0] + 1 # this is because diff(a)[n] = a[n+1] - a[n]
        self.ends   = np.where(diff_win<-1.+1.e-10)[0] + 1 # this is to make the point inclusive a[n:m] ranges from a[m] to a[n-1]
        # check if day 0 is an observation window
        #print 'starts', self.starts, self.obs_cut[0]
        #print 'ends', self.ends
        if 0 not in self.starts and self.obs_cut[0]==True:
            self.starts = np.concatenate([[0],self.starts])
        if self.days[-1] not in self.ends and self.obs_cut[-1]==True:
            self.ends = np.concatenate([self.ends, [self.days[-1]]])
        self.durations = self.ends - self.starts
        #self.durations = 
        #print 'starts', self.starts, self.obs_cut[0]
        #print 'ends', self.ends
        start_end_array=[]
        for k in range(0,len(self.starts)):
            start_end_array.append([self.starts[k], self.ends[k]])
        #print start_end_array
        return start_end_array
    
    
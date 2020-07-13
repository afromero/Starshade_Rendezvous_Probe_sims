from pylab import *
import os.path
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline


class Kepler:
    '''
    A class to forward model Keplerian orbits
    '''
    def __init__(self, force_new_LUT = False, LUT_path='./'):
        '''
        Check if Eccentric Anomaly lookup table exists
        '''
        if LUT_path[-1]=='/': LUT_path = LUT_path[:-1] 

        if (not os.path.isfile('%s/Ecc_Anom_LUT.npz'%LUT_path)) or force_new_LUT==True:
            print ('In Kepler.py:\n\t-Look-up table not found.' )
            print ('\t-Generating look-up table.')
            print ('\t-This may take a few minutes.')
            print ('\t-It only needs to be done once.')
            self.Generate_LUT()

        # Load LUT
        fl = np.load('%s/Ecc_Anom_LUT.npz'%LUT_path)
        #print fl.keys()
        self.eccentricity_array = fl['eccentricity_array']
        self.mean_anomaly_array = fl['mean_anomaly_array']
        self.nE_grid            = fl['nE_grid']
        self.get_eccentric_anomaly = RectBivariateSpline(self.eccentricity_array, 
                                                         self.mean_anomaly_array, 
                                                         self.nE_grid, 
                                                         kx=1, ky=1) # linear grid interpolation

            
        self.AU_in_km = 1.496e+8
        self.mas_in_radians = 4.84814e-9
        self.Grav_constant = 6.67408e-11 #m^3 kg^-1 s^-2
        self.mass_of_Sun_in_kg = 1.989e30 #kg
        self.AU_in_meters     = 1.496e+11
        self.parsec_in_meters = 3.086e+16

    def Mean_anomaly(self, Eccentric_anomaly, eccentricity):
        return Eccentric_anomaly - eccentricity*np.sin(Eccentric_anomaly)
    
    def Generate_LUT(self,delta_e=1.e-3, delta_M=1.e-3):
        eccentric_anomaly_array = np.arange(0., 2.*pi + delta_M/2., delta_M) 
        eccentricity_array = np.arange(0., 1. + delta_e/2.,    delta_e)
        # make a grid of M(E,e)
        print ('Making grid of points')
        EA_grid, ec_grid = np.meshgrid(eccentric_anomaly_array, eccentricity_array)
        M_grid = self.Mean_anomaly(EA_grid, ec_grid)
        
        print ('M_grid.shape', M_grid.shape )
        print ('Interpolating to a new grid')
        # Now we want to make an evenly spaced E(M,e) grid
        points = np.array( (M_grid.flatten(), ec_grid.flatten()) ).T
        values = EA_grid.flatten()
        mean_anomaly_array = np.arange(0., 2.*pi + delta_M/2., delta_M) 
        nM_grid, ne_grid = np.meshgrid(mean_anomaly_array, eccentricity_array)
        nE_grid = griddata( points, values, (nM_grid, ne_grid), method='linear' )

        print ('nE_grid.shape', nE_grid.shape)
        
        np.savez('Ecc_Anom_LUT.npz', 
                 eccentricity_array=eccentricity_array, 
                 mean_anomaly_array=mean_anomaly_array, 
                 nE_grid = nE_grid)

        print ('Plotting figures')
        
        figure()
        contourf(EA_grid, ec_grid, M_grid, 20)
        xlabel('Eccentric Anomaly, rad')
        ylabel('eccentricity')
        colorbar()
        
        figure()
        contourf(M_grid, ec_grid, EA_grid, 20)
        xlabel('Mean Anomaly, rad')
        ylabel('eccentricity')
        colorbar()

        figure()
        contourf(nM_grid, ne_grid, nE_grid, 20)
        xlabel('Mean Anomaly, rad')
        ylabel('eccentricity')
        colorbar()
        
    def get_Orbit_plane_points(self, time_array, a_SM, t_0, Period, eccentricity):
        # Calculate the Mean Anomaly
        M_array = 2.*pi*(time_array-t_0)/Period

        # The LUT needs the mean anomaly to be between 0 and 2*pi
        M_array = np.mod(M_array, 2.*pi)

        # Calculate the Eccentric Anomaly
        E_array=[]
        # this silly loop is because the RectBivariateSpline function expects values to be strictly increasing
        # We expect few time points for starshade so this is not expected to be a bottleneck.
        for k in range(0,len(M_array)):
            E_val = self.get_eccentric_anomaly(eccentricity, M_array[k])
            E_array.append(E_val[0])
        E_array = np.array(E_array)
        #print E_array
        
        # Calculate the positions in the orbit plane
        x_p = a_SM * (np.cos(E_array)-eccentricity)
        y_p = a_SM * np.sqrt(1.-eccentricity**2) * np.sin(E_array)
        return np.array([x_p, y_p])
    
    def get_Orbit_vectors(self, time_array, a_SM, t_0, Period, eccentricity, LAN, inclination, arg_peri):
        x_p, y_p = self.get_Orbit_plane_points(time_array, a_SM, t_0, Period, eccentricity)
        # Compute rotation factor
        cos_LAN = np.cos(LAN) # Omega
        sin_LAN = np.sin(LAN)
        cos_inc = np.cos(inclination) # i
        sin_inc = np.sin(inclination)
        cos_per = np.cos(arg_peri) # omega
        sin_per = np.sin(arg_peri)
        
        A = +cos_LAN*cos_per - sin_LAN*cos_inc*sin_per
        B = +sin_LAN*cos_per + cos_LAN*cos_inc*sin_per
        C = -cos_LAN*sin_per - sin_LAN*cos_inc*cos_per
        D = -sin_LAN*sin_per + cos_LAN*cos_inc*cos_per
        U = +sin_LAN*sin_inc
        V = -cos_LAN*sin_inc
        return np.array( [ (A * x_p) + (B * y_p), (C * x_p) + (D * y_p), (U * x_p) + (V * y_p) ] )

    def get_Orbit_points(self, time_array, a_SM, t_0, Period, eccentricity, LAN, inclination, arg_peri):
        x,y,z = self.get_Orbit_vectors(time_array, a_SM, t_0, Period, eccentricity, LAN, inclination, arg_peri)
        '''
        x_p, y_p = self.get_Orbit_plane_points(time_array, a_SM, t_0, Period, eccentricity)
        # Compute rotation factor
        cos_LAN = np.cos(LAN) # Omega
        sin_LAN = np.sin(LAN)
        cos_inc = np.cos(inclination) # i
        sin_inc = np.sin(inclination)
        cos_per = np.cos(arg_peri) # omega
        sin_per = np.sin(arg_peri)
        
        A = +cos_LAN*cos_per - sin_LAN*cos_inc*sin_per
        B = +sin_LAN*cos_per + cos_LAN*cos_inc*sin_per
        C = -cos_LAN*sin_per - sin_LAN*cos_inc*cos_per
        D = -sin_LAN*sin_per + cos_LAN*cos_inc*cos_per
        '''
        return np.array( [ x, y ] )

    def ravel_Orbit_points(self, time_array, a_SM, t_0, Period, eccentricity, LAN, inclination, arg_peri):
        return np.ravel(self.get_Orbit_points(time_array, a_SM, t_0, Period, eccentricity, LAN, inclination, arg_peri))
    

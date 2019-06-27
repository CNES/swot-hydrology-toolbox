import mathematical_function as math_fct
import numpy as np
import lib.height_model as height_model
import lib.true_height_model as true_height_model
import scipy
import time
import pyproj
import utm

from lib.my_variables import COEFF_X2, COEFF_Y2, COEFF_X, COEFF_Y, COEFF_XY, COEFF_CST, RAD2DEG, DEG2RAD


class Lac:

    def __init__(self, num):
        
        self.num = num
        self.seed = int(str(time.time()).split('.')[1])
        self.hmean = None
        
    def compute_pixels_in_given_lac(self, OUT_ind_lac_data):
        self.pixels = np.where(OUT_ind_lac_data == self.num)        
    def set_hmean(self, hmean):
        self.hmean = hmean
        
class Constant_Lac(Lac):
    
    def __init__(self, num, IN_attributes, lat, IN_cycle_number):
        Lac.__init__(self, num)
        self.height_model_a = IN_attributes.height_model_a
        self.lat_init = IN_attributes.lat_init[1:-1]
        self.cycle_number = IN_cycle_number
        self.cycle_duration = IN_attributes.cycle_duration
        self.height_model_t0 = IN_attributes.height_model_t0
        self.height_model_period = IN_attributes.height_model_period
        self.orbit_time = IN_attributes.orbit_time
        
        
        self.time = math_fct.linear_extrap(lat, self.lat_init, IN_attributes.orbit_time)
        self.mode = 'orbit_time'
        
    def compute_h(self, lat, lon):
    
        if self.mode == 'az':
            self.time = math_fct.linear_extrap(lat, self.lat_init, self.orbit_time)
    
        if self.mode == 'orbit_time':
            self.mode = 'az' 
        return self.height_model_a + self.height_model_a * np.sin(2*np.pi * (self.time + self.cycle_number * self.cycle_duration) - self.height_model_t0) / self.height_model_period


class Reference_height_Lac(Lac):
    

    def __init__(self, num, polygon_index, layer, IN_attributes):
        Lac.__init__(self, num)
        self.height_name = IN_attributes.height_name
        self.height = 0.
        for i in range(layer.GetFieldCount()):
        # Test 'HEIGHT' parameter in input shapefile fields
            if layer.GetFieldDefn(i).GetName() == self.height_name:
                if polygon_index.GetField(str(self.height_name)) is not None :
                    self.height = np.float(polygon_index.GetField(str(self.height_name)))
            
    def compute_h(self, lat, lon):
        return self.height
        
class Gaussian_Lac(Lac):
    
    def __init__(self, num, IN_attributes, lat, lon, IN_cycle_number):
        Lac.__init__(self, num)
        self.height_model_a = IN_attributes.height_model_a
        self.lat_init = IN_attributes.lat_init[1:-1]
        self.cycle_number = IN_cycle_number
        self.cycle_duration = IN_attributes.cycle_duration
        self.height_model_t0 = IN_attributes.height_model_t0
        self.height_model_period = IN_attributes.height_model_period
        self.orbit_time = IN_attributes.orbit_time
        
        self.time = math_fct.linear_extrap(lat, self.lat_init, self.orbit_time)
        self.mode = 'orbit_time'
        
        self.height_model_stdv = IN_attributes.height_model_stdv

        self.dlon =  10e-7
        self.dlat =  10e-7
        
                
        lonmin, lonmax, latmin, latmax = lon.min(), lon.max(), lat.min(), lat.max()

        taille_lon, taille_lat = np.int((lonmax-lonmin)/self.dlon), np.int((latmax-latmin)/self.dlat)
        
        self.height = height_model.generate_2d_profile_gaussian(self.dlat, latmin, latmax, self.dlon, lonmin, lonmax, self.height_model_stdv, seed = self.seed)
        print("gaussian min height",np.min(self.height))
        print("gaussian max height",np.max(self.height))
        
        self.h_interp = scipy.interpolate.RectBivariateSpline(latmin + self.dlat*np.arange(taille_lat),lonmin + self.dlon*np.arange(taille_lon),  self.height)
        
    def compute_h(self, lat, lon):
    
        if self.mode == 'az':
            self.time = math_fct.linear_extrap(lat, self.lat_init, self.orbit_time)
    
        if self.mode == 'orbit_time':
            self.mode = 'az' 
            
        h0 =  np.mean(self.height_model_a * np.sin(2*np.pi * (self.time + self.cycle_number * self.cycle_duration) - self.height_model_t0) / self.height_model_period)

        return h0 + self.h_interp.ev(lat,lon)
                
class Polynomial_Lac(Lac):
    
    def __init__(self, num, IN_attributes, lat , lon, IN_cycle_number):
        
        Lac.__init__(self, num)

        self.height_model_a = IN_attributes.height_model_a
        self.lat_init = IN_attributes.lat_init[1:-1]
        self.cycle_number = IN_cycle_number
        self.cycle_duration = IN_attributes.cycle_duration
        self.height_model_t0 = IN_attributes.height_model_t0
        self.height_model_period = IN_attributes.height_model_period
        self.orbit_time = IN_attributes.orbit_time
        
        self.time = math_fct.linear_extrap(lat, self.lat_init, IN_attributes.orbit_time)
        self.mode = 'orbit_time'
        
        x_c, y_c, zone_number, zone_letter = utm.from_latlon(lat[0], lon[0])
        # Convert pixel cloud to UTM (zone of the centroid)
        self.latlon = pyproj.Proj(init="epsg:4326")
        self.utm_proj = pyproj.Proj("+proj=utm +zone={}{} +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(zone_number, zone_letter))
        X, Y = pyproj.transform(self.latlon, self.utm_proj, lon, lat)                    
        k0 = np.random.randint(len(lat))
        self.X0, self.Y0 = X[k0], Y[k0]
        self.COEFF_X2 = np.random.random(1)[0] * COEFF_X2
        self.COEFF_Y2 = np.random.random(1)[0] * COEFF_Y2
        self.COEFF_X = np.random.random(1)[0] * COEFF_X
        self.COEFF_Y = np.random.random(1)[0] * COEFF_Y 
        self.COEFF_XY = np.random.random(1)[0] * COEFF_XY
        self.COEFF_CST = np.random.random(1)[0] * COEFF_CST      
        
    def compute_h(self, lat, lon):
    
        if self.mode == 'az':
            self.time = math_fct.linear_extrap(lat, self.lat_init, self.orbit_time)
    
        if self.mode == 'orbit_time':
            self.mode = 'az' 
                        
        h0 = np.mean(self.height_model_a * np.sin(2*np.pi * (self.time + self.cycle_number * self.cycle_duration) - self.height_model_t0) / self.height_model_period)
        X, Y = pyproj.transform(self.latlon, self.utm_proj, lon, lat)
        
        
        height_water = height_model.generate_2d_profile_2nd_order_list(self.X0, self.Y0, X, Y, self.COEFF_X2, self.COEFF_Y2, self.COEFF_X, self.COEFF_Y, self.COEFF_XY, self.COEFF_CST)
        
        return h0 + height_water

class Height_in_file_Lac(Lac):
    # Process true height model from Kevin Larnier
    # TDB : Add security for lat lon boundaries
    # TBD : Add specific model for 1D model (river)
    def __init__(self, num, IN_attributes):
        Lac.__init__(self, num)
        self.height_name = IN_attributes.height_name
        self.trueheight_file = IN_attributes.trueheight_file
            
    def compute_h(self, lat, lon):
        true_height_model_inst = true_height_model.TrueHeightModel(self.trueheight_file, lat, lon, verbose=True)
        true_height_model_inst.apply_model()
        return true_height_model_inst.final_height
                   
                    

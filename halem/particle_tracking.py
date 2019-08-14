import numpy as np

class particle_tracking():
    """
    Particle tracking function
    EF is euler foreward scheme
    RK4 is  the Runge-Kutta scheme after: 
    https://nl.wikipedia.org/wiki/Runge-Kuttamethode
    for the interpolation between gridpints the nearest-point method is used.
    
    This method assumes:
    - rectangular or curvilinear grid
    - WGS84 coordinates
    
    the WGS84 transformation is done with the method presented in:
    https://en.wikipedia.org/wiki/Geographic_coordinate_system
    
    the class flow should contain:
    - x,y,t as a numpy meshgrid (y, x, t = np.meshgrid(y, x, t))
    - u, v
    - u.shape = v.shape = x.shape = y.shape = t.shape
    """
    
    def __init__(self, flow):
        """Initialis the class with the input data"""
        self.flow = flow


        
    def EF(self, x0, y0, T):
        """This method finds the path with the Euler foreward scheme"""
        self.x0 = x0
        self.y0 = y0
        self.T = T
        self.dt = self.T[1] - self.T[0]
        
        lat = y0 * 2* np.pi / 360
        self.mdeg_lat = 111132.92\
                        - 559.82 * np.cos(2*lat)\
                        + 1.175 * np.cos(4*lat)\
                        - 0.0023 * np.cos( 6 * lat)
        self.mdeg_lon = 111412.84 * np.cos(lat)\
                        - 93.5 * np.cos(3*lat)\
                        + 0.118 * np.cos(5*lat)
        x = [self.x0]
        y = [self.y0]
        for t in self.T:
            if t > self.flow.t[-1] or t < self.flow.t[0]:
                print('time out of range of the hydrodynamic model')
                break
            else:
                x0, y0 = self.numeric_scheme_EF(x0,y0,t)
                x.append(x0)
                y.append(y0)
        return x, y
    
    def RK4(self, x0, y0, T):
        """ This method finds the path with the RK4 scheme"""
        self.x0 = x0
        self.y0 = y0
        self.T = T
        self.dt = self.T[1] - self.T[0]
        
        lat = y0 * 2* np.pi / 360
        self.mdeg_lat = 111132.92\
                        - 559.82 * np.cos(2*lat)\
                        + 1.175 * np.cos(4*lat)\
                        - 0.0023 * np.cos( 6 * lat)
        self.mdeg_lon = 111412.84 * np.cos(lat)\
                        - 93.5 * np.cos(3*lat)\
                        + 0.118 * np.cos(5*lat)
        
        x = [self.x0]
        y = [self.y0]
        for t in self.T:
            if t > self.flow.t[-1] or t < self.flow.t[0]:
                print('time out of range of the hydrodynamic model')
                break
            else:
                x0, y0 = self.numeric_scheme_RK4(x0,y0,t)
                x.append(x0)
                y.append(y0)
        return x, y      
        
    def numeric_scheme_EF(self, x, y, t):
        """This method contains the euler foreward sheme"""
        u, v = self.interpol(x, y, t)
        
        x = x + self.dt * u  / self.mdeg_lon
        y = y + self.dt * v  / self.mdeg_lat
        
        return x,y
    
    def numeric_scheme_RK4(self, x, y, t):
        """This method contains the RK4 sheme"""

        u0, v0 = self.interpol(x,
                               y,
                               t,
                              )
        K1_x =  self.dt * u0 / self.mdeg_lon
        K1_y =  self.dt * v0 / self.mdeg_lat
        
        u1, v1 = self.interpol(x+K1_x/2,
                               y+K1_y/2,
                               t + self.dt/2, 
                              )
        K2_x = self.dt * u1 / self.mdeg_lon
        K2_y = self.dt * v1 / self.mdeg_lat
        
        u2, v2 = self.interpol(x+K2_x/2,
                               y+K2_y/2,
                               t + self.dt/2, 
                              )
        K3_x = self.dt * u2 / self.mdeg_lon
        K3_y = self.dt * v2 / self.mdeg_lat
        
        u3, v3 = self.interpol(x+K3_x,
                               y+K3_y,
                               t + self.dt, 
                              )
        K4_x = self.dt * u3 / self.mdeg_lon
        K4_y = self.dt * v3 / self.mdeg_lat
 
        x = x + (1/6)*(K1_x + K2_x + K3_x + K4_x)
        y = y + (1/6)*(K1_y + K2_y + K3_y + K4_y)
        
        return x,y
    
    def interpol(self, x, y, t):
        """ this function returns the specific flow condtions 
        for a given point in space and time"""
        DD = (self.flow.x[:,:] - x) ** 2 + (self.flow.y[:,:] - y)** 2 
        i = np.argwhere(DD == DD.min())[0]
        TT = (self.flow.t - t)**2
        ii = np.argwhere(TT == TT.min())[0,0]
       
        return self.flow.u[ii, i[0], i[1]], self.flow.v[ii, i[0], i[1]]
from collections import defaultdict
import math
import numpy as np
from numpy import ma

def Squat_in_sea(h,T,V_max, WWL,ukc):
    ghv2 = 9.81*h/(V_max**2)
    VhV1 = 1.002 + 0.005*(np.sqrt(WWL*T)/h)-0.1159*(np.sqrt(WWL*T)/h)**2+0.0191*(np.sqrt(WWL*T)/h)**3
    V1Vinf = ((np.exp(ghv2) - np.exp(-ghv2))/(np.exp(ghv2)+np.exp(-ghv2)))**0.5
    V = V_max * V1Vinf * VhV1
    
    return V

def haversine(coord1, coord2):
    R = 6372800                                      # https://janakiev.com/blog/gps-points-distance-python/
    lat1, lon1 = coord1
    lat2, lon2 = coord2                              # use the Haversine function to determine the distance between two points in the WGS84 coordinate system
    
    phi1, phi2 = math.radians(lat1), math.radians(lat2) 
    dphi       = math.radians(lat2 - lat1)
    dlambda    = math.radians(lon2 - lon1)
    
    a = math.sin(dphi/2)**2 + \
        math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2
    
    return 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))

def costfunction_timeseries(edge,V_max, WD_min ,flow):
    xfrom = flow.nodes[edge[0]][1]
    yfrom = flow.nodes[edge[0]][0]
    xto = flow.nodes[edge[1]][1]
    yto = flow.nodes[edge[1]][0]

    vship = Squat_in_sea(flow.WD[edge[0]], WD_min, V_max, flow.WWL, flow.ukc)

    v_w = (flow.v[edge[0]] + flow.v[edge[1]])/2
    u_w = (flow.u[edge[0]] + flow.u[edge[1]])/2
    U_w = (u_w**2 + v_w**2)**0.5
    
    alpha1 = np.arctan2((yto - yfrom),(xto - xfrom))
    alpha2 = np.arctan2(v_w , u_w) - alpha1
    
    s_t = (U_w * np.cos(alpha2)) + (vship ** 2 -  (U_w * np.sin(alpha2))**2) ** 0.5     # velocity in the direction on the course

    u_t = np.cos(alpha1)*( s_t)
    v_t = np.sin(alpha1)*( s_t)
    
    L = haversine((yfrom,xfrom), (yto,xto))

    U_t = (u_t**2 + v_t**2)**0.5
    t = L/U_t
    
    t[flow.WD[edge[0]] + flow.ukc < WD_min] = np.inf
    t[flow.WD[edge[1]] + flow.ukc < WD_min] = np.inf
    t[(U_w * np.sin(alpha2))**2 > vship ** 2] = np.inf
    t[s_t < 0] = np.inf  
    return np.array(t)

def costfunction_spaceseries(edge,V_max, WD_min ,flow):
    xfrom = flow.nodes[edge[0]][1]
    yfrom = flow.nodes[edge[0]][0]
    xto = flow.nodes[edge[1]][1]
    yto = flow.nodes[edge[1]][0]

    vship = Squat_in_sea(flow.WD[edge[0]], WD_min, V_max, flow.WWL, flow.ukc)

    v_w = (flow.v[edge[0]] + flow.v[edge[1]])/2
    u_w = (flow.u[edge[0]] + flow.u[edge[1]])/2
    U_w = (u_w**2 + v_w**2)**0.5
    
    alpha1 = np.arctan2((yto - yfrom),(xto - xfrom))
    alpha2 = np.arctan2(v_w , u_w) - alpha1
    
    s_t = (U_w * np.cos(alpha2)) + (vship ** 2 -  (U_w * np.sin(alpha2))**2) ** 0.5     # velocity in the direction on the course

    u_t = np.cos(alpha1)*( s_t)
    v_t = np.sin(alpha1)*( s_t)
    
    L = haversine((yfrom,xfrom), (yto,xto))

    U_t = (u_t**2 + v_t**2)**0.5
    t = L/U_t
    
    t[flow.WD[edge[0]] + flow.ukc < WD_min] = np.inf
    t[flow.WD[edge[1]] + flow.ukc < WD_min] = np.inf
    t[(U_w * np.sin(alpha2))**2 > vship ** 2] = np.inf
    t[s_t < 0] = np.inf
    t[t!= np.inf] = L 

    return np.array(t)
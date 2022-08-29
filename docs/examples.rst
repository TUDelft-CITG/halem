========
Examples
========

This page contains a few examples to help understand the software.

Explanation example
-------------------

The First example is a fictive Explanation example.In this example 
the shortest, fastest, cheapest, and least pollutant route is 
determind for a simple rotating flow with a random bathymetry. 


The First step is to load the dependencies.

.. code:: python3

    from collections import defaultdict
    import math
    import numpy as np
    from numpy import ma
    import datetime, time
    from datetime import datetime
    import netCDF4
    from netCDF4 import Dataset, num2date
    from scipy.interpolate import griddata
    from scipy.spatial import Delaunay


    import halem
    import halem.Mesh_maker as Mesh_maker
    import halem.Functions as Functions
    import halem.path_finder as path_finder
    import halem.Flow_class as Flow_class


    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    % matplotlib inline

    import os
    from cartopy import config
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

The second step is to determine the hydrodynamic conditions. Since
this is a fictive case the hydrodynamic conditions are made up. 
For real use cases the output of a CFD model can be used. 

.. code:: python3

    class flow_potentiaalveld():
        def __init__(self, name):
            d = datetime.strptime("23/03/2019 00:00:00", "%d/%m/%Y %H:%M:%S")
            t0 = d.timestamp() 
            x = np.arange(0,1,0.04)
            y = np.arange(0,1,0.04)
            t = np.arange(t0, (t0+2*30*60*30),30*60)

            y, x = np.meshgrid(y, x)
            y = y.reshape(y.size)
            x = x.reshape(y.size)

            nodes = np.zeros((len(x),2))
            nodes[:,0] = y
            nodes[:,1] = x

            u = []
            v = []
            for node in nodes:
                ut = 0*t + 2 * np.cos(np.pi*(node[0]))
                vt = 0*t - 2 * np.cos(np.pi*(node[1]))
                u.append(ut)
                v.append(vt)

            self.v = np.transpose(np.array(v))
            self.u = np.transpose(np.array(u))
            self.WD = np.random.rand(self.u.shape[0],self.u.shape[1])*12 + 8
            #self.WD = self.u * 0 + 20
            self.t = t
            self.nodes = nodes
            self.tria = Delaunay(nodes)

    def compute_cost(week_rate, fuel_rate):
        second_rate = week_rate/7/24/60/60
        return lambda travel_time, speed: (travel_time*second_rate + fuel_rate*travel_time * speed**3)

    def compute_co2(fuel_rate):
        return lambda travel_time, speed: (fuel_rate*travel_time * speed**3)


The third step is to describe the project parameters, and the vessel parameters. 

.. code:: python3

    nl = (3,2.5)
    dx_min = 0.01
    blend = 1
    vship = np.array([[1,4,8],
                      [1,5,10]
                     ])
    WD_min = np.array([8,7])
    WPVI = np.array([4000, 10000])
    ukc = 1.5
    WWL = 40
    name_textfile_flow = 'maaktnietuit'
    Load_flow = flow_potentiaalveld
    number_of_neighbor_layers = 3

    Q_cost = compute_cost(700_000, 0.008)
    Q_co2 = compute_co2(1)


The Fourth step is to generate the Roadmap. This step is the pre-processing step.

.. code:: python3

    start_time = time.time()

    Roadmap_t = Mesh_maker.GraphFlowModel(name_textfile_flow, 
                                            dx_min, 
                                            blend, 
                                            nl, 
                                            number_of_neighbor_layers, 
                                            vship, 
                                            Load_flow, 
                                            WD_min,
                                            WPVI,
                                            WWL = WWL,
                                            ukc = ukc,
                                            compute_cost = Q_cost,
                                            compute_co2 = Q_co2,
                                           )

    stop_time = time.time()
    computation_time = stop_time - start_time
    print("the computational time is:", round(computation_time,2), "sec")


4/4
the computational time is: 448.21 sec



The Fift step is the calculate the optimal shipping route. In this step the four routes are calculated. 

.. code:: python3

    t0 = "23/03/2019 03:00:00"          #Departure time

    start = (0.2, 0.2)
    stop = (0.8, 0.8)
    v_max = 7

    start_time = time.time()

    path_t, time_t, dist_t = halem.HALEM_time(start, stop, t0, v_max, Roadmap_t)
    path_s, time_s, dist_s = halem.HALEM_space(start, stop, t0, v_max, Roadmap_t)
    path_c, time_c, dist_c = halem.HALEM_cost(start, stop, t0, v_max, Roadmap_t)
    path_p, time_p, dist_p = halem.HALEM_co2(start, stop, t0, v_max, Roadmap_t)

    stop_time = time.time()
    computation_time = stop_time - start_time
    print("the computational time is:", round(computation_time,2), "sec")

the computational time is: 14.31 sec


The last step is to plot the results. The results are according to the expectations. 

.. code:: python3

    fig = plt.figure(figsize=(23, 17))

    ax = plt.subplot(3,2,1)
    plt.axis('square')
    a = 1

    x_r = np.arange(0,1, 0.075)
    y_r = np.arange(0,1, 0.075)
    y_r, x_r = np.meshgrid(y_r,x_r)

    WD_r = griddata((Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0]), Roadmap_t.WD[:,0], (x_r, y_r), method= 'linear')
    u_r = griddata((Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0]), Roadmap_t.u[:,0], (x_r, y_r), method= 'linear')
    v_r = griddata((Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0]), Roadmap_t.v[:,0], (x_r, y_r), method= 'linear')

    cval = np.arange(0,21,0.5)
    im = plt.contourf(x_r,y_r,WD_r, cval)
    fig.colorbar(im, ax=ax, label = 'Waterdepth in meters')

    plt.quiver(x_r[::a,::a],
               y_r[::a,::a],
               u_r[::a,::a],
               v_r[::a,::a],
               label = 'flow directions'
              )
    plt.plot(path_t[:,0], path_t[:,1],'r', label = 'fastest Route', linewidth = 3)
    plt.plot(path_s[:,0], path_s[:,1],'b', label = 'shortest Route', linewidth = 3)
    plt.plot(path_c[:,0], path_c[:,1],'g', label = 'cheapest Route', linewidth = 3)
    plt.plot(path_p[:,0], path_p[:,1],'m', label = 'cleanest Route', linewidth = 3)

    plt.plot(Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0], 'k.', label = 'Nodes of the graph', markersize = 1)
    plt.plot(path_t[0, 0], path_t[0, 1], 'go', label = 'start')
    plt.plot(path_t[-1, 0 ], path_t[-1, 1], 'ro', label = 'target')


    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title(r'x/y diagram of the simple rotating flow')
    ax.legend(loc='upper center', bbox_to_anchor=(-0.5, 0.5),
              ncol=1, fancybox=True, shadow=True)
    plt.xlabel('lat')
    plt.ylabel('lon')

    plt.subplot(3,2,2)
    halem.plot_timeseries2(path_t, time_t, Roadmap_t, 'r')
    plt.title('s/t diagram of the fastest route')


    plt.subplot(3,2,3)
    halem.plot_timeseries2(path_s, time_s, Roadmap_t, 'b')
    plt.title('s/t diagram of the shortest route')


    plt.subplot(3,2,4)
    halem.plot_timeseries2(path_c, time_c, Roadmap_t, 'g')
    plt.title('s/t diagram of the cheapest route')

    plt.subplot(3,2,5)
    halem.plot_timeseries2(path_p, time_p, Roadmap_t, 'm')
    plt.title('s/t diagram of the cleanest route')

    plt.show()

.. image:: Results1.PNG


Validation based on trival test cases
-------------------------------------

In this example the python code for the validation cases is displayed. 

Import the required libaries

.. code:: python3

    from collections import defaultdict
    import math
    import numpy as np
    from numpy import ma
    import datetime, time
    from datetime import datetime
    import netCDF4
    from netCDF4 import Dataset, num2date
    from scipy.interpolate import griddata
    from scipy.spatial import Delaunay


    import halem.Base_functions as halem
    import halem.Mesh_maker as Mesh_maker
    import halem.Functions as Functions
    import halem.path_finder as path_finder
    import halem.Flow_class as Flow_class


    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    % matplotlib inline

    import os
    from cartopy import config
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

Define the hydrodynamic conditions

.. code:: python3

    class flow_potentiaalveld():
        def __init__(self, name):
            d = datetime.strptime("23/03/2019 00:00:00", "%d/%m/%Y %H:%M:%S")
            t0 = d.timestamp() 
            x = np.arange(0,1,0.025)
            y = np.arange(0,1,0.025)
            t = np.arange(t0, (t0+30*60*30),30*60)

            y, x = np.meshgrid(y, x)
            y = y.reshape(y.size)
            x = x.reshape(y.size)

            nodes = np.zeros((len(x),2))
            nodes[:,0] = y
            nodes[:,1] = x
            tria = Delaunay(nodes)

            u = []
            v = []
            for node in nodes:
                ut = 0*t + 2 * np.cos(np.pi*(node[0]))
                vt = 0*t - 2 * np.cos(np.pi*(node[1]))
                u.append(ut)
                v.append(vt)

            self.v = np.transpose(np.array(v))
            self.u = np.transpose(np.array(u))
            self.WD = self.u*0 + 20
            self.t = t
            self.nodes = nodes
            self.tria = Delaunay(nodes)

    class flow_dyn_cur():
        def __init__(self, name):
            d = datetime.strptime("23/03/2019 00:00:00", "%d/%m/%Y %H:%M:%S")
            t0 = d.timestamp() 
            x = np.arange(0,1,0.01)
            y = np.arange(0,1,0.01)
            t = np.arange(t0, (t0+30*60*30),30*60)

            y, x = np.meshgrid(y, x)
            y = y.reshape(y.size)
            x = x.reshape(y.size)

            nodes = np.zeros((len(x),2))
            nodes[:,0] = y
            nodes[:,1] = x
            tria = Delaunay(nodes)

            u = []
            v = []
            for node in nodes:
                ut =   2 * np.cos(np.pi*(node[0])) * np.cos(2*np.pi * (t - t0)/24000)
                vt = - 2 * np.cos(np.pi*(node[1]))* np.cos(2*np.pi * (t - t0)/24000)
                u.append(ut)
                v.append(vt)

            self.v = np.transpose(np.array(v))
            self.u = np.transpose(np.array(u))
            self.WD = self.u*0 + 20
            self.t = t
            self.nodes = nodes
            self.tria = Delaunay(nodes)


    class flow_wait():
        def __init__(self, name):
            d = datetime.strptime("23/03/2019 00:00:00", "%d/%m/%Y %H:%M:%S") 
            t0 = d.timestamp() 
            x = np.arange(0,1,0.025)
            y = np.arange(0,1,0.025)
            t = np.arange(t0, (t0+30*60*30),10*60)

            y, x = np.meshgrid(y, x)
            y = y.reshape(y.size)
            x = x.reshape(y.size)

            nodes = np.zeros((len(x),2))
            nodes[:,0] = y
            nodes[:,1] = x
            tria = Delaunay(nodes)

            u = []
            v = []
            for node in nodes:
                ut = 0*t + 0
                vt = 0*t + 0
                u.append(ut)
                v.append(vt)

            self.v = np.transpose(np.array(v))
            self.u = np.transpose(np.array(u))
            self.t = t
            self.nodes = nodes
            self.tria = Delaunay(nodes)

            WD = self.u*0 + 20
            for i in range(len(nodes)):
                node = nodes[i]
                if node[0] > 0.4 and node[0] < 0.6:
                    if node[1] > 0.1 and node[1] < 0.9:
                        WD[5:50,i] = 0
                        WD[75:-1,i] = 0
            self.WD = WD

Define the project parameters.

.. code:: python3

    nl = (3,2.5)
    dx_min = 0.01
    blend = 1
    vship = np.array([[4]])
    WD_min = np.array([[5]])
    ukc = 1.5
    WWL = 40
    WVPI = [10000]
    name_textfile_flow = 'maaktnietuit'

    def compute_cost(week_rate, fuel_rate):
        second_rate = week_rate/7/24/60/60
        return lambda travel_time, speed: (travel_time*second_rate + fuel_rate*travel_time * speed**3)

    QQ = compute_cost(700_000, 0.0008)

Define the Roadmaps

.. code:: python3

    Load_flow = flow_potentiaalveld
    number_of_neighbor_layers = 3

    Roadmap_t = Mesh_maker.GraphFlowModel(name_textfile_flow, 
                                            dx_min, 
                                            blend, 
                                            nl, 
                                            number_of_neighbor_layers, 
                                            vship, 
                                            Load_flow, 
                                            WD_min,
                                            WVPI,
                                            WWL = WWL,
                                            ukc = ukc,
                                            optimization_type = ['time']
                                       )
    Load_flow = flow_dyn_cur
    number_of_neighbor_layers = 3
    Roadmap_d = Mesh_maker.GraphFlowModel(name_textfile_flow, 
                                            dx_min, 
                                            blend, 
                                            nl, 
                                            number_of_neighbor_layers, 
                                            vship, 
                                            Load_flow, 
                                            WD_min,
                                            WVPI,
                                            WWL = WWL,
                                            ukc = ukc,
                                            optimization_type = ['time']
                                       )

    Load_flow = flow_wait
    number_of_neighbor_layers = 1
    Roadmap = Mesh_maker.GraphFlowModel(name_textfile_flow, 
                                            dx_min, 
                                            blend, 
                                            nl, 
                                            number_of_neighbor_layers, 
                                            vship, 
                                            Load_flow, 
                                            WD_min,
                                            WVPI,
                                            WWL = WWL,
                                            ukc = ukc,
                                            optimization_type = ['time']
                                       )


Calulate the optimized routes

.. code:: python3

	path_t, time_t, dist_t = halem.HALEM_time(start, stop, t0, 4, Roadmap_t)
	path_d, time_d, dist_d = halem.HALEM_time(start, stop, t0, 4, Roadmap_d)
	path_w, time_w, dist_w = halem.HALEM_time(start, stop, t0, 4, Roadmap)

Plot the results:

.. code:: python3

    fig = plt.figure(figsize=(23, 17))

    ax = plt.subplot(2,3,1)
    plt.axis('square')
    a = 1

    x_r = np.arange(0,1, 0.075)
    y_r = np.arange(0,1, 0.075)
    y_r, x_r = np.meshgrid(y_r,x_r)

    WD_r = griddata((Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0]), Roadmap_t.WD[:,0], (x_r, y_r), method= 'linear')
    u_r = griddata((Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0]), Roadmap_t.u[:,0], (x_r, y_r), method= 'linear')
    v_r = griddata((Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0]), Roadmap_t.v[:,0], (x_r, y_r), method= 'linear')

    cval = np.arange(0,21,0.5)
    im = plt.contourf(x_r,y_r,WD_r, cval)
    fig.colorbar(im, ax=ax, label = 'Waterdepth in meters')

    plt.quiver(x_r[::a,::a],
               y_r[::a,::a],
               u_r[::a,::a],
               v_r[::a,::a],
               label = 'flow directions'
              )
    plt.plot(path_t[:,0], path_t[:,1], 'm', label = 'Route')
    plt.plot(Roadmap_t.nodes[:,1], Roadmap_t.nodes[:,0], 'k.', label = 'Nodes of the graph', markersize = 1)
    plt.plot(path_t[0, 0], path_t[0, 1], 'go', label = 'start')
    plt.plot(path_t[-1, 0 ], path_t[-1, 1], 'ro', label = 'target')


    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title(r'x/y diagram of the simple rotating flow')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              ncol=2, fancybox=True, shadow=True)
    plt.xlabel('lat')
    plt.ylabel('lon')

    plt.subplot(2,3,4)
    halem.plot_timeseries2(path_t, time_t, Roadmap_t)
    plt.title('s/t diagram of the simple rotating flow')

    ax = plt.subplot(2,3,2)
    plt.axis('square')
    a = 1

    x_r = np.arange(0,1, 0.075)
    y_r = np.arange(0,1, 0.075)
    y_r, x_r = np.meshgrid(y_r,x_r)

    WD_r = griddata((Roadmap_d.nodes[:,1], Roadmap_d.nodes[:,0]), Roadmap_d.WD[:,0], (x_r, y_r), method= 'linear')
    u_r = griddata((Roadmap_d.nodes[:,1], Roadmap_d.nodes[:,0]), Roadmap_d.u[:,0], (x_r, y_r), method= 'linear')
    v_r = griddata((Roadmap_d.nodes[:,1], Roadmap_d.nodes[:,0]), Roadmap_d.v[:,0], (x_r, y_r), method= 'linear')


    im = plt.contourf(x_r,y_r,WD_r, cval)
    fig.colorbar(im, ax=ax, label = 'Waterdepth in meters')

    plt.quiver(x_r[::a,::a],
               y_r[::a,::a],
               u_r[::a,::a],
               v_r[::a,::a],
               label = 'flow directions'
              )
    plt.plot(path_d[:,0], path_d[:,1], 'm', label = 'Route')
    plt.plot(Roadmap_d.nodes[:,1], Roadmap_d.nodes[:,0], 'k.', label = 'Nodes of the graph', markersize = 0.8)
    plt.plot(path_d[0, 0], path_d[0, 1], 'go', label = 'start')
    plt.plot(path_d[-1, 0 ], path_d[-1, 1], 'ro', label = 'target')


    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title(r'x/y diagram of the time dependent rotating flow')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              ncol=2, fancybox=True, shadow=True)
    plt.xlabel('lat')
    plt.ylabel('lon')

    plt.subplot(2,3,5)
    halem.plot_timeseries2(path_d, time_d, Roadmap_d)
    plt.title('s/t diagram of the time dependent rotating flow')

    ax = plt.subplot(2,3,3)
    plt.axis('square')
    a = 1

    x_r = np.arange(0,1, 0.075)
    y_r = np.arange(0,1, 0.075)
    y_r, x_r = np.meshgrid(y_r,x_r)

    WD_r = griddata((Roadmap.nodes[:,1], Roadmap.nodes[:,0]), Roadmap.WD[:,10], (x_r, y_r), method= 'linear')
    u_r = griddata((Roadmap.nodes[:,1], Roadmap.nodes[:,0]), Roadmap.u[:,0], (x_r, y_r), method= 'linear')
    v_r = griddata((Roadmap.nodes[:,1], Roadmap.nodes[:,0]), Roadmap.v[:,0], (x_r, y_r), method= 'linear')

    im = plt.contourf(x_r,y_r,WD_r, cval)
    fig.colorbar(im, ax=ax, label = 'Waterdepth in meters')


    plt.quiver(x_r[::a,::a],
               y_r[::a,::a],
               u_r[::a,::a],
               v_r[::a,::a],
               label = 'flow directions'
              )
    plt.plot(path_w[:,0], path_w[:,1], 'm', label = 'Route')
    plt.plot(Roadmap.nodes[:,1], Roadmap.nodes[:,0], 'k.', label = 'Nodes of the graph', markersize = 1)
    plt.plot(path_w[0, 0], path_w[0, 1], 'go', label = 'start')
    plt.plot(path_w[-1, 0 ], path_w[-1, 1], 'ro', label = 'target')


    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.title(r"x/y diagram of the flooding and drying flow")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              ncol=2, fancybox=True, shadow=True)
    plt.xlabel('lat')
    plt.ylabel('lon')

    ax = plt.subplot(2,3,6)
    halem.plot_timeseries2(path_w, time_w, Roadmap)
    plt.title('s/t diagram of the flooding and drying flow')


    plt.savefig('D:testcases', dpi = 200)
    plt.show()


.. image:: Results2.png



Flow class for dd Zuno-v4 hirlam-kf (DCSMv6 zunov4)
--------------------------------------------------

Data available at

'http://noos.matroos.rws.nl'_ 

.. _http://noos.matroos.rws.nl: http://noos.matroos.rws.nl

With this class real hydrodynamic data can be used to optimize real shipping routes. 

.. code:: python3

    class flow_NOOS():
        def __init__(self, name):
            nc = Dataset(name)
            x_domain = (250,400)                        # general-waddden sea (250, 380)
            y_domain = (450,760)                        # (530,760)
            # x_domain = (300,390)                      # Texel-case
            # y_domain = (650,760)

            v = nc.variables['VELV'][:,:,:]
            u = nc.variables['VELU'][:,:,:]
            d = nc.variables['SEP'][:,:,:]
            x = nc.variables['x'][:,:]
            y = nc.variables['y'][:,:]
            t = nc.variables['time'][:]
            t = t *60
            x = x[x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
            y = y[x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
            u = u[:,x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
            v = v[:,x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
            d = d[:,x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]

            x_temp = ma.array(x.reshape(x.size))
            y_temp = ma.array(y.reshape(x.size))

            nodes = np.zeros((y_temp[y_temp.mask == False].size,2))
            nodes[:,0] = y_temp[y_temp.mask == False]
            nodes[:,1] = x_temp[y_temp.mask == False]
            print('1/3')

            bat, nodesb = self.bat()
            Db_new = griddata((nodesb[:,1],nodesb[:,0]), bat, (x,y), method='cubic')

            WD = d * 0
            for i in range(d.shape[0]):
                WD[i,:,:] = d[i,:,:] - Db_new

            print('2/3')

            u_n = []
            v_n = []
            d_n = []

            for node in nodes:
                xloc = np.argwhere(x == node[1])[0,1]
                yloc = np.argwhere(y == node[0])[0,0]
                u_n.append(u[:,yloc,xloc])
                v_n.append(v[:,yloc,xloc])
                d_n.append(WD[:,yloc,xloc])

            d_n = np.array(d_n)
            d_n[d_n < -600] = 0
            v_n = np.array(v_n)
            v_n[v_n < -600] = 0
            u_n = np.array(u_n)
            u_n[u_n < -600] = 0
            
            self.nodes = nodes
            self.u = np.transpose(u_n)
            self.v = np.transpose(v_n)
            self.WD = np.transpose(d_n)
            self.tria = Delaunay(nodes)
            self.t = t

	
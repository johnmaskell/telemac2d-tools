from tidal.tide import VALIDOBJECT

print(VALIDOBJECT.__doc__)

optionals = {
    "respath" : "C:/North_Sea_Model_v1/",#results file location
    "startdate" : "2007/01/01 00:00:00",#model start date
    "obspath" : "C:/Observations/sample/subsample/",#observations file locations
    "locations" : ['all'],#plot all observation location in model domain
    "variable" : "FREE SURFACE",#variable to plot
    "sampledist" : 5000.0,#find fit within distance from tide gauge (metres)
    "timestep" : 30.0,#model time step seconds
    "spinup" : 72.0,#ignore period after cold start (hours)
    "lprintout" : 120,#model output period (no. time steps)    
    "meshproj" : "epsg:32630"#mesh projection
    }
kwargs = optionals
resfile = "North_Sea_Tide_fric77_utm_newsource_jan01tojan152007.slf"#Telemac results file (selafin)


newvalobj = VALIDOBJECT(resfile,**kwargs)
newvalobj.plotTS()

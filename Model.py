import transit
import numpy as np

def lightcurve(star,planet,t):
    
    """
    Generate a light curve for the given values.
    Star    - [mass,radius] in stellar radii
    Planet  - [radius,mass,period,t0]
    t       - the array for time points from the data
    """    
    
    # Build the transiting system.
    star = transit.Central(star[0],star[1])
    
    s = transit.System(star)
    
    #body = transit.Body(r=0.11/0.74,m=0.002375/0.74,period=1.337,to=0.99,b=0.2,e=0.0167)
    body = transit.Body(r=planet[0],mass=planet[1], period=planet[2], t0=planet[3], b=0.029)
    s.add_body(body)
    
    # Compute the light curve for given values of t
    f = s.light_curve(t)
    
    return f
    
#lightcurve([0.74,0.713],[0.1,0.001,1.33712,1], np.arange(0,2,1e-4))
#    
#star = transit.Central(0.74,0.713)
#s = transit.System(star)
#body = transit.Body(r=0.1,mass=0.01, period=1.33712, t0=0, b=0)
#s.add_body(body)
#
#t = np.arange(0, 2, 1e-4)
#
#f = s.light_curve(t)
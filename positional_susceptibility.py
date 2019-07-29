#!/usr/bin/env python3

#from Jadrich code via. Truskett group webpage
'''TODO:

* fix the imports to use the style "import numpy as np" in order to 
  match the philosophy "Explicit is better than implicit"
  subsequent uses of the functions would need to be flagged np.function

* document this code more completely

* Jadrich uses a dict to format his data structures.  the coords are hiding in there.  We want this to work on our basic format.  

id type      x     y     z  
0     1    0.0   0.0   0.0
1     1    1.0   0.0   0.0
2     1    0.0   1.0   0.0
  

note that in a subroutine, Jadrich creates frames with the following:

frames = []

#static quantities
box = traj[0].configuration.box
D = traj[0].configuration.dimensions
N = len(traj[0].particles.position)
    
#loop over the configurations and stor in new format
for snap in traj:
    #dynamic quantities
    diameters = snap.particles.diameter
    coords = snap.particles.position
     
        
    #store the three lengths in a vector
    L = box[0:D]
        
    #get the particle types
    possible_types = snap.particles.types
    types = array([possible_types[type_id] for type_id in snap.particles.typeid])
        
    #replace with random positions if randomize is selected (for comparing to randomized PCA result and useful information content)
    if randomize:
        coords = L*rand(N, D) - L/2.0
            
    #remove a component from the trajectory
    for removal_type in remove_types:
        coords = coords[types == removal_type]
        diameters = diameters[types == removal_type]
        types = types[types == removal_type]
        
    #create our new data structure and shift to upper right quadrant
    frames.append({'coords': (coords[:,0:D]+L/2.0), 'diameters': diameters, 'types': types, 'L': L, 'D': D})

'''

from scipy.special import j0
from numpy import histogram, pi, power, rint
from numpy.linalg import norm
from numpy import trapz

#function for calculating the rdf
def RDF2D(frames, dr):
    M = len(frames)
    N = float(len(frames[0]['coords']))
    L = frames[0]['L']
    rho = float((N-1))/power(L, 2)
    r_edg = arange(0.0, L/2.0, dr)
    r = (r_edg[:-1] + r_edg[1:])/2.0
    hist = 0.0*r
    
    #loop over the frames
    for frame in frames:
        coords = frame['coords']

        #loop over each particle
        for coord in coords:
            #nearest neighbor coordinate wrapping
            Rpj = coord - coords
            Rpj = Rpj - rint(Rpj/L)*L
            Rpj = norm(Rpj, axis=1)
            
            #calculate the histogram
            hist = hist + histogram(Rpj, bins=r_edg)[0]
    
    #normalize out the number of frames and 
    hist = hist/float(M*(N-1))
    gr = hist/((2.0*pi*r*dr)*rho)
    
    return r, gr

#function for calculating the integral of h(r)
def PositionalSuceptibility2D(frames, dr, k):
    N = float(len(frames[0]['coords']))
    L = frames[0]['L']
    rho = float(N)/power(L, 2)
    
    d = float(L)/rint(L)
    k = 2.0*pi/d
    
    #compute the rdfs and calculate S(k0)
    r, gr = RDF2D(frames, dr)
    #Sk0 = 1.0 + 2.0*pi*rho*trapz(r*(gr-1.0)*j0(k*r))
    Sk0 = 1.0 + 2.0*pi*rho*trapz(r*abs(gr-1.0))
    
    return (Sk0, r, gr)


if __name__ == "__main__":

    #frames is a specially formatted dict that we aren't using
    
    x, r, gr = PositionalSuceptibility2D(frames[0::100], 0.025, 2.0*pi)

from visual import *
import random
from numpy import *
import time


x=0.
y=0.
z=0.

m=[(x,y,z)]

# curve is part of the 'visual' module above: see: http://vpython.org/contents/docs/visual/curve.html
polymer=curve(color=color.cyan)

# update the position of the curve (m is a list of 3D points)
polymer.pos=m


# a random walk loop in 3D:

for i in range(10000):
    
    # sleep for a little at each update so we can watch each step of the polymer growth
    time.sleep(.01)
    
    # generate either a -1, or a +1, randomly:
    r=2*random.randint(0,2)-1
    #print r
    # update the x coordinate
    x=x+r

    r=2*random.randint(0,2)-1
    #print r
    y=y+r

    r=2*random.randint(0,2)-1
    #print r
    z=z+r

    # add the coordinate to the position list
    m.append((x,y,z))
    #print m
    
    # update the visualization:
    polymer.pos=m
 
    
   
    



    


#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

k=9.e9 #Nwm^2cou^-2 Aire en la Tierra

# u position vector
# q charge
# v position vector of q
def ElectrostaticField(u,q,X,Y):
    q1=-1.0
    x0 = u[0]
    y0 = u[1]
    #v-u = v[0]-u[0], v[1]-u[1]
    r_vec_x= X - x0
    r_vec_y = Y - y0
    r = np.hypot(r_vec_x,r_vec_y)
    #r = np.sqrt(r_vec.dot(r_vec)) #(r_vec)
    E0 = [k*(q/r**3) * r_vec_x , k*(q/r**3) * r_vec_y]

    q1=1.0
    x0 = u[0] + 1.0
    y0 = u[1]
    #v-u = v[0]-u[0], v[1]-u[1]
    r_vec_x= X - x0
    r_vec_y = Y - y0
    r = np.hypot(r_vec_x,r_vec_y)
    #r = np.sqrt(r_vec.dot(r_vec)) #(r_vec)
    E = [  E0[0] + (k*(q/r**3) * r_vec_x) , E0[1]+ (k*(q/r**3) * r_vec_y)]
    return E

u=np.array([-0.5,0])
q=1.0


nx, ny = 64, 64
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

#Ex=[]
#Ey=[]
#for a in x:
#    for b in y:

[Ex,Ey] = ElectrostaticField(u,q,X,Y) 

#Ex.append(E[0])
#Ey.append(E[1])
        
fig = plt.figure()
ax = fig.add_subplot(111)

ax.streamplot(x, y, Ex, Ey, linewidth=1, cmap=plt.cm.inferno,
              density=2, arrowstyle='->', arrowsize=1.5)

plt.show()
    
    

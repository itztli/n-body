#!/usr/bin/env python3
#
# n-body.py Solve the n-body problem using Newton
# 
# Copyright (C) 2019  Victor De la Luz (vdelaluz@enesmorelia.unam.mx)
#                      
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import math

G=6.674e-11         #m^3kg^-1s^-2

class Particle:
    
    def __init__(self, p, v, m, dt=1):
        self.p = p #position
        self.v = v #velocity
        self.m = m #mass
        self.dt = dt

    def setdt(self,dt):
        self.dt = dt

    def computeR(self,p1):
        r = math.sqrt( (p1[0]-self.p[0])**2 + (p1[1]-self.p[1])**2 + (p1[2]-self.p[2])**2)
        return r

    def computeU(self,p1):
        u=[0,0,0]
        i=0
        for a,b in zip(self.p,p1):
            u[i] = b - a
            i+=1
        return u
    
    #def integrate(self,dt,p1,m1):
    def integrate(self,B):
        r = self.computeR(B.p)
        u = self.computeU(B.p)

        Vx=(G*B.m*self.dt/(r**3))*u[0]
        Vy=(G*B.m*self.dt/(r**3))*u[1]
        Vz=(G*B.m*self.dt/(r**3))*u[2]

        self.v[0] += Vx
        self.v[1] += Vy
        self.v[2] += Vz
        
        self.p = [self.p[0]+ (self.v[0]) *dt,self.p[1]+ (self.v[1])*dt,self.p[2]+ (self.v[2])*dt]

    def getPosition(self):
        return self.p

    def getKineticEnergy(self):
        k= (1/2)*self.m*(math.sqrt( self.v[0]^2 +self.v[1]^2+self.v[2]^2))
        return k    
    
p0=[1.0, 0.0, 0.0]  #km
v0=[0.0, 0.0, 0.0]  #km/s
m=1.0               #kg

p1=[0.0, 0.0, 0.0]  #km
v1=[1.0, 0.0, 0.0]  #km/s
m1=1.0               #kg

dt=0.01              #sec

A = Particle(p0,v0,m)
B = Particle(p1,v1,m1)

B.setdt(dt)

for t in range(100):
    B.integrate(A)
    print(B.getPosition())






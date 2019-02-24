#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: HKÁ, AES og MBH
"""

import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

from numpy import sin, cos, pi, sqrt, sign, log10


"""
The following function finds the value of f given theta.
"""
def f(theta):
    A2 = L3*cos(theta)-x1
    B2 = L3*sin(theta)
    A3 = L2*(cos(theta)*cos(gamma) - sin(theta)*sin(gamma))-x2
    B3 = L2*(cos(theta)*sin(gamma) + sin(theta)*cos(gamma))-y2

    N1 = B3*(p2**2-p1**2-A2**2-B2**2) - B2*(p3**2-p1**2-A3**2-B3**2)
    N2 = -A3*(p2**2-p1**2-A2**2-B2**2) + A2*(p3**2-p1**2-A3**2-B3**2)

    D=2*(A2*B3-B2*A3)

    return N1**2 + N2**2 - (p1**2)*(D**2)

"""
The following function finds the corresponding values of x and y given theta.
"""
def xy(theta):
    A2 = L3*cos(theta)-x1
    B2 = L3*sin(theta)
    A3 = L2*(cos(theta)*cos(gamma) - sin(theta)*sin(gamma))-x2
    B3 = L2*(cos(theta)*sin(gamma) + sin(theta)*cos(gamma))-y2

    N1 = B3*(p2**2-p1**2-A2**2-B2**2) - B2*(p3**2-p1**2-A3**2-B3**2)
    N2 = -A3*(p2**2-p1**2-A2**2-B2**2) + A2*(p3**2-p1**2-A3**2-B3**2)

    D = 2*(A2*B3-B2*A3)

    x = N1/D
    y = N2/D

    return x,y

"""
The following function draws the function f.
"""
def drawf():
    x = np.linspace(-pi, pi, 1000)
    y = f(x)
    y0 = 0*x
    plt.plot(x, y, label='f')
    plt.plot(x, y0, label='0')
    plt.legend()
    plt.show()
    return

"""
The following function draws a logarithmic version of the function f.
"""
def drawfZeroes():
    x = np.linspace(-pi, pi, 1000)
    y = f(x)
    t = sign(y)*log10(1+abs(y)/10**3)
    y0 = 0*x
    plt.plot(x, t, label='f')
    plt.plot(x, y0, label='0')
    plt.legend()
    plt.show()
    return

"""
The following function draws a diagram of the placement of the Stewart
for a given value of theta.
"""
def drawStewart(theta):
     gamma=np.pi/2
     x,y = xy(theta)

     """ X stores the vertices of the Stewart platform while
      Y stores the colors which will become vertices"""
     X = np.array([[0,0], [x1,0], [x2, y2], [x, y], [x+L2*cos(theta+gamma), y+L2*sin(theta+gamma)], [x+L3*cos(theta), y+L3*sin(theta)]])
     Y = ['blue', 'blue', 'blue', 'red', 'red', 'red']
     plt.figure()

     """ We draw the platform in pink for clarity. """
     t1 = plt.Polygon(X[3:6,:], color='pink')
     plt.gca().add_patch(t1)

     plt.scatter(X[:, 0], X[:, 1], s = 10, color = Y[:])

     """ We draw p1, p2 and p3. """
     plt.plot([0,x],[0,y], color='blue')
     plt.plot([x1,x+L3*cos(theta)],[0,y+L3*sin(theta)],color='blue')
     plt.plot([x2,x+L2*cos(theta+gamma)],[y2,y+L2*sin(theta+gamma)],color='blue')

     """ We draw L1, L2 and L3 """
     plt.plot([x+L3*cos(theta),x],[y+L3*sin(theta),y], color='red')
     plt.plot([x+L3*cos(theta),x+L2*cos(theta+gamma)],[y+L3*sin(theta),y+L2*sin(theta+gamma)], color='red')
     plt.plot([x+L2*cos(theta+gamma),x],[y+L2*sin(theta+gamma),y], color='red')
     plt.show()

     return

"""
The following function finds the end points of the intervals which contain the roots of f.
"""
def findIntervals():
    x = np.linspace(-pi, pi, 1000)
    y = f(x)
    endPoints = []
    fSign = sign(y[0])
    for k in range(len(x)):
        if(sign(y[k]) != fSign):
            endPoints.append(x[k-1])
            endPoints.append(x[k])
            fSign = sign(y[k])
    return endPoints

"""
The following function finds the number of roots of the function f given
the vector of endpoints found by the method findIntervals.
"""
def numRoots(endPoints):
    return len(endPoints)/2


"""
Given an interval with endpoints a and b the following function finds
the root of f contained in that interval by using the secant method.
"""
def secant(a,b,tol):
    x = a
    x1 = b
    x0 = a
    while abs(f(x)) > tol:
        x = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        x0 = x1
        x1 = x
    return x

"""
Given an interval with endpoints a and b the following function finds
the root of f contained in that interval by using binary search.
"""
def binary(a,b,tol):
    while (b-a)/2 > tol:
        c = (a+b)/2
        if f(c) == 0:
            return c
        if f(a)*f(c) < 0:
            b = c
        else:
            a = c
    return c

"""
The following function finds the roots of f up to a tolerance tol
"""
def fRoots(tol):
    x = np.linspace(-pi,pi, 1000)
    y = f(x)
    roots = []
    fSign = sign(y[0])
    for k in range(len(x)):
        if sign(y[k]) != fSign:
            roots.append(binary(x[k-1],x[k],tol))
            fSign = sign(y[k])
    return roots

"""
The following function verifies if the calculated values
of (x,y,theta) correspond to the strut lengths (p1,p2,p3) up to a tolerance tol.
"""
def verifySolution(tol):
    bool = True
    A2 = L3*cos(theta)-x1
    B2 = L3*sin(theta)
    A3 = L2*(cos(theta)*cos(gamma) - sin(theta)*sin(gamma))-x2
    B3 = L2*(cos(theta)*sin(gamma) + sin(theta)*cos(gamma))-y2
    if abs(p1^2-x^2-y^2) > tol:
        bool = False
    if abs(p2^2-(x+A2)^2-(y+B2)^2) > tol:
        bool = False
    if abs(p3^2 - (x+A3)^2-(y+B3)^2) > tol:
        bool = False
    return bool

"""
The following function gives the endpoints of the intervals of p2
for which there are 0, 2, 4 and 6 poses respectivly.
"""
def findp2Intervals():
    p = np.linspace(0,14,1000)
    y = np.zeros(len(p))
    for i in range(len(p)):
        global p2
        p2 = p[i]
        y[i] = numRoots(findIntervals())
    plt.plot(p,y)
    plt.show()
    return

L1=2
L2=sqrt(2)
L3=sqrt(2)
x1=4
x2=0
y2=4
gamma=pi/2
theta=pi/4

p1=sqrt(5)
p2=sqrt(5)
p3=sqrt(5)

print(fRoots(0.5e-10))
drawf()

drawStewart(theta)

bil = findIntervals()
print(bil)
print(secant(bil[0], bil[1], 1e-10))
print(secant(bil[2], bil[3], 1e-10))
print(numRoots(findIntervals()))
findp2Intervals()

L1=3
L2=3*sqrt(2)
L3=3
x1=5
x2=0
y2=6
gamma=pi/4

p1=5
p2=5
p3=3

findp2Intervals()

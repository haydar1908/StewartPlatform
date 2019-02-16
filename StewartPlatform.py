#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 23:41:47 2019

@authors: HKÁ, AES og MBH
"""
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

from numpy import sin, cos, pi, sqrt


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
The following function finds the corresponding value of x and y given theta.
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
The following function finds the neighbourhoods of roots of f
"""
def finnaBil():
    x = np.linspace(-pi, pi, 10000)
    y = f(x)
    endaPunktarBila = []
    formerki = np.sign(y[0])
    for k in xrange(len(x)):
        if(np.sign(y[k]) != formerki):
            endaPunktarBila.append(x[k-1])
            endaPunktarBila.append(x[k])
            formerki = np.sign(y[k])
    return endaPunktarBila

"""
The following function finds the roots of f by using the secant method.
"""
def Secant(end1,end2,tol):
    x = end2
    x1 = end2
    x0 = end1
    while abs(f(x)) > tol:
        x = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        x0 = x1
        x1 = x
    return x

"""
The following function draws the corresponding StewartPlatform
"""
def drawStewartPlatform():
     x_number_list = np.array([0, x1, x2, x, x+L2*cos(theta+gamma), x+L3*cos(theta)])
     y_number_list = np.array([0, 0, y2, y, y+L2*sin(theta+gamma), y+L3*sin(theta)])
     plt.plot(x_number_list, y_number_list, s=1000)
     plt.show()
     return

def drawStewart(theta):
     gamma=np.pi/2
     x,y = xy(theta)
     X = np.array([[0,0], [x1,0], [x2, y2], [x, y], [x+L2*cos(theta+gamma), y+L2*sin(theta+gamma)], [x+L3*cos(theta), y+L3*sin(theta)]])
     Y = ['blue', 'blue', 'blue', 'red', 'red', 'red']
     plt.figure()
     t1 = plt.Polygon(X[3:6,:], color='pink')
     plt.scatter(X[:, 0], X[:, 1], s = 10, color = Y[:])
     plt.gca().add_patch(t1)

     #lögum línuvigurinn til að teikna rétt


     #plt.plot(X)
     plt.show()
     #plt.plot(x_number_list, y_number_list, 'ro')

     return

L1=2
L2=np.sqrt(2)
L3=np.sqrt(2)
x1=4
x2=0
y2=4
gamma=np.pi/2
theta=-pi/4

p1=sqrt(5)
p2=sqrt(5)
p3=sqrt(5)

print(L1)
print(p1)
print(f(-np.pi/4))
print(f(0))
print(f(np.pi))
drawf()
#drawStewartPlatform()

drawStewart(theta)

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = range(9)
y_points = range(9)
p = ax.plot(x_points, y_points, 'b')
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()
bil = finnaBil()
print(bil)
print(Secant(bil[0], bil[1], 1e-10))
print(Secant(bil[2], bil[3], 1e-10))

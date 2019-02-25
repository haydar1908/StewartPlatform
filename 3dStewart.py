#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@authors: HKÁ, AES og MBH
"""

import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

from numpy import sin, cos, pi, sqrt, sign, log10

def fi(tx, ty, tz, psi, theta, phi, i):
    t = np.array([tx,ty,tz])
    R = RMatrix(psi,theta,phi)
    return lin.norm(t+np.matmul(R,a[i])-c[i]) - l[i]

def F(tx,ty,tz,psi,theta,phi):
    F = np.array(np.zeros(6))
    for i in range(6):
        F[i] = fi(tx,ty,tz,psi,theta,phi,i)
    return F

def Broyden1(k,x0,x1):
    n = len(x0)
    A = np.identity(n)
    for n in range(k):
        delta = x1-x0
        dell = F(x1[0], x1[1],x1[2], x1[3],x1[4], x1[5]) - F(x0[0], x0[1],x0[2], x0[3],x0[4], x0[5])
        deltadot = np.dot(delta,delta)
        A = A + np.matmul((dell - np.matmul(A,delta)),delta)/deltadot
        x = x1 - np.matmul(np.linalg.inv(A),F(x1[0], x1[1],x1[2], x1[3],x1[4], x1[5]));
        x0 = x1
        x1 = x;
    return x

def RMatrix(psi, theta, phi):
    R = np.array([[cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi), cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi), sin(psi)*sin(theta)], [-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi), -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi), cos(psi)*sin(theta)], [sin(theta)*sin(phi), -sin(theta)*cos(phi), cos(theta)]])
    return R

l = np.array([1, 1, 1, 1, 1, 1])
c = np.array([[3,0,0],[1,2,0],[-1,2,0],[-3,0,0],[-1,-2,0],[1,-2,0]])
a = np.array([[2,0,0],[1,1,0],[-1,1,0],[-2,0,0],[-1,-1,0],[1,-1,0]])

k = 100
x0 = np.array([0,0,0.8,0.5,0.5,0.5])
x1 = np.array([0.2,0.2,1,0.3,0.3,0.3])
print(Broyden1(k, x0, x1))

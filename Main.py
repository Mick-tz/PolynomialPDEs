#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 04:54:06 2018

Problem was stated as follows:
    A company must reduce prices by 10%. Profit will be modelled as the charge
    y minus a penalty P(t,dotY) (where dotY is the first derivative of y wr to t
P is given by:
        P(t,dotY) = alpha*dotY**2 + beta*(t**2 - 1)*dotY**3
    We find y that maximizes profif by minimizing:
        S(y) = integrate(P-y,0,1) w.r. to t
        
    The lagrangian for this problem is:
        L(t,y,dotY) = P(t,dotY) - y                         ...(1)
        
    with boundary conditions:
        y(0) = 1 and y(1) = 0.9                              ...(2)
    
for S to be minimal, we calculate Euler-Lagrange equation (e/l-eq) and equalize
it to 0:
        
    2*alpha*doubleDotY + 6*beta*(t*dotY**2 + (t**2-1)*dotY*doubleDotY) = 0
        
        where doubleDotY refers to y's second derivative w.r. to t
        
We explain Differential Transform Method as follows:
    objetive is to approximate y(t) using a it's taylor polynomials.
    
    for a function w(t), we define for every k>=0, W(k) as the kth derivative of
    w(t) w.r. to t evaluated in t0=0 and devided by k!
    
    We'll call W(k) the k-th differential transform of w(t)
    
    it follows that:
        w(t) is approximated by W(0) + W(1)*t + W(2)*t**2 + ... + W(k)*t**k
    
    using (e/l)-eq, we can derive a recursive method to obtain Y(k) in terms
    of Y(1) (which we'll call gamma)
    This is implemented within DTM file. Lagrangian and euler-lagrange methods
    for this particular set of problems are implemented in PoliticalPricing file.

For this model to approximate an optimal solution. Function y is assumed to be analytic within 0 and 1
    
@author: Amir
"""

import matplotlib.pyplot as plt
import numpy as np
import DTM
import PoliticalPricing as pp


domain = np.linspace(0,1,30)
k = 10       #Degree of polynomial


#case1 alpha = beta = 0
alpha = 5.0
beta = 5.0


def line(t):
    """
    A straight line connecting (0,1) and (1,0.9)
    """
    return 1.0 - 0.1*t

#fits the paramters of the taylor polynomial
gamma = DTM.findGammaNumerical(k,alpha,beta)
Ys = DTM.diffTrans(k,alpha,beta,False,gamma)

def myModel(t,k=k, Ys=Ys):
    return DTM.model(t,k=k,Ys=Ys)

plt.figure(1)
plt.plot(domain,myModel(domain))
plt.plot(domain,line(domain))
plt.title("Fitted Model vs Straight Line")

plt.figure(2)
plt.plot(domain, pp.eulerEquation(domain,myModel,alpha,beta), label ="My Model")
plt.plot(domain, pp.eulerEquation(domain,line,alpha,beta), label = "Straight line")
plt.plot(domain,[0 for i in domain])
plt.title("Euler-Lagrange equation when alpha=beta=5")
plt.legend(loc = "best")
plt.legend()

sModel = pp.integLag(myModel,alpha,beta)[0]
sLine = pp.integLag(line,alpha,beta)[0]

print("When alpha and beta are both equal to 5.")
print("The integral of the lagrangian in our model is: ")
print(sModel)
print("the integral of the lagrangian using the straight line is: " )
print(sLine)
print("\n")

#case2: alpha != beta

alpha = 7.0/4.0
beta = 5.0

#fits the paramters of the taylor polynomial
#for k in [5,10,15,20,50]:
for k in [5,7]:
    gamma = DTM.findGammaNumerical(k,alpha,beta)
    #gamma = -0.061787462699583774
#    print(gamma)
    Ys = DTM.diffTrans(k,alpha,beta,False,gamma)
#    print(Ys)

    def myModel(t,Ys=Ys):
        """
        The model fitted
        """
        return DTM.model(t, Ys=Ys)
    
#    print(myModel(1))
    
    plt.figure(k+2)
    plt.plot(domain,myModel(domain))
    plt.plot(domain,line(domain))
    plt.title("Fitted Model vs Straight Line when alpha != beta in k=" + str(k))

    plt.figure(k+3)
    plt.plot(domain, pp.eulerEquation(domain,myModel,alpha,beta), label = "My Model")
    plt.plot(domain, pp.eulerEquation(domain,line,alpha,beta), label = "Straight line")
    plt.plot(domain,[0 for i in domain])
    plt.title("Euler-Lagrange equation when alpha=7/4 and beta=5")
    plt.legend()
    
    sModel = pp.integLag(myModel,alpha,beta)[0]
    sLine = pp.integLag(line,alpha,beta)[0]

    print("When alpha equals 7/4 and beta equals 5.")
    print("The integral of the lagrangian in our model is: ")
    print(sModel)
    print("the integral of the lagrangian using the straight line is: " )
    print(sLine)
    print("\n")

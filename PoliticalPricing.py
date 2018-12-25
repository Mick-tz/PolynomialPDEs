#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 19:38:34 2018

This sheet is thought to be an example on how to use DTM for a particular case. 

Problem set is proposed in pictures within the repositories

@author: Chino
"""

from scipy.misc import derivative
from scipy.integrate import quad

def P(t, y, alpha, beta):
    """
    penalty function associated with costumer's tendancy to
    wait for the product's price to go down
    Args:
        t: time (a float between 0 and 1)
        y: profit function
        alpha,beta: a float
    """
    dotY = derivative(y, t , dx = 0.01)
    return alpha*(dotY**2) + beta*(t**2-1)*(dotY**3)

def Lag(t, y, alpha, beta):
    """
    Lagrangian for the given problem set. Associates production with the penalty
    function.
    Args:
        t: a float
        y: profit function (of t)
    """
    return P(t, y, alpha, beta) - y(t)

def eulerEquation(t, y, alpha, beta):
    """
    Euler-Lagrange equation (0 at stationary points) for the given problem set 
    (derivatives have been calculated by hand)
    Args:
        t: time (a float between 0 and 1)
        y: profit function
        alpha,beta: a float
    """
    dotY = derivative(y,t, dx = 0.01)
    doubleDotY = derivative(y, t, dx = 0.01, n = 2)    #second derivative of y wrt time
    return 2*alpha*doubleDotY + 6*beta*(t*dotY**2+(t**2 - 1)*dotY*doubleDotY) + 1

def integLag(y, alpha, beta):
    """
    numerical integration from 0 to 1 of the lagrangian associated with function
    y. alpha anf beta are to be passed to the lagrangian implementation.
    Args:
        y:profit function
        alpha,beta: a float
    """
    return quad(Lag,0,1,args=(y,alpha,beta))

def recursiveY(k, gamma, alpha, beta):
    """
    
    """
    if (0 == k):
        return [1]            #Y(0)=y(0)=1
    elif (1 == k):
        return [1, gamma]
    elif (2 == k):
        return [1, gamma, 1/(4*(3*beta*gamma-alpha))]
    elif (2 < k):
        result = 0
        r = k-2     #r stated in (3)
        previous = recursiveY(k-1, alpha, beta, gamma)
        delta = 2*(r+1)*(r+2)*(3*beta*gamma-alpha)
        addTerm = r*previous[r]*gamma
        addTerm = addTerm - 2*r*(r+1)*previous[r+1]*previous[2]
        summatory = 0
        for m in range(r-1):        #range goes up to one number before
            summatory += (m+1)*(m+2)*previous[m+2]*(r-m-1)*previous[r-m-1]
            summatory += (m+1)*(r-m)*previous[m+1]*previous[r-m]
            summatory -= (m+1)*(m+2)*previous[m+2]*(r-m+1)*previous[r-m+1]
        result += 6*beta*(summatory + addTerm)
        result /= delta
        previous.append(result)
        return previous

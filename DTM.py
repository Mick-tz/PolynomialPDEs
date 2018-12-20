#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 21:41:40 2018

This implementation of the differential transform method (DTM) is thought to comply 
with the following problem:
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
    
    for S to be minimal, we require Euler-Lagrange equation (e/l-eq) is given by:
        2*alpha*doubleDotY + 6*beta*(t*dotY**2 + (t**2-1)*dotY*doubleDotY) = 0
        
        where doubleDotY refers to y's second derivative w.r. to t
        
    In order to use DTM, we have calculated the e/l-eq differential r-th
    transform as follows:
        
      0=E(r) = delta_0(r) + 2*(r+1)*(r+2)*Y(r+2)*(alpha-3*beta*Y(1)) +
         6*beta(r*Y(r)*Y(1)-2*r*(r+1)*Y(r+1)*Y(2) +
                sum^{r-2}_0((m+1)*(m+2)*Y(m+2)*(r-m-1)*Y(r-m-1) +
                    (m+1)*Y(m+1)*(r-m)*Y(r-m) -
                        (m+1)*(m+2)*Y(m+2)*(r-m+1)*Y(r-m+1)))       ...(3)

    where delta_0(r) is equal to 1 when r=0 and 0 everywhere else
    
    using equation (3) we are able to find Y(r+2) for every r>=0

    E(r) equals 1/r!*(d^r e(0)/dt^r) where e(t) refers to e/l-eq
    for more information on how to derive (3) we refer the reader to:
    https://pdfs.semanticscholar.org/1836/5d9d7295bb58a7efbfe8d6bb3e4e2d7b3554.pdf
    
    throughout this code, gamma refers to a numeric approximation of Y(1).
    alpha and beta refer to given values in P(t,dotY)
@author: Amir
"""

#import numpy as np
#import scipy as sp
from sympy.solvers import solve
from sympy import Symbol
from scipy import optimize



def diffTrans(k, alpha, beta, symbolic = True, gama = 0.0):
    """
    kth differential transform of the function y(t) satisfy (3). This helper
    function is though to work only when t0 = 0 (as in (3))
    Args:
        k: an int
        alpha: a float
        beta: a float
        previous: a list with values of the previous terms
        symbolic: if true, transform will return a symbolic value, returns a 
        numeric value otherwise
        gama: a float (the numerica value of gamma)
    Return:
        Y(k) in terms of gamma = Y(1)
    """
    if (0 == k):
        return [1]            #Y(0)=y(0)=1
    elif (1 == k):
        if symbolic:
            return [1, Symbol("gamma")]      #gamma will be the fit of our model
        else:
            return [1, gama]
    elif (2 == k):
        if (symbolic):
            return [1 , Symbol("gamma"), (1/(4*(3*beta*Symbol("gamma")-alpha)))]  #calculated by hand
        else:
            return [1, gama, 1/(4*(3*beta*gama-alpha))]
    elif (2 < k):
        result = 0      
        if symbolic:
            gamma = Symbol("gamma")
        else:
            gamma = gama
        r = k-2     #r stated in (3)
        previous = diffTrans(k-1, alpha, beta, symbolic, gama)
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
    
    
def findGamma(k, alpha, beta):
    """
    Obtains the first k differential transforms from the function y(t).
    Uses this and the fact that $0.9=y(1)\sim\sum_k Y(k)*1**k$ to find Y(1)
    Args:
        k: an int greater than 1 (greater than 6 for a good approx)
        alpha: a float
        beta: a float
    Return:
        approximate value of gamma
    """
    equation = -0.9      #0.9= 1 + sum_{k>1}Y(k), to use sympy solver equation = 0
    for Y in diffTrans(k,alpha,beta):
        equation += Y
    gammas = solve(equation,Symbol("gamma")) #normally more than 1 gamma is find
    return float(gammas[0])    #first value works best
    

def evalGammaNumerical(gamma, k, alpha, beta):
    """
    Args:
        k: an int greater than 1 (greater than 6 for a good approx)
        alpha: a float
        beta: a float
        gamma:
    Return:
        sum over the first k e/l-eq's differential transforms evaluated in gamma
    """
    equation = -0.9 #0.9= 1 + sum_{k>1}Y(k), to use sympy solver equation = 0
    equation += sum(diffTrans(k,alpha,beta,False,gamma))
    return equation
    

def findGammaNumerical(k, alpha, beta, x0 = -0.5):
    """
    Args:
        k: an int greater than 1 (greater than 6 for a good approx)
        alpha: a float
        beta: a float
        gamma:
    Return:
        sum over the first k e/l-eq's differential transforms evaluated in gamma
    """
    root = optimize.newton(evalGammaNumerical, x0, args = (k,alpha,beta), maxiter = 1500)
    return root
    
def fitModel(k, alpha, beta):
    """
    Fits a numeric model using taylor polynoms
    Args:
        k: number of terms (an int greater than 1)
        alpha: a float
        beta: a float
    return:
        a list with Y(m) for every m<=k
    """
    gamma = findGamma(k, alpha, beta)
    return diffTrans(k, alpha, beta, False, gamma)

def model(t, alpha = 5.0, beta = 5.0, k = 6, Ys = []):
    """
    Uses fitModel to find the correspondant Y(m)'s for m<= k then evaluates the
    model at time 0<=t<=1
    Args:
        t: time to evaluate the numeric approx (a float between 0 and 1 inclusive)
        alpha: a float
        beta: a float
        k: number of terms 
    """
    if len(Ys) == 0:
        Ys = fitModel(k, alpha, beta)
    result = 0
    for r in range(len(Ys)):
        result = result + Ys[r]*(t**r)
    return result

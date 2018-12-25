#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 21:41:40 2018

This implementation of the differential transform method (DTM)

    for more information on how to derive (3) we refer the reader to:
    https://pdfs.semanticscholar.org/1836/5d9d7295bb58a7efbfe8d6bb3e4e2d7b3554.pdf
    
    throughout this code, gamma refers to a numeric approximation of Y(1).
    alpha and beta refer to given values in P(t,dotY)
@author: Chino
"""

from sympy.solvers import solve
from sympy import Symbol
from scipy import optimize


def differential_transform(k, Y, args, symbolic = True, gama = 0.0):
    """
    kth differential transform of the function y(t). 
    This helper function is though to work only when t0 = 0
    Relies on the fact that Y is a recursive function that returns recursively
    the specific k-th diff trasnsform. 
    Y must receive values as follows, Y(k, gamma, args) where gamma is an aprox
    to Y(1) and should return a list with every Y(m) for m <= k.
    Depending on the value of symbolic, second argument will for Y will be either
    symbolic (using sympy) or nummeric
    
    Args:
        k: an int
        Y: a function that returns Y(k) when (obtained using DTM over euler-lagrange equation)
        args: a tuple containing any extra arguments for Y ordered as they must be introduced
        symbolic: if true, transform will return a symbolic value (using symbol "gamma" from sympy)
            a numeric value otherwise
        gama: a float (the numerical value of gamma)
    Return:
        a list with every Y(j) in terms of gamma = Y(1) for j<=k
    """
    if symbolic:
        gamma = Symbol("gamma")
        return Y(k, gamma, *args)
    else:
        gamma = gama
        return Y(k, gamma, *args)
    
    
def polynomial_model(t, Ys):
    """
    
    Args:
        t: parameter for polynomial
        Ys: a list with the terms for the polynomial in order
    """
    n = len(Ys)
    result = 0
    for j in range(n):
        result += Ys[j]*t**n
    return result


def findGammas(k, Y, args, y1, t1=1):
    """
    Obtains and approximation to gamma using sympy solver

    Args:
        k: an int greater than 1 (greater than 6 for a good approx)
        y1: actual value at second boundary
        Y: a function that returns Y(k) when (obtained using DTM over euler-lagrange equation)
        args: a tuple containing any extra arguments for Y ordered as they must be introduced
        t1: parameter of second boundary
    Return:
        approximate value of gamma
    """
    equation = -y1
    Ys = differential_transform(k, Y, args)
    equation += polynomial_model(t1, Ys)
    gammas = solve(equation,Symbol("gamma")) #normally more than 1 gamma is found
    return gammas    


def evalGammaNumerical(k, Y, args, gamma, y1, t1=1):
    """
    
        k: an int greater than 1 (greater than 6 for a good approx)
        y1: actual value at second boundary
        Y: a function that returns Y(k) when (obtained using DTM over euler-lagrange equation)
        args: a tuple containing any extra arguments for Y ordered as they must be introduced
        t1: parameter of second boundary
    Return:
        sum over the first k e/l-eq's differential transforms evaluated in gamma
    """
    equation = -y1 #0.9= 1 + sum_{k>1}Y(k), to use sympy solver equation = 0
    Ys = differential_transform(k, Y, args, False, gamma)
    equation += polynomial_model(t1, Ys)
    return equation
    

def findGammaNumerical(k, Y, args, y1, t1=1, gamma0 = -0.5):
    """
    Args:
        k: an int greater than 1 (greater than 6 for a good approx)
        y1: actual value at second boundary
        Y: a function that returns Y(k) when (obtained using DTM over euler-lagrange equation)
        args: a tuple containing any extra arguments for Y ordered as they must be introduced
        t1: parameter of second boundary
        gamma0:
    Return:
        sum over the first k e/l-eq's differential transforms evaluated in gamma
    """
    root = optimize.newton(evalGammaNumerical, gamma0, args = (k, Y, args, y1, t1),
                           maxiter = 1500)
    return root
    

def fitModel(k, Y, args, y1, t1=1):
    """
    Fits a numeric model using taylor polynoms
    Args:
        k: an int greater than 1 (greater than 6 for a good approx)
        y1: actual value at second boundary
        Y: a function that returns Y(k) when (obtained using DTM over euler-lagrange equation)
        args: a tuple containing any extra arguments for Y ordered as they must be introduced
        t1: parameter of second boundary
    return:
        a list with the fitted Y(m) for every m<=k
    """
    gammas = findGammas(k, Y, args, y1, t1)
    return differential_transform(k, Y, args, False, gammas[0])

def integrate_lagrangian(L, args, t1=1):
    """
    numerical integreation from lagrangian L between 0 and t1
    Args:
        L: the lagrangian (as a function of t and some args)
        args: extra arguments for L in the order they appear
        t1: parameter at second boundary
    """
    return quad(L, 0, t1, args)



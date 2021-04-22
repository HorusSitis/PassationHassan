#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:57 2020

@author: amoreau
"""
# import dolfin

import numpy as np


def scalar_prod(u, v):
    """
    return the scalar product of u and v
    """
    # if pars['scalar_product'] == 'L2':
    #     form = inner(u, v)*dx
    # elif pars['scalar_product'] == 'H1':
    #     form = (dot(u, v) + inner(grad(u), grad(v)))*dx
    # else:
    #     raise NotImplementedError()
    # form = (dot(u, v) + inner(grad(u), grad(v)))*dx
    sc = assemble((dot(u, v) + inner(grad(u), grad(v)))*dx)

    return sc


def sq_dist(u, v):
    """
    return the squared distance betweeen u and v according to scalar_prod
    """
    sq_norm_u = scalar_prod(u, u)
    sq_norm_v = scalar_prod(v, v)

    sq_dist_uv = sq_norm_u + sq_norm_v - 2*scalar_prod(u, v)

    return sq_dist_uv

def var_rel_func(u, v):
    """
    returns the relative variance of the sample of two functions u and v
    """
    sq_norm_u = scalar_prod(u, u)
    sq_norm_v = scalar_prod(v, v)

    sq_dist_uv = sq_dist(u, v)

    sc_prod_uv = scalar_prod(u, v)

    var_rel = (2*sc_prod_uv + sq_dist_uv)/(0.5*(sq_norm_u + sq_norm_v) + sc_prod_uv) - 1

    return(var_rel)

def var_rel_scal(nu, mu):
    """
    returns the relative variance of the sample of two functions u and v
    """
    var_rel_scal = 2*(nu**2 + mu**2)/((nu + mu)**2) - 1

    return(var_rel_scal)

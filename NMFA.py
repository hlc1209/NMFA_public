'''
Deciphering the Generating Rules and Functionalities of Complex Networks
Author X. Xiao, H. Chen, P. Bogdan
Oct. 31 2020
The complete code will be made public after the article is published
'''

from scipy.stats import linregress
import numpy as np
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool
from tqdm import tqdm
import time

#The parts of the box-growing method and calculation of the sum function and tau will be made public after the article is published

def is_isomorphic(G1,G2,Q = [q/100 for q in range(-2000,2001,10)]):
    ntauls1 = np.array(nfd_nk(G1,Q))
    ntauls2 = np.array(nfd_nk(G2,Q))
    if np.array_equal(ntauls1, ntauls2):
        return True
    else:
        return False

def nspectrum(tau_list,q_list=[q/100 for q in range(-2000,2001,10)],name=None,linewidth=4):
    al_list = []
    fal_list = []
    for i in range(1,len(q_list)):
        al=(tau_list[i]-tau_list[i-1])/(q_list[i]-q_list[i-1])
        al_list.append(al)
    for j in range(len(q_list)-1):
        fal=q_list[j]*al_list[j]-tau_list[j]
        fal_list.append(fal)
    if name is not None:
        plt.plot(al_list,fal_list,label=name,linewidth=linewidth)
        plt.legend()
    else:
        plt.plot(al_list,fal_list,linewidth=linewidth)
    plt.xlabel('Lipschitz-Holder exponent, 'r'$\alpha$')
    plt.ylabel('Multi-fractal spectrum, 'r'$f(\alpha)$')
    alpha_0 = np.argmax(fal_list)
    asymmetry_ratio = np.log((al_list[-1] - al_list[alpha_0])/(al_list[alpha_0]-al_list[0]))
    return al_list,fal_list, asymmetry_ratio


def ndimension(tau_list,q_list=[q/100 for q in range(-2000,2001,10)],linewidth=4,name=None):
    dim_list = []
    qd_list = []
    for i in range(len(q_list)):
        if q_list[i] != 0:
            dim = tau_list[i]/q_list[i]
            dim_list.append(dim)
            qd_list.append(q_list[i])
    if name is not None:
        plt.plot(qd_list,dim_list,label=name,linewidth=linewidth)
        plt.legend()
    else:
        plt.plot(qd_list,dim_list,linewidth=linewidth)
    plt.xlabel('Distortion exponent, 'r'$q$')
    plt.ylabel('Generalized fractal dimension, 'r'$D(q)$')
    return dim_list, qd_list

def distance(ndim_list):
    dn = len(ndim_list)
    dist = np.zeros((dn,dn))
    for i in range(dn):
        for j in range(dn):
            d1 = np.array(ndim_list[i])
            d2 = np.array(ndim_list[j])
            dist[i][j] = np.sqrt(np.sum((d1-d2)**2)/len(ndim_list[i]))
    return dist


def nheat(tau_list, q_list=[q/100 for q in range(-2000,2001,10)]):
    heat_list = []
    for i in range(2,len(tau_list)):
        heat = -100 * (tau_list[i] - 2*tau_list[i-1] + tau_list[i-2])
        heat_list.append(heat)
    return heat_list, q_list[2:]

#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!pip install igraph
import igraph as ig
import pylab
import random
import matplotlib.pyplot as plt
import math
import statistics
from statistics import mean
from statistics import median
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import random
#shape,loc,scale=scipy.stats.lognorm.fit(samples,floc=0)
clr="#EFEFEF"
import time
import pandas as pd
from random import shuffle
import numpy as np
import csv
import numpy as np
from numpy.linalg import eig


# In[6]:


def epi_nbd_index(G_1,G_2,in_v_r_lst):
    v_r_lst = in_v_r_lst
    f_prob_s_low = [] #final low sus probability lst G1
    f_prob_v_low = [] #final low vac probability lst G1
    f_prob_v_med = [] #final med vac probability lst G1
    f_prob_v_high = [] #final high vac probability lst G1
    f_prob_low = [] #final low probability lst G2
    f_l_lst_nbd = []  #final low nbs lst
    f_m_lst_nbd = []  #final medium nbs lst
    f_h_lst_nbd = []  #final high nbs lst
    f_s_lst_nbd = []  #final sus nbs lst
    f_v_lst_nbd = []  #final vac nbs lst
    f_t_lst_nbd = []  #final total nbs lst
    f_t_lst_nbd_g2 = []  #final total nbs lst G2
    i_n_t_v=[]    #temp number of total neigbors in each node
    i_n_t_v_g2 =[]    #temp number of total neigbors in each node G2
    for v in G_1.vs:
        i_l_s_v=[]    #temp list of index of sus in each node
        i_l_v_v=[]    #temp list of index of vac in each node
        f_s=list(G_1.vs[G_1.neighbors(v)](name='S'))
        for i in range(len(f_s)):
            i_l_s=f_s[i].index
            i_l_s_v.append(i_l_s)
        f_v=list(G_1.vs[G_1.neighbors(v)](name='V'))
        for i in range(len(f_v)):
            i_l_v=f_v[i].index
            i_l_v_v.append(i_l_v)
        f_total=len(list(G_1.vs[G_1.neighbors(v)]))
        i_n_t_v.append(f_total)
        t_l_lst_nbd = []  #list of index of low in each node
        t_m_lst_nbd = []  #list of index of medi in each node
        t_h_lst_nbd = []  #list of index of high in each node
        neis = G_2.vs[G_2.neighbors(v)]
        f_total_g2 = len(list(G_2.vs[G_2.neighbors(v)]))
        i_n_t_v_g2.append(f_total_g2)
        for i in range(len(list(neis))):
            s_x=neis[i]
            i_s=s_x.index    #index of ith nbs of v
            if(v_r_lst[i_s] >= 0 and v_r_lst[i_s] < 0.25):
                t_l_lst_nbd.append(i_s)
            elif(v_r_lst[i_s] >= 0.25 and v_r_lst[i_s] < 0.75):
                t_m_lst_nbd.append(i_s)
            elif(v_r_lst[i_s] >= 0.75 and v_r_lst[i_s] < 1):
                t_h_lst_nbd.append(i_s)
        f_l_lst_nbd.append(t_l_lst_nbd)
        f_m_lst_nbd.append(t_m_lst_nbd)
        f_h_lst_nbd.append(t_h_lst_nbd)
        f_s_lst_nbd.append(i_l_s_v)
        f_v_lst_nbd.append(i_l_v_v)
    f_s_lst_nbd = f_s_lst_nbd   #final sus nbs lst in G1
    f_v_lst_nbd = f_v_lst_nbd   #final vac nbs lst in G1
    f_t_lst_nbd = i_n_t_v   #final total nbs lst in G1
    f_l_lst_nbd = f_l_lst_nbd   #final low nbs lst in G2
    f_m_lst_nbd = f_m_lst_nbd   #final med nbs lst in G2
    f_h_lst_nbd = f_h_lst_nbd   #final high nbs lst in G2
    f_t_lst_nbd_g2 = i_n_t_v_g2 #final total nbs lst in G2
    for h in range(len(f_s_lst_nbd)):
        s_low = intersection(f_s_lst_nbd[h], f_l_lst_nbd[h])
        n_s_low = len(s_low)
        p_s_low = n_s_low/f_t_lst_nbd[h]
        p_s_low = round(p_s_low,4)
        f_prob_s_low.append(p_s_low)
        v_low = intersection(f_v_lst_nbd[h], f_l_lst_nbd[h])
        n_v_low = len(v_low)
        p_v_low = n_v_low/f_t_lst_nbd[h]
        p_v_low = round(p_v_low,4)
        f_prob_v_low.append(p_v_low)
        v_med = intersection(f_v_lst_nbd[h], f_m_lst_nbd[h])
        n_v_med = len(v_med)
        p_v_med = n_v_med/f_t_lst_nbd[h]
        p_v_med = round(p_v_med,4)
        f_prob_v_med.append(p_v_med)
        v_high = intersection(f_v_lst_nbd[h], f_h_lst_nbd[h])
        n_v_high = len(v_high)
        p_v_high = n_v_high/f_t_lst_nbd[h]
        p_v_high = round(p_v_high,4)
        f_prob_v_high.append(p_v_high)
        if(f_t_lst_nbd_g2[h] == 0):
            f_prob_low.append(0)
        else:
            p_low = len(f_l_lst_nbd[h])/f_t_lst_nbd_g2[h]
            p_low = round(p_low,4)
            f_prob_low.append(p_low)
    f_prob_s_low = f_prob_s_low
    f_prob_v_low = f_prob_v_low
    f_prob_v_med = f_prob_v_med
    f_prob_v_high = f_prob_v_high
    f_prob_low = f_prob_low
    return f_prob_s_low, f_prob_v_low, f_prob_v_med, f_prob_v_high,f_prob_low


# In[12]:


def epi_fermi_direc_function(G_1,i_v_c_lst,i_i_c_lst):
    epi_prob_vac_lst=[]
    epi_o=delta_pay_off(G_1,i_v_c_lst,i_i_c_lst)
    for t in range(len(epi_o)):
        if epi_o[t] < 0:
            epi_prob_vac_lst.append(0)
        elif epi_o[t] == 0:
            epi_prob_vac_lst.append(0.5)
        elif epi_o[t] > 0:
            j_i=0.1
            epi_f_f_d = math.exp((-1)*(j_i)*(epi_o[t]))
            epi_f_f_d = 1+epi_f_f_d
            epi_f_f_d = 1/epi_f_f_d
            epi_f_f_d = round(epi_f_f_d,4)
            epi_prob_vac_lst.append(epi_f_f_d)
    epi_prob_vac_lst = epi_prob_vac_lst
    return epi_prob_vac_lst


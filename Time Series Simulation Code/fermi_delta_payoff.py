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
import sys
import pandas as pd
import numpy as np
import csv
import time


# In[9]:


def delta_pay_off(G_1,v_r_lst,i_r_lst,jii):
    v_r_lst=v_r_lst
    i_r_lst=i_r_lst
    pay_off_lst=[]
    for v in G_1.vs:
        i=v.index
        t_nbs = len(G_1.neighbors(v.index))
        vac_nbs = len(list(G_1.vs[G_1.neighbors(v)](name='V'))) # total vaccinated nbs of v
        nvac_nbs = t_nbs-vac_nbs       # total non-vacc nbs of v
        c =1          
        h_p=math.exp((-1)*c*nvac_nbs)
        h_p=1-h_p
        theta_i=jii*h_p
        v_p_i=v_r_lst[i]
        i_p_i=i_r_lst[i]
        p_off_i= ((-1)*(v_p_i))+((i_p_i)*(theta_i))
        pay_off_lst.append(p_off_i)
    return pay_off_lst


# In[10]:


def fermi_direc_function(G_1,i_v_c_lst,i_i_c_lst,jii):
    prob_vac_lst=[]
    o=delta_pay_off(G_1,i_v_c_lst,i_i_c_lst,jii)
    for t in range(len(o)):
        if o[t] < 0:
            prob_vac_lst.append(o[t])
        elif o[t] == 0:
            prob_vac_lst.append(o[t])
        elif o[t] > 0:
            j_i=0.1
            f_f_d = math.exp((-1)*(j_i)*(o[t]))
            f_f_d = 1+f_f_d
            f_f_d = 1/f_f_d
            f_f_d = round(f_f_d,4)
            prob_vac_lst.append(f_f_d)
    prob_vac_lst = prob_vac_lst
    return prob_vac_lst


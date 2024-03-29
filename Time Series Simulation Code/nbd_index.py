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


# In[5]:


def nbd_index(G_2,in_v_r_lst):
    s=[]
    total_nbs_lst = []
    l_lst = []
    m_lst = []
    h_lst = []
    v_r_lst = in_v_r_lst
    f_l_lst_nbd = []  #final low nbs lst
    f_m_lst_nbd = []  #final medium nbs lst
    f_h_lst_nbd = []  #final high nbs lst
    for v in G_2.vs:
        l=[]
        n_l = 0
        n_m = 0
        n_h = 0
        t_l_lst_nbd = []  #temp low nbs lst
        t_m_lst_nbd = []  #temp medium nbs lst
        t_h_lst_nbd = []  #temp high nbs lst
        neis = G_2.vs[G_2.neighbors(v)]
        t_nbs = len(list(neis))
        total_nbs_lst.append(t_nbs)
        for i in range(len(list(neis))):
            s_x=neis[i]
            i_s=s_x.index
            l.append(i_s)
            if(v_r_lst[i_s] >= 0 and v_r_lst[i_s] < 0.25):
                n_l += 1
                t_l_lst_nbd.append(i_s)
            elif(v_r_lst[i_s] >= 0.25 and v_r_lst[i_s] < 0.75):
                n_m += 1
                t_m_lst_nbd.append(i_s)
            elif(v_r_lst[i_s] >= 0.75 and v_r_lst[i_s] <= 1):
                n_h += 1
                t_h_lst_nbd.append(i_s)
        h=l
        s.append(h)
        l_lst.append(n_l)
        m_lst.append(n_m)
        h_lst.append(n_h)
        f_l_lst_nbd.append(t_l_lst_nbd)
        f_m_lst_nbd.append(t_m_lst_nbd)
        f_h_lst_nbd.append(t_h_lst_nbd)
    return s,l_lst,m_lst,h_lst,total_nbs_lst,f_l_lst_nbd,f_m_lst_nbd,f_h_lst_nbd 


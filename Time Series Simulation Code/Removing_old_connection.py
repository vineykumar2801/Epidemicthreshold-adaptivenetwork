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


# In[11]:


def remov_connection(G_2,v_r_lst,n):
    v_r_lst=v_r_lst
    nbd_lst=nbd_index(G_2,v_r_lst)[0]
    #n_remove_conn = 0
    for i in range(n):
        v_nbs=nbd_lst[i]
        for j in range(len(v_nbs)):
            ed_we=1
            dell_c=v_r_lst[i]-v_r_lst[v_nbs[j]]
            dell_c=abs(dell_c)
            dell_c=ed_we*dell_c
            h_p=math.exp((-1)*(dell_c))
            delta=1
            prob_delta=delta*(1-h_p)
            m_r_n=random.uniform(0,10)
            m_r_n=m_r_n*1.5
            if(m_r_n<prob_delta):
                if(G_2.are_connected(i,v_nbs[j])==True):
                    G_2.delete_edges([(i,v_nbs[j])])
                    #n_remove_conn = n_remove_conn+1
                    continue
                else:
                    continue
            else:
                continue
    return G_2


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


# In[12]:


def new_connection(G_2,v_r_lst,n):
    v_r_lst=v_r_lst
    start=0
    end=n
    n_list = [item for item in range(start, end)]
    nbd_lst=nbd_index(G_2,v_r_lst)[0]
    #n_new_connection = 0
    for i in range(n):
        v_nbs=nbd_lst[i]
        non_nbs_list=list(set(n_list) - set(v_nbs))
        for j in range(len(non_nbs_list)):
            if(i != non_nbs_list[j]):
                dell_c=v_r_lst[i]-v_r_lst[non_nbs_list[j]]
                dell_c=abs(dell_c)
                eps=0.0000000001
                #eps = eps*eps
                #eps= sys.float_info.epsilon
                if(dell_c < eps):
                    G_2.add_edges([(i,non_nbs_list[j])])
                    #n_new_connection = n_new_connection+1
                    continue
                else:
                    continue
            else:
                continue
    return G_2


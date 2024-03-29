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


# In[3]:


def set_inf(G_1,lst_1,lst_2,lst_3,lst_4,lst_5):
    G_1.vs["name"]="S"
    G_1.vs[lst_1]["name"]= "I"     #infected nodes
    G_1.vs[lst_2]["name"]="R"
    G_1.vs[lst_3]["name"]="V"
    G_1.vs[lst_4]["name"]="i"
    G_1.vs[lst_5]["name"]="v"
    return G_1


# In[4]:


def visualige(G_1,i):
    color_dict={}
    color_dict["S"]='green'
    color_dict["I"]='red'
    color_dict["R"]='orange'
    color_dict["V"]='blue'
    G_1.vs['color']=[color_dict[c] for c in G_1.vs["name"]]
    out= ig.plot(G_1,vertex_label=[v.index for v in G_1.vs]) 
    out.save("LT_"+str(i)+".png")


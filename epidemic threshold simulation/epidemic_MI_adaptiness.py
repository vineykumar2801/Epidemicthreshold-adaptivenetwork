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


# In[6]:


def epi_nbd_index(G_1,G_2,in_v_r_lst):
    v_r_lst = in_v_r_lst
    f_prob_s_low = [] #final low sus probability lst G1
    f_prob_v_low = [] #final low vac probability lst G1
    f_prob_v_med = [] #final med vac probability lst G1
    f_prob_v_high = [] #final high vac probability lst G1
    f_prob_low = [] #final low probability lst G2
    f_prob_med = [] #final medium probability lst G2
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
        if(f_t_lst_nbd_g2[h]==0):
            f_prob_low.append(0)
            f_prob_med.append(0)
        else:
            p_low = len(f_l_lst_nbd[h])/f_t_lst_nbd_g2[h]
            p_low = round(p_low,4)
            f_prob_low.append(p_low)
            p_med = len(f_m_lst_nbd[h])/f_t_lst_nbd_g2[h]
            p_med = round(p_med,4)
            f_prob_med.append(p_med)
    f_prob_s_low = f_prob_s_low
    f_prob_v_low = f_prob_v_low
    f_prob_v_med = f_prob_v_med
    f_prob_v_high = f_prob_v_high
    f_prob_low = f_prob_low
    f_prob_med = f_prob_med
    return f_prob_s_low, f_prob_v_low, f_prob_v_med, f_prob_v_high,f_prob_low, f_prob_med


# In[7]:


data = pd.read_csv("I_rec.csv")
arr = data["Initial_recovery"].to_numpy()
rec_lst = arr.tolist()
data = pd.read_csv("I_inf.csv")
arr = data["Initial_inf"].to_numpy()
inf_lst = arr.tolist()
data = pd.read_csv("I_vac.csv")
arr = data["Initial_vacc"].to_numpy()
vac_lst = arr.tolist()
data = pd.read_csv("I_v_c.csv")
arr = data["I_v_c"].to_numpy()
i_v_c_lst = arr.tolist()
data = pd.read_csv("I_i_c.csv")
arr = data["I_i_c"].to_numpy()
i_i_c_lst = arr.tolist()


# In[8]:


def divide_lists(list1, list2, list3, list4):
    result1 = []
    result2 = []
    result3 = []
    a = 0 
    for i in range(len(list1)):
        if list4[i] == 0:
            result1.append(a)
            result2.append(a)
            result3.append(a)
        else:
            result1.append(list1[i] / list4[i])
            result2.append(list2[i] / list4[i])
            result3.append(list3[i] / list4[i])
    return result1,result2,result3


# In[9]:


def av_vac_risk_list(G_2,n,in_v_r_lst):       # n=no of nodes in the graph
    q_1_v=in_v_r_lst
    new_r_v_lst=[]
    op_change = nbd_index(G_2,in_v_r_lst)
    index_op_change = op_change[0]      # index of the neighbourhood of all the nodes 
    l_lst_op_change = op_change[1]      # number of low neighbourhood of all the nodes  
    m_lst_op_change = op_change[2]    # number of medium neighbourhood of all the nodes    
    h_lst_op_change = op_change[3]    # number of high neighbourhood of all the nodes
    total_lst_op_change = op_change[4]   # total number of neighbourhood of all the nodes
    f_l_lst_nbd = op_change[5]
    f_m_lst_nbd = op_change[6]
    f_h_lst_nbd = op_change[7]
    prob_lst = divide_lists(l_lst_op_change,m_lst_op_change,h_lst_op_change,total_lst_op_change)  #prob of low,medium,high
    prob_l_lst = prob_lst[0]
    prob_m_lst = prob_lst[1]
    prob_h_lst = prob_lst[2]
    for i in range(n):
        t_l_v_risk = 0
        t_m_v_risk = 0
        t_h_v_risk = 0
        if(total_lst_op_change[i] == 0):
            r_avg=q_1_v[i]
            new_r_v_lst.append(r_avg)
            continue
        elif prob_l_lst[i] >= prob_m_lst[i] and prob_l_lst[i] >= prob_h_lst[i]:
            i_th_node_l_nbs = f_l_lst_nbd[i]
            to_l_nbs = len(i_th_node_l_nbs)
            for l_nbs in range(len(i_th_node_l_nbs)):
                l_nbs_index = i_th_node_l_nbs[l_nbs]
                l_v_risk = q_1_v[l_nbs_index]
                t_l_v_risk = t_l_v_risk + l_v_risk
            new_av = t_l_v_risk/to_l_nbs
            new_r_v_lst.append(new_av)
        elif prob_m_lst[i] >= prob_l_lst[i] and prob_m_lst[i] >= prob_h_lst[i]:
            i_th_node_m_nbs = f_m_lst_nbd[i]
            to_m_nbs = len(i_th_node_m_nbs)
            for m_nbs in range(len(i_th_node_m_nbs)):
                m_nbs_index = i_th_node_m_nbs[m_nbs]
                m_v_risk = q_1_v[m_nbs_index]
                t_m_v_risk = t_m_v_risk + m_v_risk
            new_av = t_m_v_risk/to_m_nbs
            new_r_v_lst.append(new_av)
        elif prob_h_lst[i] >= prob_l_lst[i] and prob_h_lst[i] >= prob_m_lst[i]:
            i_th_node_h_nbs = f_h_lst_nbd[i]
            to_h_nbs = len(i_th_node_h_nbs)
            for h_nbs in range(len(i_th_node_h_nbs)):
                h_nbs_index = i_th_node_h_nbs[h_nbs]
                h_v_risk = q_1_v[h_nbs_index]
                t_h_v_risk = t_h_v_risk + h_v_risk
            new_av = t_h_v_risk/to_h_nbs
            new_r_v_lst.append(new_av)
    updated_av_r_v_lst=new_r_v_lst
    return updated_av_r_v_lst


# In[10]:


def delta_pay_off(G_1,v_r_lst,i_r_lst):
    v_r_lst=v_r_lst
    i_r_lst=i_r_lst
    pay_off_lst=[]
    for v in G_1.vs:
        i=v.index
        t_nbs = len(G_1.neighbors(v.index))
        vac_nbs = len(list(G_1.vs[G_1.neighbors(v)](name='V'))) # total vaccinated nbs of v
        nvac_nbs = t_nbs-vac_nbs       # total non-vacc nbs of v
        jii=1
        c =1          
        h_p=math.exp((-1)*c*nvac_nbs)
        h_p=1-h_p
        theta_i=jii*h_p
        v_p_i=v_r_lst[i]
        i_p_i=i_r_lst[i]
        p_off_i= ((-1)*(v_p_i))+((i_p_i)*(theta_i))
        pay_off_lst.append(p_off_i)
    return pay_off_lst


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


# In[13]:


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


# In[14]:


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


# In[16]:


def Cumulative(lists):
    cu_list = []
    length = len(lists)
    cu_list = [sum(lists[0:x:1]) for x in range(0, length+1)]
    return cu_list[1:]


# In[17]:


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


# In[18]:


def LT(G_1,G_2,n,inf_lst,rec_lst,vac_lst,i_v_c_lst,i_i_c_lst,i):
        it=0
        days_lst_v = [0.16,0.14,0.2,0.125,0.111]
        days_lst_i = [0.14,0.125,0.111,0.10,0.09,0.08,0.077]
        gamma_lst=[0.14,0.1,0.07,0.067,0.055,0.05]
        lst_it=[]
        lst_inf=[]
        lst_new_inf=[]
        lst_sus=[]
        lst_v=[]
        lst_new_v=[]
        lst_r=[]
        beta= 0.6
        #gama= 0.14
        n_spe_inf = 20
        #a=initial_risk_list(n)
        r_v_lst=i_v_c_lst
        r_i_lst=i_i_c_lst
        main_temp_inf = []
        main_temp_vac = []
        #main_risk_list_v=[]
        #main_risk_list_i=[]
        #while(n_spe_inf!=0):
        new_conne_lst=[]
        remo_conn_lst = []
        n_inf_per_day = []
        n_vac_per_day = []
        for m in range(70):
            t_inf_temp = []
            t_vac_temp = []
            n_inf_p = 0
            n_vac_p = 0
            #print("------------------------------------------------------")
            it=it+1
            lst_it.append(it)
            n_s=len(list(G_1.vs(name='S')))
            lst_sus.append(n_s)
            n_spe_inf=len(list(G_1.vs(name='i')))
            lst_new_inf.append(n_spe_inf)
            n_i=len(list(G_1.vs(name='I')))
            lst_inf.append(n_i)
            #print(n_i)
            n_v=len(list(G_1.vs(name='V')))
            lst_v.append(n_v)
            n_spe_vac=len(list(G_1.vs(name='v')))
            lst_new_v.append(n_spe_vac)
            n_r=len(list(G_1.vs(name='R')))
            lst_r.append(n_r)
            #p_off_lst=delta_pay_off(G_1,r_v_lst,r_i_lst)           #delta p list (for vaccination)
            pay_off_lst_fermi_dirac = fermi_direc_function(G_1,r_v_lst,r_i_lst)  # fermi_dirac on positive number
            already_adopting_s=G_1.vs(name='S')
            already_adopting_i=G_1.vs(name='I')
            already_adopting_new_i = G_1.vs(name='i')
            already_adopting_r=G_1.vs(name='R')
            already_adopting_v=G_1.vs(name='V')
            already_adopting_new_v = G_1.vs(name='v')
            #main_risk_list_v.append(r_v_lst)
            for v in G_1.vs:
                if(v in already_adopting_v):
                    G_1.vs[v.index]["name"] = "V"
                    continue
                if(v in already_adopting_new_v):
                    G_1.vs[v.index]["name"] = "v"
                    continue
                if(v in already_adopting_i):
                    G_1.vs[v.index]["name"] = "I"
                    continue
                if(v in already_adopting_r):
                    G_1.vs[v.index]["name"] = "R"  
                    continue
                if(v in already_adopting_s):
                    s_1=len(list(G_1.vs[G_1.neighbors(v)](name='I')))
                    m=(1-beta)
                    m=m**s_1
                    prob_in=1-m
                    r_1=random.uniform(0,1)
                    if(r_1<prob_in):
                        G_1.vs[v.index]["name"] = "I"
                        t_inf_temp.append(v.index)
                        continue
                    else:
                        G_1.vs[v.index]["name"] = "S"
                        continue
                if(v in already_adopting_new_i):
                    prob_rec=random.uniform(0,1)
                    gama = random.choice(gamma_lst)
                    if(prob_rec<gama):
                        G_1.vs[v.index]["name"] = "R"
                        continue
                    else:
                        G_1.vs[v.index]["name"] = "i"
                        continue 
            new_already_adopting_s=G_1.vs(name='S')
            new_already_adopting_i=G_1.vs(name='I')
            new_already_adopting_r=G_1.vs(name='R')
            new_already_adopting_v=G_1.vs(name='V')
            already_adopting_new_v = G_1.vs(name='v')
            already_adopting_new_i = G_1.vs(name='i')
            for v in G_1.vs:
                if(v in new_already_adopting_v):
                    G_1.vs[v.index]["name"] = "V"
                    continue
                if(v in new_already_adopting_r):
                    G_1.vs[v.index]["name"] = "R" 
                    continue
                if(v in already_adopting_new_v):
                    G_1.vs[v.index]["name"] = "v"
                    continue
                if(v in already_adopting_new_i):
                    G_1.vs[v.index]["name"] = "i"
                    continue
                if(v in new_already_adopting_i):
                    G_1.vs[v.index]["name"] = "I"  
                    continue
                if(v in new_already_adopting_s):
                    h_i = v.index
                    if(pay_off_lst_fermi_dirac[h_i] > 0):
                        v_r_n=random.uniform(0,1)
                        v_r_n=v_r_n
                        #print(v_r_n)
                        if(v_r_n < pay_off_lst_fermi_dirac[h_i]):
                            #print(phi_p_off)
                            G_1.vs[v.index]["name"] = "V"
                            t_vac_temp.append(v.index)
                            continue
                        else:
                            G_1.vs[v.index]["name"] = "S"
                            continue
                    elif(pay_off_lst_fermi_dirac[h_i] < 0):
                        G_1.vs[v.index]["name"] = "S"
                        continue
                    elif(pay_off_lst_fermi_dirac[h_i] == 0):
                        v_z_r_n=random.uniform(0,1)
                        if(v_z_r_n>0.5):
                            G_1.vs[v.index]["name"] = "V"
                            t_vac_temp.append(v.index)
                            continue
                        else:
                            G_1.vs[v.index]["name"] = "S"
                            continue
            new_already_adopting_s=G_1.vs(name='S')
            new_already_adopting_i=G_1.vs(name='I')
            new_already_adopting_r=G_1.vs(name='R')
            new_already_adopting_v=G_1.vs(name='V')
            already_adopting_new_v = G_1.vs(name='v')
            already_adopting_new_i = G_1.vs(name='i')
            for v in G_1.vs:
                if(v in already_adopting_new_v):
                    G_1.vs[v.index]["name"] = "v"
                    continue
                if(v in new_already_adopting_s):
                    G_1.vs[v.index]["name"] = "S"
                    continue
                if(v in new_already_adopting_r):
                    G_1.vs[v.index]["name"] = "R" 
                    continue
                if(v in already_adopting_new_i):
                    G_1.vs[v.index]["name"] = "i"  
                    continue
                if(v in new_already_adopting_i):
                    i_r = random.uniform(0,1)
                    i_r = i_r
                    i_days = random.choice(days_lst_i)
                    if(i_r < i_days):
                        G_1.vs[v.index]["name"] = "i"
                        n_inf_p = n_inf_p +1
                        continue
                    else:
                        G_1.vs[v.index]["name"] = "I"  
                        continue
                if(v in new_already_adopting_v):
                    i_v = random.uniform(0,1)
                    i_days = random.choice(days_lst_v)
                    if(i_v < i_days):
                        G_1.vs[v.index]["name"] = "v"
                        n_vac_p = n_vac_p+1
                        continue
                    else:
                        G_1.vs[v.index]["name"] = "V"  
                        continue
            edg_1 = G_2.ecount()
            #print(edg_1)
            G_2=remov_connection(G_2,r_v_lst,n)
            edg_2 = G_2.ecount()
            remo_1 = edg_1 - edg_2
            remo_conn_lst.append(remo_1)
            G_2 = new_connection(G_2,r_v_lst,n)
            edg_3 = G_2.ecount()
            #print(edg_3)
            new_c_1 = edg_3 - edg_2
            new_conne_lst.append(new_c_1)
            #print(len(r_v_lst))
            av_v_l=av_vac_risk_list(G_2,n,r_v_lst)            #update risk with weighted average risk
            r_v_lst=av_v_l
            main_temp_inf.append(t_inf_temp)
            main_temp_vac.append(t_vac_temp)
            n_inf_per_day.append(n_inf_p)
            n_vac_per_day.append(n_vac_p)
        i_v_c_lst = r_v_lst
        i_i_c_lst = r_i_lst
        G_1 = G_1
        G_2 = G_2
        list_it=lst_it
        #main_risk_list_v=main_risk_list_v
        list_sus=lst_sus
        #print(list_sus)
        #list_inf=lst_inf
        list_new_inf = lst_new_inf
        list_n_inf_per_day = n_inf_per_day
        #print(list_n_inf_per_day)
        #cu_i = Cumulative(list_new_inf)
        #list_v=lst_v
        list_new_vac = lst_new_v
        list_n_vac_per_day  = n_vac_per_day
        #print(list_n_vac_per_day)
        #cu_v = Cumulative(list_new_vac)
        list_r=lst_r
        #cu_r = Cumulative(list_r)
        list_new_connection = new_conne_lst
        final_phi_lst = epi_fermi_direc_function(G_1,i_v_c_lst,i_i_c_lst)
        final_all_lst = epi_nbd_index(G_1,G_2,i_v_c_lst)
        f_prob_s_low = final_all_lst[0]
        f_prob_v_low = final_all_lst[1]
        f_prob_v_med = final_all_lst[2]
        f_prob_v_high = final_all_lst[3]
        f_prob_low = final_all_lst[4]
        f_prob_med = final_all_lst[5]
        list_remove_connection = remo_conn_lst
        kii=[]
        for i in range(len(f_prob_s_low)):
            e_1 = (1-final_phi_lst[i])
            e_2_1 = (1-f_prob_low[i])
            e_2_2 = e_2_1*e_1
            e_2 = e_2_2*f_prob_s_low[i]
            e_3_1 = 1-f_prob_v_low[i]-f_prob_v_med[i]-f_prob_v_high[i]
            e_3 = e_2 + e_3_1
            e_final = e_1 * f_prob_med[i] * e_3
            ki=e_final
            kii.append(ki)
        kii=kii
        f_cji=[]
        A_G_1=G_1.get_adjacency()
        for i in range(n):
            f_rji=[]
            for j in range(n):
                u=A_G_1[i][j]
                if(u==1):
                    f_rji.append(kii[i]*1)
                else:
                    f_rji.append(kii[i]*0)
            f_cji.append(f_rji)
        f_cji=f_cji
        f_cji=np.array(f_cji)
        w,v=eig(f_cji)
        w=w.real
        w.sort()
        lambda_H=w[-1]
        e_thresold=1/(lambda_H)
        return e_thresold


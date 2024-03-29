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


# In[8]:


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


# In[13]:


def random_op_change(r_v_lst,m):
    my_list = r_v_lst
    elements_to_remove = random.sample(my_list, m)
    for element in elements_to_remove:
        my_list.remove(element)
    my_list = my_list
    n_c_v_lst=[]
    for x in range(m):
        n_c=random.uniform(0,1)
        n_c_v_lst.append(n_c)
    n_c_v_lst = n_c_v_lst
    f_c_v = my_list + n_c_v_lst
    random.shuffle(f_c_v)
    return f_c_v


# In[14]:


def Cumulative(lists):
    cu_list = []
    length = len(lists)
    cu_list = [sum(lists[0:x:1]) for x in range(0, length+1)]
    return cu_list[1:]


# In[15]:


def LT(G_1,G_2,n,inf_lst,rec_lst,vac_lst,i_v_c_lst,i_i_c_lst,i,jii):
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
            pay_off_lst_fermi_dirac = fermi_direc_function(G_1,r_v_lst,r_i_lst,jii)  # fermi_dirac on positive number
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
                    r_1=r_1*1.5
                    if(r_1<prob_in):
                        G_1.vs[v.index]["name"] = "I"
                        t_inf_temp.append(v.index)
                        continue
                    else:
                        G_1.vs[v.index]["name"] = "S"
                        continue
                if(v in already_adopting_new_i):
                    prob_rec=random.uniform(0,1.5)
                    prob_rec=prob_rec*1.8
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
                        v_r_n=random.uniform(0.2,1)
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
            op_change_r_lst = random_op_change(r_v_lst,3500)
            r_v_lst = op_change_r_lst
            main_temp_inf.append(t_inf_temp)
            main_temp_vac.append(t_vac_temp)
            n_inf_per_day.append(n_inf_p)
            n_vac_per_day.append(n_vac_p)
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
        list_remove_connection = remo_conn_lst
        FIL = {'iteration' : list_it, 'new_conn': list_new_connection, 'remo_conn': list_remove_connection,
             'recover':list_r,'infection': list_new_inf, 'new_inf_per_day': list_n_inf_per_day, 
               'vaccination': list_new_vac, 'new_vac_per_day': list_n_vac_per_day, 'susceptible': list_sus}
        df = pd.DataFrame(FIL)
        df.to_excel(writer, sheet_name="s_"+str(i)+"_t", index=False)


# In[16]:


st = time.time()


# In[17]:


writer = pd.ExcelWriter('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', engine='xlsxwriter')
for i in range(0,5):
    results = []
    with open("G_uw_b_csv_"+str(i)+".csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)
    A_S=results
    G_1 = ig.Graph.Adjacency(A_S,"undirected")
    results_1 = []
    with open("G_2_csv_"+str(i)+".csv") as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
        for row in reader: # each row is a list
            results_1.append(row)
    A_S_1=results_1
    G_2 = ig.Graph.Adjacency(A_S_1,"undirected")
    n_G_1=set_inf(G_1,inf_lst,rec_lst,vac_lst,[10],[500])
    G_1=n_G_1
    m=LT(G_1,G_2,5000,inf_lst,rec_lst,vac_lst,i_v_c_lst,i_i_c_lst,i,z_n)
writer.save()    


# In[18]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'vaccination' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofvac=[i['vaccination'].tolist() for i in dictoframe[0] if('vaccination' in i)]
finaldic=dict(zip(sheets,listofvac))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_vac_"+str(f_o_n)+".csv")


# In[19]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'infection' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofinf=[i['infection'].tolist() for i in dictoframe[0] if('infection' in i)]
finaldic=dict(zip(sheets,listofinf))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_inf_"+str(f_o_n)+".csv")


# In[20]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'recover' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofinf=[i['recover'].tolist() for i in dictoframe[0] if('recover' in i)]
finaldic=dict(zip(sheets,listofinf))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_rec_"+str(f_o_n)+".csv")


# In[21]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'new_inf_per_day' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofcuinf=[i['new_inf_per_day'].tolist() for i in dictoframe[0] if('new_inf_per_day' in i)]
finaldic=dict(zip(sheets,listofcuinf))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_new_inf_per_day_"+str(f_o_n)+".csv")


# In[22]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'new_vac_per_day' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofcuvac=[i['new_vac_per_day'].tolist() for i in dictoframe[0] if('new_vac_per_day' in i)]
finaldic=dict(zip(sheets,listofcuvac))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_new_vac_per_day_"+str(f_o_n)+".csv")


# In[23]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'susceptible' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofcuinf=[i['susceptible'].tolist() for i in dictoframe[0] if('susceptible' in i)]
finaldic=dict(zip(sheets,listofcuinf))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_susceptible_"+str(f_o_n)+".csv")


# In[24]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'new_conn' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofcuinf=[i['new_conn'].tolist() for i in dictoframe[0] if('new_conn' in i)]
finaldic=dict(zip(sheets,listofcuinf))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_new_conn_"+str(f_o_n)+".csv")


# In[25]:


datadic = pd.read_excel('z_'+str(z_n)+'_'+str(f_o_n)+'.xlsx', sheet_name=None)
sheets=datadic.keys()
sheets=[i+"_"+'remo_conn' for i in sheets]
dictoframe=pd.Series(datadic).to_frame()
listofcuinf=[i['remo_conn'].tolist() for i in dictoframe[0] if('remo_conn' in i)]
finaldic=dict(zip(sheets,listofcuinf))
df=pd.DataFrame.from_dict(finaldic,orient='index').transpose()
df.to_csv("z_"+str(z_n)+"_remo_conn_"+str(f_o_n)+".csv")


# In[26]:


et = time.time()
s_t= et - st
m_t= s_t/60
print('Execution time:', m_t, 'minutes')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





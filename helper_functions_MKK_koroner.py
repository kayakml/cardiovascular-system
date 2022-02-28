#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:39:50 2021

@author: kemal
"""
import numpy as np

def creat_param_list(source_list):
    return source_list
    
    
def creat_params_of_journal():
    params_of_journal = 48*[0]
    Ts1 = 0.3
    params_of_journal[0] = Ts1 
    Ts2 = 0.45
    params_of_journal[1] = Ts2
    Tpwb = 0.92
    params_of_journal[2] = Tpwb
    Tpww = 0.09
    params_of_journal[3] = Tpww
    T = 1.0
    params_of_journal[4] = T
    Elamin = 0.15                             #++
    params_of_journal[5] = Elamin
    Elamax = 0.25                             #++
    params_of_journal[6] = Elamax
    Elvs = 2.5                             #++
    params_of_journal[7] = Elvs
    Elvd = 0.1                             #++
    params_of_journal[8] = Elvd
    Eramin = 0.15                             #++
    params_of_journal[9] = Eramin
    Eramax = 0.25                             #++
    params_of_journal[10] = Eramax
    Ervs = 1.15                             #++
    params_of_journal[11] = Ervs
    Ervd = 0.1                             #++
    params_of_journal[12] = Ervd
    Pla0 = 1.0                             #++
    params_of_journal[13] = Pla0
    Plv0 = 1.0                             #++
    params_of_journal[14] = Plv0
    Vla0 = 4.0                             #++
    params_of_journal[15] = Vla0
    Vlv0 = 5.0                             #5.0?
    params_of_journal[16] = Vlv0
    Pra0 = 1.0                             #++
    params_of_journal[17] = Pra0
    Prv0 = 1.0                             #++
    params_of_journal[18] = Prv0
    Vra0 = 4.0                             #++
    params_of_journal[19] = Vra0
    Vrv0 = 10.0                             #10.0?
    params_of_journal[20] = Vrv0
    CQmi = 400.                             #++
    params_of_journal[21] = CQmi
    CQao = 350.                             #++
    params_of_journal[22] = CQao
    CQti = 400.                             #++
    params_of_journal[23] = CQti
    CQpo = 350.                             #++
    params_of_journal[24] = CQpo
    Rpvn = 0.006                             #++
    params_of_journal[25] = Rpvn
    Rsvn = 0.075                             #++
    params_of_journal[26] = Rsvn
    Rsar = 0.5                             #++
    params_of_journal[27] = Rsar
    Rscp = 0.52                             #++
    params_of_journal[28] = Rscp
    Cpas = 0.18                             #++
    params_of_journal[29] = Cpas
    Lpas = 0.000052                             #++
    params_of_journal[30] = Lpas
    Rpas = 0.002                             #++
    params_of_journal[31] = Rpas
    Cpat = 3.8                             #++
    params_of_journal[32] = Cpat
    Lpat = 0.0017                             #++
    params_of_journal[33] = Lpat
    Rpat = 0.01                             #++
    params_of_journal[34] = Rpat
    Rpar = 0.05                             #++
    params_of_journal[35] = Rpar
    Rpcp = 0.25                             #++
    params_of_journal[36] = Rpcp
    Cpvn = 20.5                             #++
    params_of_journal[37] = Cpvn
    Csvn = 20.5                             #++
    params_of_journal[38] = Csvn
    Csas = 0.08                             #++
    params_of_journal[39] = Csas
    Rsas = 0.003                             #++
    params_of_journal[40] = Rsas
    Lsas = 0.000062                             #++
    params_of_journal[41] = Lsas
    Csat = 1.6                             #++
    params_of_journal[42] = Csat
    Lsat = 0.0017                             #++
    params_of_journal[43] = Lsat
    Rsat = 0.05                             #++
    params_of_journal[44] = Rsat
    R1 = 1
    params_of_journal[45] = R1
    R2 = 1
    params_of_journal[46] = R2
    C1 = 1
    params_of_journal[47] =  C1
    return params_of_journal

def vent_act (t,param):
    tt =t%param[4]
    if ((0<=tt)*(tt<=param[0])):
        vent_act = 1.0-np.cos((tt/param[0])*np.pi)
    elif ((param[0]<tt)*(tt<=param[1])):
        vent_act = 1.0+np.cos(((tt-param[0])/(param[1]-param[0]))*np.pi)
    else:
        vent_act = 0
    return (vent_act)


def atr_act (t,param):
    tt =t%param[4]
    if (0<=tt)*(tt<param[2]):
        atr_act = 0.0
    elif (param[2]<=tt)*(tt<(param[2]+param[3])):
        atr_act = 1.0-np.cos(((tt-param[2])/(param[3]))*2*np.pi)
    else:
        atr_act = 0.0
    return (atr_act)

def diff_atr_act(t,param):
    tt = t%param[4]
    if (0<=tt)*(tt<param[2]):
        diff_atr_act = 0.0
    elif (param[2]<=tt)*(tt<(param[2]+param[3])):
        diff_atr_act = (2*np.pi/param[3]) *np.sin( ( (tt*2*np.pi/param[3]) - (2*np.pi*param[2]/param[3]) ) )
    else:
        diff_atr_act = 0.0
    return (diff_atr_act)

def elv(t,param):
    t = t%param[4]
    return param[8] + ((param[7]-param[8])/2.0)*vent_act(t,param)

def ela(t,param):
    t = t%param[4]
    return param[5] + ((param[6]-param[5])/2.0)*atr_act(t,param)

def erv(t,param):
    t = t%param[4]
    return param[12] + ((param[11]-param[12])/2.0)*vent_act(t,param)

def era(t,param):
    t = t%param[4]
    return param[9] + ((param[10]-param[9])/2.0)*atr_act(t,param)

#######################################################################################################
def solve_mycvs_alg(t,solf,param):
    Plv_t = []
    Prv_t = []
    Pla_t = []
    Pra_t = []
    Qsvn_t = []
    Qpvn_t = []
    Qao_t = []
    Qmi_t = []
    Qpo_t = []
    Qti_t = []
    Pc_t = []
    Qc_t = []
    Vlv_t = solf.y[0]
    Vrv_t = solf.y[1]
    Vla_t = solf.y[2]
    Vra_t = solf.y[3]
    Psas_t = solf.y[4]
    Qsas_t = solf.y[5]
    Psat_t = solf.y[6]
    Qsat_t = solf.y[7]
    Psvn_t = solf.y[8]
    Ppas_t = solf.y[9]
    Qpas_t = solf.y[10]
    Ppat_t = solf.y[11]
    Qpat_t = solf.y[12]
    Ppvn_t = solf.y[13]
    Pc_t = solf.y[14]
    for i in range(len(t)):
        Plv_t.append(param[14] + elv(t[i],param)*(Vlv_t[i] - param[16]))
        Prv_t.append(param[18] + erv(t[i],param)*(Vrv_t[i] - param[20]))
        Pla_t.append(param[13] + ela(t[i],param)*(Vla_t[i] - param[15]))
        Pra_t.append(param[17] + era(t[i],param)*(Vra_t[i] - param[19]))
        Qsvn_t.append((1/param[26])*(Psvn_t[i] - Pra_t[i]))
        Qpvn_t.append((1/param[25])*(Ppvn_t[i] - Pla_t[i]))
        Qao_t.append((param[22]*(np.sqrt((np.greater_equal(Plv_t[i],Psas_t[i]))*(Plv_t[i]-Psas_t[i])))))
        Qmi_t.append((param[21]*(np.sqrt((np.greater_equal(Pla_t[i],Plv_t[i]))*(Pla_t[i]-Plv_t[i])))))
        Qpo_t.append((param[24]*(np.sqrt((np.greater_equal(Prv_t[i],Ppas_t[i]))*(Prv_t[i]-Ppas_t[i])))))
        Qti_t.append((param[23]*(np.sqrt((np.greater_equal(Pra_t[i],Prv_t[i]))*(Pra_t[i]-Prv_t[i])))))  
        Qc_t.append((Psas_t[i]-Pc_t[i])/param[45])
    Plv_t = np.array(Plv_t)
    Prv_t = np.array(Prv_t)
    Pla_t = np.array(Pla_t)
    Pra_t = np.array(Pra_t)
    Qsvn_t = np.array(Qsvn_t)
    Qpvn_t = np.array(Qpvn_t)    
    Qao_t = np.array(Qao_t)
    Qmi_t = np.array(Qmi_t)
    Qpo_t = np.array(Qpo_t)
    Qti_t = np.array(Qti_t)
    Pc_t = np.array(Pc_t)
    Qc_t = np.array(Qc_t)
    

    soltn_alg = [Vlv_t,Plv_t,Vrv_t,Prv_t,Vla_t,Pla_t,Vra_t,Pra_t,Psas_t,Qsas_t,Psat_t,Qsat_t,Psvn_t,Qsvn_t,\
                 Ppas_t,Qpas_t,Ppat_t,Qpat_t,Ppvn_t,Qpvn_t,Qao_t,Qmi_t,Qpo_t,Qti_t,Pc_t,Qc_t]
    return soltn_alg



"""
Equations
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import helper_functions_MKK as fu


# INITIAL VALUES
#     [Vlv_0,Vrv_0,Vla_0,Vra_0,Psas_0,Qsas_0,Psat_0,Qsat_0,Psvn_0,Ppas_0,Qpas_0,Ppat_0,Qpat_0,Ppvn_0]            
U_0 = [500,  500,  20,   20,   100,    0,    100,    0,     0,     30,     0,     30,    0,     0]

t_span = (0,10)
t_eval=np.linspace(0,10,1000)


def solve_cvs(param):

    def cvs(t,u):
    
        Vlv,Vrv,Vla,Vra,Psas,Qsas,Psat,Qsat,Psvn,Ppas,Qpas,Ppat,Qpat,Ppvn = u
    
        def pla(t):
            return (param[13] + fu.ela(t,param)*(Vla-param[15]))
        def plv(t):
            return (param[14] + fu.elv(t,param)*(Vlv-param[16]))
        def pra(t):
            return (param[17] + fu.era(t,param)*(Vra-param[19]))
        def prv(t):
            return (param[18] + fu.erv(t,param)*(Vrv-param[20]))
        def Qsvn(t):
            return (1/param[26])*(Psvn-pra(t))
        def Qpvn(t):
            return (1/param[25])*(Ppvn-pla(t))
        def Qmi(t):
            if (pla(t)>=plv(t)):
                return (param[21]*np.sqrt(pla(t)-plv(t)))
            else:
                return (0.0)
        def Qao(t):
            if (plv(t)>=Psas):
                return (param[22]*np.sqrt(plv(t)-Psas))
            else:
                return(0.0)
        def Qti(t):
            if (pra(t)>=prv(t)):
                return (param[23]*np.sqrt(pra(t)-prv(t)))
            else:
                return (0.0)
        def Qpo(t):
            if (prv(t)>=Ppas):
                return (param[24]*np.sqrt(prv(t)-Ppas))
            else:
                return(0.0)  
        dVlv = Qmi(t) - Qao(t) #dVlv u1
        dVrv = Qti(t) - Qpo(t) #dVrv u2
        dVla = Qpvn(t) - Qmi(t) #dVla u3
        dVra = Qsvn(t) - Qti(t) #dVra u4 
        dPsas = (1/param[39])*(Qao(t)-Qsas) #dPsas u5
        dQsas = (1/param[41])*(Psas-Psat-param[40]*Qsas) #dQsas u6 
        dPsat = (1/param[42])*(Qsas-Qsat) #dPsat u7
        dQsat = (1/param[43])*(Psat-Psvn-(param[44]+param[27]+param[28])*Qsat) #dQsat u8
        dPsvn = (1/param[38])*(Qsat-Qsvn(t)) #dPsvn u9
        dPpas = (1/param[29])*(Qpo(t)-Qpas) #dPpas u10
        dQpas = (1/param[30])*(Ppas-Ppat-param[31]*Qpas) #dQpas u11
        dPpat = (1/param[32])*(Qpas-Qpat) #dPpat u12
        dQpat = (1/param[33])*(Ppat-Ppvn-(param[34]+param[35]+param[36])*Qpat) #dQpat u13
        dPpvn = (1/param[37])*(Qpat-Qpvn(t))  #dPpvn u14
    
        return [dVlv,dVrv,dVla,dVra,dPsas,dQsas,dPsat,dQsat,dPsvn,dPpas,dQpas,dPpat,dQpat,dPpvn]
        
    return solve_ivp(cvs, t_span, y0=U_0, method = 'LSODA', t_eval=t_eval)
    

param = fu.creat_params_of_journal()
solf = solve_cvs(param)

[Vlv_t,Plv_t,Vrv_t,Prv_t,Vla_t,Pla_t,Vra_t,Pra_t,Psas_t,Qsas_t,Psat_t,Qsat_t,Psvn_t,Qsvn_t,\
                 Ppas_t,Qpas_t,Ppat_t,Qpat_t,Ppvn_t,Qpvn_t,Qao_t,Qmi_t,Qpo_t,Qti_t] = fu.solve_mycvs_alg(t_eval,solf,param)

#!/usr/bin/env python3

import numpy as np
import sys

from simulators import *


_,rule,proposal,m_lg_dt,M_lg_dt=sys.argv

simulator= simMALA if proposal=="EM" else simMALA_HMC
def sin_potential(q):
    return np.sin(2*np.pi*q)

def grad_sin_potential(q):
    return 2*np.pi*np.cos(2*np.pi*q)

lg_dts=np.linspace(int(m_lg_dt),int(M_lg_dt),20)

T_corr=0.5
N_iter=1000000

M = 64

for i in range(N_iter):

    q=np.zeros(M)

    #equilibriate
    simMALA(q,sin_potential,grad_sin_potential,1e-3,1.0,10000)

    for lg_dt in lg_dts:
        dt=10**lg_dt
        N_corr=int(np.floor(T_corr/dt))
        sd_coords=np.zeros_like(q)
        grad_V_hist=np.zeros((N_corr,M))
        n_steps=10*N_corr
        C_sum=np.zeros((N_corr,M))
        simulator(q,sin_potential,grad_sin_potential,dt,1.0,n_steps,rule=rule,sd_coords=sd_coords,hist=grad_V_hist,C_sum=C_sum)

        C_hat=np.mean(C_sum/n_steps,axis=1)
        msd=np.mean(sd_coords**2)
        f=open(f"{rule}_{proposal}_{lg_dt}.out","a")
        print(C_hat[0],dt*sum(C_hat),msd/(2*n_steps*dt),file=f)
        f.close()


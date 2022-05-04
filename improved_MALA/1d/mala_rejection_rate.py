import matplotlib.pyplot as plt
import numpy as np


def potential(q):
    return np.sin(2*np.pi*q)


def grad_potential(q):
    return 2*np.pi*np.cos(2*np.pi*q)


def simMALA(q, potential, grad_potential, dt, beta, nsteps, histogram=None):
    M = q.shape[0]
    V = potential(q)
    grad_V = grad_potential(q)
    sigma = np.sqrt(2*dt)

    n_accepted=0
    for i in range(nsteps):
        G = np.random.standard_normal(M)
        qtilde = q-beta*dt*grad_V+sigma*G
        displacement = qtilde-q
        qtilde=qtilde % 1
        Vtilde = potential(qtilde)
        grad_Vtilde = grad_potential(qtilde)
        alpha=beta*(Vtilde-V)+((beta*dt*grad_Vtilde-displacement)**2)/(4*dt)-((displacement+beta*dt*grad_V)**2)/(4*dt)
        U=np.random.uniform(size=M)
        accepted=np.log(U)<-alpha

        q[accepted]=qtilde[accepted]
        V[accepted]=Vtilde[accepted]
        grad_V[accepted]=grad_Vtilde[accepted]
        n_accepted+=np.sum(accepted)
    
    return (n_accepted,nsteps*M)


M = 1000
N_bins = 200
log_dts = np.linspace(-4, -2, 40)
dts = 10 ** log_dts
hist = np.zeros_like(dts)
ref_dts = 0

for dt in dts:
    q=np.zeros(M)
    (nacc,ntot)=simMALA(q,potential,grad_potential,dt,1.0,100000)
    print(dt,nacc,ntot,1-nacc/ntot)
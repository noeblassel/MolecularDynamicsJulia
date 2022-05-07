import numpy as np


def simMALA(q, potential, grad_potential, dt, beta, nsteps, rule="metropolis",sd_coords=None,hist=None,C_sum=None):
    M = q.shape[0]
    V = potential(q)
    F = -grad_potential(q)
    sigma = np.sqrt(2*dt)
    lambd = (beta*sigma)/2

    n_accepted = 0
    for i in range(nsteps):

        ## update autocorrelation estimators 
        if hist is not None:
            hist[1:,:]=hist[:-1,:]
            hist[0,:]=F
            C_sum+=hist*hist[0,:]

        G = np.random.standard_normal(M)
        qtilde = q+beta*dt*F+sigma*G

        #record displacement for self-diffusion
        if sd_coords is not None: disp=qtilde-q

        qtilde = qtilde % 1
        Vtilde = potential(qtilde)
        Ftilde = -grad_potential(qtilde)
        alpha = beta*(Vtilde-V)+((G+lambd * (F+Ftilde))**2)/2-(G**2)/2
        U = np.random.uniform(size=M)
        accepted = (np.log(U) < -alpha) if rule == "metropolis" else (U <
                                                                      np.exp(-alpha)/(1+np.exp(-alpha)))

        q[accepted] = qtilde[accepted]
        V[accepted] = Vtilde[accepted]
        F[accepted] = Ftilde[accepted]

        #update self-diffusion

        if sd_coords is not None: sd_coords[accepted]+=disp[accepted]
        n_accepted += np.sum(accepted)

    return (n_accepted, nsteps*M)


def simMALA_HMC(q, potential, grad_potential, dt, beta, nsteps,rule="metropolis",sd_coords=None,hist=None,C_sum=None ):
    M = q.shape[0]
    n_accepted = 0

    sigma = 1/np.sqrt(beta)
    h = np.sqrt(2*beta*dt)

    p = sigma*np.random.standard_normal(M)
    V=potential(q)
    

    for i in range(nsteps):
        H = p**2/2+V

        ## update autocorrelation estimators 
        if hist is not None:
            hist[1:,:]=hist[:-1,:]
            hist[0,:]=grad_potential(q)
            C_sum+=hist*hist[0,:]

        # position verlet (ABA) with timestep h from random momentum
        q_tilde = q+h*p/2
        # only need one gradient computation !!!
        grad = grad_potential(q_tilde)
        p -= h*grad
        q_tilde += h*p/2
        if sd_coords is not None: disp=q_tilde-q
        q_tilde = q_tilde % 1

        # the reverse transition term comes from the time symmetry property of the ABA scheme
        V_tilde=potential(q_tilde)
        Htilde = p**2/2+V_tilde
        alpha = beta*(Htilde-H)  # rewrite in
        U = np.random.uniform(size=M)

        accepted = (np.log(U) < -alpha) if rule == "metropolis" else (U <
                                                                      np.exp(-alpha)/(1+np.exp(-alpha)))

        q[accepted] = q_tilde[accepted]
        V[accepted] = V_tilde[accepted]
        n_accepted += np.sum(accepted)

        if sd_coords is not None: sd_coords[accepted]+=disp[accepted]

        p = sigma*np.random.standard_normal(M)

    return (n_accepted, nsteps*M)


def simMALA_implicit(q, potential, grad_potential, dt, beta, nsteps, histogram, rule="metropolis", fp_iter=5, tol=1e-3):
    M = q.shape[0]
    V = potential(q)
    sigma = np.sqrt(2*dt)
    N_bins = histogram.shape[0]

    n_accepted = 0

    for i in range(nsteps):

        if N_bins:
            ixs = np.floor(q*N_bins).astype('int32')
            ixs[ixs>N_bins-1]=N_bins-1
            histogram[ixs] += 1

        n_iter = 0

        q_tilde = np.copy(q)
        next_q_tilde = np.zeros_like(q_tilde)

        G = np.random.standard_normal(M)

        for i in range(fp_iter):
            grad_V_tilde = grad_potential((q+q_tilde)/2)
            next_q_tilde = q-beta*dt*grad_V_tilde+sigma*G
            converged=(np.abs(q_tilde-next_q_tilde) < tol)
            q_tilde=np.copy(next_q_tilde)

            if np.all(converged):
                break

        q_backward = np.copy(q_tilde)

        for i in range(fp_iter):
            grad_V_backward = grad_potential((q_backward+q_tilde)/2)
            next_q_backward = q_tilde+beta*dt*grad_V_backward-sigma*G
            backward_converged=np.abs(q_backward-next_q_backward)<tol
            q_backward=np.copy(next_q_backward)

            if np.all(backward_converged):
                break

        q_tilde = q_tilde % 1
        V_tilde = potential(q_tilde)
        alpha = beta*(V_tilde-V-grad_V_tilde*(q_tilde-q))
        U = np.random.uniform(size=M)

        metropolis_accepted = (np.log(
            U) < -alpha) if rule == "metropolis" else (U < np.exp(-alpha)/(1+np.exp(-alpha)))
        accepted = converged  & metropolis_accepted & backward_converged
        #print(np.sum(accepted)/M)
        q[accepted] = q_tilde[accepted]
        V[accepted] = V_tilde[accepted]
        n_accepted += np.sum(accepted)

    return (n_accepted, nsteps*M)

import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import io

import parametrization_walks as PW

def limit_distr_classical(corr: float, samples):
    """
    Classical Limiting distribution
    Args:
        corr: correlation
        x: range
    Returns:
        p(x)
    """
    x = np.linspace(-10,10,samples)
    sigma = np.sqrt((1+corr)/(1-corr))
    return [x,1/np.sqrt(2*np.pi*sigma**2)*np.exp(-x**2/2/sigma**2)]

def lambda_parameter(unitary: np.ndarray, psi: np.ndarray):
    """
    Lambda Parameter for arbitrary unitary and state
    Args:
        unitary: unitary
        state: state
    Returns:
        lambda
    """
    a,b,c,d = unitary.reshape(4).tolist()
    alpha,beta = psi.tolist()
    return np.real(np.abs(a)**2 - np.abs(b)**2 + (alpha*a*np.conj(b*beta) + np.conj(a*alpha)*b*beta)/np.abs(a)**2)

def lambda_parameter_coord(theta, phi, xi):
    """
    Lambda Parameter for arbitrary unitary and state
    Args:
        theta: theta
        phi: phi
        xi: xi
    Returns:
        lambda
    """
    return np.cos(phi)**2 - np.sin(phi)**2 + 2*np.tan(theta)*np.cos(phi)*np.sin(phi)*np.cos(xi)

def limit_distr_quantum(unitary: np.ndarray, psi: np.ndarray, samples = 100):
    """
    Limiting distribution for arbitrary unitary and state
    Args:
        unitary: unitary
        state: state
        samples: number of samples in x range
    Returns:
        p(x)
    """
    a = np.abs(unitary[0][0])
    x = np.linspace(-a+0.00001,a-0.00001,samples)
    lamb = lambda_parameter(unitary,psi)
    return [x,(np.sqrt(1-a**2)*(1-lamb*x))/np.pi*(1-x**2)/np.sqrt(a**2-x**2)]

def limit_distr_quantum_param(theta, phi, xi, samples = 100):
    """
    Limiting distribution for arbitrary unitary and state
    Args:
        theta: theta
        phi: phi
        xi: xi
        samples: number of samples in x range
    Returns:
        p(x)
    """
    a = np.cos(theta)
    x = np.linspace(-a+0.00001,a-0.00001,samples)
    lamb = lambda_parameter_coord(theta,phi,xi)
    return [x,(np.sqrt(1-a**2)*(1-lamb*x))/np.pi*(1-x**2)/np.sqrt(a**2-x**2)]

def tab_limit_classical():
    samples = 1000
    corrs = np.linspace(-1,1,10)
    xs = [limit_distr_classical(corr,samples)[0] for corr in corrs][0]
    limit_dists = [limit_distr_classical(corr,samples)[1] for corr in corrs]
    with io.open(f"Tabs/Classical_Distribution.dat","w") as file:
        for i,x in enumerate(xs):
            file.write(str(x)+"\t"+"\t".join([str(limit_dist[i]) for limit_dist in limit_dists])+"\n")

def tab_limit_quant():
    """
    Table for plot of non-symmetric limiting distributions
    """
    num_phis = 10
    samples = 1000
    theta = 0.4*np.pi
    phiss = [np.linspace(theta/2,np.pi/4+theta/2,num_phis),np.linspace(np.pi/4-theta/2,(np.pi - theta)/2,num_phis)]
    xis = [0,np.pi]
    limit_dists = [[limit_distr_quantum_param(theta,phi,xi,samples) for phi in phis] for xi,phis in zip(xis,phiss)]
    with io.open(f"Tabs/Quantum_Distribution.dat","w") as file:
        file.write(str(limit_dists[0][0][0][0]-0.000001) + "\t" + "\t".join(["0.001" for _ in range(2*num_phis+2)])+"\n")
        for i in range(samples):
            file.write(str(limit_dists[0][0][0][i])+"\t"+"\t".join(["\t".join([str(limit_dist[1][i]) for limit_dist in limit_distss]) for limit_distss in limit_dists])+"\n")
        file.write(str(limit_dists[0][0][0][-1]+0.000001) + "\t" + "\t".join(["0.001" for _ in range(2*num_phis+2)])+"\n")

def plot_distributions():
    """
    Plot of the unique parametrization of the limiting distributions
    """
    num_thetas = 5
    num_phis = 10
    samples = 100
    thetas = np.linspace(0,np.pi/2,num_thetas)
    xis = [0,np.pi]
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_ylim(8*10**(-2),100)
    for i,theta in enumerate(thetas):
        for xi in xis:
            if xi == 0:
                phis = np.linspace(0,theta/2,num_phis)
            elif xi == np.pi:
                phis = np.linspace(0,np.pi/2-theta/2,num_phis)
            for phi in phis:
                unitary = PW.unitary_theta(theta)
                state = PW.psi_param(phi,xi)
                x,distr = limit_distr_quantum_param(theta,phi,xi,samples)
                ax.plot(x,distr)
    plt.savefig("Plots/Limiting_All.png",dpi=150)

def plot_parameters():
    """
    Plot of the unique parametrization of parameters
    """
    num_thetas = 10
    samples = 100
    thetas = np.linspace(0,np.pi/2,num_thetas)
    xis = [0,np.pi]
    fig,ax = plt.subplots()
    ax.set_ylim(-5,5)
    for theta in thetas:
        for xi in xis:
            unitary = PW.unitary_param(theta,0,0)
            if xi == 0:
                phis = np.linspace(theta/2,np.pi/4+theta/2,samples)
                states = [PW.psi_param(phi,0) for phi in phis]
                parameters = [lambda_parameter_coord(theta,phi,xi) for phi in phis]
                index_max = np.argmax(parameters)
                parameter_max = parameters[index_max]
                phi_max = phis[index_max]
                ax.scatter([phi_max],[parameter_max])
            elif xi == np.pi:
                phis = np.linspace(np.pi/4-theta/2,(np.pi - theta)/2,samples)
                states = [PW.psi_param(phi,np.pi) for phi in phis]
                parameters = [lambda_parameter_coord(theta,phi,xi) for phi in phis]
                index_min = np.argmin(parameters)
                parameter_min = parameters[index_min]
                phi_min = phis[index_min]
                ax.scatter([phi_min],[parameter_min])
            # parameters = [lambda_parameter(unitary,state) for state in states]
            ax.plot(phis,parameters)
    ax.plot(np.linspace(0,np.pi/4,samples),np.ones(samples))
    plt.savefig("Plots/Plot_Parameter_Limiting.png",dpi=150)


if __name__ == "__main__":
    plot_distributions()
    plot_parameters()
    tab_limit_classical()
    tab_limit_quant()


import numpy as np
import io
import matplotlib.pyplot as plt
from numpy import linalg as la
from qiskit import quantum_info as qi

import limit_theorem as LT
import parametrization_walks as PW
import simulation_walks as SW


def walk_variance(theta,time):
    return SW.simulation_symmetric_quantum_walk(PW.unitary_symmetric(theta),time)["variance"]

# Quantum
def plot_variance_quantum(angles: np.ndarray, time: int = 20, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("t")
    ax.set_ylabel("Var")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([10**(-1),time])
    t = np.arange(1,time+1,1)
    for angle in angles:
        variances = walk_variance(angle,time)
        ax.plot(t,variances[1:],c=str(0.7*np.abs(np.cos(angle))))
        ax.plot(t,np.sqrt(1-np.sin(angle))*t,c=str(0.7*np.abs(np.cos(angle))))
    plt.savefig(f"Plots/Walk_Variance{name}.png",dpi=150)

def tab_variance_quantum(angles: np.ndarray, time: int = 20, name: str = ""):
    ts = np.arange(1,time+1,1)
    variances = [walk_variance(angle,time) for angle in angles]
    with io.open(f"Tabs/Variances_Quantum_{name}.dat","w") as file:
        for i,t in enumerate(ts):
            file.write(str(t)+"\t"+"\t".join([str(variance[i]) for variance in variances])+"\n")

def correlated_variance(corr,bias,time):
    return SW.simulation_classically_correlated_walk(corr,bias,time)["variance"]

# Classical
def plot_variance_classical(corrs: float, bias: float = 0.5, time: int = 20, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("t")
    ax.set_ylabel("Var")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([10**(-1),time])
    t = np.arange(1,time+1,1)
    for corr in corrs:
        variances = correlated_variance(corr,bias,time)
        ax.plot(t,variances[1:],c=str(0.7*np.abs(corr)))
        ax.plot(t,np.sqrt((1+corr)/(1-corr))*np.sqrt(t),c=str(0.7*np.abs(corr)))
    plt.savefig(f"Plots/Correlated_Variance_{name}.png",dpi=150)
    plt.close()

def tab_variance_classical(corrs: list, bias: float = 0.5, time: int = 20, name: str = ""):
    ts = np.arange(1,time+1,1)
    variances = [correlated_variance(corr,bias,time) for corr in corrs]
    with io.open(f"Tabs/Correlated_Variance_{name}.dat","w") as file:
        for i,t in enumerate(ts):
            file.write(str(t)+"\t"+"\t".join([str(variance[i]+10**(-10)) for variance in variances])+"\n")
    with io.open(f"Tabs/Correlated_Variance_{name}_ext.dat","w") as file:
        for i,t in enumerate(ts):
            file.write(str(t)+"\t"+"\t".join([str(np.sqrt((1+corr)/(1-corr))*np.sqrt(t)) for corr in corrs])+"\n")

def plot_probs(t: int = 10, val: int = 1, theta: float = 0.1*np.pi, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    thgamma = np.linspace(0,0.25*np.pi,100)
    phi = np.linspace(0,0.5*np.pi,100)
    TH,PHI = np.meshgrid(thgamma,phi)
    ax.set_xlabel("phi")
    ax.set_ylabel("th")
    ax.set_zlabel("prob")
    ax.set_ylim([0,0.25*np.pi])
    ax.set_xlim([0,0.5*np.pi])
    prob = np.array([[np.cos(theta)**(-2*(t-1))*LT.probs_a_new_type(PW.unitary_theta(theta),PW.psi_param(thgamma,phi),t)[val] for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
    ax.scatter(PHI,TH,prob)
    plt.savefig(f"Plots/plot_probs_{name}.png",dpi=150)

def plot_entropies(t: int = 10, theta: float = 0.1*np.pi, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    thgamma = np.linspace(0,0.25*np.pi,100)
    phi = np.linspace(0,0.5*np.pi,100)
    TH,PHI = np.meshgrid(thgamma,phi)
    ax.set_xlabel("phi")
    ax.set_ylabel("th")
    ax.set_zlabel("entropy")
    ax.set_ylim([0,0.25*np.pi])
    ax.set_xlim([0,0.5*np.pi])
    def entropy(probs):
        probs = probs[probs > 0]
        return - np.sum(np.log(probs)*np.array(probs))
    entropies = np.array([[entropy(LT.probs_a_new_type(PW.unitary_theta(theta),PW.psi_param(thgamma,phi),t)) for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
    ax.scatter(PHI,TH,entropies)
    plt.savefig(f"Plots/plot_entropies_{name}.png",dpi=150)

def plot_probs(t: int = 10, val: int = 1, theta: float = 0.1*np.pi, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    thgamma = np.linspace(0,0.25*np.pi,100)
    phi = np.linspace(0,0.5*np.pi,100)
    TH,PHI = np.meshgrid(thgamma,phi)
    ax.set_xlabel("phi")
    ax.set_ylabel("th")
    ax.set_zlabel("prob")
    ax.set_ylim([0,0.25*np.pi])
    ax.set_xlim([0,0.5*np.pi])
    prob = np.array([[np.cos(theta)**(-2*(t-1))*LT.probs_a_new_type(PW.unitary_theta(theta),PW.psi_param(thgamma,phi),t)[val] for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
    ax.scatter(PHI,TH,prob)
    plt.savefig(f"Plots/plot_probs_{name}.png",dpi=150)

def plot_entropies(t: int = 10, theta: float = 0.1*np.pi, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    thgamma = np.linspace(0,0.25*np.pi,100)
    phi = np.linspace(0,0.5*np.pi,100)
    TH,PHI = np.meshgrid(thgamma,phi)
    ax.set_xlabel("phi")
    ax.set_ylabel("th")
    ax.set_zlabel("entropy")
    ax.set_ylim([0,0.25*np.pi])
    ax.set_xlim([0,0.5*np.pi])
    def entropy(probs):
        probs = probs[probs > 0]
        return - np.sum(np.log(probs)*np.array(probs))
    entropies = np.array([[entropy(LT.probs_a_new_type(PW.unitary_theta(theta),PW.psi_param(thgamma,phi),t)) for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
    ax.scatter(PHI,TH,entropies)
    plt.savefig(f"Plots/plot_entropies_{name}.png",dpi=150)


if __name__ == "__main__":
    t = 150
    angles = np.pi/2*(-np.exp(-3*np.linspace(0,1,20))+(1+np.exp(-3)))
    plot_variance_quantum(angles,t)
    tab_variance_quantum(angles,t)
    
    n = 20
    bias = 0.5
    correlations = np.linspace(-1,1,n)
    plot_variance_classical(correlations,bias,t)
    tab_variance_classical(correlations,bias,t)

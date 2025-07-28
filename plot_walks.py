import numpy as np
from numpy import linalg as la
from qiskit import quantum_info as qi

# Quantum
def plot_variance(angles: np.ndarray, n: int = 50, time: int = 20, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("t")
    ax.set_ylabel("Var")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([10**(-1),time])
    t = np.arange(1,time+1,1)
    for angle in angles:
        variances = walk_variance(angle,n,time)
        ax.plot(t,variances[1:],c=str(0.7*np.abs(np.cos(angle))))
        ax.plot(t,np.sqrt(1-np.sin(angle))*t,c=str(0.7*np.abs(np.cos(angle))))
    plt.savefig(f"Plots/Walk_Variance{name}.png",dpi=150)

def tab_variance(angles: np.ndarray, n: int = 50, time: int = 20, name: str = ""):
    ts = np.arange(1,time+1,1)
    variances = [walk_variance(angle,n,time) for angle in angles]
    with io.open(f"Tabs/Variances_Quantum_{name}.dat","w") as file:
        for i,t in enumerate(ts):
            file.write(str(t)+"\t"+"\t".join([str(variance[i]) for variance in variances])+"\n")

# Classical
def plot_variance(corrs: float, bias: float = 0.5, n: int = 50, time: int = 20, name: str = ""):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("t")
    ax.set_ylabel("Var")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([10**(-1),time])
    t = np.arange(1,time+1,1)
    for corr in corrs:
        variances = correlated_variance(corr,bias,n,time)
        ax.plot(t,variances[1:],c=str(0.7*np.abs(corr)))
        ax.plot(t,np.sqrt((1+corr)/(1-corr))*np.sqrt(t),c=str(0.7*np.abs(corr)))
    plt.savefig(f"Plots/Correlated_Variance_{name}.png",dpi=150)
    plt.close()

def tab_variance(corrs: list, bias: float = 0.5, n: int = 50, time: int = 20, name: str = ""):
    ts = np.arange(1,time+1,1)
    variances = [correlated_variance(corr,bias,n,time) for corr in corrs]
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
    prob = np.array([[np.cos(theta)**(-2*(t-1))*probs_a_new_type(unitary(theta),psi_param(thgamma,phi),t)[val] for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
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
    entropies = np.array([[entropy(probs_a_new_type(unitary(theta),psi_param(thgamma,phi),t)) for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
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
    prob = np.array([[np.cos(theta)**(-2*(t-1))*probs_a_new_type(unitary(theta),psi_param(thgamma,phi),t)[val] for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
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
    entropies = np.array([[entropy(probs_a_new_type(unitary(theta),psi_param(thgamma,phi),t)) for thgamma,phi in zip(thgammas,phis)] for thgammas,phis in zip(TH,PHI)])
    ax.scatter(PHI,TH,entropies)
    plt.savefig(f"Plots/plot_entropies_{name}.png",dpi=150)


if __name__ == "__main__":
    # Classical
    n = 300
    t = 150
    angles = np.linspace(0,np.pi/2,20)
    plot_variance(angles,n,t)
    tab_variance(angles,n,t)
    
    n = 20
    correlations = np.linspace(-1,1,n)
    correlations = -np.exp(-np.linspace(1,10,n))+1
    plot_variance(correlations,0.5,300,150)
    tab_variance(correlations,0.5,300,150)

    t = 10
    plot_entropies(t)
    t = 3
    n = 9
    for k in range(0,t,2):
        plot_probs(t,val=k,name=f"t{t}_k{k}")

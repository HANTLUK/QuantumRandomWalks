import numpy as np
import matplotlib.pyplot as plt

def prob(thetaW: np.ndarray, thetaS: np.ndarray, phi: np.ndarray):
    return np.cos(thetaW)**2*np.cos(thetaS)**2 + np.sin(thetaW)**2*np.sin(thetaS)**2 - 2*np.cos(thetaW)*np.sin(thetaW)*np.cos(thetaS)*np.sin(thetaS)*np.cos(phi)

def plot_parameters():
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel("thetaS")
    ax.set_ylabel("phi")
    ax.set_zlabel("prob")
    thetaW = np.linspace(0,np.pi/2,10)
    thetas = np.linspace(0,np.pi,100)
    phis = np.linspace(0,np.pi,100)
    TS,PHI = np.meshgrid(thetas,phis)
    for thetaw in thetaW:
        probs = prob(thetaw,TS,PHI)
        ax.plot_surface(TS,PHI,probs)
    plt.savefig(f"Plots/plot_parameter.png",dpi=150)

if __name__ == "__main__":
    plot_parameters()

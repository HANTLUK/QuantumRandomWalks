import numpy as np
import matplotlib.pyplot as plt

def parameters(alpha: float, kz: float, theta: float):
    """
    Calculates kx and ky for U = cos(alpha) 1 + sum_i i sin(alpha) k_i sigma_i
    Args:
        alpha: overall phase
        kz: z component of unitary
        theta: phase of the state?
    Returns:
        U_1, U_2: [alpha,[kx,ky,kz]]
    """
    absXY = 1. - kz**2
    phaseA = np.arctan2(np.cos(alpha),kz*np.sin(alpha))
    phaseC1 = (np.pi/2 + phaseA + theta) % (2*np.pi)
    phaseC2 = (3*np.pi/2 + phaseA + theta) % (2*np.pi)
    X1 = absXY*np.sin(phaseC1)
    X2 = absXY*np.sin(phaseC2)
    Y1 = -absXY*np.cos(phaseC1)
    Y2 = -absXY*np.cos(phaseC2)
    return [[alpha,np.array([X1,Y1,kz])],[alpha,np.array([X2,Y2,kz])]]

def plot_symmetry_parameters():
    alphas = np.linspace(0,np.pi/2,10)
    thetas = np.linspace(0,np.pi,10)
    alpha = 0.0
    theta = 0.3
    kzs = np.linspace(-1,1,100)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    for theta in thetas:
        params = [parameters(alpha, kz, theta) for kz in kzs]
        x1 = [param[0][1][0] for param in params]
        y1 = [param[0][1][1] for param in params]
        x2 = [param[1][1][0] for param in params]
        y2 = [param[1][1][1] for param in params]
        ax.plot(x1,y1,kzs, color='#F33333')
        ax.plot(x2,y2,kzs, color='#33F333')
    plt.savefig(f"Plots/Condition_theta.png",dpi=150)    

if __name__ == "__main__":
    plot_symmetry_parameters()

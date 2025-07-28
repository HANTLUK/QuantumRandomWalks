import numpy as np
from numpy import linalg as la
from qiskit import quantum_info as qi
import copy
from functools import reduce
import io

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

import parametrization_walks as PW
import simulation_walks as SW
import limiting_distribution as LD

def animation_walk(walk_output, name: str = ""):
    """
    Animation of the walk
    Args:
        walk_output: dict from simulation of the walk
    Returns:
        None
    """
    prob_list = walk_output["probs"]
    n = max([prob.shape[0] for prob in prob_list])
    time = len(prob_list)
    x = range(-(n // 2), n // 2 + 1)

    fig, ax = plt.subplots()
    ax.set_xlim(-(n // 2), n // 2)
    ax.set_ylim(0, 1)
    bars = copy.copy(ax.bar(x, [0] * len(x), alpha=0.7))

    def init():
        for bar in bars:
            bar.set_height(0)
        for rectangle in bars.patches:
            yield rectangle

    def update(frame):
        prob = prob_list[frame]
        for bar, height in zip(bars, prob):
            bar.set_height(height)
        for rectangle in bars.patches:
            yield rectangle

    ani = FuncAnimation(
        fig,
        update,
        frames=time,
        init_func=init,
        blit=True,
        interval=1000
    )
    ani.save(f"Animations/{name}.gif", writer=PillowWriter(fps=10))

# Test
if __name__ == "__main__":
    # Classical: Animate same bias, different correlation
   #  bias = 0.99
   #  num_corrs = 5
   #  t = 20
   #  correlations = np.linspace(-1,1,num_corrs)
   #  for i,corr in enumerate(correlations):
   #      walk = SW.simulation_classically_correlated_walk(corr,bias,t,animation=True)
   #      animation_walk(walk,f"Animation_CC_{i}")


    # Quantum: Animate symmetric walks with different theta
   #  num_thetas = 5
   #  t = 20
   #  thetas = np.linspace(0,np.pi/2,num_thetas)
   #  for i,theta in enumerate(thetas):
   #      unitary = PW.unitary_symmetric(theta)
   #      walk = SW.simulation_symmetric_quantum_walk(unitary,t,animation=True)
   #      animation_walk(walk,f"Animation_SQ_{i}")

    # Quantum: Test of asymmetric walks
    t = 20
    theta = 1.2
    unitary = PW.unitary_param(theta,0,0)
    phi_1, xi_1 = [0.2,0]
    lambda_1 = LD.lambda_parameter_coord(theta,phi_1,xi_1)
    psi_1 = PW.psi_param(phi_1,xi_1)
    phi_2, xi_2 = [1.0,0]
    psi_2 = PW.psi_param(phi_2,xi_2)
    lambda_2 = LD.lambda_parameter_coord(theta,phi_2,xi_2)
    print(lambda_1,lambda_2)
    psis = [psi_1,psi_2]
    for i,psi in enumerate(psis):
        walk = SW.simulation_walk(unitary,psi,t,animation=True)
        animation_walk(walk,f"Animation_Asymmetric_{i}")

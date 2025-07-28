import numpy as np
from numpy import linalg as la
from qiskit import quantum_info as qi
import copy
from functools import reduce

import io

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

import parametrization_walks as PW


# projections
ZERO_PROJ: np.ndarray = np.array([[1, 0], [0, 0]])
ONE_PROJ: np.ndarray = np.array([[0, 0], [0, 1]])

def can_vector(dim: int, i: int) -> np.ndarray:
    result = np.zeros(dim)
    result[i] = 1
    return result
def can_matrix(dim: int, i: int, j: int) -> np.ndarray:
    result = np.zeros((dim, dim))
    result[i, j] = 1
    return result

def probs_classical(vec: np.ndarray) -> list[float]:
    main_dim = len(vec) // 2
    return np.array([np.abs(vec[i]) + np.abs(vec[main_dim + i]) for i in range(main_dim)])

def probs_quantum(vec: np.ndarray) -> list[float]:
    main_dim = len(vec) // 2
    return np.array([np.abs(vec[i]) ** 2 + np.abs(vec[main_dim + i])**2 for i in range(main_dim)])

def states(vec: np.ndarray) -> list[float]:
    main_dim = len(vec) // 2
    return np.array([[vec[i],vec[main_dim+i]] for i in range(main_dim)])

def simulation_symmetric_quantum_walk(unitary,time,animation=False):
    """
    Simulates a symmetric walk with unitary U
    Args:
        unitary: Unitary
        time: time
    """
    psi_1,psi_2 = PW.symmetric_initial_state(unitary)
    return simulation_walk(unitary,psi_1,time,animation=animation)

def simulation_classically_correlated_walk(correlation,bias,time,animation=False):
    """
    Simulates a correlated classical walk
    Args:
        correlation: correlation parameter
        bias: bias parameter
        time: time parameter
    Returns:
        walk: state, prob, variance, entropy
    """
    initial_state = PW.classically_correlated_state(bias)
    walk_matrix = PW.classically_correlated_matrix(correlation)
    return simulation_walk(walk_matrix,initial_state,time,probs_classical,animation=animation)

def simulation_walk(walk_matrix, initial_state, time, probabilities = probs_quantum, animation = False):
    """
    Simulates the walk
    Args:
        walk_matrix: matrix
        initial_state: initial_state
        time: time
        probabilities: [quantum] or classical
    """
    n = 2*time + 1

    state_list = []
    probs_list = []
    var_list = []

    def state_probs(state_vec,prob):
        if not animation:
            prob = np.array(prob[prob>0])
            state_vec = np.array(state_vec[prob > 0])
        else:
            prob = np.array(prob)
            state_vec = np.array(state_vec)
        return [state_vec,prob]

    walker = reduce(np.add,[np.kron(ZERO_PROJ, can_matrix(n,(i-1)%n,i)) + np.kron(ONE_PROJ,can_matrix(n,(i+1)%n,i)) for i in range(n)]) @ np.kron(walk_matrix,np.eye(n))
    state = np.kron(initial_state,can_vector(n,n // 2))

    state_vec = states(state)
    prob = probabilities(state)
    variance = np.sqrt(np.sum(prob*np.square(np.arange(0,n,1)-n//2)))
    state_vec,prob = state_probs(state_vec,prob)

    probs_list.append(prob)
    state_list.append(state_vec)
    var_list.append(variance)

    for step in range(time):
        state = walker @ state

        state_vec = states(state)
        prob = probabilities(state)
        variance = np.sqrt(np.sum(prob*np.square(np.arange(0,n,1)-n//2)))
        state_vec,prob = state_probs(state_vec,prob)

        probs_list.append(prob)
        state_list.append(state_vec)
        var_list.append(variance)

    var_list = np.array(var_list)
    return {"states":state_list,"probs":probs_list,"variance":var_list}

if __name__ == "__main__":
    t = 5
    
    # Quantum Walk
    theta, phi_1, phi_2 = [0.2*np.pi,0.4*np.pi,0]
    phi,xi = [0.25*np.pi,0]
    unitary = PW.unitary_param(theta,phi_1,phi_2)
    psi = PW.psi_param(phi,xi)
    walk_quantum = simulation_walk(unitary,psi,t)
    print("Variance Quantum: \n",walk_quantum["variance"])

    # Symmetric Walk
    unitary = PW.unitary_symmetric(theta)
    walk_symmetric = simulation_symmetric_quantum_walk(unitary,t)
    print("States Symmetric: \n",walk_symmetric["states"])

    # Classical Walk
    bias,corr = [0.5,0.7]
    walk_classical = simulation_classically_correlated_walk(corr,bias,t)
    print("Probs Classical: \n",walk_classical["probs"])

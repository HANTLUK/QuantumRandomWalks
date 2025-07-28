import numpy as np
from numpy import linalg as la
from qiskit import quantum_info as qi


def classically_correlated_matrix(correlation_parameter: float) -> np.ndarray:
    """
    Matrix for the classically correlated_ walk
    Args:
        correlation_parameter: parametrization
    Returns:
        matrix: parametrized by corr
    """
    corr = correlation_parameter
    return np.array([[(1+corr)/2,(1-corr)/2],[(1-corr)/2,(1+corr)/2]])

def classically_correlated_state(bias_parameter: float) -> np.ndarray:
    """
    Initial state for the classically correlated_ walk
    Args:
        bias_parameter: parametrization
    Returns:
        state: parametrized by corr
    """
    bias = bias_parameter
    return np.array([bias,1.-bias])

def unitary_param(theta: float, phi_1: float, phi_2: float):
    """
    Unitary in the Hopf Parametrization
    Args:
        theta: theta
        phi_1: phi_1
        phi_2: phi_2
    Returns:
        unitary: parametrized by Hopf parameters
    """
    return np.array([[np.cos(theta)*np.exp(1.j*(phi_1+phi_2)/2),np.sin(theta)*np.exp(1.j*(phi_1 - phi_2)/2)],[-np.sin(theta)*np.exp(-1.j*(phi_1 - phi_2)/2),np.cos(theta)*np.exp(-1.j*(phi_1+phi_2)/2)]])

def unitary_theta(theta: float):
    """
    Unitary in the Hopf Parametrization
    Args:
        theta: theta
    Returns:
        unitary: parametrized by Hopf parameter theta
    """
    return np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])

def unitary_symmetric(theta: float):
    """
    Unitary for symmetric walks
    Args:
        theta: angle for amplitudes
    Returns:
        unitary: (cos(theta), i sin(theta), i sin(theta), cos(theta))
    """
    return np.array([[np.cos(theta),1.j*np.sin(theta)],[1.j*np.sin(theta),np.cos(theta)]])

def psi_param(phi: float, xi: float):
    """
    Initial State
    Args:
        phi: angle for amplitudes
        xi: relative phase
    Returns:
        psi = cos(phi) up + sin(phi) exp(i xi) down
    """
    return np.array([np.cos(phi),np.sin(phi)*np.exp(1.j*xi)])

def symmetric_initial_state(U: np.ndarray):
    """
    Calculates the states that give a symmetric walk for that unitary.
    Args:
        U: Unitary
    Returns:
        psi_1, psi_2: States
    """
    U_dagger = U.conj().T
    ARGA = np.angle(U[0][0])
    ARGC = np.angle(U[1][0])
    theta1 = (np.pi/2 + ARGC - ARGA) % (2*np.pi)
    theta2 = (3*np.pi/2 + ARGC - ARGA) % (2*np.pi)
    psi_1 = 1/np.sqrt(2) * np.array([1,np.exp(1.j*theta1)])
    psi_2 = 1/np.sqrt(2) * np.array([1,np.exp(1.j*theta2)])
    psi_1 = np.dot(U_dagger,psi_1)
    psi_2 = np.dot(U_dagger,psi_2)
    return [psi_1, psi_2]

def check_symmetric_initial_state(U: np.ndarray):
    """
    Checks whether the initial state(s) and the unitary are compatible
    """
    U_dagger = U.conj().T
    Id = U_dagger @ U
    Id_T = U @ U_dagger
    psi_1,psi_2 = symmetric_initial_state(U)
    psi_1 = np.dot(U,psi_1)
    psi_2 = np.dot(U,psi_2)
    for psi in [psi_1,psi_2]:
        theta = np.angle(psi[1]/psi[0])
        ARGAC = np.angle(U[1][0]) - np.angle(U[0][0]) + np.pi/2
        print(f"Check for Args: {theta}, {ARGAC}, {ARGAC + np.pi}")
    return True

# Tests
if __name__ == "__main__":
    theta, phi_1, phi_2 = [0.5*np.pi,0.4*np.pi,0]
    U = unitary_param(theta,phi_1,phi_2)
    print("Unitary: \n",U)
    phi,xi = [0.25*np.pi,0]
    psi = psi_param(phi,xi)
    print("State: \n",psi)
    U_symm = unitary_symmetric(theta)
    psi_symm_1,psi_symm_2 = symmetric_initial_state(U_symm)
    print("Symmetric :\n",U_symm,psi_symm_1,psi_symm_2)
    check_symmetric_initial_state(U_symm)
    bias,corr = [0.5,0.7]
    matrix = classically_correlated_matrix(corr)
    print("Matrix: \n",matrix)
    state = classically_correlated_state(bias)
    print("State: \n",state)

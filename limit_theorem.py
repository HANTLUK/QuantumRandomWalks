import numpy as np
from numpy import linalg as la
from qiskit import quantum_info as qi
import copy
from functools import reduce
import math
import matplotlib.pyplot as plt

import parametrization_walks as PW
import simulation_walks as SW

def kappa(gamma: int, delta: int, n: int, k: int):
    """
    Combinatorial Prefactor kappa from "A new type of limit theorem"
    Args:
        gamma, delta, n, k
    Returns:
        Expression
    """
    return math.comb(k-1,gamma-1)*math.comb(k-1,delta-1)*math.comb(n-k-1,gamma-1)*math.comb(n-k-1,delta-1)

def probs_a_new_type(unitary:np.ndarray, psi:np.ndarray, t: int):
    """
    Probabilities from "A new type o limit theorem"
    Args:
        unitary: walk matrix
        psi: initial_state
        t: time
    Returns:
        probs: probability distribution
    """
    n = t
    a,b,c,d = unitary.reshape(4).tolist()
    alpha, beta = psi.tolist()
    pre = (-np.abs(b)**2/np.abs(a)**2)
    probs = np.zeros(2*n+1-n%2)
    aabcbc = a*alpha*np.conj(b)*np.conj(beta)
    acacbb = np.conj(a)*np.conj(alpha)*b*beta
    va = np.abs(a)**2
    vb = np.abs(b)**2
    vbeta = np.abs(beta)**2
    valpha = np.abs(alpha)**2
    for k in range(1,n//2+1):
        prob = 0.
        probneg = 0.
        for gamma in range(1,k+1):
            for delta in range(1,k+1):
                probsum = 0.
                probsum += (k**2*vb+(n-k)**2*va-(gamma+delta)*(n-k))*valpha
                probsum += (k**2*va+(n-k)**2*vb-(gamma+delta)*k)*vbeta
                probsum += 1/vb*(((n-k)*gamma-k*delta+n*(2*k-n)*vb)*aabcbc + (-k*gamma+(n-k)*delta+n*(2*k-n)*vb)*acacbb + gamma*delta)
                probsum *= pre**(gamma+delta)*kappa(gamma,delta,n,k)/(gamma*delta)
                prob += probsum

                probsum = 0.
                probsum += (k**2*vb+(n-k)**2*va-(gamma+delta)*k)*valpha
                probsum += (k**2*va+(n-k)**2*vb-(gamma+delta)*(n-k))*vbeta
                probsum += 1/vb*((k*gamma-(n-k)*delta-n*(2*k-n)*vb)*aabcbc + (-(n-k)*gamma+k*delta-n*(2*k-n)*vb)*acacbb + gamma*delta)
                probsum *= pre**(gamma+delta)*kappa(gamma,delta,n,k)/(gamma*delta)
                probneg += probsum

        probs[-2*k-1] = np.real(va**(n-1)*prob)
        probs[2*k] = np.real(va**(n-1)*probneg)
    probs[-1] = np.real(va**(n-1)*(vb*valpha+va*vbeta-(aabcbc+acacbb)))
    probs[0] = np.real(va**(n-1)*(va*valpha+vb*vbeta+(aabcbc+acacbb)))
    return probs[probs > 0]

def test_random(t):
    """
    Tests Probabilities for random unitary
    Args:
        t: time
    Returns:
        None
    """
    unitary = qi.random_unitary(2).data
    unitary /= np.sqrt(la.det(unitary))
    psi,p = PW.symmetric_initial_state(unitary)
    probs_simulation = SW.simulation_symmetric_quantum_walk(unitary,t)["probs"][-1]
    probs_limit = probs_a_new_type(unitary,psi,t)
    print(probs_simulation)
    print(probs_limit)

if __name__ == "__main__":
    t = 4
    test_random(t)

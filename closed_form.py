import numpy as np
from numpy import linalg as la
from qiskit import quantum_info as qi
import math

from scipy.special import gamma

import simulation_walks as SW
import parametrization_walks as PW

def fact(n: float):
    return math.factorial(n)

def closed_form_classical(corr: float, bias: float, t: int):
    """
    Calculates the closed form analog for the classical walk (Our calculation)
    Args:
        corr: correlation_parameter
        bias: bias parameter
        t: time
    Returns:
        probs: probabilities
        state: amplitudes
    """
    alphat = []
    betat = []
    probt = []
    beta0,alpha0 = [bias,1. - bias]
    c = corr
    for x in range(-t,t+1,2):
        alphaxt = 0.
        betaxt = 0.
        for h in range(0,t+1):
            if (h+t)%2 == 0 and h+x >= 0 and t+h >= 0 and h-x >= 0 and t-h >= 0:
                alphaxt += (-c)**((t-h)//2)*fact((t+h)//2)/(fact((t-h)//2)*fact((h-x)//2)*fact((h+x)//2))*((c+1)/2)**h*alpha0
                betaxt += (-c)**((t-h)//2)*fact((t+h)//2)/(fact((t-h)//2)*fact((h-x)//2)*fact((h+x)//2))*((c+1)/2)**h*beta0
        for h in range(0,t):
            if (t-1+h)%2 == 0 and t-1-h >= 0 and t-1+h >= 0 and h+1+x >= 0 and h-1-x >= 0:
                alphaxt -= (-c)**((t-1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h+1+x)//2)*fact((h-1-x)//2))*((c+1)/2)**(h+1)*alpha0
                betaxt += (-c)**((t-1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h-1-x)//2)*fact((h+1+x)//2))*((c+1)/2)**h*((1-c)/2)*alpha0
            if (t-1+h)%2 == 0 and t-1+h >= 0 and t-1-h >= 0 and h-1+x >= 0 and h+1-x >= 0:
                alphaxt += (-c)**((t-1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h-1+x)//2)*fact((h+1-x)//2))*((c+1)/2)**h*((1-c)/2)*beta0
                betaxt -= (-c)**((t-1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h+1-x)//2)*fact((h-1+x)//2))*((c+1)/2)**(h+1)*beta0
        betat.append(betaxt)
        alphat.append(alphaxt)
        probt.append(np.abs(betaxt)+np.abs(alphaxt))
    state = np.array([[beta,alpha] for alpha,beta in zip(alphat,betat)])
    return {"probs":probt,"state":state}

def closed_form(unitary: np.ndarray, psi: np.ndarray, t: int):
    """
    Calculates the expressions from "closed-form expressions" for comparison
    Args: 
        unitary: unitary
        psi: state
        t: time
    Returns:
        probs: probabilities
        state: amplitudes
    """
    a,b,c,d = unitary.reshape(4).tolist()
    exphi12 = - d/a
    exphi1 = -np.conj(c*d/(a*b))*exphi12
    exphi2 = exphi12*np.conj(exphi1)
    costheta = np.abs(a)
    sintheta = np.abs(b)
    beta0, alpha0 = psi.tolist()
    alphat = []
    betat = []
    probt = []
    for x in range(-t,t+1,2):
        alphaxt = 0.
        betaxt = 0.
        for h in range(0,t+1):
            if (h+t)%2 == 0 and h+x >= 0 and t+h >= 0 and h-x >= 0 and t-h >= 0:
                alphaxt += (-1)**((h+x)//2)*fact((t+h)//2)/(fact((t-h)//2)*fact((h-x)//2)*fact((h+x)//2))*costheta**h*exphi12**((t+x)//2)*alpha0
                betaxt += (-1)**((h+x)//2)*fact((t+h)//2)/(fact((t-h)//2)*fact((h-x)//2)*fact((h+x)//2))*costheta**h*exphi12**((t+x)//2)*beta0
        for h in range(0,t):
            if (t-1+h)%2 == 0 and t-1-h >= 0 and t-1+h >= 0 and h+1-x >= 0 and h-1+x >= 0:
                alphaxt += (-1)**((h-1+x)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h+1-x)//2)*fact((h-1+x)//2))*costheta**(h+1)*exphi12**((t+x)//2)*alpha0
                betaxt += (-1)**((h-1+x)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h+1-x)//2)*fact((h-1+x)//2))*costheta**h*sintheta*exphi12**((t+x)//2)*np.conj(exphi1)*alpha0
            if (t-1+h)%2 == 0 and t-1+h >= 0 and t-1-h >= 0 and h-1-x >= 0 and h+1+x >= 0:
                alphaxt += (-1)**((h+1+x)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h-1-x)//2)*fact((h+1+x)//2))*costheta**h*sintheta*exphi12**((t+x)//2)*exphi1*beta0
                betaxt -= (-1)**((h+1+x)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h-1-x)//2)*fact((h+1+x)//2))*costheta**(h+1)*exphi12**((t+x)//2)*beta0
        betat.append(betaxt)
        alphat.append(alphaxt)
        probt.append(np.abs(betaxt)**2+np.abs(alphaxt)**2)
    state = np.array([[alpha,beta] for alpha,beta in zip(alphat,betat)])
    return {"probs":probt,"state":state}

def corrected_closed_form(unitary: np.ndarray, psi: np.ndarray, t: float):
    """
    Calculates the corrected closed form expression (Our Work)
    Args:
        unitary: unitary
        psi: state
        t: time
    Returns:
        state: final state
        probs: probabilities
    """
    a,b,c,d = unitary.reshape(4).tolist()
    costheta = np.abs(a)
    sintheta = np.abs(b)
    exphis12 = a/costheta
    exphid12 = b/sintheta
    exmphis12 = d/costheta
    exmphid12 = np.conj(exphid12)
    alpha0,beta0 = psi.tolist()
    alphat = []
    betat = []
    probt = []
    for x in range(-t,t+1,2):
        alphaxt = 0.
        betaxt = 0.
        for h in range(0,t+1):
            if (h+t)%2 == 0 and h+x >= 0 and t+h >= 0 and h-x >= 0 and t-h >= 0:
                alphaxt += (-1)**((t-h)//2)*fact((t+h)//2)/(fact((t-h)//2)*fact((h-x)//2)*fact((h+x)//2))*costheta**h*exmphis12**(-x)*alpha0
                betaxt += (-1)**((t-h)//2)*fact((t+h)//2)/(fact((t-h)//2)*fact((h-x)//2)*fact((h+x)//2))*costheta**h*exmphis12**(-x)*beta0
        for h in range(0,t):
            if (t-1+h)%2 == 0 and t-1-h >= 0 and t-1+h >= 0 and h+1+x >= 0 and h-1-x >= 0:
                alphaxt += (-1)**((t+1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h+1+x)//2)*fact((h-1-x)//2))*costheta**(h+1)*exmphis12**(-x)*alpha0
                betaxt += (-1)**((t-1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h+1+x)//2)*fact((h-1-x)//2))*costheta**h*sintheta*exmphis12**(-x-1)*exmphid12*alpha0
            if (t-1+h)%2 == 0 and t-1+h >= 0 and t-1-h >= 0 and h-1+x >= 0 and h+1-x >= 0:
                alphaxt += (-1)**((t-1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h-1+x)//2)*fact((h+1-x)//2))*costheta**h*sintheta*exmphis12**(-x+1)*exphid12*beta0
                betaxt += (-1)**((t+1-h)//2)*fact((t-1+h)//2)/(fact((t-1-h)//2)*fact((h-1+x)//2)*fact((h+1-x)//2))*costheta**(h+1)*exmphis12**(-x)*beta0
        betat.append(betaxt)
        alphat.append(alphaxt)
        probt.append(np.abs(betaxt)**2+np.abs(alphaxt)**2)
    state = np.array([[alpha,beta] for alpha,beta in zip(alphat,betat)])
    return {"probs":probt,"state":state}

if __name__ == "__main__":
    # Quantum
    t = 2
    unitary = qi.random_unitary(2).data
    unitary /= np.sqrt(la.det(unitary))
    psi,p = PW.symmetric_initial_state(unitary)
    probs_closed_form = closed_form(unitary,psi,t)["probs"]
    probs_closed_form_corr = corrected_closed_form(unitary,psi,t)["probs"]
    probs_simulation = SW.simulation_symmetric_quantum_walk(unitary,t)["probs"][-1]
    print("Comparison Symmetric Quantum: ")
    print(probs_closed_form)
    print(probs_closed_form_corr)
    print(probs_simulation)

    # Classical
    t = 5
    corr = 0.9
    bias = 0.8
    probs_closed_form = closed_form_classical(corr,bias,t)["probs"]
    probs_simulation = SW.simulation_classically_correlated_walk(corr,bias,t)["probs"][-1]
    print("Classical: \n",probs_closed_form,"\n",probs_simulation)

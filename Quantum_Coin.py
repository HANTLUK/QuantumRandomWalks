import numpy as np
from qiskit import quantum_info as qi

import parametrization_walks as PW


def coin_correlation(U: np.ndarray, state_ind: int = 0, time: int = 20):
    """
    Simulates the correlation of the quantum coin for (0,1) or (1,0) m
    """
    psi_0 = np.array([1,0])
    psi_1 = np.array([0,1])
    states = [psi_0,psi_1]
    state = states[state_ind]
    prob_list = []
    corr_list = []
    prob_list.append(probs(state))
    for step in range(time):
        state = U @ state
        prob_vec = probs(state)
        prob_list.append(prob_vec)
        corr_list.append((prob_vec[state_ind] - prob_vec[state_ind-1]))
    return [prob_list,corr_list]

if __name__ == "__main__":
    for i in range(1):
        U = qi.random_unitary(2).data
        string = "\n".join([str(state) for state in symmetric_coin(U,10)])
        print(f"{string}")
        prob_coin, corr_coin = coin_correlation(U,0)
        string = "\n".join([str(state) for state in prob_coin])
        print(f"{string}")
        string = "\n".join([str(corr) for corr in corr_coin])
        print(f"{string}")
        prob_coin, corr_coin = coin_correlation(U,1)
        string = "\n".join([str(state) for state in prob_coin])
        print(f"{string}")
        string = "\n".join([str(corr) for corr in corr_coin])
        print(f"{string}")

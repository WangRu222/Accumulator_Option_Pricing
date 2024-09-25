import numpy as np
from scipy.sparse import diags
import time
import pandas as pd

def acc_cn(S0, X, H, r, sigma, settlement, NT, Ns, Mesht):
    """
    Python version of the MATLAB function `acc_cn`.
    
    Parameters:
    S0 (float): Initial price
    X (float): Strike price
    H (float): Barrier level
    r (float): Risk-free rate
    sigma (float): Volatility
    settlement (list): Array of settlement days
    NT (int): Total number of days per year
    Ns (int): Number of steps in S mesh
    Mesht (int): Number of time steps per day
    
    Returns:
    float: Price
    """
    
    
    Smax = 1.25 * H
    NT = settlement[-1]  # Last element
    AccValue = np.zeros(Ns - 1)
    
    start_time = time.time()
    
    if len(settlement) > 1:
        for i in range(len(settlement) - 1, 0, -1):
            AccValue = acc_cn_month(S0, X, H, r, sigma, NT, settlement[i], settlement[i] - settlement[i-1], AccValue, Smax, Ns, Mesht)
    
    
    AccValue = acc_cn_month(S0, X, H, r, sigma, NT, settlement[0], settlement[0], AccValue, Smax, Ns, Mesht)
    
    price = AccValue[int(np.floor(Ns * S0 / Smax))]
    
    print(f"Elapsed time: {time.time() - start_time} seconds")
    return price

def acc_cn_month(S0, X, H, r, sigma, NT, Nt, Nday, AccValue, Smax, Ns, Mesht):
    h = Smax / Ns
    dt = 1 / (NT * Mesht)
    indX = int(np.floor(X / h))
    indH = int(np.floor(H / h))
    S = np.arange(1, Ns) * h
    
    payoff = np.outer(S - X, np.arange(2 * Nday + 1)) + np.tile(AccValue, (2 * Nday + 1, 1)).T
    
    u = payoff
    
    
    u0 = np.outer((-X) * np.exp(-r * np.arange(Nday * Mesht, -1, -1) * dt) , np.arange(2 * Nday + 1)) + np.outer((-X) * np.exp(-r * np.arange(Nday * Mesht, -1, -1) * dt) * np.floor(np.arange(Nday * Mesht, -1, -1) / Mesht), np.ones(2 * Nday + 1))
    
    u1 = np.outer((H - X) * np.exp(0 * np.arange(Nday * Mesht, -1, -1) * dt) , np.arange(2 * Nday + 1))
    
    i = np.arange(1, Ns)
    alpha = sigma**2 / 2 * (i**2)
    beta = (r / 2) * i
    
    t1=alpha - beta
    t2=t1[1:]+t1[:1]
    t3=-alpha - beta
    t4=t3[1:]+t3[:1]

    B = diags([alpha + beta, -2 * alpha - r + 2 / dt, t2], [-1, 0, 1], shape=(Ns-1, Ns-1)).toarray().T
    A = diags([-alpha + beta, 2 * alpha + r + 2 / dt, t4], [-1, 0, 1], shape=(Ns-1, Ns-1)).toarray().T
    B1=B
    np.save("D:/vswork/B1.npy", B1)
    
    uu0 = (alpha[0] - beta[0]) * (u0[0:Nday*Mesht, :] + u0[1:Nday*Mesht + 1, :])
    uu1 = (alpha[-1] + beta[-1]) * (u1[0:Nday*Mesht, :] + u1[1:Nday*Mesht + 1, :])
    np.save("D:/vswork/t1.npy", t1)
    np.save("D:/vswork/t2.npy", t2)
    np.save("D:/vswork/t3.npy", t3)


    for n in range(1, Nday+1):
        u[:indX, (Nday - n):(2 * (Nday - n)+1)] = u[:indX, (Nday - n + 2):(2 * (Nday - n) + 3)]
        u[indX:indH, (Nday - n):(2 * (Nday - n)+1)] = u[indX:indH, (Nday - n + 1):(2 * (Nday - n) + 2)]
        u1=u
        np.save("D:/vswork/u1.npy", u1)
    
        if n > 1:
            u[indH:Ns-1, (Nday - n):(2 * (Nday - n))] = (S[indH:Ns-1] - X).reshape(-1, 1) @ np.arange(Nday - n, 2 * (Nday - n)).reshape(1, -1)
        for m in range(Mesht, 0, -1):
            w = B@u  
            w[0, :] += uu0[n * Mesht-m, :]
            w[Ns - 2, :] += uu1[n * Mesht - m, :]
            w1=uu1[n * Mesht - m, :]
            np.save("D:/vswork/w1.npy", w1)
            u = np.linalg.solve(A, w)
        
    return u[:, 0]



# The functions are now converted to Python.
# You can call `acc_cn` with appropriate parameters to test it.
S0 = 100  # Initial price
X =80  # Forward price (80% of initial price)
H = 120  # Knock-out price (120% of initial price)
r = 0  # Risk-free rate
sigma = 0.3  # Volatility
settlement = [20 ,39 ,62 ,80 ,101 ,122 ,142 ,164 ,187 ,208 ,229 ,250]  # Monthly settlement days
T = 252  # Number of trading days in 1 year
M = 100  # Number of simulations
ngrid = 32  # Observations per day
epsilon=0.1

price = acc_cn(S0, X, H, r, sigma, settlement, T, M, ngrid)
print(f"Estimated Price with Gearing Ratio: {price}")

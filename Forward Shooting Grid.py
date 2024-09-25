import numpy as np

def acc_FSGM(S0, X, H, r, sigma, settlement, T, mesh):
    settlement = np.array([0] + settlement)
    NT = T
    AccValue = np.zeros(2 * NT * mesh + 1)
    
    for i in range(len(settlement) - 1, 1, -1):
        AccValue = fsAccMonthPartmesh(S0, X, H, r, sigma, NT, settlement[i], 
                                       settlement[i] - settlement[i - 1], 
                                       AccValue, mesh)
    return AccValue

def fsAccMonthPartmesh(S0, X, H, r, sigma, NT, Nt, Nday, AccValue, mesh):
    dt = 1 / (NT * mesh)
    discount = np.exp(-r * dt)
    u = np.exp(sigma * np.sqrt(dt))
    d = 1 / u
    p = (np.exp(r * dt) - d) / (u - d)
    S = S0 * np.exp(sigma * np.sqrt(dt) * np.arange(-Nt * mesh, Nt * mesh + 1))
    payoff = (S - X)[:, None] @ (np.arange(0, 2 * Nday + 1)[None, :]) + AccValue[:, None]
    V = payoff.copy()
    Vtemp = payoff.copy()
    
    indS = Nt * mesh
    indj = 0
    
    for j in range(Nday - 1, -1, -1):
        # Last time step of each day, daily fixings
        for i in range((Nday - Nt - j - 1) * mesh + 1 + indS, (Nt - Nday + j + 1) * mesh - 1 + indS):
            if S[i - 1] > X:
                Vu = V[i + 1, j + indj:2 * j + 1 + indj]
                Vd = V[i - 1, j + indj:2 * j + 1 + indj]
                if S[i + 1] >= H:
                    Vu = (S[i] - X) * (np.arange(j + indj, 2 * j + indj + 1))
                if S[i - 1] > H:
                    Vd = (S[i] - X) * (np.arange(j + indj, 2 * j + indj + 1))
            elif S[i + 1] > X:
                Vu = V[i + 1, j + 1 + indj:2 * j + 1 + indj]
                Vd = V[i - 1, j + 2 + indj:2 * j + 2 + indj]
            else:
                Vu = V[i + 1, j + 2 + indj:2 * j + 2 + indj]
                Vd = V[i - 1, j + 2 + indj:2 * j + 2 + indj]
            
            Vtemp[i, j + indj:2 * j + indj + 1] = discount * (p * Vu + (1 - p) * Vd)

        V = Vtemp.copy()

        # Other time steps for each day, trivial binomial expectation
        if mesh > 1:
            for m in range(2, mesh + 1):
                for i in range((Nday - Nt - j - 1) * mesh + m + indS, (Nt - Nday + j + 1) * mesh - m + indS):
                    Vtemp[i, j + indj:2 * j + indj + 1] = discount * (p * V[i + 1, j + indj:2 * j + indj + 1] +
                                                                     (1 - p) * V[i - 1, j + indj:2 * j + indj + 1])
                V = Vtemp.copy()

    return V[(Nday - Nt) * mesh + indS:(Nt - Nday) * mesh + indS, indj].T





# Given assumptions
S0 = 100  # Initial price
X =80  # Forward price (80% of initial price)
H = 120  # Knock-out price (120% of initial price)
r = 0.02  # Risk-free rate
sigma = 0.3  # Volatility
settlement =  [20 ,39 ,62 ,80 ,101 ,122 ,142 ,164 ,187 ,208 ,229 ,250]  # Monthly settlement days
T = 252  # Number of trading days in 1 year
M = 1000  # Number of simulations
ngrid = 1  # Observations per day
epsilon = 0.1

price = acc_FSGM(S0, X, H, r, sigma, settlement, T, M)
print(f"Estimated Price with Gearing Ratio: {price}")

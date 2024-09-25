import numpy as np

def stockpath(S, H, r, sigma, T, ngrid):
    """
    Generates random stock paths and checks for knockout events.

    Parameters:
    S (float): Initial price
    H (float): Barrier level
    r (float): Risk-free rate
    sigma (float): Volatility
    T (int): Number of trading days per year
    ngrid (int): Number of barrier observations per day

    Returns:
    tuple: (stock path array, knockout day, discount factors array)
    """
    
    
    dt = 1 / (T * ngrid)  # Time increment
    discount = np.exp(-r * np.arange(1, T + 1) / T)  # Discount factor for each day
    path = np.zeros(T)  # Stock price path
    path[0] = S

    knockout = 0  # Initialize knockout day to 0

    for i in range(1, T):
        # Generate the next price in the path using geometric Brownian motion
        Z = np.random.randn()  # Random normal variable
        #path[i] = path[i-1] * np.exp((r - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * Z)
        path[i] = path[i-1] * np.exp((r - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt)*Z)

        # Check for knockout
        if path[i] >= H and knockout == 0:
            knockout = i

    return path, knockout, discount

def acc_mc(S, X, H, r, sigma, settlement, T, M, ngrid):
    """
    Monte Carlo simulation for SEMBCORP accumulator.


    Returns:
    tuple: (Estimated value from Monte Carlo simulation, Standard error of the estimate)
    """
    dt = 1 / T
    npath = 1000  # Number of paths per sample
    nsample = M // npath  # Number of samples
    prices = np.zeros(nsample)
    strikepath = np.full(T, X)

    for i in range(nsample):
        p = 0
        for j in range(npath):
            settle = [0] + settlement
            accumulation = np.ones(T)
            spath, knockout, discount = stockpath(S, H, r, sigma, T, ngrid)
            if knockout > 0:
                accumulation[knockout:T] = 0  # If knockout occurs, set accumulation to 0
                # Adjust settlement days if knockout occurs
                for ix in range(len(settle) - 1):
                    if settle[ix] >= knockout:
                        settle[ix] = knockout
                        settle = settle[:ix + 1]
                        break

            # Compute accumulation amounts for each day: {0, 1, 2}
            accumulation = (1 * (spath > strikepath) + 2 * (spath <= strikepath)) * accumulation

            # Sum up the PV of payoff for each settlement day
            for k in range(len(settle) - 1):
                p += np.sum(accumulation[(settle[k] + 1):settle[k + 1]]) * (spath[settle[k + 1]] - X) * discount[settle[k + 1]]

        prices[i] = p / npath

    value = np.mean(prices)  # Estimated value
    se = np.sqrt(np.var(prices) / nsample)  # Standard error of the estimate

    return value





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

price = acc_mc(S0, X, H, r, sigma, settlement, T, M, ngrid)
print(f"Estimated Price with Gearing Ratio: {price}")

# Accumulator-Option-Pricing
Compute the fair value of the accumulator using
- Forward Shooting Grid method;
- Crank-Nicolson Scheme Code;
- Monte Carlo simulation.
Based on the previous experiences by the students of an earlier class, the numerical value of
the accumulator is highly dependent on the choices of the time steps and stepwidth in the
numerical calculations.

## The option
The “accumulator” or “accumulative forward” is a daily accumulated and knock-out structured
product linked to the performance of an underlying asset. It can be considered as a portfolio
of forward contracts with the “occupation time” feature. The accumulated amount of assets
depends on the total excursion time of the asset price below the strike price. This leads to an
enhanced downside loss. The upside gain is limited by the knock-out feature with an upside
barrier.

- A typical equality-linked accumulator contract obligates an investor to buy a preset amount
of underlying stocks at the strike price X, if the closing stock price on a trading day is higher
than X. However, when the stock closes lower than X, the investor has to buy twice the
amount of stocks at X.
- Normally, the strike price X is set at a discount of the original spot
price S0. This explains why the accumulator is also called “discounted stock” among public.
- On the other hand, the profit from an accumulator contract is capped by an knock-out barrier
H which is set higher than S0.

## Replication of barrier-type derivatives
To derive the analytic formula based on discrete settlement of stock transaction on each business day, we assume continuous monitoring of the knock-out barrier $H$.
With immediate transaction of the stock, the payoff at date $t_i$ is given by

$$ 0  { if } \max_{0 \leq \tau \leq t_i} S_\tau \geq H $$

$$ S_{t_i}-X { if } \max_{0 \leq \tau \leq t_i} S_\tau < H { and } S_{t_i} \geq X  $$

$$ 2(S_{t_i}-X)   { if } \max_{0 \leq \tau t_i} S_\tau<H  { and } S_{t_i}<X $$

Assuming that there are $n$ business days, the fair value of the accumulator is given by

$$
V=\sum_{i=1}^n c_{u o}\left(t_i, X, H\right)-2 p_{u o}\left(t_i, X, H\right)
$$

where $c_{u o}$ and $p_{u o}$ denote the price function of the up-and-out call and up-and-out put, respectively.



### Monte Carlo simulation

Assume that the underlying asset price is known
distribution function of the grid, and then divide the validity period of the option into several
A small time interval, with the help of a computer, can be allocated from
Randomly sample from the sample to simulate changes in stock prices at each time interval
movement and stock price, so that we can calculate
Get the best end value. This result can be estimated as all
Possibly collect the final value in a random sample, using another of that sample
One path leads to another random sample. 

We assume the usual Black-Scholes framework where

$$
\frac{d S_t}{S_t}=r d t+\sigma d Z_t
$$


The stock prices at successive time step $t_j$ can be generated as

$$
S_j=S_{j-1} e^{\left(r-\frac{\sigma^2}{2}\right) \Delta t+\sigma \epsilon \sqrt{\Delta t}}
$$

where $\epsilon \sim N(0,1)$. We take $\Delta t=\frac{T}{N}$, where $N$ is the number of trading days until maturity date $T$.

For each trading day $i$, the following computation is implemented:
1. If $S_i \geq H, V_i=V_{i-1}+(1+K) e^{-i r \delta t}\left(S_i-e^{-3 r \delta t} X\right)$. Terminate the loop and return $V_N=V_i$; the accumulator is terminated.
2. If $X \leq S_i<H, K=K+1, V_i=V_{i-1}$; the number of stock sold is one.
3. If $S_i<X, K=K+2, V_i=V_{i-1}$; the number of stock sold is two.
4. If day $i$ is a settlement day, $V_i=V_{i-1}+(1+K) e^{-i r \delta t}\left(S_i-e^{-3 r \delta t} X\right), K=0$; after settlement, the count of stocks accumulated is set to be zero.
5. Go to next day.

Generate $M$ stock price paths do the above computations to obtain $\{V_{N, j}\}_{j=1}^M$. The price of the accumulator is:

$$
\bar{V}=\frac{1}{M} \sum_{j=1}^M V_{N, j}
$$


The standard error of the estimate is

$$
S E_{\bar{V}}=\sqrt{\frac{var(\{V_{N, j}\}_{j=1}^M)}{N}}
$$

where $var(\{V_{N, j}\}_{j=1}^M)$ is the variance of the sample result set. The standard error of the estimate is expected to be inversely proportional to the square root of the number of samples $M$.
### Forward Shooting Grid method

Finite difference methods mainly include intrinsic finite difference methods and
The basic idea of ​​the extrapolation finite difference method is to use numerical methods to
Solve the differential equation that a derivative asset satisfies to value a derivative asset
value, after transforming the differential equation into a series of difference equations,
Then solve these difference equations by iterative method.

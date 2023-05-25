import math
from scipy.stats import norm
import numpy as np
from GeneralFunctions import surface_plot


def blackScholesMertonCallPrice(time: float, underlying: float, strike: float, interest: float, sigma: float) -> float:
    if time == 0:
        time = 0.001
    if underlying == 0:
        underlying = 0.001
    if sigma == 0:
        sigma = 0.01
    d1 = (math.log(underlying / strike, math.e) + (interest + (sigma * sigma) / 2) * time) / (sigma * math.sqrt(time))
    d2 = d1 - sigma * math.sqrt(time)
    nd1 = norm.cdf(d1)
    nd2 = norm.cdf(d2)
    callPrice = underlying * nd1 - strike * math.exp(-interest * time) * nd2
    return callPrice


def matrixBlackScholes(t_grid: list, s_grid: list, strike: float, interest: float, sigma: float) -> list[list]:
    # grid = np.zeros(shape=(len(s_grid), len(t_grid)))

    grid = [[0 for j in range(len(t_grid))] for i in range(len(s_grid))]
    for i in range(len(s_grid)):
        for j in range(len(t_grid)):
            grid[i][j] = round(blackScholesMertonCallPrice(t_grid[j], s_grid[i], strike, interest, sigma), 10)
    print(len(grid))
    return grid


if __name__ == '__main__':
    np.set_printoptions(precision=3)
    # Parameters
    sigma = 0.2
    r = 0.03
    strike = 100
    T = 1

    # Create grid - array of different underlying price options
    nus = 10  # number of underlying price steps
    ds = 2 * strike / nus
    s_grid = [j * ds for j in range(nus)]
    print(*s_grid, sep=' ')

    # Create grid - array of different time steps
    nts = 10
    dt = float(T) / float(nts - 1)
    t_grid = [round(n * dt, 3) for n in range(nts)]
    print(*t_grid, sep=' ')

    # print(*matrixBlackScholes(t_grid, s_grid, strike, r, sigma), sep='\n')
    t_grid_exp = [1 - i for i in t_grid]
    print(t_grid_exp)
    surface_plot(matrixBlackScholes(t_grid, s_grid, strike, r, sigma), "Underlying price", "Time to expiry",
                 "Option price", t_grid_exp, s_grid, "Black-Scholes option price")

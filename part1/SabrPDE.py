import math
import numpy as np
from BlackScholesPDE import blackScholesMertonCallPrice
from GeneralFunctions import surface_plot


def zeta(vega, alpha, f, K, beta):
    t1 = vega / alpha
    t2 = pow(f * K, (1 - beta) * 0.5)
    t3 = math.log(f / K)
    return t1 * t2 * t3


def x_zeta_function(z, rho):
    denominator = 1 - rho
    numerator = math.sqrt(1 - 2 * rho * z + z ** 2) + z - rho
    res = math.log(numerator / denominator)
    return res


def term1(f, K, beta, alpha):
    denominator_term1 = pow(f * K, (1 - beta) * 0.5) * (
            1 + (pow(1 - beta, 2) / 24) * pow(math.log(f / K), 2) + (pow(1 - beta, 4) / 1920) * pow(math.log(f / K),
                                                                                                    4))
    t1 = alpha / denominator_term1
    return t1


def term2(z, rho):
    denominator = x_zeta_function(z, rho)
    if denominator == 0:
        denominator = 0.0001
    t2 = z / denominator
    return t2


def term3(alpha, beta, rho, vega, f, K, timeExp):
    t1 = (pow(1 - beta, 2) / 24) * (alpha ** 2 / pow(f * K, 1 - beta))
    t2 = 0.25 * (rho * beta * vega * alpha) / (pow(f * K, (1 - beta) * 0.5))
    t3 = 1 / 24 * (2 - 3 * rho) ** 2 * vega ** 2
    res = 1 + timeExp * (1 + t1 + t2 + t3)
    return res


def sabrImpliedVolatility(alpha, beta, rho, vega, strike, time, underlying_value):
    #     strike = int(strike)
    if underlying_value == 0:
        underlying_value = 1
    if time == 0:
        time = 0.05
    t1 = term1(underlying_value, strike, beta, alpha)
    t2 = term2(zeta(vega, alpha, underlying_value, strike, beta), rho)
    t3 = term3(alpha, beta, rho, vega, underlying_value, strike, time)
    implied_volatility = t1 * t2 * t3
    return implied_volatility


def matrix_SABR(t_grid, s_grid, strike, r, alpha, beta, rho, vega):
    matrix = [[0 for j in range(len(s_grid))] for i in range(len(t_grid))]
    # sigmaPrev = 0.01
    # matrix = list()
    for i in range(len(s_grid)):
        # row = np.array([0 for i in range(len(s_grid))])
        for j in range(len(t_grid)):
            sigma = sabrImpliedVolatility(alpha, beta, rho, vega, strike, t_grid[j], s_grid[i])
            # if sigma == 0:
            #     sigma = sigmaPrev
            # ans = blackScholesMertonCallPrice(t_grid[i], s_grid[j], strike, r, sigma)
            # if math.isnan(ans):
            #     ans = 0
            matrix[i][j] = sigma
            # sigmaPrev = sigma
            # row[j] = ans
        # matrix.append(row)
    # print(matrix, sep)
    return matrix


def matrix_SABR_BSM(t_grid, s_grid, strike, r, alpha, beta, rho, vega):
    matrix = [[0 for j in range(len(s_grid))] for i in range(len(t_grid))]
    # sigmaPrev = 0.01
    # matrix = list()
    for i in range(len(s_grid)):
        # row = np.array([0 for i in range(len(s_grid))])
        for j in range(len(t_grid)):
            sigma = sabrImpliedVolatility(alpha, beta, rho, vega, strike, t_grid[j], s_grid[i])
            # if sigma == 0:
            #     sigma = sigmaPrev
            ans = blackScholesMertonCallPrice(t_grid[j], s_grid[i], strike, r, sigma)
            # if math.isnan(ans):
            #     print(sigma)
            #     ans = 0
            matrix[i][j] = ans
            if matrix[i][j] == 0:
                print(sigma)
        # sigmaPrev = sigma
        # row[j] = ans
        # matrix.append(row)
    return matrix


if __name__ == '__main__':
    alpha = 0.18
    beta = 0.5002
    underlying = 101  # Current value of the underlying
    strike = 100  # Strike Price
    vega = 0.51
    rho = -0.02
    time = 1
    interest = 0.05

    nus = 10  # number of underlying price steps
    ds = 2 * strike / nus
    s_grid = [j * ds for j in range(nus)]
    print(*s_grid, sep=' ')

    # Create grid - array of different time steps
    nts = 10
    dt = float(time) / float(nts - 1)
    t_grid = [round(n * dt, 3) for n in range(nts)]
    print(*t_grid, sep=' ')

    print(sabrImpliedVolatility(alpha, beta, rho, vega, strike, time, underlying))
    # print(*matrix_SABR_BSM(t_grid, s_grid, strike, interest, alpha, beta, rho, vega), sep="\n")
    t_grid_exp = [1 - i for i in t_grid]
    print(t_grid_exp)
    surface_plot(matrix_SABR_BSM(t_grid, s_grid, strike, interest, alpha, beta, rho, vega), "Time to expiry",
                 "Underlying price",
                 "Option price", s_grid, t_grid_exp, "Black-Scholes combined with SABR option price")

    surface_plot(matrix_SABR(t_grid, s_grid, strike, interest, alpha, beta, rho, vega), "Time to expiry",
                 "Underlying price",
                 "Volatility", t_grid_exp, s_grid, "SABR volatility plot")
    # print(*matrix_SABR(t_grid, s_grid, strike, interest, alpha, beta, rho, vega), sep='\n')

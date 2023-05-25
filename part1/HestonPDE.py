import numpy as np
from GeneralFunctions import surface_plot

i = complex(0, 1)


def product(rho: float, sigma: float, phi: complex) -> complex:
    return rho * sigma * i * phi


def d_f(pr: complex, kappa: float, phi: complex, sigma: float) -> complex:
    term1 = (pr - kappa) ** 2
    term2 = (sigma ** 2) * (i * phi + phi ** 2)
    return np.sqrt(term1 + term2)


def g_f(kappa: float, pr: complex, d: complex) -> complex:
    numerator = kappa - pr - d
    denominator = kappa - pr + d
    return numerator / denominator


def charFterm1(phi: complex, g: complex, d: complex, sigma: float, kappa: float, theta: float, interest: float,
               time: float,
               underlying: float) -> float:
    p1 = underlying ** (i * phi) * np.exp(i * phi * interest * time)
    numerator = 1 - g * np.exp(-d * time)
    denominator = 1 - g
    return p1 * (numerator / denominator) ** (-2 * theta * kappa / (sigma ** 2))


def charFterm2(theta: float, kappa: float, time: float, sigma: float, pr: complex, d: complex, sigmasigma: float,
               g: complex) -> float:
    p1 = theta * kappa * time / (sigma ** 2)
    p2 = kappa - pr - d
    p3 = sigmasigma / (sigma ** 2)
    p4 = (1 - np.exp(-d * time)) / (1 - g * np.exp(-d * time))
    return np.exp((p1 * p2) + (p3 * p2 * p4))


def hestonCharacteristicFunction(phi: complex, underlying: float, strike: float, interest: float, time: float,
                                 sigma: float, kappa: float, theta: float, sigmasigma: float, rho: float) -> complex:
    pr = product(rho, sigma, phi)
    d = d_f(pr, kappa, phi, sigma)
    g = g_f(kappa, pr, d)
    term1 = charFterm1(phi, g, d, sigma, kappa, theta, interest, time, underlying)
    term2 = charFterm2(theta, kappa, time, sigma, pr, d, sigmasigma, g)
    return (term1 * term2)


def hestonPrice(underlying: float, strike: float, interest: float, time: float, sigma: float, kappa: float,
                theta: float, sigmasigma: float, rho: float) -> float:
    price, iter, maxN = 0, 1000, 100
    dphi = maxN / iter
    if underlying == 0:
        underlying = 0.001
    if time == 0:
        time = 0.001

    # Calculate the complex integral
    # Using j instead of i to avoid confusion
    for j in range(1, iter):
        phi1 = dphi * (2 * j + 1) / 2
        phi2 = phi1 - i

        numerator1 = hestonCharacteristicFunction(phi2, underlying, strike, interest, time, sigma, kappa, theta,
                                                  sigmasigma, rho)
        numerator2 = strike * hestonCharacteristicFunction(phi1, underlying, strike, interest, time, sigma, kappa,
                                                           theta, sigmasigma, rho)
        denominator = np.exp(np.log(strike) * i * phi1) * i * phi1

        price = price + dphi * (numerator1 - numerator2) / denominator

    term1 = 0.5 * (underlying - strike * np.exp(-interest * time))
    term2 = price / np.pi

    return np.real((term1 + term2))


def matrixHeston(t_grid: list, s_grid: list, strike: float, interest: float, sigma: float, kappa: float, theta: float,
                 sigmasigma: float, rho: float) -> list[list]:
    # grid = np.zeros(shape=(len(s_grid), len(t_grid)))

    grid = [[0 for j in range(len(t_grid))] for i in range(len(s_grid))]
    for i in range(len(s_grid)):
        for j in range(len(t_grid)):
            grid[i][j] = round(
                hestonPrice(s_grid[i], strike, interest, t_grid[j], sigma, kappa, theta, sigmasigma, rho), 2)
            if (grid[i][j] < 0):
                grid[i][j] = 0
    return grid


if __name__ == '__main__':
    underlying = 100.  # asset price
    strike = K = 100.  # strike
    sigma = v0 = 0.1  # initial volatility
    interest = r = 0.03  # risk free rate
    kappa = 1.5768  # rate of mean reversion of variance process
    theta = 0.0398  # long-term mean variance
    sigmasigma = 0.3  # volatility of volatility
    lambd = 0.575  # risk premium of variance
    rho = 0.5711  # correlation between variance and stock process
    time = tau = 1  # time to maturity

    nus = 10  # number of underlying price steps
    ds = 2 * strike / nus
    s_grid = [j * ds for j in range(nus)]
    print(*s_grid, sep=' ')

    # Create grid - array of different time steps
    nts = 10
    dt = float(time) / float(nts - 1)
    t_grid = [round(n * dt, 3) for n in range(nts)]
    print(*t_grid, sep=' ')

    call_price = hestonPrice(underlying, strike, time, interest, sigma, kappa, rho, theta, sigmasigma)
    print(call_price)

    # print(*matrixHeston(t_grid, s_grid, strike, r, sigma, kappa, theta, sigmasigma, rho), sep='\n')
    t_grid_exp = [1 - i for i in t_grid]
    print(t_grid_exp)
    surface_plot(matrixHeston(t_grid, s_grid, strike, r, sigma, kappa, theta, sigmasigma, rho), "Time to expiry",
                 "Underlying price",
                 "Option price", t_grid_exp, s_grid, "Heston option price")

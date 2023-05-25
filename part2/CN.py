from EFD import setBoundary
from GeneralFunctions import surface_plot
import numpy as np
import scipy.linalg as linalg
import time
from part1.BlackScholesPDE import matrixBlackScholes
from GeneralFunctions import calcErrorPlot


def setCoefficients(interest, dt, i_values, sigma, M):
    a = [0.5 * (interest * dt * i_values[i] - sigma ** 2 * dt * i_values[i] ** 2) for i in range(len(i_values))]
    b = [1 + sigma ** 2 * dt * i_values[i] ** 2 + interest * dt for i in range(len(i_values))]
    c = [-0.5 * (interest * dt * i_values[i] + sigma ** 2 * dt * i_values[i] ** 2) for i in range(len(i_values))]
    coeffs = np.diag(a[2:M], -1) + \
             np.diag(b[1:M]) + \
             np.diag(c[1:M - 1], 1)
    return coeffs, a


def setMatrix(coeffs, M, N, grid, a):
    P, L, U = linalg.lu(coeffs)
    aux = [0 for i in range(M - 1)]

    for j in reversed(range(N)):
        aux[0] = np.dot(-a[1], grid[-1, j])
        x1 = linalg.solve(L, grid[1:M, j + 1] + aux)
        x2 = linalg.solve(U, x1)
        grid[1:M, j] = x2
    return grid


def crankNicholson(strike, interest, time, sigma, maxUnderlying, N, M, dt, ds):
    underlyingS = nus = M
    timeS = nts = N

    i_values = [ds for i in range(M)]
    j_values = [dt for i in range(N)]
    matrix = np.zeros(shape=(M + 1, N + 1))
    print(len(matrix))

    boundaryConds = [30 * ds * i for i in range(M + 1)]
    setBoundary(matrix, boundaryConds, strike, maxUnderlying, interest, dt, N, j_values)
    coeffs, a = setCoefficients(interest, dt, i_values, sigma, M)
    matrix = setMatrix(coeffs, M, N, matrix, a)
    return matrix, i_values, j_values


if __name__ == '__main__':
    # Parameters
    strike = 50;
    interest = 0.05;
    mytime = 1.;
    sigma = 0.1;
    maxUnderlying = 200
    M = 50
    N = 100
    ds = maxUnderlying / M
    dt = mytime / N

    # Plotting CN
    CN_matrix, i_val, j_val = crankNicholson(strike, interest, mytime, sigma, maxUnderlying, M, N, dt, ds)
    u_val = [i * ds for i in range(len(CN_matrix))]
    t_val = [i * dt for i in range(len(CN_matrix[0]))]

    # Recording time taken to complete plot
    st = time.time()
    surface_plot(CN_matrix, "time", "underlying", "option", t_val, u_val, "Crank-Nicholson")
    et = time.time()

    # Calculate error produced
    BSM_matrix = matrixBlackScholes(t_val, u_val, strike, interest, sigma)
    print(len(BSM_matrix))
    print(len(CN_matrix))
    surface_plot(BSM_matrix, "Time (t)", "Asset Price (S)", "Price of option (V)", t_val, u_val,
                 "BSM model of an option price")
    diff = calcErrorPlot(BSM_matrix, CN_matrix, t_val, u_val, "Error")
    print("Total time taken: ", et - st)
    print("Total difference is : {}".format(diff))

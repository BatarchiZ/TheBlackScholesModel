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


# def setMatrix(coeffs, M, N, grid, a):
#     P, L, U = linalg.lu(coeffs)
#     aux = [0 for i in range(M - 1)]
#
#     for j in reversed(range(N)):
#         aux[0] = np.dot(-a[1], grid[-1,j])
#         x1 = linalg.solve(L, grid[1:M, j + 1] + aux)
#         x2 = linalg.solve(U, x1)
#         grid[1:M, j] = x2
#     return grid
# def truncMatrix(array):
#     for i in range(len(array)):
#         for j in range(len(array[i])):
#             array[i][j] = round(array[i][j],2)
#     return array
def setMatrix(coeffs, M, N, grid, a):
    for j in reversed(range(N)):
        x1 = linalg.solve(coeffs, grid[1:M, j + 1])
        grid[1:M, j] = x1
        # grid = truncMatrix(grid)
        # print(*grid, sep= '\n')
        # print()
    return grid


def IFD(strike, interest, time, sigma, maxUnderlying, M, N):
    underlyingS = nus = M
    timeS = nts = N

    i_values = [i for i in range(M)]
    j_values = [i for i in range(N)]
    matrix = np.zeros(shape=(M + 1, N + 1))
    ds = maxUnderlying / M
    dt = time / N

    boundaryConds = [ds * i for i in range(M + 1)]  # Poidet v konex i sverhu nujno ubrat sdes; strike
    # print(boundaryConds)

    setBoundary(matrix, boundaryConds, strike, maxUnderlying, interest, dt, N, j_values)
    # print(*matrix, sep = '\n')
    coeffs, a = setCoefficients(interest, dt, i_values, sigma, M)
    matrix = setMatrix(coeffs, M, N, matrix, a)
    return matrix, i_values, j_values


if __name__ == '__main__':
    # Parameters
    strike = 10;
    interest = 0.05;
    mytime = 0.1;  # Stranno chto ona shoditsya k dt esli umen'shat' no esli menyat cherez n to vse ploho, I wonder why... Maybe smthng to do with corner
    sigma = 0.1;
    maxUnderlying = 200
    M = 50
    N = 25
    ds = maxUnderlying / M
    dt = mytime / N

    # Plotting CN
    IFD_matrix, i_val, j_val = IFD(strike, interest, mytime, sigma, maxUnderlying, M, N)
    u_val = [i * ds for i in range(len(IFD_matrix))]
    t_val = [i * dt for i in range(len(IFD_matrix[0]))]

    # Recording time taken to complete plot
    st = time.time()
    surface_plot(IFD_matrix, "time", "underlying", "option", t_val, u_val, "IFD")
    et = time.time()

    # Calculate error produced
    BSM_matrix = matrixBlackScholes(t_val, u_val, strike, interest, sigma)
    print(len(BSM_matrix))
    print(len(IFD_matrix))
    surface_plot(BSM_matrix, "Time (t)", "Asset Price (S)", "Price of option (V)", t_val, u_val,
                 "BSM model of an option price")
    diff = calcErrorPlot(BSM_matrix, IFD_matrix, t_val, u_val, "Error")
    print("Total time taken: ", et - st)
    print("Total difference is : {}".format(diff))

import numpy as np
from GeneralFunctions import surface_plot


def setBoundary(grid, boundaryConds, strike, maxUnderlying, interest, dt, N, j_values):
    for i in range(len(grid)):
        grid[i][-1] = np.maximum(0, boundaryConds[i] - strike)
    print(*grid, sep='\n')
    print()
    for j in range(len(grid[0]) - 1):
        grid[-1][j] = (maxUnderlying - strike) * np.exp(-interest * dt * (N - j_values[j]))
    print(*grid, sep='\n')
    return grid


def setCoefficients(dt, sigma, i_values, interest):
    a = [0.5 * dt * (sigma ** 2 * i_values[i] ** 2 - interest * i_values[i]) for i in range(len(i_values))]
    b = [1 - dt * (sigma ** 2 * i_values[i] ** 2 + interest) for i in range(len(i_values))]
    c = [0.5 * dt * (sigma ** 2 * i_values[i] ** 2 + interest * i_values[i]) for i in range(len(i_values))]
    return a, b, c


def setMatrix(a, b, c, j_values, M, grid):
    for j in reversed(j_values):
        for i in range(M)[2:]:
            grid[i][j] = a[i] * grid[i - 1][j + 1] + b[i] * grid[i][j + 1] + c[i] * grid[i + 1][j + 1]
    return grid


def EFD(underlying, strike, interest, time, sigma, maxUnderlying, M, N):
    underlyingS = nus = M
    timeS = nts = N
    i_values = [i for i in range(M)]
    j_values = [i for i in range(N)]
    matrix = [[0 for j in range(N + 1)] for i in range(M + 1)]
    ds = maxUnderlying / M
    dt = time / N
    boundaryConds = [ds * i for i in range(M + 1)]

    setBoundary(matrix, boundaryConds, strike, maxUnderlying, interest, dt, N, j_values)
    a, b, c = setCoefficients(dt, sigma, i_values, interest)
    matrix = setMatrix(a, b, c, j_values, M, matrix)
    return matrix, i_values, j_values


if __name__ == '__main__':
    # # Stable
    matrixStable, i_val, j_val = EFD(50, 10, 0.1, 5. / 12., 0.4, 100, 5, 5)

    # # Unstable

    matrixUnstable, i2_val, j2_val = EFD(50, 10, 0.1, 5. / 12., 0.4, 100, 100, 10)

    u_val = [i for i in range(len(matrixStable))]
    t_val = [i for i in range(len(matrixStable[0]))]
    surface_plot(matrixUnstable, "time", "underlying", "option", t_val, u_val, "EFD")
    surface_plot(matrixStable, "time", "underlying", "option", t_val, u_val, "EFD")

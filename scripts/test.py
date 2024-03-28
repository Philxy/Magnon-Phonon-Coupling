import numpy as np


N = 10
for nk in range(0, N):
    for nq in range(0, N):
        for nkp in range(0, N):

            m = (2*(nk+nkp+nq) - 3*N)/(2*N)

            # print(m)
            # nkp = (2*N*m+3*N)/(2) - nk - nq


for nk in range(0, N):
    for nq in range(0, N):
        print('\n')

        for m in range(-1, 2):
            nkp = (2 * N * m + 3 * N) / 2 - nk - nq
            print(m, nkp)
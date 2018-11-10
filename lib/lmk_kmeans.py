from sklearn.cluster import KMeans
import scipy.sparse.linalg as sla
import mosek
import numpy as np
import time
import mosek
import sys
from mosek.fusion import *


def lmkkmeans_train(Km, iteration_count=2, cluster_count=10):
    '''
    Localized Multiple Kernel k-Means Clustering
    From: https://github.com/mehmetgonen/lmkkmeans
    Parameters
    ----------
    Km : numpy.array
        `Km` is a similarity matrix of PxNxN, N = # of items and P = # of kernels
    iteration_count : int
        `iteration_count` sets how many iterations will be run
    cluster_count : int
        `cluster_count` sets how many clusters will be produced by kmeans
    Returns
    -------
    dict
        dictionary with following
        clustering: K-Means clustering result
        objective: objective matrix
        parameters: run parameters
        Theta: resultant theta
        elapsed: elapsed time for calculation
    '''

    start_time = time.time()
    print(Km.shape)
    N = Km.shape[1]
    P = Km.shape[0]
    if Km.shape[1] != Km.shape[2]:
        raise ValueError('Km dimensions must be of form PxNxN')

    Theta = np.ones((N, P)) / P
    K_Theta = calculate_localized_kernel_theta(Km, Theta)

    objective = np.zeros(iteration_count)
    for i in range(iteration_count):
        print('Running iteration %2d...' % i)
        H, _ = sla.eigs(K_Theta, cluster_count, which='LR')
        HHT = H @ H.transpose()

        Q = np.zeros((N * P, N * P))
        for m in range(P):
            st = m * N
            ed = st + N
            Q[st:ed, st:ed] = np.eye(N, N) * Km[m, :, :] - HHT * Km[m, :, :]

        res = mosek(Q, np.zeros(N * P), np.tile(np.eye(N, N), (1, P)), np.ones((N, 1)), np.ones((N, 1)), np.zeros((N * P, 1)), np.ones((N * P, 1)), [], 'minimize echo(0)')
        Theta = np.reshape(res.sol.itr.xx, N, P)
        K_Theta = calculate_localized_kernel_theta(Km, Theta)

        objective[i] = np.trace(H.transpose() * K_Theta * H) - np.trace(K_Theta)

    H_normalized = H / np.tile(sqrt(np.sum(H*H, 1)), (1, cluster_count))

    clustering = KMeans(n_clusters=cluster_count, max_iter=1000, n_jobs=10).fit(H_normalized)
    print('Calculation took %d ms', time.time() - start_time)
    return clustering

def calculate_localized_kernel_theta(Km, Theta):
    K_Theta = np.zeros((Km.shape[1], Km.shape[2]))
    for m in range(Km.shape[0]):
        K_Theta = K_Theta + (Theta[:, m] @ Theta[:, m].transpose()) * Km[m, :, :]
    return K_Theta

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

if __name__ == '__main__':
    inf = 0.0 # symbollic value since it will be ignored
    with mosek.Env() as env:
        # Attach a printer to the environment
        env.set_Stream(mosek.streamtype.log, streamprinter)
        # Create a task
        with env.Task() as task:
            task.set_Stream(mosek.streamtype.log, streamprinter)
            # Set up and input bounds and linear coefficients
            bkc = [mosek.boundkey.lo]
            blc = [1.0]
            buc = [inf]
            numvar = 3
            bkx = [mosek.boundkey.lo] * numvar
            blx = [0.0] * numvar
            bux = [inf] * numvar
            c = [0.0, -1.0, 0.0]
            asub = [[0], [0], [0]]
            aval = [[1.0], [1.0], [1.0]]

            numvar = len(bkx)
            numcon = len(bkc)

            # Append 'numcon' empty constraints.
            # The constraints will initially have no bounds.
            task.appendcons(numcon)

            # Append 'numvar' variables.
            # The variables will initially be fixed at zero (x=0).
            task.appendvars(numvar)

            for j in range(numvar):
                # Set the linear term c_j in the objective.
                task.putcj(j, c[j])
                # Set the bounds on variable j
                # blx[j] <= x_j <= bux[j]
                task.putbound(mosek.accmode.var, j, bkx[j], blx[j], bux[j])
                # Input column j of A
                task.putacol(j,                  # Variable (column) index.
                             # Row index of non-zeros in column j.
                             asub[j],
                             aval[j])            # Non-zero Values of column j.
            for i in range(numcon):
                task.putbound(mosek.accmode.con, i, bkc[i], blc[i], buc[i])

            # Set up and input quadratic objective
            qsubi = [0, 1, 2, 2]
            qsubj = [0, 1, 0, 2]
            qval = [2.0, 0.2, -1.0, 2.0]

            task.putqobj(qsubi, qsubj, qval)

            # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)

            # Optimize
            task.optimize()
            # Print a summary containing information
            # about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)

            prosta = task.getprosta(mosek.soltype.itr)
            solsta = task.getsolsta(mosek.soltype.itr)

            # Output a solution
            xx = [0.] * numvar
            task.getxx(mosek.soltype.itr, xx)

            if solsta == mosek.solsta.optimal or solsta == mosek.solsta.near_optimal:
                print("Optimal solution: %s" % xx)
            elif solsta == mosek.solsta.dual_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif solsta == mosek.solsta.prim_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif solsta == mosek.solsta.near_dual_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif solsta == mosek.solsta.near_prim_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif mosek.solsta.unknown:
                print("Unknown solution status")
            else:
                print("Other solution status")
    # lmkkmeans_train(np.array([[[1, 2], [2, 3]], [[1, 3], [2, 3]], [[1, 1], [0, 0 ]]]))

from .. import config
import heapq
import mosek
import sys

import numpy as np
import numpy.matlib as npML
from sklearn.cluster import KMeans


def lmkkmeans_train(Km, iteration_count=2, cluster_count=10):
    if iteration_count < 1:
        raise ValueError(f'iteracation_count cannot be non positive give: {iteration_count}')
    N = Km.shape[1]
    P = Km.shape[0]
    Theta = np.zeros((N, P))
    for i in range(N):
        for j in range(P):
            Theta[i][j] = 1 / P
    K_Theta = calculate_localized_kernel_theta(Km, Theta)
    objective = np.zeros((iteration_count, 1))

    for it in range(iteration_count):
        eig, H = np.linalg.eig(K_Theta)
        maxIndexes = heapq.nlargest(cluster_count, range(eig.shape[0]), eig.take)
        H = H[:, maxIndexes]
        H_con = np.conjugate(H.transpose())
        HHT = H @ H_con

        qsubi = []
        qsubj = []
        qval = []
        for m in range(P):
            start_index = m * N
            tmpMatrix = np.eye(N) * Km[m] - HHT * Km[m]
            for i in range(len(tmpMatrix)):
                for j in range(i + 1):
                    if tmpMatrix[i, j] != 0:
                        qsubi.append(start_index + i)
                        qsubj.append(start_index + j)
                        qval.append(tmpMatrix[i, j])
        resxx = call_mosek(qsubi, qsubj, qval, N, P)
        Theta = np.array(resxx).reshape((P, N)).T
        K_Theta = calculate_localized_kernel_theta(Km, Theta)
        objective[it] = np.trace(H_con @ K_Theta @ H) - np.trace(K_Theta)
        print()

    tempH = (npML.repmat(np.sqrt((H ** 2).sum(axis=1)), cluster_count, 1)).transpose()  #
    H_normalized = np.divide(H, tempH, out=np.zeros_like(H), where=tempH != 0)
    clustering = KMeans(n_clusters=cluster_count, max_iter=1000).fit_predict(H_normalized)
    return clustering, H_normalized


def calculate_localized_kernel_theta(K, Theta):
    N = K.shape[1]
    K_Theta = np.zeros((N, N))
    for i in range(K.shape[0]):
        K_Theta += Theta[:, i] * Theta[:, [i]] * K[i]
    return K_Theta


def call_mosek(qsubi, qsubj, qval, N, P):
    with mosek.Env() as env:
        # Attach a printer to the environment
        env.set_Stream(mosek.streamtype.log, streamprinter)
        # Create a task
        with env.Task() as task:
            task.set_Stream(mosek.streamtype.log, streamprinter)
            # Set up and input bounds and linear coefficients
            bkc = [mosek.boundkey.ra] * N
            blc = np.ones((N, 1))
            buc = np.ones((N, 1))
            bkx = [mosek.boundkey.ra] * N * P
            blx = np.zeros((N * P, 1))
            bux = np.ones((N * P, 1))
            c = np.zeros((N * P, 1))
            aval = np.ones((N * P, 1)).astype(int)
            asub = np.zeros((N * P, 1)).astype(int)
            for i in range(P):
                for j in range(N):
                    asub[j + N * i] = j

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
                task.putvarbound(j, bkx[j], blx[j], bux[j])
                # Input column j of A
                task.putacol(j,  # Variable (column) index.
                             # Row index of non-zeros in column j.
                             asub[j],
                             aval[j])  # Non-zero Values of column j.
            for i in range(numcon):
                task.putconbound(i, bkc[i], blc[i], buc[i])
            # Set up and input quadratic objective
            task.putqobj(qsubi, qsubj, qval)
            # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)
            # Optimize
            task.optimize()
            # Print a summary containing information
            # about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)
            # prosta = task.getprosta(mosek.soltype.itr)
            # solsta = task.getsolsta(mosek.soltype.itr)
            # Output a solution
            xx = [0.] * numvar
            task.getxx(mosek.soltype.itr, xx)
            return xx


def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()


def main():
    try:
        print(config.MOSEK_LICENCE_FILE_PATH)
        v1, v2, v3 = views.getLinearKernel()
        Kmm = np.stack((v1, v2, v3))

        lmkkmeans_train(Kmm, cluster_count=3, iteration_count=10)
    except mosek.MosekException as e:
        print("ERROR: %s" % str(e.errno))
        if e.msg is not None:
            import traceback

            traceback.print_exc()
            print("\t%s" % e.msg)
        sys.exit(1)


if __name__ == '__main__':
    main()

import sys, mosek
import numpy as np
import heapq
from sklearn.cluster import KMeans


def lmkkmeans_train(Km, iteration_count=2, cluster_count=10):
    N = Km.shape[1]
    P = Km.shape[0]
    Theta = np.zeros((N, P))
    for i in range(N):
        for j in range(P):
            Theta[i][j] = 1 / P
    K_Theta = calculate_localized_kernel_theta(Km, Theta)
    objective = np.zeros((iteration_count, 1))

    for iter in range(iteration_count):
        # print("running iteration "+str(iter)+"...\n" )
        temp, H = np.linalg.eig(K_Theta)
        maxIndexes = heapq.nlargest(cluster_count, range(len(temp)), temp.take)
        H = H[:, maxIndexes]
        HHT = np.matmul(H, np.conjugate(H.transpose()))

        Q = np.zeros((N * P, N * P))
        for m in range(P):
            start_index = m * N
            end_index = (m + 1) * N
            Q[start_index:end_index, start_index:end_index] = np.eye(N) * Km[m] - HHT * Km[m]
        resxx = call_mosek(Q)
        Theta = np.array(resxx).reshape((P, N)).T
        K_Theta = calculate_localized_kernel_theta(Km, Theta)
        np.matmul(np.conjugate(H.transpose()), K_Theta)
        objective[iter] = np.trace(np.matmul(np.matmul(np.conjugate(H.transpose()), K_Theta), H)) - np.trace(K_Theta)
        print()

    tempH = (np.matlib.repmat(np.sqrt((np.power(H, 2)).sum(axis=1)), cluster_count, 1)).transpose()
    H_normalized = np.divide(H, tempH)
    clustering = KMeans(n_clusters=cluster_count, max_iter=1000, ).fit(H_normalized)
    return clustering, objective, Theta


def calculate_localized_kernel_theta(K, Theta):
    N = len(K[0])
    P = len(K)
    K_Theta = np.zeros((N, N))
    for i in range(0, P):
        K_Theta = K_Theta + (Theta[:, i] * np.array([Theta[:, i]]).transpose()) * (K[i])
    return K_Theta


def call_mosek(Q):
    N = 120
    P = 3
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
            numvar = 3
            bkx = [mosek.boundkey.ra] * N * P
            blx = np.zeros((N * P, 1))
            bux = np.ones((N * P, 1))
            c = np.zeros((N * P, 1))
            A = np.eye(N)
            aval = np.ones((N * P, 1)).astype(int)
            asub = np.zeros((N * P, 1)).astype(int)
            for i in range(P):
                for j in range(N):
                    asub[j + N * i] = j

            qsubi = []
            qsubj = []
            qval = []
            for i in range(len(Q[1, :])):
                for j in range(i + 1):
                    if (Q[i, j].real) != 0:
                        qsubi.append(i)
                        qsubj.append(j)
                        qval.append(Q[i, j])

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
                task.putacol(j,  # Variable (column) index.
                             # Row index of non-zeros in column j.
                             asub[j],
                             aval[j])  # Non-zero Values of column j.
            for i in range(numcon):
                task.putbound(mosek.accmode.con, i, bkc[i], blc[i], buc[i])
            # Set up and input quadratic objective
            task.putqobj(qsubi, qsubj, qval)
            # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)
            # Optimize
            task.optimize()
            # Print a summary containing information
            # about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)
            #prosta = task.getprosta(mosek.soltype.itr)
            #solsta = task.getsolsta(mosek.soltype.itr)
            # Output a solution
            xx = [0.] * numvar
            task.getxx(mosek.soltype.itr, xx)
            return xx
            print("Optimal solution: %s" % xx)


def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()


if __name__ == '__main__':
    try:
        main()
    except mosek.MosekException as e:
        print("ERROR: %s" % str(e.errno))
        if e.msg is not None:
            import traceback
            traceback.print_exc()
            print("\t%s" % e.msg)
        sys.exit(1)

def main():
    v1, v2, v3 = views.getLinearKernel()
    Kmm = np.stack((v1, v2, v3))

    lmkkmeans_train(Kmm, cluster_count=3, iteration_count=10)

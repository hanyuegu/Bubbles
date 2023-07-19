from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os
import multiprocessing as mp

def solver(N):
    # ---------------------------------
    # Define subdomains
    # ---------------------------------
    fname = "data_stellar/{}".format(N)
    mesh = Mesh(fname + ".xml")
    subdomains = MeshFunction("size_t", mesh, fname + "_physical_region.xml")
    # File("subdomains.pvd") << subdomains
    # print(subdomains.array())
    # for cell in cells(mesh):
    #     print(cell.midpoint()[0],cell.midpoint()[1])

    # ---------------------------------
    # fill cell value in k
    # ---------------------------------
    k_values = [2,8]
    V0 = FunctionSpace(mesh,"DG",0)
    k = Function(V0)
    kvalues = k_values # values of k in the two subdomains
    for cell in range(len(subdomains.array())):
        subdomain_cell = subdomains.array()[cell]
        index = int(subdomain_cell-1)
        k.vector()[cell] = kvalues[index]
    # print('k degree of freedoms:', k.vector().get_local())

    # ---------------------------------
    # class used to define the periodic boundary map
    # ---------------------------------
    vertices = np.array([[0, 0],[2, 0],[2, 2],[0, 2]])

    class PeriodicBoundary(SubDomain):
        def __init__(self, vertices, tolerance=DOLFIN_EPS):
            """ vertices stores the coordinates of the 4 unit cell corners"""
            SubDomain.__init__(self, tolerance)
            self.tol = tolerance
            self.vv = vertices
            self.a1 = self.vv[1,:]-self.vv[0,:] # first vector generating periodicity
            self.a2 = self.vv[3,:]-self.vv[0,:] # second vector generating periodicity
            # check if UC vertices form indeed a parallelogram
            assert np.linalg.norm(self.vv[2, :]-self.vv[3, :] - self.a1) <= self.tol
            assert np.linalg.norm(self.vv[2, :]-self.vv[1, :] - self.a2) <= self.tol

        def inside(self, x, on_boundary):
            # return True if on left or bottom boundary AND NOT on one of the
            # bottom-right or top-left vertices
            return bool((near(x[0], self.vv[0,0] + x[1]*self.a2[0]/self.vv[3,1], self.tol) or
                        near(x[1], self.vv[0,1] + x[0]*self.a1[1]/self.vv[1,0], self.tol)) and
                        (not ((near(x[0], self.vv[1,0], self.tol) and near(x[1], self.vv[1,1], self.tol)) or
                        (near(x[0], self.vv[3,0], self.tol) and near(x[1], self.vv[3,1], self.tol)))) and on_boundary)

        def map(self, x, y):
            if near(x[0], self.vv[2,0], self.tol) and near(x[1], self.vv[2,1], self.tol): # if on top-right corner
                y[0] = x[0] - (self.a1[0]+self.a2[0])
                y[1] = x[1] - (self.a1[1]+self.a2[1])
            elif near(x[0], self.vv[1,0] + x[1]*self.a2[0]/self.vv[2,1], self.tol): # if on right boundary
                y[0] = x[0] - self.a1[0]
                y[1] = x[1] - self.a1[1]
            else:   # should be on top boundary
                y[0] = x[0] - self.a2[0]
                y[1] = x[1] - self.a2[1]

    # ---------------------------------
    # define corrector equation
    # ---------------------------------
    Ve = FiniteElement('Lagrange', mesh.ufl_cell(), 2)
    Re = FiniteElement('Real', mesh.ufl_cell(), 0)
    X = FunctionSpace(mesh, MixedElement([Ve, Re]), constrained_domain = PeriodicBoundary(vertices, tolerance=1e-10))
    u, m = TrialFunctions(X)
    v, r = TestFunctions(X)
    F1 = inner(nabla_grad(u) + Constant((1,0)), k*nabla_grad(v))*dx + u*r*dx + m*v*dx
    F2 = inner(nabla_grad(u) + Constant((0,1)), k*nabla_grad(v))*dx + u*r*dx + m*v*dx

    # ---------------------------------
    # Compute first order corrector
    # ---------------------------------
    x = Function(X)
    a1 = lhs(F1)
    L1 = rhs(F1)
    solve(a1 == L1, x)
    chi_x, m_x = x.split()    # chi_1
    # File("chi_x.pvd") << chi_x

    y = Function(X)
    a2 = lhs(F2)
    L2 = rhs(F2)
    solve(a2 == L2, y)
    chi_y, m_y = y.split()
    # File("chi_y.pvd") << chi_y

    A_bar_11 = assemble((k*nabla_grad(chi_x)[0] + k)*dx)/4
    A_bar_21 = assemble((k*nabla_grad(chi_x)[1])*dx)/4
    A_bar_12 = assemble((k*nabla_grad(chi_y)[0])*dx)/4
    A_bar_22 = assemble((k*nabla_grad(chi_y)[1] + k)*dx)/4
    A_bar = np.array([[A_bar_11, A_bar_12], [A_bar_21, A_bar_22]])
    return A_bar

if __name__ == '__main__':

    cpus = mp.cpu_count()

    from time import time
    t1 = time()
    with mp.Pool(processes=50) as pool:
        A_bar = pool.map(solver, np.arange(1,1001,1))
        A_bar = np.array(A_bar)

    np.save("A_bar_1000.npy", A_bar)
    print(f"Time elapsed: {time() - t1}")

    A_bar_11 = A_bar[:,0,0]
    A_bar_12 = A_bar[:,0,1]
    A_bar_21 = A_bar[:,1,0]
    A_bar_22 = A_bar[:,1,1]

    plt.subplot(2,2,1)
    plt.hist(A_bar_11,label="The mean of A_bar_11: %f\n" %(np.array(A_bar_11).mean()) +
             "The variance of A_bar_11: %f" %(np.array(A_bar_11).var()))
    plt.legend()
    plt.title("results of A_bar_11")

    plt.subplot(2,2,2)
    plt.hist(A_bar_12,label="The mean of A_bar_12: %f\n" %(np.array(A_bar_12).mean()) +
             "The variance of A_bar_12: %f" %(np.array(A_bar_12).var()))
    plt.legend()
    plt.title("results of A_bar_12")

    plt.subplot(2,2,3)
    plt.hist(A_bar_21,label="The mean of A_bar_21: %f\n" %(np.array(A_bar_21).mean()) +
             "The variance of A_bar_21: %f" %(np.array(A_bar_21).var()))
    plt.legend()
    plt.title("results of A_bar_21")

    plt.subplot(2,2,4)
    plt.hist(A_bar_22,label="The mean of A_bar_22: %f\n" %(np.array(A_bar_22).mean()) +
             "The variance of A_bar_22: %f" %(np.array(A_bar_22).var()))
    plt.legend()
    plt.title("results of A_bar_22")
    plt.savefig('A_bar_1000.png')
    plt.close()

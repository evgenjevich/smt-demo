import fipy as fp
import numpy as np

mesh = fp.Grid2D(nx=400, ny=400, dx=0.5, dy=0.5)

c_alpha = 0.05
c_beta = 0.95
A = 2.0
kappa = 2.0
c_m = (c_alpha + c_beta) / 2.
B = A / (c_alpha - c_m)**2
D = D_alpha = D_beta = 2. / (c_beta - c_alpha)
c_0 = 0.45
q = np.sqrt((2., 3.))
epsilon = 0.01

c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

r = np.array((mesh.x, mesh.y))
c_var[:] = c_0 + epsilon * np.cos((q[:, None] * r).sum(0))

f_0_var = -A + 3*B*(c_var - c_m)**2 + 3*c_alpha*(c_var - c_alpha)**2 + 3*c_beta*(c_var - c_beta)**2

eqn = fp.TransientTerm(coeff=1.) == fp.DiffusionTerm(D * f_0_var) - fp.DiffusionTerm((D, kappa))

elapsed = 0.0
steps = 0
dt = 0.01
total_sweeps = 2
tolerance = 1e-1
total_steps = 300

c_var[:] = c_0 + epsilon * np.cos((q[:, None] * r).sum(0))

c_var.updateOld()

from fipy.solvers.pysparse import LinearLUSolver as Solver

solver = Solver()

viewer = fp.Viewer(c_var)
while steps < total_steps:
    res0 = eqn.sweep(c_var, dt=dt, solver=solver)

    for sweeps in range(total_sweeps):
        res = eqn.sweep(c_var, dt=dt, solver=solver)

        print ' '
        print 'steps',steps
        print 'res',res
        print 'sweeps',sweeps
        print 'dt',dt


    if res < res0 * tolerance:
        steps += 1
        elapsed += dt
        dt *= 1.1
        if steps % 1 == 0:
             viewer.plot('image{0}.png'.format(steps))
        c_var.updateOld()
    else:
        dt *= 0.8
        c_var[:] = c_var.old

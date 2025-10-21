from two_body.two_body_solver import TwoBodySolver
import numpy as np

two_body_solver = TwoBodySolver(5.972e24, 200)

r1 = [0, 0, 0]
r2 = [6700, 0, 0]
v1 = [0, 0, 0]
v2 = [0, 9.7, 0]

output = two_body_solver.solve(r1, r2, v1, v2, time_span=350*60)

y = output.y
t = output.t

R1_t = np.empty((0, 5))
for T, X0, Y0, Z0, X1, Y1, Z1 in zip(t, y[0], y[1], y[2], y[3], y[4], y[5]):
    a = np.array([[T, X0, Y0, X1, Y1]])
    R1_t = np.append(R1_t, a, axis=0)

np.savetxt("my-file.txt", R1_t, delimiter=",")

exit()

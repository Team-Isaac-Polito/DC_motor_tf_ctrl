import csv
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import fmin

def read_csv(filename):
    with open(filename) as source_file:
        reader = csv.reader(source_file)
        rows = list(reader)
    return rows

def pred_error(theta,y,u):
    a1 = theta[0]
    a2 = theta[1]
    b1 = theta[2]

    N = len(u)
    yp = N*[0]
    e = N*[0]
    for t in range(2,N):
        # strictly proper tf of the 2nd order
        yp[t]=-a1*yp[t-1]-a2*yp[t-2]+b1*u[t-1]
        e[t] = yp[t]-y[t]
    error_norm = np.linalg.norm(e)
    return error_norm

def sim(theta,u):
    a1 = theta[0]
    a2 = theta[1]
    b1 = theta[2]

    N = len(u)
    yp = N*[0]
    for t in range(2,N):
        # strictly proper tf of the 2nd order
        yp[t]=-a1*yp[t-1]-a2*yp[t-2]+b1*u[t-1]
    return yp

rows = read_csv('step1.csv')
r1 = [float(num) for num in rows[0]]
y1 = [float(num) for num in rows[1]]
r2 = [float(num) for num in rows[2]]
y2 = [float(num) for num in rows[3]]
N = len(r1)
Ts = 1/100
t1 = np.linspace(0, Ts*(N-1), num=N)

i_start = 50
i_end = 200
u = r1[i_start:i_end]
y = y1[i_start:i_end]
t = t1[i_start:i_end]

# SYS ID
xopt = fmin(pred_error, [0,0,0], args=(y,u), xtol=1e-8)
xopt = list(xopt)

zero_G = 0
poles_G = np.roots([1]+xopt[0:2])
kg = xopt[2]

yp = sim(xopt,u)

# PLOT
plt.plot(t,u)
plt.plot(t,y)
plt.plot(t,yp)
plt.show()
from FuncDesigner import *
from openopt import SNLE
import numpy as np
import matplotlib.pyplot as plt
from math import *
from numpy import genfromtxt
from DerApproximator import *
import csv


def der(f, x, y, point):
	S = len(f)
	T = len(y) #y.size x,y-> tuple
	U = len(x)
	# df = f.D(point)
	df = []
	for s in range(S):
		 df.append(f[s].D(point))
	return df


# Set parameters
alpha = 0.36
beta = 1.011 ** 16
tau_l = 0.2241
tau_k = 0.2921
Lambda = alpha * tau_k + (1 - alpha) * tau_l
delta = 0.179
A = 3.11

data = genfromtxt('cal_data_cs_22.csv', delimiter=',')

# # 4 States ##########################################
# data = genfromtxt('cal_data_cs_4.csv', delimiter=',')
# tau_l = 0.20
# tau_k = 0.25
# #####################################################

results = np.empty((len(data),3))
i = -1

for row in data:
# row = data[0,:]
	i = i + 1
	nu = row[1].tolist()
	Y_K_ratio = row[2].tolist()

	phi_k_0 = row[3].tolist()
	phi_l = row[4].tolist()
	gamma_0 = row[6].tolist()
	theta_0 = row[7].tolist()
	N = row[8].tolist()

	eta_0 = (Y_K_ratio / A)**(-1/(1-alpha))
	mu_0 = row[5].tolist()
	L_0 = N * (mu_0 * (1-mu_0/2) + theta_0 * (1-mu_0))

	x1, x2, x3 = gamma_0, phi_k_0, theta_0
	y1, y2, y3 = eta_0, mu_0,  L_0

	# Set equations
	gamma, phi_k, theta = oovars('gamma', 'phi_k', 'theta')
	eta, mu, L = oovars('eta', 'mu', 'L')

	point = {gamma:x1, phi_k:x2, theta:x3, eta:y1, mu:y2, L:y3}

	f1 = mu * phi_k * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) + (1-mu)*phi_l* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

	c_yk = 1/(1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

	c_ok = beta * (1-tau_k) * (1+alpha*A*eta**(alpha-1)*L**(1-alpha) - delta) / (1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

	c_yl = (1-tau_l)* theta * (1-alpha)* A* eta**(alpha-1)* L**(-alpha)

	c_ol = 1/(N*(1-mu)) * (tau_l*(1-alpha)*A*eta**(alpha-1)*L**(1-alpha) + tau_k* (alpha*A*eta**(alpha-1)*L**(1-alpha) +N*mu*gamma) - eta**(-1))

	U_k = c_yk * c_ok**beta
	U_l = c_yl * c_ol**beta

	f2 = U_k-U_l

	f3 = L - N * (mu * (1-mu/2) + theta * (1-mu))

	# f1 = x1 * y1 + x2 * y2**2 + x3 * y3**3 
	# f2 = x2 * y1 + x3 * y2**2 + x1 * y3**3 

	x = [gamma ,phi_k, theta]
	y = [eta, mu, L] 
	f = [f1, f2, f3]


	# Calculate dydx
	df1, df2, df3 = der(f,x,y,point)
	# print(df1, df2, df3)

	dfdy = np.empty((3,3))
	dfdx = np.empty((3,3))

	dfdy[0,0] = df1[eta]
	dfdy[0,1] = df1[mu]
	dfdy[0,2] = df1[L]

	dfdy[1,0] = df2[eta]
	dfdy[1,1] = df2[mu]
	dfdy[1,2] = df2[L]

	dfdy[2,0] = 0
	dfdy[2,1] = df3[mu]
	dfdy[2,2] = df3[L]

	dfdx[0,0] = df1[gamma]
	dfdx[0,1] = df1[phi_k]
	dfdx[0,2] = 0

	dfdx[1,0] = df2[gamma]
	dfdx[1,1] = 0
	dfdx[1,2] = df2[theta]

	dfdx[2,0] = 0
	dfdx[2,1] = 0
	dfdx[2,2] = df3[theta]

	inv_dfdy = np.linalg.inv(dfdy)

	dydx = - np.dot(inv_dfdy, dfdx)
	# print(dydx)

	# Calculate dvdx
	dnudy1 = - N * mu_0 * beta * (1-tau_l) * (1-mu_0/2) * (1-alpha)**2 * L_0**(-alpha) * A * eta_0**(alpha-2) / (1+beta)

	dnudy2 =  N * beta * ( (1-tau_l) * (1-alpha) * L_0**(-alpha) * A * eta_0**(alpha-1) * (1-mu_0) - gamma_0 ) / (1+beta)

	dnudy3 = - N * mu_0 * beta * (1-tau_l) * (1-mu_0/2) * alpha * (1-alpha) * L_0**(-alpha-1) * eta_0**(alpha-1) / (1+beta)

	dnudx1 = - N * mu_0 * beta / (1+beta)

	tot_diff_nu_x1 = dnudy1 * dydx[0,0] + dnudy2 * dydx[1,0] + dnudy3 * dydx[2,0] + dnudx1

	tot_diff_nu_x2 = dnudy1 * dydx[0,1] + dnudy2 * dydx[1,1] + dnudy3 * dydx[2,1]

	tot_diff_nu_x3 = dnudy1 * dydx[0,2] + dnudy2 * dydx[1,2] + dnudy3 * dydx[2,2]

	print(tot_diff_nu_x1, tot_diff_nu_x2, tot_diff_nu_x3)

	results[i,:] = tot_diff_nu_x1, tot_diff_nu_x2, tot_diff_nu_x3 

# csv_out = open('total_diff_A.csv', 'wb')
# mywriter = csv.writer(csv_out)
# mywriter.writerows(results)
# csv_out.close()

csv_out = open('total_diff_B.csv', 'wb')
mywriter = csv.writer(csv_out)
mywriter.writerows(results)
csv_out.close()
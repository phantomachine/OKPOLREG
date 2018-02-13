from FuncDesigner import *
from openopt import SNLE
import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt

def solve_eqn(nu, eta, phi_k, phi_l): #y={eta_fixed, mu_fixed, nu_fixed}

	theta, mu, gamma = oovars(3)
	start_point = {theta:0.5, mu:0.8, gamma:0.2}	

	global alpha, beta, tau_l, tau_k, Lambda, delta, A

	# Set equations
	L = 1
	N = L / (mu*(1-mu/2)+(1-mu)*theta)

	f_1 = -nu - delta + N * mu * beta / (1 + beta) * ((1 - tau_l) * (1 - mu / 2) * (1 - alpha) * A * L**(-alpha) *eta ** (alpha - 1) - gamma)

	f_2 = mu * phi_k * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) + (1-mu)*phi_l* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

	c_yk = 1/(1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

	c_ok = beta * (1-tau_k) * (1+alpha*A*eta**(alpha-1)*L**(1-alpha) - delta) / (1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

	c_yl = (1-tau_l)* theta * (1-alpha)* A* eta**(alpha-1)* L**(-alpha)

	c_ol = 1/(N*(1-mu)) * (tau_l*(1-alpha)*A*eta**(alpha-1)*L**(1-alpha) + tau_k* (alpha*A*eta**(alpha-1)*L**(1-alpha) +N*mu*gamma) - eta**(-1))

	U_k = c_yk * c_ok**beta
	U_l = c_yl * c_ol**beta

	f_3 = U_k-U_l

	equations = (f_1 == 0, f_2 == 0, f_3==0)

	p = SNLE(equations, start_point)

	cons_1 = Lambda * A * eta ** alpha + tau_k * mu * gamma * eta - 1
	cons_2 = eta - ((1 - alpha) * A) ** (-1 / alpha)
	cons_3 = (1 - tau_l) * (1 - mu) * (1 - alpha) * A * eta ** (alpha - 1) - gamma
	
	p.constraints = [0 < mu, mu < 1,  cons_1 > 0, cons_2 > 0, cons_3 >=0, N>1, gamma > 0, theta>0, theta<1 ]

	r = p.solve('nssolve')

	theta_s, mu_s, gamma_s = r(theta, mu, gamma)
	print('Solution: x1=%f,   x2=%f, x3=%f'  %(theta_s, mu_s, gamma_s))

	# theta = ( 1/N - mu_s*( 1 - mu_s/2) )/(1 - mu_s)
	N = 1 / (mu_s*(1-mu_s/2)+(1-mu_s)*theta_s)

	return theta_s, mu_s, gamma_s, N, eta


# Set parameters
alpha = 0.36
beta = 1.011 ** 16

# For Group A
tau_l = 0.2241
tau_k = 0.2921

# For Group B
tau_l = 0.2 
tau_k = 0.25 

# For Group C
tau_l = 0.24
tau_k = 0.2921

# For Group D
# tau_l = 0.26
# tau_k = 0.2921

# For Group E (MD)
# tau_l = 0.28
# tau_k = 0.2921

Lambda = alpha * tau_k + (1 - alpha) * tau_l
delta = 0.179
A = 3.11

# data = genfromtxt('cal_data_A.csv', delimiter=',')
# data = genfromtxt('cal_data_B.csv', delimiter=',')
data = genfromtxt('cal_data_C.csv', delimiter=',')
# data = genfromtxt('cal_data_D.csv', delimiter=',')

results = []

for row in data:
	nu = row[1].tolist()
	Y_K_ratio = row[2].tolist()
	eta = (Y_K_ratio / A)**(-1/(1-alpha))
	phi_k = row[3].tolist()
	phi_l = row[4].tolist()
	solutions = solve_eqn(nu, eta, phi_k, phi_l) # mu gamma theta
	results.append(solutions)

# data = genfromtxt('cal_data_E.csv', delimiter=',')
# nu = data[1].tolist()
# Y_K_ratio = data[2].tolist()
# eta = (Y_K_ratio / A)**(-1/(1-alpha))
# phi_k = data[3].tolist()
# phi_l = data[4].tolist()
# solutions = solve_eqn(nu, eta, phi_k, phi_l) # mu gamma theta
# results.append(solutions)


# csv_out = open('bench_GRP_A.csv', 'wb') # For Group A
# csv_out = open('bench_GRP_B.csv', 'wb') # For Group B
csv_out = open('bench_GRP_C.csv', 'wb') # For Group C
# csv_out = open('bench_GRP_D.csv', 'wb') # For Group D
# csv_out = open('bench_GRP_E.csv', 'wb') # For Group D

mywriter = csv.writer(csv_out)
mywriter.writerows(results)
csv_out.close()
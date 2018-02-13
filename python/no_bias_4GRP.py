from FuncDesigner import *
from openopt import SNLE
import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt

# Case for non-politics 

def solve_eqn_np(theta, gamma, N): #y={eta_fixed, mu_fixed, nu_fixed}

	nu, eta, mu, = oovars(3)
	# start_point = {nu:0.08, eta:4, mu:0.8}		
	start_point = {nu:0.05, eta:3., mu:0.5}	
	global alpha, beta, tau_l, tau_k, Lambda, delta, A

	L = N*(mu*(1-mu/2)+(1-mu)*theta)

	f_1 = -nu - delta + N * mu * beta / (1 + beta) * ((1 - tau_l) * (1 - mu / 2) * (1 - alpha) * A * L**(-alpha) *eta ** (alpha - 1) - gamma)

	f_2 = mu * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) + (1-mu)* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

	c_yk = 1/(1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

	c_ok = beta * (1-tau_k) * (1+alpha*A*eta**(alpha-1)*L**(1-alpha) - delta) / (1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

	c_yl = (1-tau_l)* theta * (1-alpha)* A* eta**(alpha-1)* L**(-alpha)

	c_ol = 1/(N*(1-mu)) * (tau_l*(1-alpha)*A*eta**(alpha-1)*L**(1-alpha) + tau_k* (alpha*A*eta**(alpha-1)*L**(1-alpha) +N*mu*gamma) - eta**(-1))

	U_k = c_yk * c_ok**beta
	U_l = c_yl * c_ol**beta

	f_3 = U_k-U_l

	equations = (f_1 == 0, f_2 == 0, f_3==0)

	p = SNLE(equations, start_point)

	# cons_1 = Lambda * A * eta ** alpha + tau_k * mu * gamma * eta - 1
	# cons_2 = eta - ((1 - alpha) * A) ** (-1 / alpha)
	cons_3 = (1 - tau_l) * (1 - mu) * (1 - alpha) * A * eta ** (alpha - 1) - gamma
	
	p.constraints = [nu>0, eta>0,  0< mu, mu < 1,  cons_3 >=0]

	r = p.solve('nssolve')

	nu_s, eta_s, mu_s = r(nu, eta, mu)
	print('Solution: x1=%f,   x2=%f, x3=%f'  %(nu_s, eta_s, mu_s))

	L = N*(mu_s*(1-mu_s/2)+(1-mu_s)*theta)
	Y_K_ratio = A*eta_s**(alpha-1)*L**(1-alpha)

	return nu_s, Y_K_ratio, eta_s, mu_s, L


# Set parameters
alpha = 0.36
beta = 1.011 ** 16

# For Group A
# tau_l = 0.2241
# tau_k = 0.2921

# For Group B
# tau_l = 0.2 
# tau_k = 0.25 

# For Group C
# tau_l = 0.24
# tau_k = 0.2921

# For Group D 
# tau_l = 0.26
# tau_k = 0.2921

# For Group E(MD)
tau_l = 0.28
tau_k = 0.2921

Lambda = alpha * tau_k + (1 - alpha) * tau_l
delta = 0.179
A = 3.11

# data = genfromtxt('bench_GRP_A_result.csv', delimiter=',')
# data = genfromtxt('bench_GRP_B_result.csv', delimiter=',')
# data = genfromtxt('bench_GRP_C_result.csv', delimiter=',')
# data = genfromtxt('bench_GRP_D_result.csv', delimiter=',')
data = genfromtxt('bench_GRP_E_result.csv', delimiter=',')

results = []

for row in data:
	theta = row[6].tolist()
	gamma = row[8].tolist()
	N = row[9].tolist()

	solutions = solve_eqn_np(theta, gamma, N) # mu gamma theta
	results.append(solutions)

# csv_out = open('no_bias_GRP_A.csv', 'wb')
# csv_out = open('no_bias_GRP_B.csv', 'wb')
# csv_out = open('no_bias_GRP_C.csv', 'wb')
# csv_out = open('no_bias_GRP_D.csv', 'wb')
csv_out = open('no_bias_GRP_E.csv', 'wb')

mywriter = csv.writer(csv_out)
mywriter.writerows(results)
csv_out.close()
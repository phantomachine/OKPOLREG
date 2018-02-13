from FuncDesigner import *
from openopt import SNLE
import numpy as np
import matplotlib.pyplot as plt
from math import *
from numpy import genfromtxt

# from c_statics import *

# def example_plot(ax):
# 	f_size = 12
# 	ax.set_xlabel(fontsize=f_size)
# 	ax.set_ylabel(fontsize=f_size)
# 	ax.set_title(fontsize=f_size)

def solve_gamma_3D(gammas): #y={eta_fixed, mu_fixed, nu_fixed}
	gamma_array = []
	sol_1_array = []
	sol_2_array = []
	sol_3_array = []
	sol_4_array = []

	eta, mu, nu = oovars(3)
	start_point = {eta:3.0, mu:0.7, nu:0.20}

	for gamma in gammas:
		# gamma = gamma.tolist()
		global alpha, beta, tau_l, tau_k, Lambda, delta, A, N, phi_k, phi_l, theta

		# Set equations
		L = N * (mu * (1-mu/2) + theta * (1-mu))
		
		f_1 = -nu - delta + N * mu * beta / (1 + beta) * ((1 - tau_l) * (1 - mu / 2) * (1 - alpha) * A * L**(-alpha) *eta ** (alpha - 1) - gamma)

		f_2 = mu * phi_k * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) + (1-mu)*phi_l* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

		c_yk = 1/(1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

		c_ok = beta * (1-tau_k) * (1+alpha*A*eta**(alpha-1)*L**(1-alpha) - delta) / (1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

		c_yl = (1-tau_l)* theta * (1-alpha)* A* eta**(alpha-1)* L**(-alpha)

		c_ol = 1/(N*(1-mu)) * (tau_l*(1-alpha)*A*eta**(alpha-1)*L**(1-alpha) + tau_k* (alpha*A*eta**(alpha-1)*L**(1-alpha) +N*mu*gamma) - eta**(-1))

		U_k = c_yk * c_ok**beta
		U_l = c_yl * c_ol**beta

		f3 = U_k-U_l

		equations = (f_1 == 0, f_2 == 0, f3 == 0)

		p = SNLE(equations, start_point)

		cons_1 = Lambda * A * eta ** alpha + tau_k * mu * gamma * eta - 1
		cons_2 = eta - ((1 - alpha) * A) ** (-1 / alpha)
		cons_3 = (1 - tau_l) * (1 - mu) * (1 - alpha) * A * eta ** (alpha - 1) - gamma
		p.constraints = [0 < mu, mu < 1,  cons_1 > 0, cons_2 > 0, cons_3 >=0, nu>0]

		r = p.solve('nssolve')
		sol_1, sol_2, sol_3 = r(eta, mu, nu)
		sol_4 = N * (sol_2 * (1-sol_2/2) + theta * (1-sol_2))
	
		print('Solution: eta=%f,   mu=%f, nu=%f, L=%f' % (sol_1, sol_2, sol_3, sol_4))
		if r.stopcase>0:
			gamma_array.append(gamma)
			sol_1_array.append(sol_1)
			sol_2_array.append(sol_2)
			sol_3_array.append(sol_3)
			sol_4_array.append(sol_4)				

	return gamma_array, sol_1_array, sol_2_array, sol_3_array, sol_4_array

def solve_phi_k_3D(phi_ks): #y={eta_fixed, mu_fixed, nu_fixed}
	phi_k_array = []
	sol_1_array = []
	sol_2_array = []
	sol_3_array = []
	sol_4_array = []

	eta, mu, nu = oovars(3)
	start_point = {eta:3.0, mu:0.6, nu:0.20}

	for phi_k in phi_ks:
		global alpha, beta, tau_l, tau_k, Lambda, delta, A, N, phi_l, theta

		# Set equations
		L = N * (mu * (1-mu/2) + theta * (1-mu))
		
		f_1 = -nu - delta + N * mu * beta / (1 + beta) * ((1 - tau_l) * (1 - mu / 2) * (1 - alpha) * A * L**(-alpha) *eta ** (alpha - 1) - gamma)

		f_2 = mu * phi_k * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) + (1-mu)*phi_l* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

		c_yk = 1/(1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

		c_ok = beta * (1-tau_k) * (1+alpha*A*eta**(alpha-1)*L**(1-alpha) - delta) / (1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

		c_yl = (1-tau_l)* theta * (1-alpha)* A* eta**(alpha-1)* L**(-alpha)

		c_ol = 1/(N*(1-mu)) * (tau_l*(1-alpha)*A*eta**(alpha-1)*L**(1-alpha) + tau_k* (alpha*A*eta**(alpha-1)*L**(1-alpha) +N*mu*gamma) - eta**(-1))

		U_k = c_yk * c_ok**beta
		U_l = c_yl * c_ol**beta

		f3 = U_k-U_l

		equations = (f_1 == 0, f_2 == 0, f3 == 0)

		p = SNLE(equations, start_point)

		cons_1 = Lambda * A * eta ** alpha + tau_k * mu * gamma * eta - 1
		cons_2 = eta - ((1 - alpha) * A) ** (-1 / alpha)
		cons_3 = (1 - tau_l) * (1 - mu) * (1 - alpha) * A * eta ** (alpha - 1) - gamma
		p.constraints = [0 < mu, mu < 1,  cons_1 > 0, cons_2 > 0, cons_3 >=0, nu>0]

		r = p.solve('nssolve')
		sol_1, sol_2, sol_3 = r(eta, mu, nu)
		phi_k = phi_k / phi_l
		sol_4 = N * (sol_2 * (1-sol_2/2) + theta * (1-sol_2))	
		print('Solution: eta=%f,   mu=%f, nu=%f' % (sol_1, sol_2, sol_3))
		if r.stopcase>0:
			phi_k_array.append(phi_k)
			sol_1_array.append(sol_1)
			sol_2_array.append(sol_2)
			sol_3_array.append(sol_3)			
			sol_4_array.append(sol_4)				

	return phi_k_array, sol_1_array, sol_2_array, sol_3_array, sol_4_array

def solve_theta_3D(thetas): #y={eta_fixed, mu_fixed, nu_fixed}
	theta_array = []
	sol_1_array = []
	sol_2_array = []
	sol_3_array = []
	sol_4_array = []

	eta, mu, nu = oovars(3)
	start_point = {eta:4.0, mu:0.7, nu:0.20}

	for theta in thetas:
		global alpha, beta, tau_l, tau_k, Lambda, delta, A, N, phi_l, phi_k

		# Set equations
		L = N * (mu * (1-mu/2) + theta * (1-mu))
		
		f_1 = -nu - delta + N * mu * beta / (1 + beta) * ((1 - tau_l) * (1 - mu / 2) * (1 - alpha) * A * L**(-alpha) *eta ** (alpha - 1) - gamma)

		f_2 = mu * phi_k * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) + (1-mu)*phi_l* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

		c_yk = 1/(1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

		c_ok = beta * (1-tau_k) * (1+alpha*A*eta**(alpha-1)*L**(1-alpha) - delta) / (1+beta) * ( (1-tau_l)* (1-mu)* (1-alpha)* A* eta**(alpha-1)* L**(-alpha) - gamma)

		c_yl = (1-tau_l)* theta * (1-alpha)* A* eta**(alpha-1)* L**(-alpha)

		c_ol = 1/(N*(1-mu)) * (tau_l*(1-alpha)*A*eta**(alpha-1)*L**(1-alpha) + tau_k* (alpha*A*eta**(alpha-1)*L**(1-alpha) +N*mu*gamma) - eta**(-1))

		U_k = c_yk * c_ok**beta
		U_l = c_yl * c_ol**beta

		f3 = U_k-U_l

		equations = (f_1 == 0, f_2 == 0, f3 == 0)

		p = SNLE(equations, start_point)

		cons_1 = Lambda * A * eta ** alpha + tau_k * mu * gamma * eta - 1
		cons_2 = eta - ((1 - alpha) * A) ** (-1 / alpha)
		cons_3 = (1 - tau_l) * (1 - mu) * (1 - alpha) * A * eta ** (alpha - 1) - gamma
		p.constraints = [0 < mu, mu < 1,  cons_1 > 0, cons_2 > 0, cons_3 >=0, nu>0]

		r = p.solve('nssolve')
		sol_1, sol_2, sol_3 = r(eta, mu, nu)
		sol_4 = N * (sol_2 * (1-sol_2/2) + theta * (1-sol_2))
		print('Solution: eta=%f,   mu=%f, nu=%f' % (sol_1, sol_2, sol_3))
		if r.stopcase>0:
			theta_array.append(theta)
			sol_1_array.append(sol_1)
			sol_2_array.append(sol_2)
			sol_3_array.append(sol_3)			
			sol_4_array.append(sol_4)				

	return theta_array, sol_1_array, sol_2_array, sol_3_array, sol_4_array


# Set parameters
alpha = 0.36
beta = 1.011 ** 16
delta = 0.179
A = 3.11

# data = genfromtxt('cal_data_cs_22.csv', delimiter=',')

# 4 States ##########################################
# data = genfromtxt('cal_data_cs_4.csv', delimiter=',')
# tau_l = 0.20
# tau_k = 0.25
#####################################################

# Alabama
tau_l = 0.24
tau_k = 0.2921
Lambda = alpha * tau_k + (1 - alpha) * tau_l
nu = 0.08
Y_K_ratio = 1.431975901
eta = (Y_K_ratio / A)**(-1/(1-alpha))
phi_k = 0.29
phi_l = 0.319
mu = 0.542882302
gamma = 0.102273664
theta = 0.144793787
N = 2.165863997
eta = (Y_K_ratio / A)**(-1/(1-alpha))

# row = data[0,:]
# nu = row[1].tolist()
# Y_K_ratio = row[2].tolist()
# eta = (Y_K_ratio / A)**(-1/(1-alpha))
# phi_k = row[3].tolist()
# phi_l = row[4].tolist()
# mu = row[5].tolist()
# gamma = row[6].tolist()
# theta = row[7].tolist()
# N = row[8].tolist()
# # fixed_values = (eta, mu, nu)
# eta = (Y_K_ratio / A)**(-1/(1-alpha))

gammas = np.linspace(0.005, 0.150, num=20)
phi_ks = np.linspace(0.5* phi_l, 2.0*phi_l, num=20)
thetas = np.linspace(0.005, 0.3, num=20)

gamma_x, gamma_eta, gamma_mu, gamma_nu, gamma_L = solve_gamma_3D(gammas)
phi_x, phi_eta, phi_mu, phi_nu, phi_L  = solve_phi_k_3D(phi_ks)
theta_x, theta_eta, theta_mu, theta_nu, theta_L = solve_theta_3D(thetas)

# AL
# gammas = np.linspace(0.005, 0.200, num=20)
# phi_ks = np.linspace(0.5* phi_l, 2.0*phi_l, num=20)
# thetas = np.linspace(0.005, 0.3, num=20)

# CA
# gammas = np.linspace(0.005, 0.200, num=20)
# phi_ks = np.linspace(0.7* phi_l, 2.0*phi_l, num=20)
# thetas = np.linspace(0.005, 0.2, num=20)

# CT
# gammas = np.linspace(0.005, 0.200, num=20)
# phi_ks = np.linspace(0.7* phi_l, 2.0*phi_l, num=20)
# thetas = np.linspace(0.05, 0.2, num=20)

# FL
# gammas = np.linspace(0.005, 0.03, num=20)
# phi_ks = np.linspace(0.5* phi_l, 2.0*phi_l, num=20)
# thetas = np.linspace(0.05, 0.3, num=20)

# GA
# gammas = np.linspace(0.005, 0.03, num=20)
# phi_ks = np.linspace(0.5* phi_l, 2.0*phi_l, num=20)
# thetas = np.linspace(0.05, 0.3, num=20)

# IL
# gammas = np.linspace(0.005, 0.10, num=20)
# phi_ks = np.linspace(0.6* phi_l, 1.0*phi_l, num=20)
# thetas = np.linspace(0.05, 0.4, num=20)

# IN
# gammas = np.linspace(0.005, 0.1, num=10)
# phi_ks = np.linspace(0.6* phi_l, 1.2*phi_l, num=10)
# thetas = np.linspace(0.05, 0.1, num=10)

# IA
# gammas = np.linspace(0.05, 0.2, num=10)
# phi_ks = np.linspace(0.6* phi_l, 0.95*phi_l, num=10)
# thetas = np.linspace(0.0001, 0.001, num=10)

# MD
# MA
# MI
# MN
# MS
# NJ
# NC
# OH
# OR
# TN
# VA
# WA
# WI
# CO
# gammas = np.linspace(0.005, 0.100, num=10)
# phi_ks = np.linspace(0.9* phi_l, 1.9*phi_l, num=10)
# thetas = np.linspace(0.09, 0.3, num=10)

# NY
# gammas = np.linspace(0.005, 0.100, num=10)
# phi_ks = np.linspace(0.9* phi_l, 1.9*phi_l, num=10)
# thetas = np.linspace(0.2, 0.4, num=10)

# PA
# TX

fig, axes = plt.subplots(nrows=4, ncols=3)

eta_min = 2.0
eta_max = 4.5

nu_max = 0.25
nu_min = 0.0

mu_min = 0.0
mu_max = 1.0

L_min = 0.9
L_max = 1.2

gamma_min = 0.0
gamma_max = 0.15

phi_min = 0.50
phi_max = 2.00

theta_min = 0.0
theta_max = 0.15


f_size = 12

l1, = axes[0,0].plot(gamma_x, gamma_eta, 'r', linewidth=2)
axes[0,0].set_title('(A1)', fontsize=f_size)
axes[0,0].set_ylim([eta_min, eta_max])
axes[0,0].set_xlim([gamma_min, gamma_max])
axes[0,0].set_ylabel(r'$\eta^* $', fontsize=f_size*1.5)
axes[0,0].set_xticks(np.linspace(gamma_min,gamma_max,4, endpoint=True))

l2, = axes[1,0].plot(gamma_x, gamma_mu, 'b', linewidth=2)
axes[1,0].set_title('(A2)', fontsize=f_size)
axes[1,0].set_ylim([mu_min, mu_max])
axes[1,0].set_xlim([gamma_min, gamma_max])
axes[1,0].set_ylabel(r'$\mu^* $', fontsize=f_size*1.5)
axes[1,0].set_xticks(np.linspace(gamma_min,gamma_max,4, endpoint=True))

l4, = axes[2,0].plot(gamma_x, gamma_L, 'm', linewidth=2)
axes[2,0].set_title('(A3)', fontsize=f_size)
axes[2,0].set_ylim([L_min, L_max])
axes[2,0].set_xlim([gamma_min, gamma_max])
axes[2,0].set_ylabel(r'$L^* $', fontsize=f_size*1.5)
axes[2,0].set_yticks(np.linspace(L_min,L_max, 4, endpoint=True))
axes[2,0].set_xticks(np.linspace(gamma_min,gamma_max,4, endpoint=True))

l3, = axes[3,0].plot(gamma_x, gamma_nu, 'g', linewidth=2)
axes[3,0].set_title('(A4)', fontsize=f_size)
axes[3,0].set_ylim([nu_min, nu_max])
axes[3,0].set_xlim([gamma_min, gamma_max])
axes[3,0].set_xlabel(r'$\rho$' '\n(Market Institution)', fontsize=f_size*1.5)
axes[3,0].set_ylabel(r'$\nu$', fontsize=f_size*1.5)
axes[3,0].set_xticks(np.linspace(gamma_min,gamma_max,4, endpoint=True))

axes[0,1].plot(phi_x, phi_eta, 'r', linewidth=2)
axes[0,1].set_title('(B1)', fontsize=f_size)
axes[0,1].set_ylim([eta_min, eta_max])
axes[0,1].set_xlim([phi_min, phi_max])
axes[0,1].set_xticks(np.linspace(phi_min,phi_max,4, endpoint=True))

axes[1,1].plot(phi_x, phi_mu, 'b', linewidth=2)
axes[1,1].set_title('(B2)', fontsize=f_size)
axes[1,1].set_ylim([mu_min, mu_max])
axes[1,1].set_xlim([phi_min, phi_max])
axes[1,1].set_xticks(np.linspace(phi_min,phi_max,4, endpoint=True))

axes[2,1].plot(phi_x, phi_L, 'm', linewidth=2)
axes[2,1].set_title('(B3)', fontsize=f_size)
axes[2,1].set_ylim([L_min, L_max])
axes[2,1].set_xlim([phi_min, phi_max])
axes[2,1].set_xticks(np.linspace(phi_min,phi_max,4, endpoint=True))
axes[2,1].set_yticks(np.linspace(L_min,L_max, 4, endpoint=True))

axes[3,1].plot(phi_x, phi_nu, 'g', linewidth=2)
axes[3,1].set_title('(B4)', fontsize=f_size)
axes[3,1].set_ylim([nu_min, nu_max])
axes[3,1].set_xlim([phi_min, phi_max])
axes[3,1].set_xticks(np.linspace(phi_min,phi_max,4, endpoint=True))
axes[3,1].set_xlabel(r'$\phi^{*k}(0) / \phi^{*\ell}(0) $' '\n (Political Institution)', fontsize=f_size*1.5)

axes[0,2].plot(theta_x, theta_eta, 'r', linewidth=2)
axes[0,2].set_title('(C1)', fontsize=f_size)
axes[0,2].set_ylim([eta_min, eta_max])
axes[0,2].set_xlim([theta_min, theta_max])
axes[0,2].set_xticks(np.linspace(theta_min,theta_max,4, endpoint=True))

axes[1,2].plot(theta_x, theta_mu, 'b', linewidth=2)
axes[1,2].set_title('(C2)', fontsize=f_size)
axes[1,2].set_ylim([mu_min, mu_max])
axes[1,2].set_xlim([theta_min, theta_max])
axes[1,2].set_xticks(np.linspace(theta_min,theta_max,4, endpoint=True))

axes[2,2].plot(theta_x, theta_L, 'm', linewidth=2)
axes[2,2].set_title('(C3)', fontsize=f_size)
axes[2,2].set_ylim([L_min, L_max])
axes[2,2].set_xlim([theta_min, theta_max])
axes[2,2].set_xticks(np.linspace(theta_min,theta_max,4, endpoint=True))
axes[2,2].set_yticks(np.linspace(L_min,L_max, 4, endpoint=True))

axes[3,2].plot(theta_x, theta_nu, 'g', linewidth=2)
axes[3,2].set_title('(C4)', fontsize=f_size)
axes[3,2].set_ylim([nu_min, nu_max])
axes[3,2].set_xlim([theta_min, theta_max])
axes[3,2].set_xlabel(r'$\theta$' '\n (Inherent Potential)', fontsize=f_size*1.5)
axes[3,2].set_xticks(np.linspace(theta_min,theta_max,4, endpoint=True))

# fig.legend((l1, l2, l3), (r'$\eta$', r'$\mu^*$', r'$\nu$'), loc=8, ncol=3, shadow=True, fancybox=True)
fig.tight_layout()
plt.show()
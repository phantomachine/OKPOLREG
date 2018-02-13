from FuncDesigner import *
from openopt import SNLE
import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt

# Case for non-politics 

def solve_eqn_np(theta, gamma, N, init_nu, init_eta, init_mu): #y={eta_fixed, mu_fixed, nu_fixed}

	nu, eta, mu, = oovars(3)
	# start_point = {nu:0.2, eta:2.25, mu:0.6}		
	start_point = {nu:init_nu, eta:init_eta, mu:init_mu}		

	global alpha, beta, tau_l, tau_k, Lambda, delta, A

	L = N*(mu*(1-mu/2)+(1-mu)*theta)

	f_1 = -nu - delta + N * mu * beta / (1 + beta) * ((1 - tau_l) * (1 - mu / 2) * (1 - alpha) * A * L**(-alpha) *eta ** (alpha - 1) - gamma)

	Log = log(((1-tau_l)*(1-mu)*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) -gamma)/ ((1-tau_l)*(1-alpha)*A*L**(-alpha)*eta**(alpha-1)) -gamma)

	dV_yk = (1-alpha)*mu*eta - gamma/((1-tau_l)*A*L**(-alpha)*eta**(alpha-2))*Log
	dV_yl =  (1-mu)*(1-alpha)*eta
	dV_ok = beta* mu * alpha * (1-alpha)*A* L**(1-alpha) *eta**alpha / (alpha*A*L**(1-alpha)*eta**(alpha-1) + N*mu*gamma) 
	dV_ol = beta * (1-mu)* (tau_l*(1-alpha)**2*A*L**(-alpha)*eta**alpha + tau_k*alpha*(1-alpha)*A*L**(1-alpha)*eta**alpha-1) / (tau_l*(1-alpha)*A*L**(-alpha)*eta**(alpha-1) + tau_k*(alpha*A*L**(1-alpha)*eta**(alpha-1)+N*mu*gamma)-eta**(-1))

	f_2 = dV_yk + dV_yl + dV_ok + dV_ol

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
	
	p.constraints = [nu>0, eta>0,  0< mu, mu < 1, cons_3 >0]
	# cons_1>=0, cons_2 >=0,
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
tau_l = 0.2241
tau_k = 0.2921

# For Group B
tau_l = 0.2 
tau_k = 0.25 

# For Group C
tau_l = 0.24
tau_k = 0.2921

# For Group D 
tau_l = 0.26
tau_k = 0.2921

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
	init_nu = row[11].tolist()
	init_eta = row[12].tolist()
	init_mu = row[13].tolist()

	solutions = solve_eqn_np(theta, gamma, N, init_nu, init_eta, init_mu) # mu gamma theta
	results.append(solutions)

# csv_out = open('social_planner_GRP_A.csv', 'wb')
# csv_out = open('social_planner_GRP_B.csv', 'wb')
# csv_out = open('social_planner_GRP_C.csv', 'wb')
# csv_out = open('social_planner_GRP_D.csv', 'wb')
csv_out = open('social_planner_GRP_E.csv', 'wb')

mywriter = csv.writer(csv_out)
mywriter.writerows(results)
csv_out.close()


# data = genfromtxt('bench_GRP_E_result.csv', delimiter=',')
# J = 1

# inp = data[J]
# theta = inp[6].tolist()
# gamma = inp[8].tolist()
# N = inp[9].tolist()

# init_nu = 0.2
# init_eta = 2.
# init_mu = 0.6
# solutions = solve_eqn_np(theta, gamma, N, init_nu, init_eta, init_mu)

# Alabama (2nd)
# theta = 0.14479378872802878
# gamma =0.1022736641778145
# N =  2.1658639929755723
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.240032
# eta=1.972062
# mu=0.659140
# start_point = {nu:0.15, eta:2.35, mu:0.5}

# Arkansas
# theta = 0.16238501
# gamma = 0.0767882942958462
# N = 2.10141175517456
# init_nu = 0.2
# init_eta = 2.25
# init_mu = 0.6
# solutions = solve_eqn_np(theta, gamma, N, init_nu, init_eta, init_mu)

# California (OK)
# theta = 0.064732172
# gamma =0.133999274808461
# N = 2.283556702
# solutions = solve_eqn_np(theta, gamma, N)

# Conneticat
# theta = 0.280948168
# gamma =0.019065051
# N = 1.880792455
# solutions = solve_eqn_np(theta, gamma, N)

# Florida (2nd)
# theta = 0.093131907463400693
# gamma = 0.14838503800618663
# N = 2.2925673146446255
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.26039457402283606
# eta=2.2757675383308298
# mu=0.66250993063860142
# start_point = {nu:0.15, eta:2.15, mu:0.6}

# Georgia
# theta = 0.090151558
# gamma = 0.105313834
# N = 2.206046757
# solutions = solve_eqn_np(theta, gamma, N)

# Illinois
# theta = 0.314872777
# gamma = 0.016658683
# N = 1.847249252
# solutions = solve_eqn_np(theta, gamma, N)
# Init 

# Indiana
# theta = 0.13274732
# gamma = 0.102842185
# N = 2.180015437
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.2, eta:2.55, mu:0.6}

# Iowa(2nd)
# theta = 0.13934857481510821
# gamma = 0.10513619533381147
# N =  2.1782243301275224
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.24240916416528943
# eta=1.9574830714101905
# mu=0.66008781785079229
# start_point = {nu:0.15, eta:1.9, mu:0.6}

# Maryland(2nd)
# theta = 0.17369235904507718
# gamma = 0.12068220489970838
# N = 2.1583687364463646
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.25350266459287696
# eta=1.6025155504212414
# mu=0.62166052384715464
# start_point = {nu:0.18, eta:1.91, mu:0.65}

# Massachusetts
# theta = 0.042295089
# gamma = 0.1171106
# N = 2.22487717
# solutions = solve_eqn_np(theta, gamma, N)
	# start_point = {nu:0.15, eta:2.5, mu:0.55}

# Michigan
# theta = 0.135045284
# gamma = 0.134714002
# N = 2.223362828
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:2.55, mu:0.55}

# Minnesota
# theta = 0.24941988
# gamma = 0.026986358
# N = 1.934062862
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:3.2, mu:0.55}

# Mississipi
# theta = 0.058338673
# gamma = 0.156713333
# N = 2.310831942
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:2.5, mu:0.55}

# New Jersey
# theta = 0.158986463
# gamma = 0.077561123
# N = 2.095195891
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:3.1, mu:0.55}

# North Carolina
# theta = 0.02130909
# gamma = 0.183301391
# N = 2.298406324
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:2.6, mu:0.5}

# Ohio
# theta = 0.127258256
# gamma = 0.115963405
# N = 2.205251503
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:2.6, mu:0.5}

# Oregon
# theta = 0.356816701
# gamma = 0.010042274
# N = 1.786573475
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.15, eta:2.6, mu:0.5}

# Tennessee
# theta = 0.038075645
# gamma = 0.158283782
# N = 2.295097846
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.8, eta:2.5, mu:0.55}

# Virginia
# theta = 0.038695539
# gamma = 0.187632081
# N = 2.303485281
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.9, eta:2.5, mu:0.55}

# Washington
# theta = 0.017434791
# gamma = 0.172122095
# N = 2.27954689
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.2, eta:2.5, mu:0.55}

# Wisconsin
# theta = 0.057273101
# gamma = 0.136933921
# N = 2.277971087
# solutions = solve_eqn_np(theta, gamma, N)
# start_point = {nu:0.2, eta:2.5, mu:0.55}

# Cololado
# theta = 0.045083748
# gamma = 0.0904642687042879
# N = 2.14837990895219
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.195618
# eta=2.734511
# mu=0.766077
# start_point = {nu:0.15, eta:3.2, mu:0.6}

# New York
# theta = 0.324649304
# gamma = 0.033029477
# N = 1.814479382
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.107547
# eta=3.576764
# mu=0.665519
# start_point = {nu:0.15, eta:3.2, mu:0.6}

# Pennsylvania
# theta = 0.176339004367691
# gamma = 0.0693073013812467
# N = 2.022219072
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.140509
# eta=3.193431
# mu=0.683428
# start_point = {nu:0.15, eta:3.5, mu:0.6}

# Texas
# theta = 0.183275453
# gamma = 0.048007881
# N = 1.995482978
# solutions = solve_eqn_np(theta, gamma, N)
# nu=0.150546
# eta=3.286514
# mu=0.708223
# start_point = {nu:0.15, eta:3.5, mu:0.6}
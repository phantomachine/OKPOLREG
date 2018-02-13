import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    s = 12
    for i, point in a.iterrows():
        ax.text(point['x'], point['y'], str(point['val']), fontsize=s)


df_A = pd.read_csv('4GRP_result_GRP_A_plot.csv', usecols=[
                 'Abb', 'nu', 'phi', 'theta', 'mu', 'gamma', 'N', 'eta'])
df_B = pd.read_csv('4GRP_result_GRP_B_plot.csv', usecols=[
                 'Abb', 'nu', 'phi', 'theta', 'mu', 'gamma', 'N', 'eta'])
df_C = pd.read_csv('4GRP_result_GRP_C_plot.csv', usecols=[
                 'Abb', 'nu', 'phi', 'theta', 'mu', 'gamma', 'N', 'eta'])
df_D = pd.read_csv('4GRP_result_GRP_D_plot.csv', usecols=[
                 'Abb', 'nu', 'phi', 'theta', 'mu', 'gamma', 'N', 'eta'])
df_E = pd.read_csv('4GRP_result_GRP_E_plot.csv', usecols=[
                 'Abb', 'nu', 'phi', 'theta', 'mu', 'gamma', 'N', 'eta'])

ms = 8
f_size = 12

gamma_min = 0.0
gamma_max = 0.2

phi_min = 0.60
phi_max = 1.8

theta_min = 0.0
theta_max = 0.4


fig2, ax2 = plt.subplots(nrows=1, ncols=2)
g_A, = ax2[0].plot(df_A.gamma, df_A.phi, color='#551A8B', marker='o', markersize=ms, linestyle='None' )
g_B, = ax2[0].plot(df_B.gamma, df_B.phi, color='#551A8B', marker='^', markersize=ms, linestyle='None' )
g_C, = ax2[0].plot(df_C.gamma, df_C.phi, color='#551A8B', marker='s', markersize=ms, linestyle='None' )
g_D, = ax2[0].plot(df_D.gamma, df_D.phi, color='#551A8B', marker='D', markersize=ms, linestyle='None' )
g_E, = ax2[0].plot(df_E.gamma, df_E.phi, color='#551A8B', marker='p', markersize=ms+2, linestyle='None' )

# g_B, = ax2[0].plot(df_4.gamma, df_4.phi, mfc='#9B30FF', marker='^', markersize=ms+2, linestyle='None' )
label_point(df_A.gamma, df_A.phi, df_A.Abb, ax2[0])
label_point(df_B.gamma, df_B.phi, df_B.Abb, ax2[0])
label_point(df_C.gamma, df_C.phi, df_C.Abb, ax2[0])
label_point(df_D.gamma, df_D.phi, df_D.Abb, ax2[0])
label_point(df_E.gamma, df_E.phi, df_E.Abb, ax2[0])

ax2[0].set_title('(1)', fontsize=f_size)
ax2[0].set_ylim([phi_min, phi_max])
ax2[0].set_xlim([gamma_min, gamma_max])
ax2[0].set_xlabel(r'$\rho $', fontsize=f_size*1.5)
ax2[0].set_ylabel(r'$\phi^{*k}(0) / \phi^{*l}(0)$', fontsize=f_size*1.2)

ax2[1].plot(df_A.gamma, df_A.theta, color='#551A8B', marker='o', markersize=ms, linestyle='None' )
ax2[1].plot(df_B.gamma, df_B.theta, color='#551A8B', marker='^', markersize=ms, linestyle='None' )
ax2[1].plot(df_C.gamma, df_C.theta, color='#551A8B', marker='s', markersize=ms, linestyle='None' )
ax2[1].plot(df_D.gamma, df_D.theta, color='#551A8B', marker='D', markersize=ms, linestyle='None' )
ax2[1].plot(df_E.gamma, df_E.theta, color='#551A8B', marker='p', markersize=ms+2, linestyle='None' )

label_point(df_A.gamma, df_A.theta, df_A.Abb, ax2[1])
label_point(df_B.gamma, df_B.theta, df_B.Abb, ax2[1])
label_point(df_C.gamma, df_C.theta, df_C.Abb, ax2[1])
label_point(df_D.gamma, df_D.theta, df_D.Abb, ax2[1])
label_point(df_E.gamma, df_E.theta, df_E.Abb, ax2[1])

ax2[1].set_title('(2)', fontsize=f_size)
ax2[1].set_ylim([theta_min, theta_max])
ax2[1].set_xlim([gamma_min, gamma_max])
ax2[1].set_xlabel(r'$\rho $', fontsize=f_size*1.5)
ax2[1].set_ylabel(r'$\theta$', fontsize=f_size*1.5)

fig2.legend((g_A, g_B, g_C, g_D, g_E), ('Group A', 'Group B', 'Group C', 'Group D', 'Group E'), loc=8, ncol=5, shadow=True, fancybox=True, numpoints = 1)

fig2.tight_layout()
plt.show()
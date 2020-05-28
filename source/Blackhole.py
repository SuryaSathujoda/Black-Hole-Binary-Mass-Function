import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_sin(date, alpha, beta, gamma):
    return alpha*np.sin(beta*date + gamma)


date, v, v_err = np.loadtxt('../data/blackhole_data.csv', unpack=True, delimiter=',',  skiprows=2)

params, covar = curve_fit(fit_sin, date, v, sigma=v_err, absolute_sigma=True)
alpha, beta, gamma = params
print(alpha, beta, gamma)
alpha_err, beta_err, gamma_err = np.diag(covar)
print(alpha_err, beta_err, gamma_err)

best_T = 2*np.pi/beta
best_T_err = np.sqrt((-2*np.pi/beta**2)**2 * beta_err**2)
print("Time Period = ", best_T)
print("Error in Time Period = ", best_T_err)

print("Max Velocity = ", np.abs(alpha))
print("Error in Max Velocity = ", np.abs(alpha_err))

date_smooth = np.arange(np.min(date), np.max(date), 0.005)
best_v = fit_sin(date, alpha, beta, gamma)
best_v_smooth = fit_sin(date_smooth, alpha, beta, gamma)

best_chi2 = np.sum((best_v - v)**2/(v_err**2))
print("Chi_squared = ", best_chi2)
print("No. of data points = ", len(date))

plt.plot(date, v, color='red', marker='o', markersize=5, linestyle='None')
plt.plot(date_smooth, best_v_smooth, color='blue', marker='None', linewidth=2, linestyle='-')
#plt.show()

plt.errorbar(date, v, v_err, color='red', linewidth=2, linestyle='None')
#plt.show()

best_T_SI = best_T * 3600 * 24
v_peak_SI = alpha * 10**3

# Calc error in f

f = (best_T_SI * np.abs(v_peak_SI)**3)/(2*np.pi*6.673*10**-11)
print("Binary Mass Function = ", f)

inc_angle = np.deg2rad(56)
inc_angle_err = np.deg2rad(4)

q = 0.067
q_err = 0.005

M = (f * 1.067**2)/(np.sin(inc_angle))**3
print("Mass of Black Hole = ", M)


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_sin(date, alpha, beta, gamma):
    return alpha*np.sin(beta*date + gamma)


date, v, v_err = np.loadtxt('../data/blackhole_data.csv', unpack=True, delimiter=',',  skiprows=2)

params, covar = curve_fit(fit_sin, date, v, sigma=v_err, absolute_sigma=True)
alpha, beta, gamma = params
print(alpha, beta, gamma)
alpha_err, beta_err, gamma_err = np.sqrt(np.diag(covar))
print(alpha_err, beta_err, gamma_err)
print("correlation = ", covar[1,2])

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
print("Reduced Chi_squared = ", best_chi2/(len(date)-3.0))

plt.xlabel("Days (mjd)")
plt.ylabel("Velocity (km/s)")
plt.plot(date, v, color='red', marker='o', markersize=5, linestyle='None')
plt.plot(date_smooth, best_v_smooth, color='blue', marker='None', linewidth=2, linestyle='-')
plt.show()

plt.xlabel("Days (mjd)")
plt.ylabel("Velocity (km/s)")
plt.errorbar(date, v, v_err, color='red', linewidth=2, linestyle='None')
plt.show()

best_T_SI = best_T * 3600 * 24
best_T_err_SI = best_T_err * 3600 * 24
v_peak_SI = alpha * 10**3
v_peak_err_SI = alpha_err * 10**3

# Calc error in f

G_denom = 2*np.pi*6.673*10**-11
f = (best_T_SI * np.abs(v_peak_SI)**3)/G_denom
f_err = np.sqrt((np.abs(v_peak_SI)**3/G_denom)**2*best_T_err_SI**2 + 
                ((3*best_T_SI*np.abs(v_peak_SI)**2)/G_denom)**2*v_peak_err_SI**2)
print("Binary Mass Function = ", f)
print("Binary Mass Function Err = ", f_err)


inc_angle = np.deg2rad(56)
inc_angle_err = np.deg2rad(4)

q = 0.067
q_err = 0.005

M = (f * 1.067**2)/(np.sin(inc_angle))**3
M_err = ((2*(1+q)*f)/np.sin(inc_angle)**3)*q_err**2
M_err = M_err + ((1+q)**2/np.sin(inc_angle)**3)*f_err**2
M_err = M_err + ((-3*(1+q)**2*f*np.cos(inc_angle))/np.sin(inc_angle)**4)**2*inc_angle_err**2
M_err = np.sqrt(M_err)
print("Mass of Black Hole = ", M)
print("Mass of Black Hole Err = ", M_err)


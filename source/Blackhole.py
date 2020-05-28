import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_sin(date, alpha, beta, gamma):
    return alpha*np.sin(beta*date + gamma)


date, v, v_err = np.loadtxt('../data/blackhole_data.csv', unpack=True, delimiter=',',  skiprows=2)

params, covar = curve_fit(fit_sin, date, v)
alpha, beta, gamma= params
date_smooth = np.arange(np.min(date), np.max(date), 0.005)
best_v = fit_sin(date, alpha, beta, gamma)
best_v_smooth = fit_sin(date_smooth, alpha, beta, gamma)

plt.plot(date, v, color='red', marker='o', markersize=5, linestyle='None')
plt.plot(date_smooth, best_v_smooth, color='blue', marker='None', linewidth=2, linestyle='-')
plt.show()

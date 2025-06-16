import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, FuncFormatter)

from scipy.optimize import curve_fit

import os
import sys

my_exp = lambda x,a,b,c: b-a*np.exp(-c*x)

def fit_curve_to_data(curve, x_data, y_data, fitting_range = [0, None], p0 = None):
	return curve_fit(curve, x_data[fitting_range[0]:fitting_range[1]], y_data[fitting_range[0]:fitting_range[1]], p0=p0)

def sma_smooth(data, shift):
	data_smooth = np.zeros(len(data))
	for i in range(shift,len(data)-shift-1):
		data_smooth[i] = np.mean(data[i-shift:i+shift+1])
	data_smooth[:shift] = data[:shift]
	data_smooth[len(data)-shift-1:] = data[len(data)-shift-1:]
	return data_smooth

fig, ax = plt.subplots()
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]/2, DefaultSize[1]/2) )
# Import config (FIXME: change the path accordingly)
configfile = 'Config.pyc'
sys.path.append(os.path.dirname(os.path.expanduser(configfile)))
from Config import (plot_config,my_formatter)
# Set the config
plot_config()
# label formatter (in Config.py)
formatter = FuncFormatter(my_formatter)
# create grid
nx=1
ny=1

char_labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
layer_labels = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]

N_gdp = len(np.transpose( np.genfromtxt('data/cylindrical/gdp_{}{}_phi.txt'.format("A", "16")) ))
N_gtp = len(np.transpose( np.genfromtxt('data/cylindrical/gtp_{}{}_phi.txt'.format("A", "16")) ))

ts = np.arange(40000)/10000
ts_gdp = np.arange(20000,20000+N_gdp)/10000
ts_gtp = np.arange(20000,20000+N_gtp)/10000
gdp_r_per_layer = {i:np.zeros(N_gdp) for i in range(1, 17)}
gtp_r_per_layer = {i:np.zeros(N_gtp) for i in range(1, 17)}
gdp_z_per_layer = {i:np.zeros(N_gdp) for i in range(1, 17)}
gtp_z_per_layer = {i:np.zeros(N_gtp) for i in range(1, 17)}
gdp_phi_per_layer = {i:np.zeros(N_gdp) for i in range(1, 17)}
gtp_phi_per_layer = {i:np.zeros(N_gtp) for i in range(1, 17)}
dt = ts_gdp[1]-ts_gdp[0]

gdp_r_per_layer_ref = {i:np.zeros(40000) for i in range(1, 17)}
gtp_r_per_layer_ref = {i:np.zeros(40000) for i in range(1, 17)}
gdp_z_per_layer_ref = {i:np.zeros(40000) for i in range(1, 17)}
gtp_z_per_layer_ref = {i:np.zeros(40000) for i in range(1, 17)}

gdp_phi_per_fil = {i:np.zeros(N_gdp) for i in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]}
gtp_phi_per_fil = {i:np.zeros(N_gtp) for i in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]}

for label in char_labels:
	for layer in layer_labels:
		phii = np.transpose( np.genfromtxt('data/cylindrical/gdp_{}{}_phi.txt'.format(label, layer)) ) - np.transpose( np.genfromtxt('data/cylindrical/gdp_{}{}_phi.txt'.format(label, "16")) )
		gdp_phi_per_fil[label]+=phii
		phii = np.transpose( np.genfromtxt('data/cylindrical/gtp_{}{}_phi.txt'.format(label, layer)) ) - np.transpose( np.genfromtxt('data/cylindrical/gtp_{}{}_phi.txt'.format(label, "16")) )
		gtp_phi_per_fil[label]+=phii
	gdp_phi_per_fil[label]/=16
	gtp_phi_per_fil[label]/=16

for layer in layer_labels:
	for label in char_labels:
		ri = np.transpose( np.genfromtxt('data/cylindrical/gdp_{}{}_r.txt'.format(label, layer)) )
		zi = np.transpose( np.genfromtxt('data/cylindrical/gdp_{}{}_z.txt'.format(label, layer)) ) - np.transpose( np.genfromtxt('data/cylindrical/gdp_{}{}_z.txt'.format(label, "16")) ) 
		gdp_r_per_layer[int(layer)]+=ri
		gdp_z_per_layer[int(layer)]+=zi
		ri = np.transpose( np.genfromtxt('data/cylindrical/gtp_{}{}_r.txt'.format(label, layer)) )
		zi = np.transpose( np.genfromtxt('data/cylindrical/gtp_{}{}_z.txt'.format(label, layer)) ) - np.transpose( np.genfromtxt('data/cylindrical/gtp_{}{}_z.txt'.format(label, "16")) ) 
		gtp_r_per_layer[int(layer)]+=ri
		gtp_z_per_layer[int(layer)]+=zi
		ri_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/data/cylindrical/{}_{}{}_r.txt'.format("gdp", label, layer)) )
		zi_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/data/cylindrical/{}_{}{}_z.txt'.format("gdp", label, layer)) )
		gdp_r_per_layer_ref[int(layer)]+=ri_ref
		gdp_z_per_layer_ref[int(layer)]+=zi_ref
		ri_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/data/cylindrical/{}_{}{}_r.txt'.format("gtp", label, layer)) )
		zi_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/data/cylindrical/{}_{}{}_z.txt'.format("gtp", label, layer)) )
		gtp_r_per_layer_ref[int(layer)]+=ri_ref
		gtp_z_per_layer_ref[int(layer)]+=zi_ref

	gdp_r_per_layer[int(layer)]/=14
	gtp_r_per_layer[int(layer)]/=14
	gdp_z_per_layer[int(layer)]/=14
	gtp_z_per_layer[int(layer)]/=14
	gdp_r_per_layer_ref[int(layer)]/=14
	gtp_r_per_layer_ref[int(layer)]/=14
	gdp_z_per_layer_ref[int(layer)]/=14
	gtp_z_per_layer_ref[int(layer)]/=14


for j, layer in enumerate(layer_labels):
	grid = plt.GridSpec(nrows=ny, ncols=nx, wspace=0, hspace=0)
	axa = plt.subplot(grid[0, 0])
	axa.set_xlabel(r'Time, $t$ ($\si{\micro\second}$)')
	axa.set_ylabel(r'Radial coord., $\rho$ ($\si{\angstrom}$)')
	shift = 5000
	axa.plot(ts, gdp_r_per_layer_ref[int(layer)], '-', color = 'blue')
	axa.plot(ts, gtp_r_per_layer_ref[int(layer)], '-', color = 'red')
	axa.plot(ts_gdp, gdp_r_per_layer[int(layer)], ':', color = 'blue', label = 'GDP layer {}'.format(layer))
	axa.plot(ts_gtp, gtp_r_per_layer[int(layer)], ':', color = 'red', label = 'GTP layer {}'.format(layer))
	axa.legend()
	plt.savefig('figures/cylindrical/r_time_{}.pdf'.format(int(layer)), bbox_inches = 'tight', dpi = 100)
	del axa

for layer in layer_labels:
	grid = plt.GridSpec(nrows=ny, ncols=nx, wspace=0, hspace=0)
	axa = plt.subplot(grid[0, 0])
	axa.set_xlabel(r'Time, $t$ ($\si{\micro\second}$)')
	axa.set_ylabel(r'Vertical coord., $z$ ($\si{\angstrom}$)')
	axa.plot(ts_gdp, gdp_z_per_layer[int(layer)], color = 'blue', label = 'GDP layer {}'.format(layer))
	axa.plot(ts_gtp, gtp_z_per_layer[int(layer)], color = 'red', label = 'GTP layer {}'.format(layer))
	axa.legend()
	plt.savefig('figures/cylindrical/z_time_{}.pdf'.format(int(layer)), bbox_inches = 'tight', dpi = 100)
	del axa

for label in char_labels:
	grid = plt.GridSpec(nrows=ny, ncols=nx, wspace=0, hspace=0)
	axa = plt.subplot(grid[0, 0])
	axa.set_xlabel(r'Time, $t$ ($\si{\micro\second}$)')
	axa.set_ylabel(r'Polar coord., $\phi$ ($\si{\degree}$)')
	axa.plot(ts_gdp, gdp_phi_per_fil[label], color = 'blue', label = 'GDP')
	axa.plot(ts_gtp, gtp_phi_per_fil[label], color = 'red', label = 'GTP')
	axa.legend()
	plt.savefig('figures/cylindrical/phi_time_{}.pdf'.format(label), bbox_inches = 'tight', dpi = 100)
	del axa

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, FuncFormatter)

import os
import sys

def prepara_canva():
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
	grid = plt.GridSpec(nrows=ny, ncols=nx, wspace=0, hspace=0)
	return plt.subplot(grid[0, 0])

def plot_histogram(title, xlabel, ylabel, xlim = None):
	axa = prepara_canva()
	gdp_hist_filename = 'data/hist/gdp_{}.dat'.format(title)
	gtp_hist_filename = 'data/hist/gtp_{}.dat'.format(title)
	bc, h_gdp, dh_gdp = np.transpose( np.genfromtxt(gdp_hist_filename) )
	_, h_gtp, dh_gtp = np.transpose( np.genfromtxt(gtp_hist_filename) )
	_, h_gdp_ref, dh_gdp_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/{}'.format(gdp_hist_filename)) )
	_, h_gtp_ref, dh_gtp_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/{}'.format(gtp_hist_filename)) )
	axa.set_xlabel(xlabel)
	axa.set_ylabel(ylabel)
	if xlim is not None:
		axa.set_xlim(xlim)
	axa.plot(bc, h_gdp_ref, color = 'blue')
	axa.fill_between(bc, h_gdp_ref+dh_gdp_ref, h_gdp_ref-dh_gdp_ref, color = 'blue', alpha = 0.4)
	axa.plot(bc, h_gtp_ref, color = 'red')
	axa.fill_between(bc, h_gtp_ref+dh_gtp_ref, h_gtp_ref-dh_gtp_ref, color = 'red', alpha = 0.4)
	axa.plot(bc, h_gdp, ':', color = 'blue')
	axa.plot(bc, h_gtp, ':', color = 'red')
	name='figures/hist/{}.pdf'.format(title)
	plt.savefig(name, bbox_inches = 'tight', dpi = 100)
	del axa

# -----------------------------------------------------------------------------
# DISTANCES
# -----------------------------------------------------------------------------

xlabel_long_intra = r'Distance, $r_{\mathrm{long,intra}}$ ($\si{\angstrom}$)'
xlabel_long_inter = r'Distance, $r_{\mathrm{long,inter}}$ ($\si{\angstrom}$)'
xlabel_alpha_lat_inter = r'Distance, $r_{\alpha\mathrm{,lat,inter}}$ ($\si{\angstrom}$)'
xlabel_beta_lat_inter = r'Distance, $r_{\beta\mathrm{,lat,inter}}$ ($\si{\angstrom}$)'
xlabel_seam_lat_inter = r'Distance, $r_{\mathrm{seam,lat,inter}}$ ($\si{\angstrom}$)'
xlabel_diag1 = r'Distance, $r_{\mathrm{diag,}1}$ ($\si{\angstrom}$)'
xlabel_diag2 = r'Distance, $r_{\mathrm{diag,}2}$ ($\si{\angstrom}$)'

ylabel_dist = r'Prob. dens., $\mathcal{P}(r)$ ($\si{\per\angstrom}$)'

plot_histogram('long_intra', xlabel_long_intra, ylabel_dist, xlim=(39.0, 51.0))
plot_histogram('long_inter', xlabel_long_inter, ylabel_dist, xlim=(39.0, 51.0))
plot_histogram('alpha_lat_inter_all', xlabel_alpha_lat_inter, ylabel_dist, xlim=(41.0, 89.0))
plot_histogram('beta_lat_inter_all', xlabel_beta_lat_inter, ylabel_dist, xlim=(41.0, 89.0))
plot_histogram('seam_lat_inter_all', xlabel_seam_lat_inter, ylabel_dist, xlim=(41.0, 129.0))
plot_histogram('diag1_all', xlabel_diag1, ylabel_dist, xlim=(59.0, 91.0))
plot_histogram('diag2_all', xlabel_diag2, ylabel_dist, xlim=(41.0, 79.0))

# -----------------------------------------------------------------------------
# ANGLES
# -----------------------------------------------------------------------------

xlabel_bab_long_angle = r'Angle, $\theta_{\beta\alpha\beta\mathrm{,long}}$ ($\si{\degree}$)'
xlabel_aba_long_angle = r'Angle, $\theta_{\alpha\beta\alpha\mathrm{,long}}$ ($\si{\degree}$)'

ylabel_angle = r'Prob. dens., $\mathcal{P}(\theta)$ ($\si{\per\degree}$)'

plot_histogram('aba_long_angle_top', xlabel_aba_long_angle, ylabel_angle, xlim=(159.0, 181.0))
plot_histogram('bab_long_angle_top', xlabel_bab_long_angle, ylabel_angle, xlim=(159.0, 181.0))

# -----------------------------------------------------------------------------
# DIHEDRALS
# -----------------------------------------------------------------------------

xlabel_long_dihe = r'Dihedral angle, $\varphi_{\mathrm{long}}$ ($\si{\degree}$)'
xlabel_out_dihe = r'Dihedral angle, $\varphi_{\mathrm{out}}$ ($\si{\degree}$)'

ylabel_dihe = r'Prob. dens., $\mathcal{P}(\varphi)$ ($\si{\per\degree}$)'

plot_histogram('long_dihe_top', xlabel_long_dihe, ylabel_dihe)
plot_histogram('out_dihe_top', xlabel_out_dihe, ylabel_dihe)

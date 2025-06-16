import numpy as np
import sys

which = sys.argv[1]

potentials = ["gdp_long_intra", "gdp_long_inter", "gdp_alpha_lat_inter_all", "gdp_beta_lat_inter_all", "gdp_seam_lat_inter_all",
			  "gtp_long_intra", "gtp_long_inter", "gtp_alpha_lat_inter_all", "gtp_beta_lat_inter_all", "gtp_seam_lat_inter_all",
			  "gdp_diag1_all", "gdp_diag2_all", "gdp_aba_long_angle_top", "gdp_bab_long_angle_top",
			  "gtp_diag1_all", "gtp_diag2_all", "gtp_aba_long_angle_top", "gtp_bab_long_angle_top",
			  "gdp_long_dihe_top", "gtp_long_dihe_top", "gdp_out_dihe_top", "gtp_out_dihe_top"]

for potential in potentials:
	data = np.transpose( np.genfromtxt('data/pot/{}.dat'.format(potential)) )
	N = len(data[0])
	if 'dihe' in potential.split('_'): dihe_type = True
	else: dihe_type = False
	with open('../bd_simulation/potentials/{}.{}.table'.format(potential, which), 'w') as table_file:
		if dihe_type: table_file.write('\nENTRY\nN {} NOF\n\n'.format(N))
		else: table_file.write('\nENTRY\nN {}\n\n'.format(N))
		for i in range(N):
			if dihe_type: table_file.write('{} {} {}\n'.format(i+1, data[1][i], data[2][i]))
			else: table_file.write('{} {} {} {}\n'.format(i+1, data[1][i], data[2][i], data[3][i]))
	with open('../bd_simulation/potentials/{}.table'.format(potential), 'w') as table_file:
		if dihe_type: table_file.write('\nENTRY\nN {} NOF\n\n'.format(N))
		else: table_file.write('\nENTRY\nN {}\n\n'.format(N))
		for i in range(N):
			if dihe_type: table_file.write('{} {} {}\n'.format(i+1, data[1][i], data[2][i]))
			else: table_file.write('{} {} {} {}\n'.format(i+1, data[1][i], data[2][i], data[3][i]))
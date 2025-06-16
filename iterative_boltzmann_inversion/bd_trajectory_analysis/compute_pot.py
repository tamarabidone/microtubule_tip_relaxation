import numpy as np
from scipy.interpolate import splrep, BSpline
import sys

energy_conv = 0.5922
convergence = float(sys.argv[1])

def compute_pot(r, h, pot_type):
	# TODO: shift it somhow in order to have 0 at dissociation? but what to do with angle/dihedral?
	if pot_type == 'bond':
		norm = r**2
		# h[h==0]+=epsilon
	elif pot_type == 'angle':
		# h[h==0]+=epsilon
		norm = np.sin(np.deg2rad(r))
		norm[0] = norm[1]
		norm[-1] = norm[-2]
	elif pot_type == 'dihe':
		# h[h==0]+=epsilon
		norm = 1
	else: return 1/0
	pot0 = -energy_conv*np.log(h/norm)
	pot0[h==0] = 2*np.max(pot0[h!=0])
	return pot0
	# return BSpline(*splrep(r, pot0, s=len(r)))(r)

def compute_force(pot, dr):
	force = -np.gradient(pot)/dr
	return force

def write_table(title, pot_type):
	gdp_hist_filename = 'data/hist/gdp_{}.dat'.format(title)
	gtp_hist_filename = 'data/hist/gtp_{}.dat'.format(title)
	bc, h_gdp, dh_gdp = np.transpose( np.genfromtxt(gdp_hist_filename) )
	_, h_gtp, dh_gtp = np.transpose( np.genfromtxt(gtp_hist_filename) )
	_, h_gdp_ref, dh_gdp_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/data/hist/gdp_{}.dat'.format(title)) )
	_, h_gtp_ref, dh_gtp_ref = np.transpose( np.genfromtxt('../../md_trajectory_analysis/data/hist/gtp_{}.dat'.format(title)) )
	dr = bc[1]-bc[0]
	pot_pmf_gdp = compute_pot(bc, h_gdp, pot_type)
	pot_pmf_gtp = compute_pot(bc, h_gtp, pot_type)
	pot_pmf_gdp_ref = compute_pot(bc, h_gdp_ref, pot_type)
	pot_pmf_gtp_ref = compute_pot(bc, h_gtp_ref, pot_type)
	try:
		_, _, pot_gdp_prev, _ = np.transpose( np.genfromtxt('data/pot/gdp_{}.dat'.format(title)) )
		_, _, pot_gtp_prev, _ = np.transpose( np.genfromtxt('data/pot/gtp_{}.dat'.format(title)) )
	except:
		pot_gdp_prev = pot_pmf_gdp_ref
		pot_gtp_prev = pot_pmf_gtp_ref
	new_potential_gdp = pot_gdp_prev + convergence*(pot_pmf_gdp_ref - pot_pmf_gdp)
	new_potential_gtp = pot_gtp_prev + convergence*(pot_pmf_gtp_ref - pot_pmf_gtp)
	new_force_gdp = compute_force(new_potential_gdp, dr)
	new_force_gtp = compute_force(new_potential_gtp, dr)
	np.savetxt('data/pot/gdp_{}.dat'.format(title), np.transpose(np.array([range(1, len(bc)+1), bc, new_potential_gdp, new_force_gdp])), fmt='%i %f %f %f')
	np.savetxt('data/pot/gtp_{}.dat'.format(title), np.transpose(np.array([range(1, len(bc)+1), bc, new_potential_gtp, new_force_gtp])), fmt='%i %f %f %f')

# # -----------------------------------------------------------------------------
# # DISTANCES
# # -----------------------------------------------------------------------------

pot_type = 'bond'

write_table('long_intra', pot_type)
write_table('long_inter', pot_type)
write_table('alpha_lat_inter_all', pot_type)
write_table('beta_lat_inter_all', pot_type)
write_table('seam_lat_inter_all', pot_type)
write_table('diag1_all', pot_type)
write_table('diag2_all', pot_type)

# # -----------------------------------------------------------------------------
# # ANGLES
# # -----------------------------------------------------------------------------

pot_type = 'angle'

write_table('aba_long_angle_top', pot_type)
write_table('bab_long_angle_top', pot_type)

# # -----------------------------------------------------------------------------
# # DIHEDRALS
# # -----------------------------------------------------------------------------

pot_type = 'dihe'

write_table('long_dihe_top', pot_type)
write_table('out_dihe_top', pot_type)
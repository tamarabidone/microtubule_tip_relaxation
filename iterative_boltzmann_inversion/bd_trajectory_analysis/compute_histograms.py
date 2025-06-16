import numpy as np

def histogram_bootstrap(data, bins, number_of_experiments=10):
	np.random.shuffle(data)
	experiment_size = len(data)//number_of_experiments
	experiments = []
	for i in range(number_of_experiments):
		histogram_experiment, _ = np.histogram(data[i*experiment_size:(i+1)*experiment_size], bins = bins, density = True)
		experiments.append(histogram_experiment)
	bin_centers = 0.5 * ( bins[1:] + bins[:-1] )
	uncertainty = np.std(experiments, axis = 0)#/np.sqrt(number_of_experiments)
	return bin_centers, np.histogram(data, bins = bins, density = True)[0], uncertainty

def compute_histograms(label_tuples, hist_bins, time_cutoff, number_of_exps = 10):
	all_data_gdp = []
	all_data_gtp = []
	nbody = len(label_tuples[0])
	if nbody == 2: type_string = "dists"
	elif nbody == 3: type_string = "angles"
	elif nbody == 4: type_string = "dihes"
	else: 1/0
	for label_tuple in label_tuples:
		label_string = ('_{}'*nbody).format(*label_tuple)
		data_gdp = np.genfromtxt( 'data/{}/gdp{}.txt'.format(type_string, label_string) )
		data_gtp = np.genfromtxt( 'data/{}/gtp{}.txt'.format(type_string, label_string) )	
		all_data_gdp += list(data_gdp[time_cutoff:])
		all_data_gtp += list(data_gtp[time_cutoff:])
	bc, h_gdp, dh_gdp = histogram_bootstrap(all_data_gdp, hist_bins, number_of_exps)
	_, h_gtp, dh_gtp = histogram_bootstrap(all_data_gtp, hist_bins, number_of_exps)
	return bc, h_gdp, dh_gdp, h_gtp, dh_gtp

def write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp):
	np.savetxt('data/hist/gdp_{}.dat'.format(title), np.transpose(np.array([bc, h_gdp, dh_gdp])))
	np.savetxt('data/hist/gtp_{}.dat'.format(title), np.transpose(np.array([bc, h_gtp, dh_gtp])))

time_cutoff = 0

# -----------------------------------------------------------------------------
# DISTANCES
# -----------------------------------------------------------------------------

title = 'long_intra'
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_pairs = [("15","16"), ("13","14"), ("11","12"), ("09","10"), ("07","08"), ("05","06"), ("03","04"), ("01","02")]
long_label_pairs = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1])) for char in chars for long_numbers in long_number_pairs]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(long_label_pairs, np.linspace(0, 100, 501), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'long_inter'
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_pairs = [("14", "15"), ("12","13"), ("10","11"), ("08","09"), ("06","07"), ("04","05"), ("02","03")]
long_label_pairs = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1])) for char in chars for long_numbers in long_number_pairs]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(long_label_pairs, np.linspace(0, 100, 501), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'alpha_lat_inter_all'
char_pairs = [("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"), ("E", "F"), ("F", "G"), ("G", "H"), ("H", "I"), ("I", "J"), ("J", "K"), ("K", "L"), ("L", "M"), ("M", "N")]
lat_numbers = ["02", "04", "06", "08", "10", "12", "14"]
lat_label_pairs = [("{}{}".format(char_pair[0],lat_number), "{}{}".format(char_pair[1],lat_number)) for char_pair in char_pairs for lat_number in lat_numbers]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(lat_label_pairs, np.linspace(0, 400, 1001), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'beta_lat_inter_all'
char_pairs = [("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"), ("E", "F"), ("F", "G"), ("G", "H"), ("H", "I"), ("I", "J"), ("J", "K"), ("K", "L"), ("L", "M"), ("M", "N")]
lat_numbers = ["01", "03", "05", "07", "09", "11", "13", "15"]
lat_label_pairs = [("{}{}".format(char_pair[0],lat_number), "{}{}".format(char_pair[1],lat_number)) for char_pair in char_pairs for lat_number in lat_numbers]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(lat_label_pairs, np.linspace(0, 400, 1001), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'seam_lat_inter_all'
char_pairs = [("A", "N")]
lat_numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
lat_label_pairs = [("{}{}".format(char_pair[0],lat_numbers[j+3]), "{}{}".format(char_pair[1],lat_numbers[j])) for char_pair in char_pairs for j in range(len(lat_numbers)-3)]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(lat_label_pairs, np.linspace(0, 400, 1001), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'diag1_all'
diag1_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
diag1_label_pairs = [("{}{}".format(diag_chars[0],numbers[j]), "{}{}".format(diag_chars[1],numbers[j+1])) for diag_chars in diag1_char_pairs for j in range(len(numbers)-1)]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(diag1_label_pairs, np.linspace(0, 400, 1001), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'diag2_all'
diag2_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
diag2_label_pairs = [("{}{}".format(diag_chars[1],numbers[j]), "{}{}".format(diag_chars[0],numbers[j+1])) for diag_chars in diag2_char_pairs for j in range(len(numbers)-1)]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(diag2_label_pairs, np.linspace(0, 400, 1001), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

# -----------------------------------------------------------------------------
# ANGLES
# -----------------------------------------------------------------------------

title = 'aba_long_angle_top'
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("02","03","04")]
long_label_triples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(long_label_triples, np.arange(-0.5, 181.5), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'bab_long_angle_top'
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("01","02","03")]
long_label_triples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(long_label_triples, np.arange(-0.5, 181.5), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

# -----------------------------------------------------------------------------
# DIHEDRALS
# -----------------------------------------------------------------------------

title = 'long_dihe_top'
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_quadruples = [("04","03","02","01")]
long_label_quadruples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2]), "{}{}".format(char,long_numbers[3])) for char in chars for long_numbers in long_number_quadruples]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(long_label_quadruples, np.arange(0.5, 360.5, 5.0), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

title = 'out_dihe_top'
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("16","04","03"), ("16","03","02"), ("16","02","01")]
long_label_quadruples = [("000", "{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]
bc, h_gdp, dh_gdp, h_gtp, dh_gtp = compute_histograms(long_label_quadruples, np.arange(0.5, 360.5, 5.0), time_cutoff, number_of_exps = 10)
write_histograms(title, bc, h_gdp, dh_gdp, h_gtp, dh_gtp)

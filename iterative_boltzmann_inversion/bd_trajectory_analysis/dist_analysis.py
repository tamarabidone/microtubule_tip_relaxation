import numpy as np
import sys

which = sys.argv[1]

xyz_filename = "../bd_simulation/microtubule_{}.xyz".format(which)

def read_xyz(xyz_filename):
	counter = 0
	labels = []
	with open(xyz_filename, 'r') as xyz_file:
		number_of_beads_per_file = int( xyz_file.readline() )
		descriptive_line = xyz_file.readline().split(' ')
	trajectories = [ [] for i in range(number_of_beads_per_file) ]
	with open(xyz_filename, 'r') as xyz_file:
		for i, line in enumerate(xyz_file):
			timeframe_index = counter // number_of_beads_per_file
			bead_index = counter % number_of_beads_per_file
			line_unimportant = ( len( line.split() ) <= 1 ) or ( line.split()[0] == 'Atoms.' )
			if line_unimportant:
				continue
			else:
				if counter < number_of_beads_per_file:
					labels.append( line.split()[0] ) # here we will swap
				coords = np.array( [ float(line.split()[x]) for x in range(1, 4) ] )
				trajectories[bead_index].append( coords )
			counter += 1
	return labels, trajectories, descriptive_line

labels, trajectories, descriptive_line = read_xyz(xyz_filename)

new_labels = ["{}01", "{}02", "{}03", "{}04",
				"{}05", "{}06", "{}07", "{}08",
				"{}09", "{}10", "{}11", "{}12",
				"{}13", "{}14", "{}15", "{}16"]
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
new_labels_1 = [ x.format(y) for y in chars for x in new_labels ]
new_labels_1 += [ "000" ]

trajectories_dict = {}

for i in range(len(trajectories)):
	trajectories_dict[new_labels_1[i]] = np.array(trajectories[i])

dists = {}

def dist(x1, x2):
	return np.sqrt(np.sum((x1-x2)**2, axis = 1))

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_pairs = [("14","15"), ("12","13"), ("10","11"), ("08","09"), ("06","07"), ("04","05"), ("02","03")]
long_label_pairs = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1])) for char in chars for long_numbers in long_number_pairs]
for label_pair in long_label_pairs:
	label_i, label_j = label_pair
	if label_i == label_j: continue
	np.savetxt( 'data/dists/{}_{}_{}.txt'.format(which, label_i, label_j), dist(trajectories_dict[label_i], trajectories_dict[label_j]) )

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_pairs = [("15","16"), ("13","14"), ("11","12"), ("09","10"), ("07","08"), ("05","06"), ("03","04"), ("01","02")]
long_label_pairs = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1])) for char in chars for long_numbers in long_number_pairs]
for label_pair in long_label_pairs:
	label_i, label_j = label_pair
	if label_i == label_j: continue
	np.savetxt( 'data/dists/{}_{}_{}.txt'.format(which, label_i, label_j), dist(trajectories_dict[label_i], trajectories_dict[label_j]) )

lat_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
numbers = ["16", "15", "14", "13", "12", "11", "10", "09", "08", "07", "06", "05", "04", "03", "02", "01"]
lat_label_pairs = [("{}{}".format(lat_chars[0],number), "{}{}".format(lat_chars[1],number)) for lat_chars in lat_char_pairs for number in numbers]
for label_pair in lat_label_pairs:
	label_i, label_j = label_pair
	if label_i == label_j: continue
	np.savetxt( 'data/dists/{}_{}_{}.txt'.format(which, label_i, label_j), dist(trajectories_dict[label_i], trajectories_dict[label_j]) )

char_pairs = [("A","N")]
lat_numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
lat_label_pairs = [("{}{}".format(char_pair[0],lat_numbers[j+3]), "{}{}".format(char_pair[1],lat_numbers[j])) for char_pair in char_pairs for j in range(len(lat_numbers)-3)]
for label_pair in lat_label_pairs:
	label_i, label_j = label_pair
	if label_i == label_j: continue
	np.savetxt( 'data/dists/{}_{}_{}.txt'.format(which, label_i, label_j), dist(trajectories_dict[label_i], trajectories_dict[label_j]) )

diag_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N"), ("N","A")]
numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
diag_label_pairs = [("{}{}".format(diag_chars[i],numbers[j]), "{}{}".format(diag_chars[1-i],numbers[j+1])) for diag_chars in diag_char_pairs for j in range(len(numbers)-1) for i in range(2)]
diag_label_pairs += [("N01", "A05"), ("N02", "A06"), ("N03", "A07"), ("N04", "A08"),
		  			 ("N05", "A09"), ("N06", "A10"), ("N07", "A11"), ("N08", "A12"),
		  			 ("N09", "A13"), ("N10", "A14"), ("N11", "A15"), ("N12", "A16")]
diag_label_pairs += [("A03", "N01"), ("A04", "N02"), ("A05", "N03"), ("A06", "N04"),
		  			 ("A07", "N05"), ("A08", "N06"), ("A09", "N07"), ("A10", "N08"),
		  			 ("A11", "N09"), ("A12", "N10"), ("A13", "N11"), ("A14", "N12"),
		  			 ("A15", "N13"), ("A16", "N14")]
for diag_label_pair in diag_label_pairs:
	label_i, label_j = diag_label_pair
	if label_i == label_j: continue
	np.savetxt( 'data/dists/{}_{}_{}.txt'.format(which, label_i, label_j), dist(trajectories_dict[label_i], trajectories_dict[label_j]) )
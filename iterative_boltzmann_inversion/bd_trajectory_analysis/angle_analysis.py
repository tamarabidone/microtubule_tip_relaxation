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

angles = {}

def angle(x1, x2, x3):
	dist12 = np.sqrt(np.sum((x1-x2)**2, axis = 1))
	dist23 = np.sqrt(np.sum((x2-x3)**2, axis = 1))
	dots = np.array([np.dot(x2[i]-x1[i],x3[i]-x2[i])/dist12[i]/dist23[i] for i in range(len(x1))])
	angles = np.rad2deg( np.arccos( -dots ) )
	return angles

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("14","15","16"), ("13","14","15"), ("12","13","14"),
					   ("11","12","13"), ("10","11","12"), ("09","10","11"),
					   ("08","09","10"), ("07","08","09"), ("06","07","08"),
					   ("05","06","07"), ("04","05","06"), ("03","04","05"),
					   ("02","03","04"), ("01","02","03")]
long_label_triples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]
for label_triple in long_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = angle(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/angles/{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )

char_triples = [("A", "B", "C"), ("B", "C", "D"), ("C", "D", "E"), ("D", "E", "F"), ("E", "F", "G"), ("F", "G", "H"), ("G", "H", "I"), ("H", "I", "J"), ("I", "J", "K"), ("J", "K", "L"), ("K", "L", "M"), ("L", "M", "N")]
lat_numbers = ["16", "15", "14", "13", "12", "11", "10", "09", "08", "07", "06", "05", "04", "03", "02", "01"]
lat_label_triples = [("{}{}".format(char_triple[0],lat_number), "{}{}".format(char_triple[1],lat_number), "{}{}".format(char_triple[2],lat_number)) for char_triple in char_triples for lat_number in lat_numbers]
lat_label_triples += [("N13", "A16", "B16"), ("M13", "N13", "A16"), 
			 		  ("N12", "A15", "B15"), ("M12", "N12", "A15"),
			 		  ("N11", "A14", "B14"), ("M11", "N11", "A14"),
			 		  ("N10", "A13", "B13"), ("M10", "N10", "A13"),
			 		  ("N09", "A12", "B12"), ("M09", "N09", "A12"),
			 		  ("N08", "A11", "B11"), ("M08", "N08", "A11"),
			 		  ("N07", "A10", "B10"), ("M07", "N07", "A10"),
			 		  ("N06", "A09", "B09"), ("M06", "N06", "A09"),
			 		  ("N05", "A08", "B08"), ("M05", "N05", "A08"),
			 		  ("N04", "A07", "B07"), ("M04", "N04", "A07"),
			 		  ("N03", "A06", "B06"), ("M03", "N03", "A06"),
			 		  ("N02", "A05", "B05"), ("M02", "N02", "A05"),
			 		  ("N01", "A04", "B04"), ("M01", "N01", "A04")]
for label_triple in lat_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = angle(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/angles/{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )

char_triples = [("A", "B", "C"), ("B", "C", "D"), ("C", "D", "E"), ("D", "E", "F"), ("E", "F", "G"), ("F", "G", "H"), ("G", "H", "I"), ("H", "I", "J"), ("I", "J", "K"), ("J", "K", "L"), ("K", "L", "M"), ("L", "M", "N")]
lat_numbers = ["16", "15", "14", "13", "12", "11", "10", "09", "08", "07", "06", "05", "04", "03", "02", "01"]
lat_label_triples = [("{}{}".format(char_triple[0],lat_numbers[j]), "{}{}".format(char_triple[1],lat_numbers[j+1]), "{}{}".format(char_triple[2],lat_numbers[j])) for char_triple in char_triples for j in range(len(lat_numbers)-1)]
for label_triple in lat_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = angle(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/angles/{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("01","16","15")]
long_label_triples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]
for label_triple in long_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = angle(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/angles/{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )

char_triples = [("A", "A", "B"), ("B", "B", "C"), ("C", "C", "D"), ("D", "D", "E"), ("E", "E", "F"), ("F", "F", "G"), ("G", "G", "H"), ("H", "H", "I"), ("I", "I", "J"), ("J", "J", "K"), ("K", "K", "L"), ("L", "L", "M"), ("M", "M", "N")]
lat_numbers = ["16", "15", "14", "13", "12", "11", "10", "09", "08", "07", "06", "05", "04", "03", "02", "01"]
lat_label_triples = [("{}{}".format(char_triple[0],lat_numbers[j+1]), "{}{}".format(char_triple[1],lat_numbers[j]), "{}{}".format(char_triple[2],lat_numbers[j])) for char_triple in char_triples for j in range(len(lat_numbers)-1)]
lat_label_triples += [("N12", "N13", "A16"), 
			 		  ("N11", "N12", "A15"),
			 		  ("N10", "N11", "A14"),
			 		  ("N09", "N10", "A13"),
			 		  ("N08", "N09", "A12"),
			 		  ("N07", "N08", "A11"),
			 		  ("N06", "N07", "A10"),
			 		  ("N05", "N06", "A09"),
			 		  ("N04", "N05", "A08"),
			 		  ("N03", "N04", "A07"),
			 		  ("N02", "N03", "A06"),
			 		  ("N01", "N02", "A05")]
for label_triple in lat_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = angle(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/angles/{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )

char_triples = [("A", "B", "B"), ("B", "C", "C"), ("C", "D", "D"), ("D", "E", "E"), ("E", "F", "F"), ("F", "G", "G"), ("G", "H", "H"), ("H", "I", "I"), ("I", "J", "J"), ("J", "K", "K"), ("K", "L", "L"), ("L", "M", "M"), ("M", "N", "N")]
lat_numbers = ["16", "15", "14", "13", "12", "11", "10", "09", "08", "07", "06", "05", "04", "03", "02", "01"]
lat_label_triples = [("{}{}".format(char_triple[0],lat_numbers[j]), "{}{}".format(char_triple[1],lat_numbers[j]), "{}{}".format(char_triple[2],lat_numbers[j+1])) for char_triple in char_triples for j in range(len(lat_numbers)-1)]
lat_label_triples += [("N13", "A16", "A15"), 
			 		  ("N12", "A15", "A14"),
			 		  ("N11", "A14", "A13"),
			 		  ("N10", "A13", "A12"),
			 		  ("N09", "A12", "A11"),
			 		  ("N08", "A11", "A10"),
			 		  ("N07", "A10", "A09"),
			 		  ("N06", "A09", "A08"),
			 		  ("N05", "A08", "A07"),
			 		  ("N04", "A07", "A06"),
			 		  ("N03", "A06", "A05"),
			 		  ("N02", "A05", "A04")]
for label_triple in lat_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = angle(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/angles/{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )
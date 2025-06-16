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

dihedrals = {}

def dihedral(x1, x2, x3, x4):
	# return np.atan2()
	dist12 = np.sqrt(np.sum((x1-x2)**2, axis = 1))
	dist23 = np.sqrt(np.sum((x2-x3)**2, axis = 1))
	dist34 = np.sqrt(np.sum((x3-x4)**2, axis = 1))
	arg1 = np.array([dist23[i]*np.dot(x2[i]-x1[i], np.cross(x3[i]-x2[i], x4[i]-x3[i])) for i in range(len(x1))])
	arg2 = np.array([np.dot(np.cross(x2[i]-x1[i], x3[i]-x2[i]), np.cross(x3[i]-x2[i], x4[i]-x3[i])) for i in range(len(x1))])
	dihedrals = np.rad2deg( np.arctan2( arg1, arg2 ) )
	dihedrals[dihedrals<0.0]+=360.0
	return dihedrals

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_quadruples = [("04","03","02","01"), ("05","04","03","02"),
						  ("06","05","04","03"), ("07","06","05","04"),
						  ("08","07","06","05"), ("09","08","07","06"),
						  ("10","09","08","07"), ("11","10","09","08"),
						  ("12","11","10","09"), ("13","12","11","10"),
						  ("14","13","12","11"), ("15","14","13","12"),
						  ("16","15","14","13")]
long_label_quadruples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2]), "{}{}".format(char,long_numbers[3])) for char in chars for long_numbers in long_number_quadruples]
for label_quadruple in long_label_quadruples:
	label_i, label_j, label_k, label_l = label_quadruple
	this_angle = dihedral(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k], trajectories_dict[label_l])
	np.savetxt( 'data/dihes/{}_{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k, label_l), this_angle )


lat_number_quadruples = [("01","02","02","02"),
						  ("02","03","03","03"),
						  ("03","04","04","04"),
						  ("04","05","05","05"),
						  ("05","06","06","06"),
						  ("06","07","07","07"),
						  ("07","08","08","08"),
						  ("08","09","09","09"),
						  ("09","10","10","10"),
						  ("10","11","11","11"),
						  ("11","12","12","12"),
						  ("12","13","13","13"),
						  ("13","14","14","14"),
						  ("14","15","15","15"),
						  ("15","16","16","16")]
lat_label_quadruples = [("{}{}".format(chars[j+1],lat_numbers[0]), "{}{}".format(chars[j],lat_numbers[1]), "{}{}".format(chars[j+2],lat_numbers[2]), "{}{}".format(chars[j+1],lat_numbers[3])) for j in range(len(chars)-2) for lat_numbers in lat_number_quadruples]
lat_label_quadruples += [("A15", "N13", "B16", "A16"), ("N12", "M13", "A16", "N13"), 
			 ("A14", "N12", "B15", "A15"), ("N11", "M12", "A15", "N12"),
			 ("A13", "N11", "B14", "A14"), ("N10", "M11", "A14", "N11"),
			 ("A12", "N10", "B13", "A13"), ("N09", "M10", "A13", "N10"),
			 ("A11", "N09", "B12", "A12"), ("N08", "M09", "A12", "N09"),
			 ("A10", "N08", "B11", "A11"), ("N07", "M08", "A11", "N08"),
			 ("A09", "N07", "B10", "A10"), ("N06", "M07", "A10", "N07"),
			 ("A08", "N06", "B09", "A09"), ("N05", "M06", "A09", "N06"),
			 ("A07", "N05", "B08", "A08"), ("N04", "M05", "A08", "N05"),
			 ("A06", "N04", "B07", "A07"), ("N03", "M04", "A07", "N04"),
			 ("A05", "N03", "B06", "A06"), ("N02", "M03", "A06", "N03"),
			 ("A04", "N02", "B05", "A05"), ("N01", "M02", "A05", "N02"),
			 ("A03", "N01", "B04", "A04")]
for label_quadruple in lat_label_quadruples:
	label_i, label_j, label_k, label_l = label_quadruple
	# if label_i == label_j or label_i == label_k or label_j == label_k: continue
	this_angle = dihedral(trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k], trajectories_dict[label_l])
	np.savetxt( 'data/dihes/{}_{}_{}_{}_{}.txt'.format(which, label_i, label_j, label_k, label_l), this_angle )

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("16","15","14"), ("16","14","13"), ("16","13","12"),
					   ("16","12","11"), ("10","11","10"), ("16","10","09"),
					   ("16","09","08"), ("16","08","07"), ("16","07","06"),
					   ("16","06","05"), ("16","05","04"), ("16","04","03"),
					   ("16","03","02"), ("16","02","01")]
long_label_triples = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]
for label_triple in long_label_triples:
	label_i, label_j, label_k = label_triple
	this_angle = dihedral(trajectories_dict["000"], trajectories_dict[label_i], trajectories_dict[label_j], trajectories_dict[label_k])
	np.savetxt( 'data/dihes/{}_000_{}_{}_{}.txt'.format(which, label_i, label_j, label_k), this_angle )
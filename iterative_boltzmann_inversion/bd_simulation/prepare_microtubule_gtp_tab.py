import numpy as np

def number_from_label(label, length = 16):
	if label == "000": return 225
	protofilament_id = ord(label[0])-64
	layer_id = int(label[1:])
	return (protofilament_id-1)*length + layer_id

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
			line_unimportant = ( len( line.split() ) <= 1 ) or ( 'frame' in line.split()[0].split('.') )
			if line_unimportant:
				continue
			else:
				if counter < number_of_beads_per_file:
					labels.append( line.split()[0] ) # here we will swap
				else:
					return labels, trajectories, descriptive_line
				coords = np.array( [ float(line.split()[x]) for x in range(1, 4) ] )
				trajectories[bead_index].append( coords )
			counter += 1

def attribute_atom_type(label):
	if label == "000": return 3
	layer_index = int(label[1:])
	if layer_index == 16: return 3
	if layer_index % 2 == 1: return 2
	else: return 1

xyz_filename = "../../md_trajectory_analysis/gtp-cacb-center-caonly-align-cg-nojump.xyz"

import MDAnalysis as mda
mt = mda.Universe(xyz_filename)
mt_traj = mt.select_atoms('all')
mt.trajectory[20000]
final_positions = mt_traj.positions.copy()

labels, trajectories, descriptive_line = read_xyz(xyz_filename)
for element in range(len(trajectories)):
	trajectories[element][0] = final_positions[element]

md_labels =    ["A01", "A02", "A03", "A04",
				"B01", "B02", "B03", "B04",
				"C01", "C02", "C03", "C04",
				"A05", "A06", "A07", "A08",
				"B05", "B06", "B07", "B08",
				"C05", "C06", "C07", "C08",
				"A09", "A10", "A11", "A12",
				"B09", "B10", "B11", "B12",
				"C09", "C10", "C11", "C12",
				"A13", "A14", "A15", "A16",
				"B13", "B14", "B15", "B16",
				"C13", "C14", "C15", "C16",
				"D01", "D02", "D03", "D04",
				"E01", "E02", "E03", "E04",
				"F01", "F02", "F03", "F04",
				"D05", "D06", "D07", "D08",
				"E05", "E06", "E07", "E08",
				"F05", "F06", "F07", "F08",
				"D09", "D10", "D11", "D12",
				"E09", "E10", "E11", "E12",
				"F09", "F10", "F11", "F12",
				"D13", "D14", "D15", "D16",
				"E13", "E14", "E15", "E16",
				"F13", "F14", "F15", "F16",
				"G01", "G02", "G03", "G04",
				"H01", "H02", "H03", "H04",
				"I01", "I02", "I03", "I04",
				"G05", "G06", "G07", "G08",
				"H05", "H06", "H07", "H08",
				"I05", "I06", "I07", "I08",
				"G09", "G10", "G11", "G12",
				"H09", "H10", "H11", "H12",
				"I09", "I10", "I11", "I12",
				"G13", "G14", "G15", "G16",
				"H13", "H14", "H15", "H16",
				"I13", "I14", "I15", "I16",
				"J01", "J02", "J03", "J04",
				"K01", "K02", "K03", "K04",
				"L01", "L02", "L03", "L04",
				"J05", "J06", "J07", "J08",
				"K05", "K06", "K07", "K08",
				"L05", "L06", "L07", "L08",
				"J09", "J10", "J11", "J12",
				"K09", "K10", "K11", "K12",
				"L09", "L10", "L11", "L12",
				"J13", "J14", "J15", "J16",
				"K13", "K14", "K15", "K16",
				"L13", "L14", "L15", "L16",
				"M01", "M02", "M03", "M04",
				"N01", "N02", "N03", "N04",
				"M05", "M06", "M07", "M08",
				"N05", "N06", "N07", "N08",
				"M09", "M10", "M11", "M12",
				"N09", "N10", "N11", "N12",
				"M13", "M14", "M15", "M16",
				"N13", "N14", "N15", "N16"]

bd_labels_temp = ["{}01", "{}02", "{}03", "{}04",
				  "{}05", "{}06", "{}07", "{}08",
				  "{}09", "{}10", "{}11", "{}12",
				  "{}13", "{}14", "{}15", "{}16"]
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
bd_labels = [ x.format(y) for y in chars for x in bd_labels_temp ]

initial_positions = {}

for i in range(len(trajectories)):
	initial_positions[md_labels[i]] = np.array(trajectories[i])

base_labels = ["A16", "B16", "C16", "D16", "E16", "F16", "G16", "H16", "I16", "J16", "K16", "L16", "M16", "N16"]
initial_positions["000"] = np.mean( [ initial_positions[label] for label in base_labels ], axis = 0 )
bd_labels.append("000")

### ----------------------------------------------------------------------

def return_header(description):
	return description+'\n'

def return_atoms_bonds_angles(number_of_atoms, number_of_bonds, number_of_angles, number_of_dihedrals):
	return '{} atoms\n{} bonds\n{} angles\n{} dihedrals\n'.format(number_of_atoms, number_of_bonds, number_of_angles, number_of_dihedrals)

def return_atom_bond_angle_types(number_of_atom_types, number_of_bond_types, number_of_angle_types, number_of_dihedral_types):
	return '{} atom types\n{} bond types\n{} angle types\n{} dihedral types\n'.format(number_of_atom_types, number_of_bond_types, number_of_angle_types, number_of_dihedral_types)

def return_box(box_length):
	return '-{} {} xlo xhi\n-{} {} ylo yhi\n-{} {} zlo zhi\n'.format(*[box_length/2]*6)

def return_masses(number_of_atom_types):
	result = 'Masses\n\n'
	for j in range(number_of_atom_types):
		result += '  {} 1.0\n'.format(j+1)
	return result

def return_pair_coeffs(number_of_atom_types, pair_coeff_lines):
	result = 'Pair Coeffs\n\n'
	for j in range(number_of_atom_types):
		result += '  {} {}\n'.format(j+1, pair_coeff_lines[j])
	return result

def return_bond_coeffs(number_of_bond_types, bond_coeff_lines):
	result = 'Bond Coeffs\n\n'
	for j in range(number_of_bond_types):
		result += '  {} {}\n'.format(j+1, bond_coeff_lines[j])
	return result

def return_angle_coeffs(number_of_angle_types, angle_coeff_lines):
	result = 'Angle Coeffs\n\n'
	for j in range(number_of_angle_types):
		result += '  {} {}\n'.format(j+1, angle_coeff_lines[j])
	return result

def return_dihedral_coeffs(number_of_dihedral_types, dihedral_coeff_lines):
	result = 'Dihedral Coeffs\n\n'
	for j in range(number_of_dihedral_types):
		result += '  {} {}\n'.format(j+1, dihedral_coeff_lines[j])
	return result

def return_atoms(bd_labels, initial_positions):
	result = 'Atoms\n\n'
	for j, bd_label in enumerate(bd_labels):
		result += '  {} 1 {} {} {} {}\n'.format(j+1, attribute_atom_type(bd_label), *initial_positions[bd_label][0])
	return result

def return_bonds(bond_definitions, number_of_bond_types):
	assert number_of_bond_types == len(bond_definitions)
	result = 'Bonds\n\n'
	counter = 0
	for j, bond_definition in enumerate(bond_definitions):
		for label_pair in bond_definition:
			label_1, label_2 = label_pair
			result += '  {} {} {} {}\n'.format(counter+1, j+1, number_from_label(label_1), number_from_label(label_2))
			counter += 1
	return result

def return_angles(angle_definitions, number_of_angle_types):
	assert number_of_angle_types == len(angle_definitions)
	result = 'Angles\n\n'
	counter = 0
	for j, angle_definition in enumerate(angle_definitions):
		for label_triple in angle_definition:
			label_1, label_2, label_3 = label_triple
			result += '  {} {} {} {} {}\n'.format(counter+1, j+1, number_from_label(label_1), number_from_label(label_2), number_from_label(label_3))
			counter += 1
	return result

def return_dihedrals(dihedral_definitions, number_of_dihedral_types):
	assert number_of_dihedral_types == len(dihedral_definitions)
	result = 'Dihedrals\n\n'
	counter = 0
	for j, dihedral_definition in enumerate(dihedral_definitions):
		for label_quadruple in dihedral_definition:
			label_1, label_2, label_3, label_4 = label_quadruple
			result += '  {} {} {} {} {} {}\n'.format(counter+1, j+1, number_from_label(label_1), number_from_label(label_2), number_from_label(label_3), number_from_label(label_4))
			counter += 1
	return result

### ----------------------------------------------------------------------

number_of_atom_types = 3
pair_coeff_lines = ["0.39 54.0 56.0"]*number_of_atom_types
number_of_atoms = len(bd_labels)

### ----------------------------------------------------------------------

long_intra_chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_intra_number_pairs = [("01","02"), ("03","04"), ("05","06"), ("07","08"), ("09","10"), ("11","12"), ("13","14"), ("15","16")]
long_intra = [("{}{}".format(char,long_intra_numbers[0]), "{}{}".format(char,long_intra_numbers[1])) for char in long_intra_chars for long_intra_numbers in long_intra_number_pairs]

long_inter_chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_inter_number_pairs = [("02","03"), ("04","05"), ("06","07"), ("08","09"), ("10","11"), ("12","13"), ("14","15")]
long_inter = [("{}{}".format(char,long_inter_numbers[0]), "{}{}".format(char,long_inter_numbers[1])) for char in long_inter_chars for long_inter_numbers in long_inter_number_pairs]

lat_beta_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
lat_beta_numbers = ["01", "03", "05", "07", "09", "11", "13", "15"]
lat_beta = [("{}{}".format(lat_beta_chars[0],lat_beta_number), "{}{}".format(lat_beta_chars[1],lat_beta_number)) for lat_beta_chars in lat_beta_char_pairs for lat_beta_number in lat_beta_numbers]

lat_alpha_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
lat_alpha_numbers = ["02", "04", "06", "08", "10", "12", "14", "16"]
lat_alpha = [("{}{}".format(lat_alpha_chars[0],lat_alpha_number), "{}{}".format(lat_alpha_chars[1],lat_alpha_number)) for lat_alpha_chars in lat_alpha_char_pairs for lat_alpha_number in lat_alpha_numbers]

lat_seam_char_pairs = [("A", "N")]
lat_seam_numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
lat_seam = [("{}{}".format(lat_seam_char_pair[0],lat_seam_numbers[j+3]), "{}{}".format(lat_seam_char_pair[1],lat_seam_numbers[j])) for lat_seam_char_pair in lat_seam_char_pairs for j in range(len(lat_seam_numbers)-3)]

diag1_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
diag1 = [("{}{}".format(diag_chars[0],numbers[j]), "{}{}".format(diag_chars[1],numbers[j+1])) for diag_chars in diag1_char_pairs for j in range(len(numbers)-1)]
diag1 += [("N01", "A05"), ("N02", "A06"), ("N03", "A07"), ("N04", "A08"),
		  ("N05", "A09"), ("N06", "A10"), ("N07", "A11"), ("N08", "A12"),
		  ("N09", "A13"), ("N10", "A14"), ("N11", "A15"), ("N12", "A16")]

diag2_char_pairs = [("A","B"), ("B","C"), ("C","D"), ("D","E"), ("E","F"), ("F","G"), ("G","H"), ("H","I"), ("I","J"), ("J","K"), ("K","L"), ("L","M"), ("M","N")]
numbers = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
diag2 = [("{}{}".format(diag_chars[1],numbers[j]), "{}{}".format(diag_chars[0],numbers[j+1])) for diag_chars in diag2_char_pairs for j in range(len(numbers)-1)]
diag2 += [("A03", "N01"), ("A04", "N02"), ("A05", "N03"), ("A06", "N04"),
		  ("A07", "N05"), ("A08", "N06"), ("A09", "N07"), ("A10", "N08"),
		  ("A11", "N09"), ("A12", "N10"), ("A13", "N11"), ("A14", "N12"),
		  ("A15", "N13"), ("A16", "N14")]

bond_definitions = [ long_intra, long_inter, lat_beta, lat_alpha, lat_seam, diag1, diag2 ]
bond_coeff_lines = ["potentials/gtp_long_intra.table ENTRY", "potentials/gtp_long_inter.table ENTRY", "potentials/gtp_beta_lat_inter_all.table ENTRY", "potentials/gtp_alpha_lat_inter_all.table ENTRY", "potentials/gtp_seam_lat_inter_all.table ENTRY", "potentials/gtp_diag1_all.table ENTRY", "potentials/gtp_diag2_all.table ENTRY"]

number_of_bond_types = len(bond_definitions)
number_of_bonds = np.sum([ len(b) for b in bond_definitions ])

### ----------------------------------------------------------------------

long_aba_chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_aba_numbers = [("02","03","04"), ("04","05","06"), ("06","07","08"),
					("08","09","10"), ("10","11","12"), ("12","13","14"),
					("14","15","16")]
long_aba = [("{}{}".format(long_aba_char,long_aba_number[0]), "{}{}".format(long_aba_char,long_aba_number[1]), "{}{}".format(long_aba_char,long_aba_number[2])) for long_aba_char in long_aba_chars for long_aba_number in long_aba_numbers]

long_bab_chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_bab_numbers = [("01","02","03"), ("03","04","05"), ("05","06","07"),
					("07","08","09"), ("09","10","11"), ("11","12","13"),
					("13","14","15")]
long_bab = [("{}{}".format(long_bab_char,long_bab_number[0]), "{}{}".format(long_bab_char,long_bab_number[1]), "{}{}".format(long_bab_char,long_bab_number[2])) for long_bab_char in long_bab_chars for long_bab_number in long_bab_numbers]

angle_definitions = [ long_aba, long_bab ]
angle_coeff_lines = [ "potentials/gtp_aba_long_angle_top.table ENTRY", "potentials/gtp_bab_long_angle_top.table ENTRY" ]

number_of_angle_types = len(angle_definitions)
number_of_angles = np.sum([ len(a) for a in angle_definitions ])

### ----------------------------------------------------------------------

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_quadruples = [("04","03","02","01"),
						  ("05","04","03","02"),
						  ("06","05","04","03"),
						  ("07","06","05","04"),
						  ("08","07","06","05"),
						  ("09","08","07","06"),
						  ("10","09","08","07"),
						  ("11","10","09","08"),
						  ("12","11","10","09"),
						  ("13","12","11","10"),
						  ("14","13","12","11"),
						  ("15","14","13","12"),
						  ("16","15","14","13")]
dihe = [("{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2]), "{}{}".format(char,long_numbers[3])) for char in chars for long_numbers in long_number_quadruples]

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
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
lat_dihe = [("{}{}".format(chars[j+1],lat_numbers[0]), "{}{}".format(chars[j],lat_numbers[1]), "{}{}".format(chars[j+2],lat_numbers[2]), "{}{}".format(chars[j+1],lat_numbers[3])) for j in range(len(chars)-2) for lat_numbers in lat_number_quadruples]
lat_dihe += [("A15", "N13", "B16", "A16"), ("N12", "M13", "A16", "N13"), 
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

chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
long_number_triples = [("16","15","14"), ("16","14","13"), ("16","13","12"),
					   ("16","12","11"), ("16","11","10"), ("16","10","09"),
					   ("16","09","08"), ("16","08","07"), ("16","07","06"),
					   ("16","06","05"), ("16","05","04"), ("16","04","03"),
					   ("16","03","02"), ("16","02","01")]
dihe_out = [("000", "{}{}".format(char,long_numbers[0]), "{}{}".format(char,long_numbers[1]), "{}{}".format(char,long_numbers[2])) for char in chars for long_numbers in long_number_triples]

dihe_definitions = [ dihe, dihe_out ]
dihe_coeff_lines = [ "aat 1.0 177.0 180.0 potentials/gtp_long_dihe_top.table ENTRY", "aat 1.0 179.0 180.0 potentials/gtp_out_dihe_top.table ENTRY" ]

number_of_dihe_types = len(dihe_definitions)
number_of_dihes = np.sum([ len(a) for a in dihe_definitions ])

### ----------------------------------------------------------------------

print(return_header('LAMMPS Description'))
print(return_atoms_bonds_angles(number_of_atoms = number_of_atoms, number_of_bonds = number_of_bonds, number_of_angles = number_of_angles, number_of_dihedrals = number_of_dihes))
print(return_atom_bond_angle_types(number_of_atom_types = number_of_atom_types, number_of_bond_types = number_of_bond_types, number_of_angle_types = number_of_angle_types, number_of_dihedral_types = number_of_dihe_types))
print(return_box(box_length = 10000.0))
print(return_masses(number_of_atom_types = number_of_atom_types))
# print(return_pair_coeffs(number_of_atom_types = number_of_atom_types, pair_coeff_lines = pair_coeff_lines))
print(return_bond_coeffs(number_of_bond_types = number_of_bond_types, bond_coeff_lines = bond_coeff_lines))
print(return_angle_coeffs(number_of_angle_types = number_of_angle_types, angle_coeff_lines = angle_coeff_lines))
print(return_dihedral_coeffs(number_of_dihedral_types = number_of_dihe_types, dihedral_coeff_lines = dihe_coeff_lines))
print(return_atoms(bd_labels = bd_labels, initial_positions = initial_positions))
print(return_bonds(bond_definitions = bond_definitions, number_of_bond_types = number_of_bond_types))
print(return_angles(angle_definitions = angle_definitions, number_of_angle_types = number_of_angle_types))
print(return_dihedrals(dihedral_definitions = dihe_definitions, number_of_dihedral_types = number_of_dihe_types))
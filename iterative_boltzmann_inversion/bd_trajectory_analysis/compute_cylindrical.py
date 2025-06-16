import numpy as np
import MDAnalysis as mda
import sys

which = sys.argv[1]

xyz_filename_template = "../bd_simulation/microtubule_{}.xyz"

new_labels = ["{}01", "{}02", "{}03", "{}04",
				"{}05", "{}06", "{}07", "{}08",
				"{}09", "{}10", "{}11", "{}12",
				"{}13", "{}14", "{}15", "{}16"]
chars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
new_labels_1 = [ x.format(y) for y in chars for x in new_labels ]
new_labels_1 += [ "000" ]
labels_dict = { new_labels_1[i]: i for i in range(len(new_labels_1)) }
char_labels = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]

for which in ["gdp", "gtp"]:
	xyz_filename = xyz_filename_template.format(which)
	uni = mda.Universe(xyz_filename)
	traj = uni.select_atoms('all')
	uni.trajectory[0]
	frame0 = traj.positions.copy()

	main_axis = np.zeros(3)
	center = np.zeros(2)
	for char_label in char_labels:
		versor = frame0[labels_dict["{}15".format(char_label)]]-frame0[labels_dict["{}16".format(char_label)]]
		versor /= np.linalg.norm(versor)
		main_axis += versor
		center += frame0[labels_dict["{}16".format(char_label)]][:-1]
	main_axis/=np.linalg.norm(main_axis)
	center/=14
	maxz = np.max(frame0[:,2])

	newz = main_axis
	newz/=np.linalg.norm(newz)
	newx = (frame0[labels_dict["A16"]]+frame0[labels_dict["N16"]])/2
	newx[:-1] -= center
	newx[2] = 0.0
	newx/=np.linalg.norm(newx)
	newx -= np.dot(newx, newz)*newz
	newx/=np.linalg.norm(newx)
	newy = np.cross(newz, newx)
	newy/=np.linalg.norm(newy)
	transformation_matrix = np.array([newx, newy, newz])

	rs = { i:np.zeros(len(uni.trajectory)) for i in new_labels_1 }
	phis = { i:np.zeros(len(uni.trajectory)) for i in new_labels_1 }
	zs = { i:np.zeros(len(uni.trajectory)) for i in new_labels_1 }

	for i in range(len(uni.trajectory)):
		uni.trajectory[i]
		framei = traj.positions.copy()

		for j, coord in enumerate(framei):
			coord[:-1] -= center
			coord[-1] -= maxz
			coord = transformation_matrix@coord
			new_label = new_labels_1[j]
			r = np.sqrt(coord[0]**2+coord[1]**2)
			if coord[1] >= 0: phi = np.rad2deg(np.arccos(coord[0]/r))
			else: phi = 360-np.rad2deg(np.arccos(coord[0]/r))
			z = coord[2]
			rs[new_label][i] = r
			phis[new_label][i] = phi
			zs[new_label][i] = z

	for label in new_labels_1:
		np.savetxt( 'data/cylindrical/{}_{}_r.txt'.format(which, label), rs[label] )
		np.savetxt( 'data/cylindrical/{}_{}_phi.txt'.format(which, label), phis[label] )
		np.savetxt( 'data/cylindrical/{}_{}_z.txt'.format(which, label), zs[label] )
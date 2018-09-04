import numpy as np

x = float(open(".x.dat", "r").read())
y = float(open(".y.dat", "r").read())
z = float(open(".toppavsylinder.dat", "r").read())

data = np.loadtxt("./data/medvann.xyz", unpack=True, skiprows=2)
number_of_atoms = len(data[0])

new_data = np.column_stack((np.arange(1, number_of_atoms+1), *data))

np.savetxt("./data/medvann.data", new_data, fmt="%d %d %g %g %g", header="""# medvann
%d atoms
3 atom types
0 %g xlo xhi
0 %g ylo yhi
0 %g zlo zhi

Atoms # atomic
""" % (number_of_atoms, x, y, z), comments="")

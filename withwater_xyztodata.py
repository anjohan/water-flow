from extract_variable import extract_variable
import numpy as np

x = extract_variable("x")
y = extract_variable("y")
z = extract_variable("z_topofcylinder")

data = np.loadtxt("./data/withwater.xyz", unpack=True, skiprows=2)
number_of_atoms = len(data[0])

new_data = np.column_stack((np.arange(1, number_of_atoms + 1), *data))

np.savetxt(
    "./data/withwater.data",
    new_data,
    fmt="%d %d %g %g %g",
    header="""# withwater
%d atoms
3 atom types
0 %g xlo xhi
0 %g ylo yhi
0 %g zlo zhi

Atoms # atomic
""" % (number_of_atoms, x, y, z),
    comments="")

from extract_variable import extract_variable
from math import pi

z_bottom = extract_variable("z_bottom")
z_topofslab = extract_variable("z_topofslab")
z_topofcylinder = extract_variable("z_topofcylinder")
x = extract_variable("x")
y = extract_variable("y")
R = extract_variable("cylinderradius")
h = extract_variable("cylinderheight")

density = 1.0  # g/cm3
massperwater = 18 * 1.67E-24  # g/molecule
volume = (x * y - pi * R**2) * h * 1E-24  # cm3
number_of_molecules = volume * density / massperwater  # molecules
number_of_molecules = int(round(number_of_molecules))

with open("./packmol_template.inp", "r") as fil:
    txt = fil.read()

replacements = {
    "NUMBEROFMOLECULES": number_of_molecules,
    "MIDx": x / 2,
    "MIDy": y / 2,
    "WATERBOTTOM": z_topofslab,
    "WATERTOP": z_topofcylinder,
    "RADIUS": R,
    "HEIGHT": h,
    "X": x,
    "Y": y
}

for word, value in replacements.items():
    txt = txt.replace(word, str(value))

with open("water.inp", "w") as fil:
    fil.write(txt)

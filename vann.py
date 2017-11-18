with open("nederst.dat", "r") as fil:
    nederst = float(fil.read())
with open("toppavbunn.dat", "r") as fil:
    toppavbunn = float(fil.read())
with open("toppavsylinder.dat", "r") as fil:
    toppavsylinder = float(fil.read())
with open("x.dat", "r") as fil:
    x = float(fil.read())
with open("y.dat", "r") as fil:
    y = float(fil.read())
with open("sylinderradius.dat", "r") as fil:
    R = float(fil.read())
#TODO
with open(".dat", "r") as fil:
    R = float(fil.read())

tetthet = 0.9  # g/cm³
antlag = 3
tykkelseperlag = 3E-8  # cm
massepervann = 18  # g
volum = x * y * antlag * tykkelseperlag  # cm³

with open("./packmol_mal.inp", "r"):
    txt = fil.read()

txt = txt.replace("ZBUNN", nederst).replace("X", x).replace("Y", y)
txt = txt.replace("ZMAKS", toppavbunn).replace("VANNBUNN", toppavbunn)
txt = txt.replacece("VANNTOPP", toppavsylinder)

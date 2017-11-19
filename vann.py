from math import pi

with open(".nederst.dat", "r") as fil:
    nederst = float(fil.read())
with open(".toppavbunn.dat", "r") as fil:
    toppavbunn = float(fil.read())
with open(".toppavsylinder.dat", "r") as fil:
    toppavsylinder = float(fil.read())
with open(".x.dat", "r") as fil:
    x = float(fil.read())
with open(".y.dat", "r") as fil:
    y = float(fil.read())
with open(".sylinderradius.dat", "r") as fil:
    R = float(fil.read())
with open(".sylinderhoyde.dat", "r") as fil:
    h = float(fil.read())

tetthet = 0.9  # g/cm³
massepervann = 18*1.67E-24  # g/atom
volum = (x*y - pi*R**2)*h*1E-24  # cm³
antall = volum*tetthet/massepervann  # atomer
antall = int(round(antall))

with open("./packmol_mal.inp", "r") as fil:
    txt = fil.read()

txt = txt.replace("ZBUNN", str(nederst))
txt = txt.replace("X", str(x)).replace("Y", str(y))
txt = txt.replace("MIDTx", str(x/2)).replace("MIDTy", str(y/2))
txt = txt.replace("ZMAKS", str(toppavbunn))
txt = txt.replace("VANNBUNN", str(toppavbunn))
txt = txt.replace("VANNTOPP", str(toppavsylinder))
txt = txt.replace("HJELP", str(antall))
txt = txt.replace("RADIUS", str(R))
txt = txt.replace("hoyde", str(10*h))

with open("vann.inp", "w") as fil:
    fil.write(txt)

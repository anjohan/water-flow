with open("./data/unpassivated.data", "r") as infile:
    txt = infile.read()
    txt = txt.replace("2 atom types", "3 atom types")
    txt = txt.replace("2 15.9994", "2 15.9994\n3 1.00794")
with open("./data/unpassivatedwithH.data", "w") as outfile:
    outfile.write(txt)

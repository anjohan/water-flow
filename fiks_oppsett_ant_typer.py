with open("./data/klartilpassivering.data", "r") as fil:
    txt = fil.read()
    txt = txt.replace("2 atom types", "3 atom types")
    txt = txt.replace("2 15.9994", "2 15.9994\n3 1.00794")
with open("./data/klartilpassivering.data", "w") as fil:
    fil.write(txt)

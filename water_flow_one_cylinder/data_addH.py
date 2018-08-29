import sys
infilename, outfilename = sys.argv[1:]

with open(infilename, "r") as infile:
    txt = infile.read()
    txt = txt.replace("2 atom types", "3 atom types")
with open(outfilename, "w") as outfile:
    outfile.write(txt)

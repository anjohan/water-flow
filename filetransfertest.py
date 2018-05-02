import sys
import os
output_commands = [
    "print", "dump", "log", "restart", "write_dump", "write_restart",
    "pair_write", "write_coeff", "ave/time", "ave/chunk", "bond_write",
    "ave/histo/weight", "saed/vtk", "ave/histo", "write_data"
]
with open(sys.argv[1], "r") as infile:
    for line in infile:
        words = line.split()
        if set(words).isdisjoint(output_commands):  # Not an output command
            for word in words:
                if os.path.isfile(word):
                    print(word)

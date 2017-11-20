from ovito.io import import_file, export_file

data = import_file("./data/medvann.xyz",
                   columns=["Particle Type", "Position.X",
                            "Position.Y", "Position.Z"])

export_file(data, "data/medvann.data", "lammps/data")

from argparse import ArgumentParser
import time
from ovito.io import import_file
from ovito.pipeline import FileSource
import ovito.modifiers as mods
import matplotlib.pyplot as plt
import numpy as np

parser = ArgumentParser()
parser.add_argument("-i", default="data/vann.in.bin", help="Datafile to read")
parser.add_argument("-d", default=False, action="store_true", help="Calc. displacements")
parser.add_argument("-o", default="./data/displ_vs_R.png", help="Plot output")
args = parser.parse_args()

pipeline = import_file(args.i, multiple_frames=True,
                       columns=["Particle Identifier", "Particle Type",
                                "Position.X", "Position.Y", "Position.Z"])

pipeline.add_to_scene()

num_frames = pipeline.source.num_frames

x = float(open(".x.dat", "r").read())
y = float(open(".y.dat", "r").read())
z = float(open(".toppavbunn.dat", "r").read())
R = float(open(".sylinderradius.dat", "r").read())
d = min(x/2, y/2)


"""
moving = mods.AffineTransformationModifier(transformation=[[1, 0, 0, -x/2],
                                                           [0, 1, 0, -y/2],
                                                           [0, 0, 1, 0]])
pipeline.modifiers.append(moving)
"""
expr = "sqrt((Position.X - %g)*(Position.X - %g)\
             + (Position.Y - %g)*(Position.Y - %g))" % (x/2, x/2, y/2, y/2)

dist = mods.ComputePropertyModifier(output_property='DistanceFromCenter',
                                    expressions=[expr])
pipeline.modifiers.append(dist)

expr = "DistanceFromCenter < %g || Position.Z < %g" % (d, z)
near_cylinder = mods.ExpressionSelectionModifier(expression=expr)
pipeline.modifiers.append(near_cylinder)

pipeline.modifiers.append(mods.InvertSelectionModifier())
pipeline.modifiers.append(mods.DeleteSelectedModifier())


def set_colors(data):
    si = data.particle_properties["Particle Type"].types[0]
    o = data.particle_properties["Particle Type"].types[1]
    h = data.particle_properties["Particle Type"].types[2]
    si.name = "Si"
    o.name = "O"
    h.name = "H"
    si.radius = 1.18
    o.radius = 0.77
    si.radius = 1.18
    o.radius = 0.77
    h.radius = 0.4
    si.color = (1, 1, 0)
    o.color = (1, 0, 0)
    h.color = (1, 1, 1)


set_colors(pipeline.source)

tracking = mods.CalculateDisplacementsModifier(use_frame_offset=True,
                                               frame_offset=-20)
tracking.reference = FileSource()
tracking.reference.load(args.i, multiple_frames=True,
                        columns=["Particle Identifier", "Particle Type",
                                 "Position.X", "Position.Y", "Position.Z"])
pipeline.modifiers.append(tracking)

expr = "DisplacementMagnitude\
        * (DistanceFromCenter > %g\
           && DistanceFromCenter < %g\
           && Position.Z > %g)" % (R, d, z)
mywater = mods.ComputePropertyModifier(output_property="WaterDisplacement",
                                       expressions=[expr])
pipeline.modifiers.append(mywater)

"""
expr = "WaterDisplacement > 0"
interesting_water = mods.ExpressionSelectionModifier(expression=expr)
pipeline.modifiers.append(interesting_water)
"""

coloring = mods.ColorCodingModifier(property="WaterDisplacement",
                                    only_selected=True)
pipeline.modifiers.append(coloring)

pipeline.dataset.anim.current_frame = num_frames-1
data = pipeline.compute()

radii = data.particle_properties["DistanceFromCenter"]
displacements = data.particle_properties["WaterDisplacement"]

coloring.start_value = min(radii)
coloring.end_value = max(radii)

if args.d:
    num_timesteps = num_frames - int(num_frames/2)

    num_bins = 200
    bin_edges = np.linspace(R, d, num_bins)
    mean_displacements = np.zeros(num_bins - 1)
    bin_mids = (bin_edges[1:] + bin_edges[:-1])/2

    for n in range(int(num_frames/2), num_frames):
        t0 = time.time()
        data = pipeline.compute(n)
        t1 = time.time()
        radii = data.particle_properties["DistanceFromCenter"]
        displacements = data.particle_properties["WaterDisplacement"]

        sorted_indices = np.argsort(radii)
        radii = radii[sorted_indices]
        displacements = displacements[sorted_indices]
        t2 = time.time()

        lower_index = 0
        for i in range(num_bins-1):
            index_step = np.searchsorted(radii[lower_index:], bin_edges[i+1])
            upper_index = lower_index + index_step
            mean_displacements[i] += displacements[lower_index:upper_index].mean()
            lower_index = upper_index

        t3 = time.time()
        num_done = n - int(num_frames/2) + 1
        txt = "%3d out of %d:    Ovito: %6.5s s" % (num_done, num_timesteps, t1-t0)
        txt += "    Indexing: %6.5s s    Displacement: %6.5s s" % (t2-t1, t3-t2)
        print(txt)

    mean_displacements /= num_timesteps

    plt.plot(bin_mids - R, mean_displacements, ".")
    plt.xlabel("Distance from edge of cylinder [Å]")
    plt.ylabel("Mean displacement [Å]")
    plt.title("Time-averaged mean displacement vs radius")
    plt.savefig(args.o)
    plt.show()

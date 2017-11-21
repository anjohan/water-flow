from argparse import ArgumentParser
import time
from ovito.io import import_file
from ovito.pipeline import FileSource
import ovito.modifiers as mods
import matplotlib.pyplot as plt
import numpy as np

if 1 < 0:
    a = FileSource()
    time.time()

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
expr = "(Position.Z > %g) * sqrt((Position.X - %g)*(Position.X - %g)\
                                  + (Position.Y - %g)*(Position.Y - %g))" % (z, x/2, x/2, y/2, y/2)

dist = mods.ComputePropertyModifier(output_property='DistanceFromCenter',
                                    expressions=[expr])
pipeline.modifiers.append(dist)
"""

expr = "DistanceFromCenter < %g || Position.Z < %g" % (d, z)
near_cylinder = mods.ExpressionSelectionModifier(expression=expr)
pipeline.modifiers.append(near_cylinder)

pipeline.modifiers.append(mods.InvertSelectionModifier())
pipeline.modifiers.append(mods.DeleteSelectedModifier())
"""


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


"""
set_colors(pipeline.source)

tracking = mods.CalculateDisplacementsModifier(use_frame_offset=True,
                                               frame_offset=-20)
pipeline.modifiers.append(tracking)
"""

"""
expr = "DisplacementMagnitude\
        * (DistanceFromCenter > %g\
           && DistanceFromCenter < %g\
           && Position.Z > %g)" % (R, d, z)
mydisplacement = mods.ComputePropertyModifier(output_property="WaterDisplacement",
                                              expressions=[expr])
pipeline.modifiers.append(mydisplacement)
"""

"""
expr = "WaterDisplacement > 0"
interesting_water = mods.ExpressionSelectionModifier(expression=expr)
pipeline.modifiers.append(interesting_water)
"""

#                    _
#  _ __ ___  ___  __| |
# | '_ ` _ \/ __|/ _` |
# | | | | | \__ \ (_| |
# |_| |_| |_|___/\__,_|


if args.d:
    dt = 0.0005
    num_timesteps = num_frames - int(num_frames/2)
    time_per_frame = dt * (pipeline.compute(num_frames-1).attributes["Timestep"]
                           - pipeline.compute(num_frames-2).attributes["Timestep"])
    total_time = time_per_frame * num_frames

    time_between_t0s = 10  # HENT FRA VACF
    frames_between_t0s = int(round(time_between_t0s/time_per_frame))

    first_t0 = 1000  # WHEN EQUILIBRIATED
    first_t0_frame = int(round(first_t0/time_per_frame))

    time_window = 50  # GJETNING
    frames_per_window = int(time_window/time_per_frame)
    frame_values = np.arange(frames_per_window)
    time_values = frame_values * dt

    t0_frames = np.arange(first_t0_frame, num_frames, frames_per_window)
    num_windows = len(t0_frames)

    num_bins = int((d - R)/10)  # Ca. 1 nm per bin
    bin_edges = np.linspace(R, d, num_bins+1)
    mean_displacements = np.zeros((num_bins, frames_per_window))
    bin_mids = (bin_edges[1:] + bin_edges[:-1])/2

    for i in range(1, frames_per_window):
        mod = mods.CalculateDisplacementsModifier(use_frame_offset=True,
                                                  frame_offset=i)
        pipeline.modifiers.append(mod)
        displacement = mods.ComputePropertyModifier(output_property="Displacement%d" % i,
                                                    expressions=["DisplacementMagnitude"])
        pipeline.modifiers.append(displacement)

    for i in range(num_windows):
        print("Time window %4d of %4d" % (i+1, num_windows))
        data = pipeline.compute(i)
        radii = data.particle_properties["DistanceFromCenter"]

        sorted_indices = np.argsort(radii)
        radii = radii[sorted_indices]

        lower_index = 0
        for j in range(num_bins-1):
            print("Bin %3d of %3d" % (j+1, num_bins))
            index_step = np.searchsorted(radii[lower_index:], bin_edges[i+1])
            upper_index = lower_index + index_step
            lower_index = upper_index
            for k in range(1, frames_per_window):
                print("Time value %3d of %3d" % (k+1, frames_per_window))
                displacements = data.particle_properties["Displacement%d" % k][sorted_indices]
                mean_displacements[j, k] += displacements[lower_index:upper_index].mean()
    mean_displacements /= num_windows

    """
    txt = "%3d out of %d:    Ovito: %6.5s s" % (i, num_windows, t1-t0)
    txt += "    Indexing: %6.5s s    Displacement: %6.5s s" % (t2-t1, t3-t2)
    print(txt)
    """

    mean_displacements /= num_timesteps

    for i in range(num_bins-1):
        plt.plot(time_values, mean_displacements[i], "-",
                 label=r"Distance from cylinder $\in [%.2f\ \mathrm{nm}, %.2f\ \mathrm{nm}]"
                       % (bin_edges[i]/10, bin_edges[i+1]/10))

    plt.xlabel("t [ps]")
    plt.ylabel("Mean square displacement [Å]")
    plt.title("Time-averaged mean displacement vs radius")
    plt.legend()
    plt.savefig(args.o)
    plt.show()

from argparse import ArgumentParser
import sys
import time
from ovito.io import import_file
from ovito.pipeline import FileSource
import ovito.modifiers as mods
import matplotlib.pyplot as plt
import numpy as np
import pdb

plt.figure(figsize=(8, 6))

if 1 < 0:
    a = FileSource()
    time.time()
    sys.exit()

x = float(open(".x.dat", "r").read())
y = float(open(".y.dat", "r").read())
z = float(open(".toppavbunn.dat", "r").read())
R = float(open(".sylinderradius.dat", "r").read())
d = min(x/2, y/2)

parser = ArgumentParser()
parser.add_argument("-i", default="data/vann.in.bin", help="Datafile to read")
parser.add_argument("-d", default=False, action="store_true", help="Calc. displacements")
parser.add_argument("-o", default="./data/displ_vs_R.png", help="Plot output")
parser.add_argument("--window_length", default=200, type=float)
parser.add_argument("--t0_distance", default=100, type=float)
parser.add_argument("--bin_radius", default=5, type=float)
parser.add_argument("--max_bin_radius", default=d-R, type=float)
parser.add_argument("--frames_per_window", default=4, type=int)
parser.add_argument("--pdb", default=False, action="store_true")
args = parser.parse_args()

pipeline = import_file(args.i, multiple_frames=True,
                       columns=["ParticleIdentifier", "ParticleType",
                                "Position.X", "Position.Y", "Position.Z"])

pipeline.add_to_scene()

num_frames = pipeline.source.num_frames

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
    time_per_frame = dt * (pipeline.compute(num_frames-1).attributes["Timestep"]
                           - pipeline.compute(num_frames-2).attributes["Timestep"])
    total_time = time_per_frame * num_frames

    time_between_t0s = args.t0_distance  # HENT FRA VACF
    frames_between_t0s = int(round(time_between_t0s/time_per_frame))

    first_t0 = 200  # WHEN EQUILIBRIATED
    first_t0_frame = int(round(first_t0/time_per_frame))

    max_frames_per_window = args.frames_per_window  # REDUCE RUNTIME
    time_window = args.window_length  # GJETNING
    frames_per_window = int(time_window/time_per_frame)
    if max_frames_per_window > frames_per_window:
        max_frames_per_window = frames_per_window
    skip_frames = int((frames_per_window-1)/(max_frames_per_window-1))
    frame_values = np.arange(0, frames_per_window, skip_frames)
    time_values = frame_values * time_per_frame

    t0_frames = np.arange(first_t0_frame, num_frames - frames_per_window, frames_between_t0s)
    num_windows = len(t0_frames)

    print("num_windows = ", num_windows)
    print("num_frames_per_window = ", frames_per_window)

    num_bins = int(round(args.max_bin_radius/args.bin_radius))  # Ca. 0.5 nm per bin
    bin_edges = np.linspace(R, R+args.max_bin_radius, num_bins+1)
    mean_displacements = np.zeros((num_bins, max_frames_per_window))
    bin_mids = (bin_edges[1:] + bin_edges[:-1])/2

    if args.pdb:
        pdb.set_trace()

    for i in range(num_windows):
        t0 = time.time()
        data = pipeline.compute(t0_frames[i])
        t1 = time.time()
        current_IDs = data.particle_properties["ParticleIdentifier"]
        current_positions = data.particle_properties["Position"]
        current_sorted_indices = np.argsort(current_IDs)
        current_positions = current_positions[current_sorted_indices]
        current_radii = data.particle_properties["DistanceFromCenter"][current_sorted_indices]
        current_types = data.particle_properties["ParticleType"][current_sorted_indices]

        water = (current_positions[:, 2] > z) * (current_radii > R)
        water *= (current_radii < args.max_bin_radius + R) * (current_types > 1)
        water_indices = np.where(water > 0)

        current_positions = current_positions[water_indices]

        for k in range(1, max_frames_per_window):
            future_data = pipeline.compute(t0_frames[i] + k*skip_frames)
            future_IDs = future_data.particle_properties["ParticleIdentifier"]
            future_positions = future_data.particle_properties["Position"]
            future_sorted_indices = np.argsort(future_IDs)
            future_positions = future_positions[future_sorted_indices][water_indices]
            displacements = np.linalg.norm(future_positions - current_positions, axis=1)

            for j in range(num_bins):
                radii = data.particle_properties["DistanceFromCenter"][current_sorted_indices]
                radii = radii[water_indices]
                sorted_indices = np.argsort(radii)
                radii = radii[sorted_indices]

                lower_index = np.searchsorted(radii, bin_edges[0])
                print("Time window %4d of %4d" % (i+1, num_windows), end="    ")
                print("Bin %3d of %3d" % (j+1, num_bins), end="    ")
                print("Time value %3d of %3d" % (k+1, max_frames_per_window))
                index_step = np.searchsorted(radii[lower_index:], bin_edges[j+1])
                upper_index = lower_index + index_step

                mean_displacements[j, k] += displacements[lower_index:upper_index].mean()
                lower_index = upper_index
    mean_displacements /= num_windows

    bin_edges -= R

    for i in range(num_bins):
        plt.plot(time_values, mean_displacements[i], "-",
                 label=r"Distance from cylinder $\in$ [%.2f nm, %.2f nm]"
                       % (bin_edges[i]/10, bin_edges[i+1]/10))

    plt.xlabel("t [ps]")
    plt.ylabel("Mean square displacement [Å]")
    plt.title("Time-averaged mean displacement vs radius")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.savefig(args.o, bbox_inches="tight")
    plt.show()

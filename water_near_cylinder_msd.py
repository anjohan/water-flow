from argparse import ArgumentParser
from tqdm import trange
import sys
import time
import ovito
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
    trange(3)

x = float(open(".x.dat", "r").read())
y = float(open(".y.dat", "r").read())
z = float(open(".toppavbunn.dat", "r").read())
R = float(open(".sylinderradius.dat", "r").read())
d = min(x/2, y/2)

parser = ArgumentParser()
parser.add_argument("-i", default=None, help="Datafile to read")
parser.add_argument("-o", default="./data/displ_vs_R_tmp.png", help="Plot output")
parser.add_argument("--window_length", default=200, type=float)
parser.add_argument("--t0_distance", default=100, type=float)
parser.add_argument("--first_t0", default=500, type=float)
parser.add_argument("--bin_radius", default=5, type=float)
parser.add_argument("--max_bin_radius", default=d-R, type=float)
parser.add_argument("--frames_per_window", default=-1, type=int)
parser.add_argument("--pdb", default=False, action="store_true")
parser.add_argument("-r", "--radial", default=False, action="store_true")
parser.add_argument("-p", "--progress", default=False, action="store_true")
parser.add_argument("--columns", nargs="+",
                    default=["Particle Identifier", "Particle Type",
                             "Position.X", "Position.Y", "Position.Z"])
args = parser.parse_args()

if not args.progress:
    def trange(*args, **kwargs):
        return range(*args)


loaded = True
pipeline = ovito.dataset.selected_pipeline
if pipeline is None:
    assert args.i is not None
    pipeline = import_file(args.i, multiple_frames=True, columns=args.columns)
    pipeline.add_to_scene()
    ovito.dataset.save(args.i + ".ovito")
    loaded = False


num_frames = pipeline.source.num_frames

expr = "(Position.Z > %g) * sqrt((Position.X - %g)*(Position.X - %g)\
                                  + (Position.Y - %g)*(Position.Y - %g))" % (z, x/2, x/2, y/2, y/2)

dist = mods.ComputePropertyModifier(output_property='DistanceFromCenter',
                                    expressions=[expr])
pipeline.modifiers.append(dist)


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


set_colors(pipeline.source)
"""

#                    _
#  _ __ ___  ___  __| |
# | '_ ` _ \/ __|/ _` |
# | | | | | \__ \ (_| |
# |_| |_| |_|___/\__,_|


dt = 0.0005
time_per_frame = dt * (pipeline.compute(num_frames-1).attributes["Timestep"]
                       - pipeline.compute(num_frames-2).attributes["Timestep"])
total_time = time_per_frame * num_frames

time_between_t0s = args.t0_distance  # vacf implies 1 ps is sufficient
frames_between_t0s = int(round(time_between_t0s/time_per_frame))

first_t0 = args.first_t0  # WHEN EQUILIBRIATED
first_t0_frame = int(round(first_t0/time_per_frame))

max_frames_per_window = args.frames_per_window  # REDUCE RUNTIME
time_window = args.window_length
frames_per_window = int(time_window/time_per_frame)
if max_frames_per_window > frames_per_window or max_frames_per_window <= 0:
    max_frames_per_window = frames_per_window
skip_frames = int((frames_per_window-1)/(max_frames_per_window-1))
frame_values = np.arange(0, frames_per_window, skip_frames)
time_values = frame_values * time_per_frame

t0_frames = np.arange(first_t0_frame, num_frames - frames_per_window, frames_between_t0s)
num_windows = len(t0_frames)
frames_per_window = len(frame_values)

print("num_windows = ", num_windows)
print("num_frames_per_window = ", frames_per_window)

num_bins = int(round(args.max_bin_radius/args.bin_radius))  # Ca. 0.5 nm per bin
bin_edges = np.linspace(R, R+args.max_bin_radius, num_bins+1)
mean_displacements = np.zeros((num_bins, frames_per_window))
bin_mids = (bin_edges[1:] + bin_edges[:-1])/2

tot_iterations = num_windows * (frames_per_window-1) * num_bins
iteration_count = 0
max_progress = 20

if args.pdb:
    pdb.set_trace()

for i in trange(num_windows, desc="%-15s" % "Time window:"):
    t0 = time.time()
    data = pipeline.compute(t0_frames[i])
    t1 = time.time()
    current_IDs = data.particle_properties["Particle Identifier"]
    current_positions = data.particle_properties["Position"]
    current_sorted_indices = np.argsort(current_IDs)
    current_positions = current_positions[current_sorted_indices]
    current_radii = data.particle_properties["DistanceFromCenter"][current_sorted_indices]
    current_types = data.particle_properties["Particle Type"][current_sorted_indices]

    water = (current_positions[:, 2] > z) * (current_radii > R)
    water *= (current_radii < args.max_bin_radius + R) * (current_types > 1)
    water_indices = np.where(water > 0)

    current_positions = current_positions[water_indices]

    radii = data.particle_properties["DistanceFromCenter"][current_sorted_indices]
    radii = radii[water_indices]
    sorted_indices = np.argsort(radii)
    radii = radii[sorted_indices]

    bin_indices = np.zeros(num_bins+1, dtype=int)
    bin_indices[-1] = len(radii)

    for j in range(1, num_bins):
        bin_indices[j] = np.searchsorted(radii, bin_edges[j])
    if args.pdb:
        pdb.set_trace()

    for k in trange(1, frames_per_window, desc="%-15s" % "Time frame:"):
        future_data = pipeline.compute(t0_frames[i] + k*skip_frames)
        future_IDs = future_data.particle_properties["Particle Identifier"]
        future_positions = future_data.particle_properties["Position"]
        future_sorted_indices = np.argsort(future_IDs)
        future_positions = future_positions[future_sorted_indices][water_indices]
        drs = (future_positions - current_positions)[sorted_indices]
        max_dim = 2 if args.radial else 3
        displacements = np.sum(drs[:, :max_dim]**2, axis=1)

        for j in trange(num_bins, desc="%-15s" % "Radial bin:"):
            mean_displacements[j, k] += displacements[bin_indices[j]:bin_indices[j+1]].mean()
mean_displacements /= num_windows

bin_edges -= R

for i in range(num_bins):
    plt.plot(time_values, mean_displacements[i], "-",
             label=r"Distance from cylinder $\in$ [%.2f nm, %.2f nm]"
                   % (bin_edges[i]/10, bin_edges[i+1]/10))

plt.xlabel("t [ps]")
plt.ylabel("Mean square displacement [Å^2]")
plt.title("Time-averaged mean square " + ("radial " if args.radial else "")
          + "displacement vs radius")
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig(args.o, bbox_inches="tight")
plt.show()

with open(args.o + ".log", "w") as outfile:
    t = time.strftime("%H:%M %d/%m-%y")
    outfile.write("Generated " + t + "\n")
    if loaded:
        outfile.write("Loaded .ovito state file.\n")
    outfile.write("Args: " + " ".join(sys.argv[1:]))

import Pkg; Pkg.activate(".")

#local libraries
include("src/main.jl")
using .SISSIBoxes
const SISSI = SISSIBoxes

#=
    This script is designed to make the user familiar with how to obtain uniform fluid boxes from the SISSI simulation snapshot data on OP3PGN.
    You can use the functions introduced in this script to build a more automated analysis (i.e. loops over snapshots, bubbles, etc.).
    Paths to data files are hard-coded into the source code. This has the advantage that it allows the user to run these scripts with minimal knowledge about the data management,
    but it should be kept in mind, if one tries to run it on a different system. That simply won't work unless the code is modifed significantly.
    Currently only Leonard Romano, Andreas Burkert, Manuel Behrendt and Mingyu Hu have access to these files on OP3PGN (and system admins).
=#

# Parameters
model_name = "N1"               # Select the model you want to analyze (select from: "N1", "N10", "N1x10", "N1x50", "N0", "N0_zoom")
dx_pc      = 10.0               # desired resolution of the output uniform grid (approximately; min. ~0.18 pc)
# Output files are written to the following directory in the form "path_output/output_00XXX.dat" where XXX is the snapshot number
# NOTE: Different SNRs need different directories! Files are not labelled by SNR. If you want to analyze multiple different SNRs make sure to provide different directories for each!
path_output = "/path/to/output" 

# This function returns the range of snapshot numbers for a given simulation run
isnapshots = get_snapshots(model_name=model_name)

# This function returns the time at each snapshot for a given simulation run
# together with the above function you can pinpoint specific snapshots to analyze
time = get_time(model_name=model_name)

# This function returns the full time-series of box centers for a given SNR (labelled 1-30) in a given model, along with the extent of the ISM box and the tracer variable.
# This can be fed to "get_data_box" to obtain boxes consistent with my analysis
# NOTE: There are no SNRs in models N0 and N0_zoom and correspondingly this function won't work with these models!
# If you plan to analyze them, use the positions of the SNRs and their respective boxes of another model instead. 
centers, extent, var_tracer = get_regions(model_name=model_name, i_SNR=23)

# This function loads the data, projects it onto a 3D uniform grid and writes the results to a file
# The function inputs are:
# isnapshot[Integer] <-> which snapshot should be analyzed
# model_name[String] <-> which simulation data?
# path_output[String] <-> where to write results?
# center[Vector{<:Real}] <-> center of the box to be analyzed (3 x kpc; relative to the galactic center)
# extent[Vector{<:Real}] <-> extent of the box to be analyzed (6 x kpc; relative to the center, i.e. [xmin, xmax, ymin, ymax, zmin, zmax])
# var_tracer[Symbol] <-> which variable is the tracer variable (:var6 or :var7)
# dx_pc[Real] <-> desired resolution in pc (only approximately -> resolution multiple of 2 of grid resolution)
# corotating_frame[Bool] <-> Should the average galactic rotation be subtracted from the velocity field? 
isnapshot = 329
get_data_box(isnapshot, model_name=model_name, path_output=path_output, center=centers[1], extent=extent, var_tracer=var_tracer, dx_pc=dx_pc, corotating_frame=true)
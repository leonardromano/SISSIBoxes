__precompile__(true)
module SISSIBoxes
    using Mera, JuliaDB
    using JLD2, CodecLz4

    export
        get_snapshots, # return list of all snapshots for a given simulation
        get_time,      # return list of all snapshot times for a given simulation
        get_regions,   # return list of all box-centers, box extents and the color variable for a given SNR
        get_data_box   # write a file containing a 3D grid of the fluid variables (rho, T, vx, vy, vz, tracer) for a given ISM region

    # constants, types & utility functions
    include("utility/constants.jl")
    include("utility/utility_functions.jl")

    # I/O functionality
    include("IO/IO_utility.jl")

    # snapshot analysis drivers
    include("snapshots/analysis_driver.jl")
end
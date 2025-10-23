"Return the path to the snapshot files on OP3PGN."
get_path_snapshot(model_name::String) = "/data1/lromano/shared/SISSI/" * model_name * "/merafiles"

"Return the path to the curated data files on OP3PGN."
get_path_data() = "/data2/lromano/shared/SISSI/data/"

"Return the path to the curated data bubble-files on OP3PGN."
get_path_bubbles(model_name::String) = get_path_data() * "Simulations/" * model_name * "/Bubbles/"

"Return the range of snapshots of a given simulation model"
function get_snapshots(;model_name::String="N1")
    if model_name == "N1"
        return 329:378
    elseif model_name == "N10"
        return 329:375
    elseif model_name == "N1x10"
        return 329:430
    elseif model_name == "N1x50"
        return 329:455
    elseif model_name == "N0"
        return 329:373
    elseif model_name == "N0_zoom"
        return 329:426
    else
        println("Model $model_name not defined!")
        exit(86)
        return 329:328
    end
end

"Return the snapshot times of a given simulation model"
function get_time(;model_name::String="N1")
    if model_name ∈ ["N0", "N0_zoom"]
        path_snapshot = get_path_snapshot(model_name)
        ioutputs = get_snapshots(model_name=model_name)
        t0 = infodata(329, path=path_snapshot, verbose=false).time
        t2Myr = infodata(329, path=path_snapshot, verbose=false).scale.Myr
        return [infodata(ioutput, path=path_snapshot, verbose=false).time - t0 for ioutput in ioutputs] * t2Myr
    else
        path_data = get_path_bubbles(model_name) * "Bubble_001/"
        ioutputs = get_snapshots(model_name=model_name)
        return [JLD2.load(path_data * "ISM_box_$(lpad(ioutput,5,"0")).jld2", "header/time") for ioutput in ioutputs]
    end
end

"Return the CoMs, the relative extent of the region covering the SNR and the corresponding tracer variable."
function get_regions(;model_name::String="N1", i_SNR::Integer=30)
    if model_name ∈ ["N1", "N10", "N1x10", "N1x50"] && 1 <= i_SNR <= 31
        path_data = get_path_bubbles(model_name) * "Bubble_$(lpad(i_SNR,3,"0"))/"
        ioutputs = get_snapshots(model_name=model_name)

        CoMs  = [JLD2.load(path_data * "ISM_box_$(lpad(ioutput,5,"0")).jld2", "header/center") for ioutput in ioutputs]
        box   = JLD2.load(path_data * "ISM_box_$(lpad(329,5,"0")).jld2", "header/box")
        color = JLD2.load(path_data * "ISM_box_$(lpad(329,5,"0")).jld2", "header/color") == "red" ? :var6 : :var7

        return CoMs, box, color
    else
        println("Error: Model $model_name not defined, or does not contain any SNRs. Make sure that i_SNR ∈ [1,...,30]!")
        exit(86)
        return [zeros(3)], zeros(6), :var6 
    end
end

"Read rotation curve from jld2 file"
function read_rotation_curve(; verbose::Bool=false)
    # Check if file exists
    filename = get_path_data() * "rotation_curve.jld2"
    if !isfile(filename)
        println("Error: File $filename does not exist!")
        exit(86)
    end

    # read data
    f = jldopen(filename)

    rbins = f["rbins"]
    vrot  = f["vrot"]

    if verbose
        println("Reading rotation curve of snapshot no. $(f["snapshot"])...")
        println("Radial range from 0.0 to $(Rmax) kpc")
        println("There are $(length(vrot)) bins of width $(rbins[2] - rbins[1]) kpc.")
    end

    close(f)

    return rbins, vrot
end
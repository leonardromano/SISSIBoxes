"Load ISM box, project it onto uniform grid and write results to file."
function get_data_box(ioutput::Integer; model_name::String="N1", path_output::String="", center::Vector{<:Real}=[0, 0, 0], extent::Vector{<:Real}=[0, 0, 0, 0, 0, 0], 
                      var_tracer::Symbol=:var6, dx_pc::Real=10.0, corotating_frame::Bool=true)
        # setting resolution
        lmax   = min(ceil(Int, log2(48000.0 / dx_pc)), 18)
        
        # status message
        println("Setting resolution level to $lmax (max: 18), corresponding to $(round(48000 / 2^lmax, digits=2)) pc.")
        println("The box will be resolved by approx. $((extent[2] - extent[1]) * 2^lmax / 48.0) x $((extent[4] - extent[3]) * 2^lmax / 48.0) x $((extent[6] - extent[5]) * 2^lmax / 48.0) cells.")

        # read in gas
        println("Loading data...")
        path_snapshot = get_path_snapshot(model_name)
        gas = loaddata(ioutput, path=path_snapshot, datatype=:hydro, xrange=extent[1:2], yrange=extent[3:4], zrange=extent[5:6], center=center.+24, range_unit=:kpc);

        # Remove rotation velocity
        if corotating_frame
            println("Removing rotation around galactic center...")
            remove_rotation!(gas)
        end

        # get time
        t0 = infodata(329, path=path_snapshot, verbose=false).time
        t  = gas.info.time
        t_Myr = (t - t0) * gas.info.scale.Myr

        # project onto coarse grid (IMPORTANT: density HAS TO come LAST!!!)
        println("Interpolating AMR grid onto uniform grid...")
        cubes = interpolate_hydro(gas, [:T, :vx, :vy, :vz, var_tracer, :rho], [:K, :km_s, :km_s, :km_s, :standard, :nH], lmax=lmax, center=center.+24, verbose=false)

        # get cube sizes
        ss = size(cubes[:rho])

        # write results to file
        println("writing results to $(path_output) /output_$(lpad(ioutput,5,"0")).dat...")
        open(path_output * "/output_$(lpad(ioutput,5,"0")).dat", "w") do f
            write(f, "t = $t_Myr Myr\n")
            write(f, "x,y,z, rho[amu/cm^3], vx[km/s], vy[km/s], vz[km/s], T[K], tracer \n")
            write(f, "x-cells=$(ss[1]), y-cells=$(ss[2]), z-cells=$(ss[3]) \n")
            write(f, "cellsize=$(48000 / 2^lmax) pc \n")
            for ijk in CartesianIndices(cubes[:rho])
                write(f, "$(ijk[1]), $(ijk[2]), $(ijk[3]), $(cubes[:rho][ijk]), $(cubes[:vx][ijk]), $(cubes[:vy][ijk]), $(cubes[:vz][ijk]), $(cubes[:T][ijk]), $(cubes[var][ijk])\n")
            end
        end
end
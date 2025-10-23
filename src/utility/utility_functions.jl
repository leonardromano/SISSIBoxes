#########################################################################################################
# Data interpolation

"Project gas data onto unfiform, 3D, cartesian mesh."
function interpolate_hydro(dataobject::HydroDataType, vars::Array{Symbol,1}, units::Array{Symbol,1}; lmax::Int=dataobject.lmax,
                           xrange::Array{<:Any,1}=[missing, missing], yrange::Array{<:Any,1}=[missing, missing], zrange::Array{<:Any,1}=[missing, missing],
                           center::Array{<:Any,1}=[0., 0., 0.], range_unit::Symbol=:standard, verbose::Bool=true)

    lmin = dataobject.lmin

    if  lmax > dataobject.lmax
        println("Current data lmax=$(dataobject.lmax) < your lmax=$lmax")
        sleep(1)
        quit()
    elseif lmax < lmin
        println("Current data lmin=$(lmin) > your lmin=$lmax")
        sleep(1)
        quit()
    end

    # get ranges
    ranges = Mera.prepranges(dataobject.info,range_unit, verbose, xrange, yrange, zrange, center, dataranges=dataobject.ranges)

    fmax = 1<<lmax
    shiftx = Int(round( ranges[1] * fmax ) )
    a_size = Int(round( ranges[2] * fmax ) ) - shiftx

    shifty = Int(round( ranges[3] * fmax ) )
    b_size = Int(round( ranges[4] * fmax) ) - shifty

    shiftz = Int(round( ranges[5] * fmax ) )
    c_size = Int(round( ranges[6] * fmax ) ) - shiftz

    # output object
    cubes = Dict{Symbol, Array{FloatType}}()

    #allocate memory
    for var in vars
        cubes[var] = zeros(FloatType, a_size,b_size,c_size)
    end

    # number of cells
    Nrows = length(dataobject.data)

    Threads.@threads for irow in 1:Nrows
        row = dataobject.data[irow]
        if row.level <= lmax
            factor = 1<<(lmax - row.level)

            a_min = max((row.cx-1) * factor - shiftx, 1)
            a_max = min(row.cx     * factor - shiftx, a_size)
            b_min = max((row.cy-1) * factor - shifty, 1)
            b_max = min(row.cy     * factor - shifty, b_size)
            c_min = max((row.cz-1) * factor - shiftz, 1)
            c_max = min(row.cz     * factor - shiftz, c_size)

            for var in vars
                cubes[var][ a_min:a_max, b_min:b_max, c_min:c_max ] .= FloatType(my_getproperty(row, var))
            end
        elseif row.level > lmax #todo: average cells together to lower level
            factor = 1<<(row.level-lmax)

            a = round(Int32, (row.cx-0.5) / factor - shiftx)
            b = round(Int32, (row.cy-0.5) / factor - shifty)
            c = round(Int32, (row.cz-0.5) / factor - shiftz)

            if (1 <= a <= a_size) && (1 <= b <= b_size) && (1 <= c <= c_size)
                for var in vars
                    cubes[var][ a, b, c ] += FloatType(my_getproperty(row, var)) / factor^3
                end
            end
        end
    end

    #units & conversions
    for var in vars
        # convert to primitve variable
        if var in [:T, :vx, :vy, :vz, :var6, :var7] cubes[var] ./= cubes[:rho] end
        
        selected_unit= getunit(dataobject, var, vars, units)
        cubes[var] *= selected_unit
    end

    return  cubes
end

"Return conserved quantities associated with given variable."
function my_getproperty(row, var::Symbol)
    if var == :T
        return getproperty(row, :p)
    elseif var in [:vx, :vy, :vz, :var6, :var7]
        return getproperty(row, :rho) * getproperty(row, var)
    else
        return getproperty(row, var)
    end
end

#########################################################################################################
# Rotation calculation

"Remove mean rotation velocity from each cell's/particle's velocity"
function remove_rotation!(dataobject::HydroPartType)
    # Load rotation curve
    rbins, vrot = read_rotation_curve()

    # convert to code units
    vrot /= dataobject.info.scale.km_s

    # remove rotation curve
    remove_rotation!(dataobject, rbins, vrot)
end

"Update velocity vector in dataobject by removing rotation velocity"
function remove_rotation!(dataobject::HydroPartType, r_bins::Vector{Float64}, v_rot::Vector{Float64})
    # Get velocities and coordinates
    x, y, z = my_getpositions(dataobject, :kpc, center=[:bc])
    r       = zeros(size(x))
    Threads.@threads for i in eachindex(r) r[i] = sqrt(x[i]^2 + y[i]^2) end
    vx      = select(dataobject.data, :vx)
    vy      = select(dataobject.data, :vy)

    # get bin width
    dr = 0.5 * (r_bins[2] - r_bins[1])

    # iterate over cells
    Threads.@threads for icell in eachindex(r)
        irot = min(round(Int, (r[icell] / dr + 1) / 2), length(r_bins))

        vx[icell] += v_rot[irot] * y[icell] / r[icell] 
        vy[icell] -= v_rot[irot] * x[icell] / r[icell] 
    end

    # Update entries in data base
    dataobject.data = transform(dataobject.data, :vx => vx)
    dataobject.data = transform(dataobject.data, :vy => vy)
end

#########################################################################################################
# Mera patches

"""
#### Get the x,y,z positions from the dataset (cells/particles/clumps/...):
```julia
my_getpositions( dataobject::DataSetType;
                    unit::Symbol=:standard,
                    direction::Symbol=:z,
                    center::Array{<:Any,1}=[0., 0., 0.],
                    center_unit::Symbol=:standard)

return x, y, z
```


#### Arguments
##### Required:
- **`dataobject`:** needs to be of type: "DataSetType"
##### Predefined/Optional Keywords:
- **`center`:** in unit given by argument `center_unit`; by default [0., 0., 0.]; the box-center can be selected by e.g. [:bc], [:boxcenter], [value, :bc, :bc], etc..
- **`center_unit`:** :standard (code units), :Mpc, :kpc, :pc, :mpc, :ly, :au , :km, :cm (of typye Symbol) ..etc. ; see for defined length-scales viewfields(info.scale)
- **`direction`:** todo
- **`unit`:** return the variables in given unit

### Defined Methods - function defined for different arguments

- my_getpositions( dataobject::DataSetType; ...) # one given dataobject
- my_getpositions( dataobject::DataSetType, unit::Symbol; ...) # one given dataobject and position unit

"""
function my_getpositions(dataobject::HydroDataType, unit::Symbol;
    center::Array{<:Any,1}=[0., 0., 0.], center_unit::Symbol=:standard)

    # convert center to standard notation
    center = Mera.center_in_standardnotation(dataobject.info, center, center_unit)

    # get box dimensions
    funit = getfield(dataobject.info.scale, unit)
    Lbox  = dataobject.info.boxlen * funit
    lmax  = dataobject.info.levelmax
    lmin  = dataobject.info.levelmin

    # get level
    Δx = lmax == lmin ? 0.5^lmax : 0.5.^select(dataobject.data, :level)

    # Get position variables
    x = ((select(dataobject.data, :cx) .- 0.5) .* Δx .- center[1]) * Lbox
    y = ((select(dataobject.data, :cy) .- 0.5) .* Δx .- center[2]) * Lbox
    z = ((select(dataobject.data, :cz) .- 0.5) .* Δx .- center[3]) * Lbox

    return x, y, z
end

function my_getpositions(dataobject::PartDataType, unit::Symbol;
    center::Array{<:Any,1}=[0., 0., 0.], center_unit::Symbol=:standard)
    return getpositions(dataobject, unit, center=center, center_unit=center_unit)
end
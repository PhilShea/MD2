module MD2 
# Two dimensional Molecular Dynamics
using Random
using LinearAlgebra
using StaticArrays
using Plots
using Dates
using JLD2
using DataFrames
using CSV
using FLoops
using Printf

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

export MDRecord, InitMDRecord!, CorrShell!, PlotMDRecord, Nscalars, RunSystem, 
        InitSystem

struct MDRecord{T <: Real}
    N::Int32 # Number of particles.
    n::Int32 #
    totalsteps::Int32 # Number of steps since init or heating/cooling.
    tau::T # Time step of simulation.
    shell::T # Radius of involved particles.
    rskn::T # The range of the pe & forc; must be <= shell.
    forceFun::AbstractString # The force function with an implicit 1/r.
    peFun::AbstractString # THe potentiol energy as a function of r.
    skip::Int32 # Number of steps to skip between updating involved particles.
    bx::T # Box size of the x axis {-bx to bx}.
    by::T # y box size.
    p::Vector{ SVector{2, T}} # particle positions.
    v::Vector{ SVector{2, T}} # particle velocities.
    function MDRecord{s}(N::Integer, n::Integer, totalsteps::Integer, tau::s, 
        shell::s, rskn::s, forceFun::AbstractString, peFun::AbstractString, 
        skip::Integer, bx::s, by::s, p::Vector{ SVector{ 2, s}},
        v::Vector{ SVector{ 2, s}}) where s <: Real
        @show size( p), size( v)
        if !( size( p) == size( v) ==  (N,))
            error("MDRecord: invalid matrix sizes: p: $(size(p)), \
            v: $(size(v)), should be ($N,)")
        end
        if !( rskn <= shell)
            error( "MDRecord: rskn ($rskn) must be <= shell ($shell)")
        end
        new(N, n, totalsteps, tau, shell, rskn, forceFun, peFun, skip, bx, by, 
            p, v)
    end
end

"""
    updateMDRecord( ss::MDRecord{T}, totalsteps::Integer) where T <: Real

Creates a new record with updated total steps.
"""
updateMDRecord( ss::MDRecord{T}, totalsteps::Integer) where T <: Real = 
    MDRecord{T}( ss.N, ss.n, totalsteps, ss.tau, ss.shell, ss.rskn, 
    ss.forceFun, ss.peFun, ss.skip, ss.bx, ss.by, ss.p, ss.v)

"""
    MDRecord{T}(N::Integer, n::Integer, tau, shell, rskn, skip::Integer, 
                     bx, by) where {T<:Real}

Creates a record with zero vectors for p and v, and creates the functions based
    on an input n.
"""
function MDRecord{T}(N::Integer, n::Integer, tau, shell, rskn, skip::Integer, 
                     bx, by) where T<:Real
    funs = rn(T, n)
    MDRecord{T}( N, n, 0, T(tau), T(shell), T(rskn), funs.force, funs.PE, 
                skip, T(bx), T(by),
                 zeros( SVector{ 2, T}, N), zeros( SVector{ 2, T}, N))
end

function InitMDRecordSimple!( ss::MDRecord{<:Real})
    for i in 1:ss.N
        ss.p[ i] .= [SVector( 2.0 * mod( i, ss.bx) - ss.bx,
            floor( i / ss.bx) - ss.bx) for i in 1:ss.N]
    end
    ss.v .= (0.1 * ss.bx) .* randn( eltype( ss.v), ss.N)
    ss
end

TDAM( p::SVector{2, <: Real}, v::SVector{2, <: Real}) = 
            p[1] * v[2] - p[2] * v[1]
    # Two Dimensional Angular momentum

AngularMomentum(  p::Vector{ SVector{2, T}}, 
                  v::Vector{ SVector{2, T}}) where {T <: Real} = 
                  sum( TDAM.( p, v))
    # the cross product of two vectors in a plane will be in the orthogonal
    # to the plane.
   
function InitMDRecord!( ss::MDRecord{T}; Temp=0.30, io=stdout) where T <: Real
        # box is centered at zero, and size 2 bx x 2 by
    density = ss.N/(4 * ss.bx * ss.by)
    c = sqrt(3) / 2
    a = sqrt( 2.0 / (density * sqrt(3))) # Hexagonal grid spacing
    Nc = round( Int, 2 * ss.bx / a) # number of rows
    Nr = round( Int, ss.N / Nc) # number of columns
    println( io, "InitMDRecord----")
    println( io, "Density: $density, Grid: $a, Columns: $Nc, Rows: $Nr")
    # place particles by row and column.
    k = 1 # particle k
    for i in 0:(Nr - 1) # 
        for j in 0:(Nc - 1)
            ss.p[ k] = SVector( 
                (a / 4) * ( 1 + 2 * mod( i, 2)) - ss.bx + j * a, 
                a * c * (i + 0.5) - ss.by)
            #@show i, j, k, ss.p[k]
            k += 1
        end
    end
    ss.v .= (sqrt( Temp)) .* randn( eltype( ss.v), ss.N)
    # each dimention has energy kT/2 = V^2/2, E[V^2]=T
    ke = T(0.5) * sum( dot.(ss.v, ss.v)) # Kinetic Energy
    println( io, "Initial Total Kinetic Energy: $ke")
    m = sum( ss.v) / ss.N
    println( io, "Linear Momentum: $m")
    ss.v .= [ x - m for x in ss.v]
    println( io, "Corrected Linear Momentum: ", sum( ss.v) / ss.N)
    L = AngularMomentum( ss.p, ss.v)
    println( io, "Angular Momentum: ", L)
    r = norm.(ss.p)
    S = sum( r)
    println( io, "Moment of Inertia: ", S)
    w = L/S
    l = w .* r
    dv = zeros( SVector{ 2, T}, ss.N)
    for i in 1:ss.N
        dv[i] = (l[i] / r[i]^2) .* SVector( -ss.p[i][2], ss.p[i][1])
    end
    ss.v .= ss.v .- dv
    println( io, "Corrected Angular Momentum: ", AngularMomentum( ss.p, ss.v))
    println( io, "Final Linear Momentum: ", sum( ss.v) / ss.N)
    ke = T(0.5) * sum( dot.(ss.v, ss.v)) # Kinetic Energy
    println( io, "Final Total Kinetic Energy: ", ke)
    return ss
end

function PlotMDRecord( pn::Vector{ SVector{2, T}}, ss::MDRecord{T}; 
            margin = 1.2, vecs = true,
             show = true, virtual = true) where T <: Real
    x  = [s[1] for s in pn] 
    y  = [s[2] for s in pn]
    if vecs 
        vx = [s[1] for s in ss.v]
        vy = [s[2] for s in ss.v]
        #@show x
        quiver( x, y, quiver=( vx, vy), show=show, legend = false,
            xlim = margin .* [ -ss.bx, ss.bx], ylim = margin .* [-ss.by, ss.by])
    else
        scatter( x, y, show=show, legend = false,
        xlim = margin .* [ -ss.bx, ss.bx], ylim = margin .* [-ss.by, ss.by])
    end
    Plots.vline!( [-ss.bx, ss.bx], linecolor = :gray30)
    Plots.hline!( [-ss.by, ss.by], linecolor = :gray30)
    tbx = 2 * ss.bx
    tby = 2 * ss.by
    if virtual
        for xl in [-tbx, 0.0, tbx]
            for yl in [-tby, 0.0, tby]
                if (xl == 0.0) && (yl == 0.0) 
                    continue 
                end
                scatter!( x .+ xl, y .+ yl, markercolor = :gray40, 
                          markersize = 1)
            end # yl
        end # xl
    end # if virtual...
end

PlotMDRecord( ss::MDRecord; margin = 1.2) = PlotMDRecord( ss.p, ss; margin)

function saveMDRecord( filename::AbstractString, SystemState::MDRecord) 
    dt = Dates.now()
    fname = filename * ".jld2"
    isfile( fname) && error( "file $fname already exists.")
    jldsave( fname; SystemState, dt)
end 

function readMDRecord( filename::AbstractString; io=stdout)
    f = jldopen( filename, "r")
    k = JLD2.keys( f)
    SystemState = read( f, k[1])
    println( io, "System State restored from ", read( f, k[2]))
    close( f)
    return SystemState
end

"""
    pbc( x::T, b::T) where T<: Real

Returns the delta such that -b <= x - delta <= b
"""
pbc( x::T, b::T) where T <: Real = (abs(x) < b ? T(0) : (sign(x) * 2 * b))

"""
    pbcDelta( v::SVector{2,T}, bx::T, by::T) where T <: Real

Returns a vector delta such that v-delta will be in the box.  
This is used to correct pnm1 points when pn moved out of the box.  The call 
sequence is:

delta = pbcDelta( pn, bx, by)
pn = pn - delta
pnm1 = pnm1 - delta

Note that the pnm1 points may not be in the box anymore.
"""
pbcDelta( v::SVector{2,T}, bx::T, by::T) where T <: Real = 
    SVector{2, T}( pbc( v[1], bx), pbc( v[2], by)) 

"""
    pbc( v::SVector{2,T}, bx::T, by::T) where T <: Real

returns v in periodic boundary conditions of +- bx & +- by
"""
pbc( v::SVector{2,T}, bx::T, by::T) where T <: Real = v - pbcDelta( v, bx, by)

struct CorrVec
    Mmax::Int32 # Length of corr. Should be greater than 
                  # N * density * pi * shell^2
                  # for density = 1/pi, N * shell^2 
    m::Vector{Int32} # Hold the number of defined items for each particle
    corr::Vector{Int32} 
    # corr[m[i]:(m[i + 1] - 1)] hold a list of interacting particles for 
    # particle i
    function CorrVec( Mmax::Int, m::Vector{Int}, corr::Vector{Int})
        if !(length( corr) >= length(m)) 
            error("CorrVec object corr must be greater then m ( $(length(m))) but is $(length( corr))")
        elseif !(length( corr) === Mmax)
            error( "CorrVec object corr must be vectors of size $Mmax but are $(length( corr))")
        end
        new( Mmax, m , corr)
    end
end

CorrVec(N::Integer, Mmax::Integer) = CorrVec( Mmax, zeros(Int, N), 
        zeros( Int, Mmax))

function CorrShell!( cv::CorrVec, ss::MDRecord{T}, 
                     pn::Vector{ SVector{2, T}}) where T <: Real
    # Find all particles within shell
    # cv should be declared in the calling program, be size (ss.N-1),
    # where M is large enough to hold the biggest list of interacting particles 
    # possible.
    # The m vector will be reset, and if the cv vectors have left over values, 
    # they will be ignored.
    if (length( cv.m) < ss.N) 
        error( "corr is too small: cv.m is length $(length(cv.m)), system is $(ss.N)") 
    end
    shell = ss.shell#^2
    mt = 0 # Init temp m
    @fastmath @inbounds for i::Int32 = 1:(ss.N - 1) 
        # The last particle will be done by default
        cv.m[ i] = mt # mt will not have changed if no particles were added
        # cv.m[1] will be zero.
        for j::Int32 in (i + 1):ss.N
            #rs = pbc( pn[i] - pn[j], ss.bx, ss.by) 
            ds::T = LinearAlgebra.BLAS.nrm2( 2, 
                        pbc( pn[i] - pn[j], ss.bx, ss.by), 1) #dot( rs, rs) 
            # @show i j pn[i] pn[j] rs ds
            if ds < shell 
                mt += 1
                if mt > cv.Mmax 
                    error( "MD2.CorrShell!: m > Mmax ($(cv.Mmax)) @ part. $i") 
                end
                cv.corr[ mt] = j
            end # if ds ...
        end # for j ...  
    end # for i...
    cv.m[ ss.N] = mt # likely cv.m[ ss.N - 1]
end # CorrShell!

struct MD2State{ T <: Real} # This contains constants for the 
    ## run, including pre-allocated Arrays
    N::Int32 # Number of particles. Copied from ss
    oo2t::T # one over 2 tau (1/(2*tau))
    density::T # calculated density
    NHOPr::T # Range for Nelson-Halperin order parm calc.
    stepbskip::Int32 # Number of CorrShell! calls, so steps=ss.skip*stepbskip
    ss::MDRecord{ T} # the initiating record
    peFun::Function # The potential energy function.
    forceFun::Function # The force function.
    Mmax:: Int32 # Length of CorrVec.cv.
    cv::CorrVec # vector of correlatted particles.
    pn::Vector{ SVector{2, T}} # current particle positions.
    pnm1::Vector{ SVector{2, T}} # last positions (may be outside the box)
    pnp1::Vector{ SVector{2, T}} # hold the new positions during the calculation
    delta::Vector{ SVector{2, T}} # hold any changes mapping particles into box
    forces::Vector{ SVector{2, T}} # forces pre-allocated vecotor
    scalars::Matrix{ Float64} # matrix to hold scalars for each step.
end

"""
    initpnm1( ss::MDRecord{T}) where T<:Real

Creates the n-1 positions from the velocity vectors in ss.
"""
function initpnm1( ss::MDRecord{T}) where T<:Real
    map( i -> (ss.p[i] - (ss.tau .* ss.v[i])), 1:ss.N)
end    

"""
    MD2State( ss::MDRecord{T}, Mmax::Integer, Steps::Integer, 
        Nscalar::Integer) where T <: Real

Inits MD2State simply. ss should be initialized with InitMDrecord! or read from
a file.  Mmax is the size of CorrVec, and must be > N*density*pi*shell^2.
Steps is the number of steps being run, and, with Nscalar (the number of
scalars) will set the size of the scalar matrix.
 Note that this initializer will set st.pn = ss.pn, 
and this means that they share memory (i.e., they are the same item).
"""
function MD2State( ss::MDRecord{T}, Mmax::Integer, Steps::Integer, 
                   Nscalar::Integer = Nscalar; io = stdout) where T <: Real
    N = ss.N # used a lot
    peFun    = @RuntimeGeneratedFunction( Meta.parse( ss.peFun))
    forceFun = @RuntimeGeneratedFunction( Meta.parse( ss.forceFun))
    # execute an eval of the new functions to force compilation
    println( io, "MD2State: Pe(rskn) = $(peFun( T(0.9) * ss.rskn)), \
    Force(rskn) = $(forceFun( T(0.9) * ss.rskn))")
    density = N / ( 4 * ss.bx * ss.by)
    a = sqrt(2/(density * sqrt(3)))
    NHOPr = sqrt(3) * a
    pefcmin = density * pi * ss.shell^2 # avereage # parts. in shell.
    minMmax = N * pefcmin
    println( io, "Density: $density, Shell: $(ss.shell), NHOPr: $NHOPr")
    println( io, "Mean # in shell: $pefcmin, Min Mmax: $minMmax")
    Mmax > minMmax || error(
        "MD2State: Mmax ($Mmax) < minimum ($minMmax)")
    stepbskip = cld( Steps, ss.skip) #Ceiling( Steps/ss.skip)
    MD2State{T}( N, 1/( 2 * ss.tau), density, NHOPr, stepbskip, ss, peFun, 
                forceFun, Mmax, CorrVec( N, Mmax), ss.p, initpnm1( ss), 
                 zeros( SVector{ 2, T}, N), zeros( SVector{ 2, T}, N), 
                 zeros( SVector{ 2, T}, N), zeros(Float64, 
                 stepbskip * ss.skip, Nscalar))
end

include( ".//GetForce.jl")
include( ".//MoveParts.jl")
include( ".//StartRun.jl")

println("MD2 included.")

end
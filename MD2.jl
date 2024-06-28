module MD2 
# Two dimensional Molecular Dynamics

export MDRecord, ttt, InitMDRecord!, CorrVec, CorrShell

struct MDRecord{T<:Real}
    N::Int
    tau::T # Time step
    shell::T # Radius of involved particles.
    forceFun::Function # The force function
    skip::Int # Number of steps to skip between updating involved particles.
    bsize::T # Box size of the x axis {-bsize to bsize}
    boxratio::T # size of the y axis: boxratio*{-bsize to bsize}, area=4*boxratio*bsize^2
    p::Matrix{T}
    v::Matrix{T}
    function MDRecord{s}(N::Int, tau::s, shell::s,  forceFun::Function, skip::Int, bsize::s, 
        boxratio::s, p::Matrix{s}, v::Matrix{s}) where s <: Real
        if !( size( p) === size( v) === (N, 2))
            error("invalid matrix sizes")
        end
        new(N, tau, shell, forceFun, skip, bsize, boxratio, p, v)
    end
end

function MDRecord( ::Type{T}, N::Integer, tau, shell, forceFun::Function, 
    skip, size, boxratio) where {T<:Real}
    MDRecord( N, T(tau), T(shell), forceFun, Int(skip), T(size), T(boxratio), 
    zeros(T,N,2), zeros(T,N,2))
end

struct CorrVec
    Mmax::Int
    m::Vector{Int} # Hold the number of defined items for each particle
    corr::Matrix{Int} # corr[i,:] hold a list of interacting particles for particle i
    function CorrVec( Mmax::Int, m::Vector{s}, corr::Matrix{s}) where s <: Real
        # Note that corr is at twice as big as it needs to be.  corr[N-1] need only be length 1
        if !( size( corr, 1) === (size(m), Mmax)) error("corr must be size (size(m), Mmax)") end
        new( Mmax, m , corr)
    end
end

CorrVec(N::Int, Mmax::Int) = CorVec( Mmax, zeros(Int, N), zeros( Int, N, Minit))

struct ttt{q<:Real}
    i::Int
    w::q
    f::Function
    x::Matrix{q}
    y::Matrix{q}
    function ttt{s}(i::Int, w::s, f::Function, x::Matrix{s}, y::Matrix{s}) where s<:Real
        if !(size(x)===size(y)===(i,2))
            error("failed")
        end
        new(i,w,f,x,y)
    end
end

function InitMDRecord!( ss::MDRecord{<:Real})
    for i in 1:ss.N
        ss.p[ i, 1] = 2.0 * mod( i, ss.bsize) - ss.bsize
        ss.p[ i, 2] = floor( i / ss.bsize) - ss.bsize
    end
    ss.v .= (0.1 * ss.bsize) .* randn( eltype( ss.v), ss.N, 2)
    ss
end

function CorrShell( ss::MDRecord{T}, pn::Matrix{T}, corr::CorrVec) where T<:Real
    # Find all particles within shell
    # corr should be declared in the calling program, be size (ss.N-1),
    # where M is large enough to hold the biggest list of interacting particles possible.
    # The m vector will be reset, and if the corr vectors have left over values, 
    # they will be ignored.
    shell = ss.shell^2 # we will use the squared distance as the threshold
    for i = 1:(ss.N - 1) # The last particle will be done by default
        m = 0 # Init m
        for j in i:ss.N
            ds = (pn[i,:] - pn[ j, :]).^2 # distance squared.
            if ds > shell 
                m += 1
                if m > corr.Mmax error( "m > Mmax") end
                corr.corr[ i, m] = j
            end
            corr.m = m # this will be zero if no additional particles were added.
        end
    end
end

end
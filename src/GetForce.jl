# This file should be included in MD2

function gfij( Pi::SVector{2,T}, Pj::SVector{2,T}, 
                st::MD2State{T}) where T <: Real
    vd = pbc( Pi - Pj, st.ss.bx, st.ss.by)  
    r = LinearAlgebra.BLAS.nrm2( 2, vd, 1) # scalar distance
    r > st.ss.rskn && return (false, T(0.0), zero( typeof(vd)))
    #forcevec = st.forceFun( r) * vd
    return (true, st.peFun(r), st.forceFun( r) * vd)
end

"""
    GetForce!( st::MD2State)

GetForce return a Vector of two element StaticArray vectors with the force 
vector on each particle.
"""
function GetForce!( st::MD2State{T}) where T <: Real
    #st.forces .= zeros( eltype(st.forces), st.ss.N) # zero the force vectors.
    fill!( st.forces, zero( st.forces[1]))
    pe::T = 0.0 # init the pe accumulator.
    psi = complex( T(0))
    psicnt::Int32 = 0
     @fastmath @inbounds for i::Int32 in range(1, stop = st.ss.N - 1) # Step through particles.
        mfirst::Int32 = st.cv.m[ i] + 1 # cv.m[1] will be zero.
        mlast::Int32  = st.cv.m[ i + 1]
        if mlast >= mfirst # high particle numbers may have 
                           # no additional interactions.
            mrange = st.cv.corr[ mfirst:mlast] # list of interacting particles.
            @floop begin #@floop 
                petmp::T = 0.0
                for j::Int32 in mrange
                    vd = pbc( st.pn[i] - st.pn[j], st.ss.bx, st.ss.by)  
                    # Vector distance to closest particle image
                    r::T = LinearAlgebra.BLAS.nrm2( 2, vd, 1) # scalar distance
                    #(peflg, peij, fij) = gfij( st.pn[i], st.pn[j], st)
                    if  r <= st.ss.rskn
                        petmp += st.peFun(r)
                        fij::SVector{2,T} = st.forceFun( r) * vd
                        st.forces[ i] += fij
                        st.forces[ j] -= fij
                        if r <= st.NHOPr
                            psi += cis( 6 * atan(vd[2], vd[1]))
                            psicnt += 1
                        end
                    end # if fij ...
            end # for j ...
        end #FLoops@floop...
            pe += petmp
        end # if mlast ...
    end # for i ..
    return (pe, psi / psicnt)
end # GetForce

"""
    rn( n::Int)

Returns a named tuple with 1/r^n PE and force functions expressed as strings:
a potentiol energy function and a force function.
The force function is the negative derivative of the potential energy function 
with an additional 1/r factor.  n must be greater than or equal to 1.
"""
function rn( T::Type, n::Int)
    n >= 1 ? 
    (PE = "r::$T -> r^(-$n)", 
    force = "r::$T -> $n * r^(-$n-2)") :
    error( "rn: n must be > 1")
end

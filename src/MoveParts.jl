# Should be included in MD2.jl

"""
    MoveOnePart( pn::SVector{2, T}, pnm1::SVector{2, T}, forces::SVector{2, T}, 
    ts::T) where T <: Real

Moves one particle. ts is Tau^2
"""
function MoveOnePart( pn::SVector{2, T}, pnm1::SVector{2, T}, forces::SVector{2, T}, 
    ts::T) where T <: Real
    (2 .* pn) .- pnm1 .+ (forces .* ts)
end

"""
    MoveParts!( pn::Vector{ SVector{2, T}}, pnm1::Vector{ SVector{2, T}}, 
    forces::Vector{ SVector{2, T}}, ss::MDRecord{T}) where T<:Real

Moves all the particles according the forces, updating pn and pnm1 in place.
returns total linear momentum and total ke.
"""
function MoveParts!( st::MD2State{T}) where T <: Real
    @floop begin
        @fastmath @inbounds for i in eachindex(st.pn)
            st.pnp1[i] = MoveOnePart( st.pn[i], st.pnm1[i], st.forces[i], st.ss.tau^2)
            st.ss.v[i] = st.oo2t .* (st.pnp1[i] - st.pnm1[i]) # distance over 2 time steps
            delta      = pbcDelta( st.pn[i], st.ss.bx, st.ss.by)
            st.pnm1[i] = st.pn[i]   - delta # pnm1 positions may be outside the box.
            st.pn[i]   = st.pnp1[i] - delta
        end
    end
    momentum = sum( st.ss.v) # total momentum (SVector)
    kt = T(0.5) * sum( dot.(st.ss.v, st.ss.v)) # Kinetic Energy
    return (momentum, kt)
end


global Nscalars::Int32 = 10

function ExecuteSteps!( st::MD2State; animate=false, anim)
    k = 0
    @inbounds for i in 1:st.stepbskip
        MD2.CorrShell!( st.cv, st.ss, st.pn)
        for j in 1:st.ss.skip
            k += 1
            (pt, psi) = MD2.GetForce!( st)
            (momentum, kt) = MD2.MoveParts!( st)
            if animate
                MD2.PlotMDRecord( st.pn, st.ss, vecs=false, margin = 1.1, show = false) 
                Plots.frame( anim)
            end
            ke = kt / st.N
            pe = pt / st.N
            pressure = st.density * ke * (1 + (st.ss.n * pt / (2 * kt)))
            st.scalars[ k, :] = [k, momentum[1], momentum[2], 
                dot( momentum, momentum), ke, pe, ke + pe, pressure, 
                real(psi), imag(psi)]
        end # for j ...
    end # for i ...
    return st.scalars
end
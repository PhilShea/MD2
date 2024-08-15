module TestMove
include( "..\\src\\MD2.jl")

using Test

using .MD2
using Random
using StaticArrays
using Plots
using DataFrames
using CSV
using Profile
using ProfileView
using InteractiveUtils

export ss, TestMoveParts, pnm1

ss = MD2.MDRecord{Float64}(9, 5, 0.01, 1.1, 1.0, 2, 1.5, 1.5)
Random.seed!(31394)
MD2.InitMDRecord!( ss)
MD2.PlotMDRecord( ss; margin=2.0)
state = MD2.MD2State( ss, 40, 1000, MD2.Nscalars)
@time MD2.MoveParts!(state)

# plot the vectors
scatter!( [s[1] for s in state.pnm1], [s[2] for s in state.pnm1])
v0 = [ss.p[i] - state.pnm1[i] for i in 1:ss.N]
@test isapprox((v0[1] ./ ss.tau), ss.v[1])

println("Press Enter to continue:")
readline() # pause to see graph.
global N = 256
global bx = sqrt(N * pi / (2 *sqrt(3)))
funs = MD2.rn(Float64, 5)

@show funs
ss = MD2.MDRecord{Float64}( N, 5, 0.01, 5.0, 3.5, 20, bx, sqrt(3) * bx / 2)
MD2.InitMDRecord!( ss)
MD2.PlotMDRecord( ss; margin=2.0)
println("Press Enter to continue:")
readline() # pause to see graph.

function TestMoveParts(ss::MDRecord{T}, steps::Int; 
                        animate=false) where T <: Real
    #cv = CorrVec( ss.N, 24*ss.N)
    #pn = ss.p
    animate || (anim = 0)
    animate && (anim = Animation()) # if animate is true...
    #step2 = cld( steps, ss.skip)
    #@show ss.skip * step2
    #scalars = zeros(Float64, ss.skip * step2, 6)
    #forces = zeros( SVector{2,T}, ss.N)
    density = ss.N/(4 * ss.bx * ss.by)
    pefcsize = round(Int, 1.5 * pi * density * ss.shell^2)
    @show pefcsize
    #pefc = zeros( SVector{3,T}, pefcsize)
    #pefcFlag = fill( false, pefcsize)
    state = MD2.MD2State( ss, pefcsize*ss.N, steps, MD2.Nscalars)
    @time MD2.ExecuteSteps!( state, animate=animate, anim=anim) 
    animate && gif(anim, "test2 gif.gif", fps=15)
    df = DataFrame( state.scalars, ["step", "xMomentum", "yMomentum", "M^2", 
                    "ke", "pe", "TE", "Pressure", "Rpsi", "Impsi"])
    MD2.writeScalars( "test", df)
    plot( state.scalars[ :, 4:end], labels=permutedims( names(df[ :, 4:end])),
          show=true)
end
@time TestMoveParts( ss, 100)

@time TestMoveParts( ss, 1000, animate=true)

#=
println("Test at float32")
ss32 = MD2.MDRecord{Float32}( N, 5, 0.01, 5.0, 3.5,funs[2], funs[1], 20, bx, sqrt(3) * bx / 2)
MD2.InitMDRecord!( ss32)

@time TestMoveParts( ss32, MD2.initpnm1( ss32), 10)

@time TestMoveParts( ss32, MD2.initpnm1( ss32), 1000)
=#

end
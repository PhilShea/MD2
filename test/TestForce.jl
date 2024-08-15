module TestForce
include( "..\\src\\MD2.jl")

using Test

using .MD2
using Random
using StaticArrays
using InteractiveUtils

const T = Float64 # set to Float 32 or 64
println( "Testing rn")
funs = MD2.rn( T, 5)
@show funs
PEfun = eval( Meta.parse( funs.PE))
forcefun = eval( Meta.parse( funs.force))
@test isapprox( PEfun( T(1.0)), 1)
@test isapprox( PEfun( T(0.5)), 32)
@test isapprox( PEfun( T(2.0)), 1/32)
println( forcefun( T( 1)))
@time( forcefun( T( 1)))
#println("Testing gfij; test 1 as $T")

ss = MD2.MDRecord{T}( 9, 5, 0.01, 1.1, 1.0, 2, 1.5, 1.5)
Random.seed!(31394)
MD2.InitMDRecord!( ss, Temp=0.0001)

MD2.PlotMDRecord( ss; margin=2.0)
state = MD2.MD2State( ss, 4*9, 10, 6)
MD2.CorrShell!( state.cv, ss, state.pn)
#f = MD2.gfij( ss.p[1], ss.p[2], state)
#@show ss.p[1], ss.p[2], f, typeof( f)

#println("gfij test 2")

#@time (f, peij, force) = MD2.gfij( state.pn[1], state.pn[2], state)
#@show state.pn[1] state.pn[2] f
#@test isapprox( force, SVector( 0.0, -5.0))
# pbc will cause the force from p3 to be in the positive y direction.
#println("gfij test 3")
#@time (f, peij, force) = MD2.gfij( state.pn[1], state.pn[3], state)
#@show state.pn[1] state.pn[3] f
#@test isapprox( force, SVector(0.0, 5.0))

println("Testing GetForce")
state.pn .= [SVector(tuple(T(i),T(j))) for i in 1.0:3.0 for j in 1.0:3.0]
MD2.PlotMDRecord( state.pn, ss, vecs=false)
MD2.CorrShell!( state.cv, ss, state.pn)
@show state.cv state.pn
(pe, psi) = MD2.GetForce!( state)
@show pe, psi
#InteractiveUtils.@code_warntype MD2.GetForce!( state)
@show state.forces
@test pe == 18
@test isapprox( state.forces, 
    SVector{2, Float32}[[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], 
                        [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], 
                        [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
end
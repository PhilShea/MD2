module TestMD2

include("..\\src\\MD2.jl")
using .MD2

using Test
using Plots
using Random   
using StaticArrays
using LinearAlgebra

x = MD2.MDRecord{Float64}( 9, 5, 0, 0.01, 1.1, 1.0, "x->x^-5", "x->6 * x^-7", 2, 9.0, 9 * sqrt(3)/2,
    zeros( SVector{2,Float64}, 9), zeros( SVector{2,Float64}, 9))
println( "Test 1: MDRecord main constructor.")
@test (size(x.p)===(9,))
println( "subtest 1: typeof x: ", typeof( x.p), "; x is correct size: ", size(x.p))
@test (eltype( x.p[1])===Float64)
println( x)

Random.seed!(31394)
MD2.InitMDRecord!( x)
s = SVector{2, Float64}[[-7.5, -5.19615], [-1.5, -5.19615], [4.5, -5.19615],
                        [-4.5, 0.], [1.5, 0.], [7.5, 0.], 
                        [-7.5, 5.19615], [-1.5, 5.19615], [4.5, 5.19615]]
println( x.p)
@test isapprox( x.p, s, atol=1e-5)
@show norm( x.p - s)

println("Testing Angular Momentum")
@test isapprox(MD2.TDAM( SVector{2, Float64}( 1, 1), SVector{2, Float64}( -1, 1)), 2)
pt = SVector{2, Float64}[ [ 1.0, 0.0], [ 0.0, 1.0], [ -1.0, 0.0], [ 0.0, -1.0]]
vt = SVector{2, Float64}[[  0.0, 1.0], [ -1.0, 0.0], [ 0.0, -1.0], [ 1.0, 0.0]]
@test isapprox(MD2.AngularMomentum( pt, vt), 4.0)

println("Testing Float32 versions")

x2 = MD2.MDRecord{Float32}(9, 5, 0.01, 1.1, 1.0, 2, 9.0, 9 * sqrt(3)/2)
@test (size( x2.p)===(9, ))
@test (eltype( x2.p[1])===Float32)

# Random.seed!(31394) # DP floats have different series of random numbers than single precision.
MD2.InitMDRecord!( x2)
println( x2.p)
println( x2.v)
@test (isapprox( x2.p, s)) # note that this is true anyway
MD2.PlotMDRecord( x2)
println("Press Enter to continue:")
readline() # pause to see graph.

x3 = MD2.MDRecord{Float32}( 9, 5, 0.01, 1.1, 1.0, 2, 1.5, 1.5)
x3.p .= [SVector(tuple(i,j)) for i in (-1.0):1.0 for j in (-1.0):1.0]
scatter( [s[1] for s in x3.p], [s[2] for s in x3.p], show=true)

# remember that the CorrVec vectors are length N-1
cv = MD2.CorrVec( 4*9, zeros( Int, 9), zeros( Int, 4*9))
cv2 = MD2.CorrVec( 9, 4*9)
@test (cv.Mmax==cv2.Mmax) && (cv.m==cv2.m) && (cv.corr==cv2.corr)

println("Testing CorrShell!")
MD2.CorrShell!( cv, x3, x3.p)
@show cv
@test (cv.m == cumsum([0, 4, 3, 2, 3, 2, 1, 2, 1]))
# this next line is easily calculated by hand.see
ca =  [ 2, 3, 4, 7, 3, 5, 8, 6, 9, 5, 6, 7, 6, 8, 9, 8, 9, 9]

io = open( "testfile.txt", "w")
MD2.InitMDRecord!( x2, io=io)
close( io)

println("Tests Finished")

end
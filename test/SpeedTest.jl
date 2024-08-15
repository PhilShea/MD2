module SpeedTest

include( "C:\\Users\\phils\\Documents\\DSP\\MD\\dev\\MD2\\src\\MD2.jl")
import .MD2
using Plots
using Random


function SpeedTest1( N::Int, density::T; tau = 0.01, boxratio = 1.0, shellM::Int=10, 
    forceFun::Function = r -> r ^ -5, skip::Int = 10) where T<:Real
    bsize = convert(T, sqrt( N / (density * boxratio)))
    shell = convert(T, sqrt( shellM / (pi * density)))
    println( "shell: ", shell, " bsize: ", bsize)

    ss = MD2.MDRecord( T, N, T(tau), T(shell), forceFun, skip, T(bsize), T(boxratio))
    Random.seed!(31394)
    MD2.InitMDRecord!( ss)
    MD2.PlotMDRecord( ss)

    corr = MD2.CorrVec( N, 5 * shellM)

    @time MD2.CorrShell!( corr, ss, ss.p)
end

end
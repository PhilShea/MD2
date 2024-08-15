# should be included in MD2.

"""
    InitSystem( filename::AbstractString; T::Type = Float64, N::Integer = 256, 
    density::Number = 1 / pi, n::Integer = 5, tau::Number = 0.01, 
    shell::Number = 5.0, rskn::Number = 3.5, skip::Integer = 50, temp = 0.30)

Creates a new saved MDRecord record for a new series of runs.
"""
function InitSystem( filename::AbstractString; T::Type = Float64, 
        N::Integer = 256, density::Number = 1 / pi, n::Integer = 5,
        tau::Number = 0.01, shell::Number = 5.0, rskn::Number = 3.5, 
        skip::Integer = 50, temp = 0.30, io=stdout)
    bx = bx = sqrt( N / (2 * density * sqrt(3)))
    by = sqrt(3) * bx / 2
    ss = MDRecord{T}( N, n, tau, shell, rskn, skip, bx, by)
    InitMDRecord!( ss, Temp = temp, io=io)
    saveMDRecord( filename, ss)
end

"""
    writeScalars( filename::AbstractString, df::DataFrame)

Since a fair amount of computation time has already been devoted, we don't 
want the write to create a duplicate filename, or otherwise fail, so we will 
look to write the file with a modified name.
"""
function writeScalars( filename::AbstractString, df::DataFrame)
    testname = filename
    fname = testname * ".scalars.csv"
    i = 0
    while isfile( fname)
        i += 1
        testname = filename * Printf.@sprintf "(%d)" i
        fname = testname * ".scalars.csv"
    end
    CSV.write( fname, df)
end

"""
    RunSystem( readfile::AbstractString, savefile::AbstractString, 
    steps::Integer; animate = false, io=dtdout)

Reads the record at "readfile.jld2" (note that extensions added and will be 
duplicated if provided), runs the number of steps, saves the scalars in file 
names "readfile scalars.jld2".  The new state is recorded in 
"savefile.jld2".  If animate is true, an animation is saved to "readfile.gif".
"""
function RunSystem( readfile::AbstractString, savefile::AbstractString, 
                    steps::Integer; animate = false, plotflag=false, io=stdout)
    println( pwd())
    ss = readMDRecord( readfile * ".jld2", io=io)
    animate || (anim = 0)
    animate && (anim = Animation()) # if animate is true...
    density = ss.N/(4 * ss.bx * ss.by)
    pefcsize = round(Int, 1.5 * pi * density * ss.shell^2)
    Mmax = pefcsize * ss.N
    state = MD2.MD2State( ss, Mmax, steps, MD2.Nscalars, io=io)
    @time MD2.ExecuteSteps!( state, animate=animate, anim=anim) 
    new = MD2.updateMDRecord( 
        state.ss, ss.totalsteps + state.stepbskip * ss.skip)
    saveMDRecord( savefile, new)
    animate && gif(anim, readfile, fps=15)
    df = DataFrame( state.scalars, ["step", "xMomentum", "yMomentum", "M^2", 
                    "ke", "pe", "TE", "Pressure", "Rpsi", "Impsi"])
    writeScalars( readfile, df)
    plotflag && plot( Matrix(df[ :, 4:end]), 
                     labels=permutedims( names(df[ :, 4:end])), show=true)
end
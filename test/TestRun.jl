module TestRun
include( "..\\src\\MD2.jl")

using .MD2

MD2.InitSystem( "test1")
MD2.RunSystem( "test1", "test2", 100)
MD2.RunSystem( "test2", "test3", 1000)

end
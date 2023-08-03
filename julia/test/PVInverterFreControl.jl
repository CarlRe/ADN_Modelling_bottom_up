using Test: @test,@testset, @test_throws
using PowerDynamics: PowerPerturbation, simulate, SwingEqLVS, SlackAlgebraic, StaticLine, PowerGrid, State,systemsize, FieldUpdateError


    nodes = [SwingEqLVS(H=1., P=-1.0, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1)]
    lines = [StaticLine(Y=-50*im, from=1, to=2)]
    grid = PowerGrid(nodes, lines)
    state = State(grid, rand(systemsize(grid)))

    pp = PowerPerturbation(node = 1, fault_power = -0.9, tspan_fault = (0.1,1))
    sol = simulate(pp, state, (0., 1.))
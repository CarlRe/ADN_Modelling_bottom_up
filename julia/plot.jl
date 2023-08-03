using PowerDynamics
using Plots

function create_plot(sol)
    generator_indices = findall(bus -> typeof(bus) == SwingEqLVS,powergrid.nodes)
    #labels = reshape(generator_indices,(1,length(generator_indices)))

    pl_v = plot(sol, generator_indices, :v, legend = (0.8, 0.7), ylabel="V [p.u.]")
    pl_p = plot(sol, generator_indices, :p,ylims=(-10.,10.) ,legend = (0.8, 0.7), ylabel="p [p.u.]" )
    pl_q = plot(sol, generator_indices, :q, ylims=(-10.,10.), legend = (0.8, 0.7), ylabel="q [p.u.]")
    pl_ω = plot(sol, generator_indices, :ω,ylims=(-5.,5.), legend = (0.8, 0.7), ylabel="ω [rad/s]")
    
    pl = plot( pl_v, pl_ω, pl_p, pl_q;
            layout=(2,2),
            size = (1000, 500),
            lw=3,
            xlabel="t[s]")
end
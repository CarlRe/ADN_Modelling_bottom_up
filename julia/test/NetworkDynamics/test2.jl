
#=
function swing_equation!(dv, v, e_s, e_d, P, H, t)
    dv[1] = v[2]
    dv[2] = P - P * v[2]
    for e in e_s
        dv[2] -= e[1]
    end
    for e in e_d
        dv[2] += e[1]
    end
end
=#
function swing_equation!(dv, v, e_s, e_d, P, H, t)
    current = total_current(edges)
    voltage = v[1] + v[2] * im
    dv[3] = sv.P - sv.D* v[3] + real(voltage * conj(current))
    dvolt = 1.0im * v[3] * voltage - (abs(voltage) - 1) * voltage
    dv[1] = real(dvolt)
    dv[2] = imag(dvolt)
end

function total_current(edges)
    # Keeping with the convention of negative sign for outging current
    current = 0.0im
    for e in edges
        current -= e[1] + e[2] * im
    end
    current
end

function SlackAlgebraic(dv, v, e_s, e_d, P, t)
    dv[1] = v[2]
    dv[2] = P[1] - P[2] * v[2]
    for e in e_s
        dv[2] -= e[1]
    end
    for e in e_d
        dv[2] += e[1]
    end
end

function powerflow!(e, v_s, v_d, K, t)
    source_voltage = v_s[1] + v_s[2] * im
    destination_voltage = v_d[1] + v_d[2] * im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = K * (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
end


using Graphs
graph1 = Graphs.SimpleGraphs.watts_strogatz(4, 2, 0.)

using NetworkDynamics

swing_vertex = ODEVertex(f = swing_equation!, dim = 2, sym=[:θ, :ω])
powerflow_edge = StaticEdge(f = powerflow!, dim = 1)

nd = network_dynamics(swing_vertex, powerflow_edge, graph1)

K = 6.0
P = [1.,-1.,1.,-1.]
p = (P,K)

u0 = find_fixpoint(nd, p, zeros(8))

using StochasticDiffEq, OrdinaryDiffEq
ode_prob = ODEProblem(nd, u0, (0.,500.), p)
ode_sol = solve(ode_prob, Tsit5())

using Plots, LaTeXStrings
plot(ode_sol, vars = syms_containing(nd, "ω"), ylabel = L"\omega", legend = false)
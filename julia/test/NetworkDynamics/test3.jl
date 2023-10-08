
using NetworkDynamics
using DifferentialEquations
using Graphs

graph2 = SimpleGraph(2)
function Slack!(dv, v, edges, U, t)
    voltage = v[1] +v[2]*im
    residual = U - voltage 
    dv[1] = real(residual)
    dv[2] = imag(residual)
end

function Swing!(dv, v, edges, P,D, t)
    current = total_current(edges)
    voltage = v[1] + v[2] * im
    dv[3] = P - D * v[3] + real(voltage * conj(current))
    dvolt = 1.0im * v[3] * voltage - (abs(voltage) - 1) * voltage
    dv[1] = real(dvolt)
    dv[2] = imag(dvolt)
end

function Admittance!(e, v_s, v_d, Y, t)
    source_voltage = v_s[1] + v_s[2] * im
    destination_voltage = v_d[1] + v_d[2] * im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = Y* (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
end

Slack_vertex = ODEVertex(f = Slack!,dim=2, sym=[:v_r, :v_i])
Swing_vertex = ODEVertex(f = Swing!,dim=3, sym = [:v_r, :v_i, :ω])
edge = StaticEdge(f=Admittance!,dim=2)
vertex = vcat(Slack_vertex, Swing_vertex)
#u0 = [1.,0.,1.,0.,0.]
nd = network_dynamics(vertex, edge, graph2)
U = 1.
P = 1.
D = 0.1
p = (U,P,D)
u0 = find_fixpoint(nd, p, [1.,0.,1.,0.,0.])
#test_prob = ODEProblem(nd,u0) 
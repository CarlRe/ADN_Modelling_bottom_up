using PowerDynamics
using NetworkDynamics
using OrderedCollections: OrderedDict
import PowerDynamics: dimension, symbolsof, construct_vertex, MassMatrix
using Plots
include("GridFollowingInverter.jl")


@DynamicNode MySwingEq(P, H, D, Ω) begin
MassMatrix(;m_u = true, m_int = [true])
end begin
@assert D > 0 "damping (D) should be >0"
@assert H > 0 "inertia (H) should be >0"
Ω_H = Ω * 2pi / H
end [[ω, dω]] begin
p = real(u * conj(i))
dϕ = ω # dϕ is only a temp variable that Julia should optimize out
du = u * im * dϕ
dω = (P - D*ω - p)*Ω_H
end

@DynamicNode MyPQAlgebraic(P,Q) begin
    MassMatrix()
end  begin
    # no prep
end [] begin
    s = u*conj(i)
    du = complex(P, Q) - s
end

@DynamicNode ZIP(P0,Q0,A,B,C,D,E,F) begin
    MassMatrix()
end begin    
    V0 = 1.
end [] begin   
    P  = real(u*conj(i))
    Q  = imag(u*conj(i))
    s = u*conj(i)
    #du = -(P0*(A*(abs(u)/V0)^2+B*(abs(u)/V0)+C)) +P -(1*im*Q0*(D*(abs(u)/V0)^2+E*(abs(u)/V0)+F))+Q
    du = s-(P0*(A*(abs(u)/V0)^2+B*(abs(u)/V0)+C))-(1*im*Q0*(D*(abs(u)/V0)^2+E*(abs(u)/V0)+F))
end

@DynamicNode ThirdOrderMachine(P,H,T,X,Xs,B) begin
MassMatrix(m_u = false)
    @assert P!= 0
    @assert T > 0 
    @assert X > 0
    @assert B > 0
end [[x_d, dx_d],[x_q, dx_q], [n,dn]] begin
    Tm = P/B
    p = real(u * conj(i))
    v_d = real(u)
    v_q = imag(u)
    dϕ = ω 
    s = ω -n
    dx_d = -(x_d+(X-Xs)*i_q)/T + s*x_q 
    dx_q = -(x_q-(X-Xs)*i_d)/T -s*x_d 
    i_d = (v_q/Xs)*(B/ω)-(x_q/Xs)
    i_q = -(v_d/Xs)*(B/ω)+(x_d/Xs)
    dn = ((1/(2*H))*(Te-Tm))/B
    Te = x_d*i_d+x_q*i_q 
    Pe = v_d*i_d+v_q*i_q
    du = Pe - real(u*conj(i))
end


#showdefinition(stdout, MySwingEq) |> println

buses = OrderedDict(
   # "bus0" => SlackAlgebraic(U=1),
   "bus0" => SwingEqLVS(H=1., P=1.0, D=1, Ω=50, Γ=20, V=1),
    #"bus1" => ZIP(P0 = 0.2, Q0 = 0.1, A=0.1,B=0.4,C=0.5,D=0.1,E=0.2,F=0.5));
    #"bus1" => MySwingEq(P=1., H=1., D=0.1, Ω=50.));
    #"bus1" => MyPQAlgebraic(P=0.1,Q=0.)
    #"bus2" => FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98,  Ω=50, X_d_dash=0.185, T_q_dash=0.4,  X_q_dash=0.36, P=0.1, H=6.54, E_f= 1),
    "bus1" => GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=40000,K_v=0.8,omega_r=50,V_r=1,P=-0.9,Q_r=-0.1));
    #"bus1" => ThirdOrderMachine(P = -0.1, H = 1, T = 1,X=2, Xs=1,B=50));

#Leitungsparameter 
R = 0.3;
X = 0.3;
#Leitungslänge als variable 
length1 = 1;

#Leitungsbeläge
#Y = 1 / (length*R + 1*im*length*X)

branches=OrderedDict(
   
    "branch1" => StaticLine(from= "bus0", to = "bus1", Y = 1 / (length1*R + 1*im*length1*X))
    #"branch2" => StaticLine(from= "bus0", to = "bus2", Y = 1 / (length1*R + 1*im*length1*X))
    );
    
#
powergrid = PowerGrid(buses, branches)
#operationpoint = find_operationpoint(powergrid,sol_method = :nlsolve)
timespan= (0.0,10.0)
states =  rand(systemsize(powergrid))
state = State(powergrid,states )
fault1 = ChangeInitialConditions(node = "bus0", var=:v, f = Inc(0.1))
fault2 = PowerPerturbation(node = "bus0",fault_power = 0.9, tspan_fault = (4.,5.))
#solution = simulate(fault1, powergrid, operationpoint, timespan)
solution = simulate(fault2, state, timespan)
plot(solution, "bus0", :v)






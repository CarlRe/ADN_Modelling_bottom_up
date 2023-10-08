using .Data
using JuMP
using Ipopt
using DifferentialEquations
using EasyFit
using Sundials

t = Data.t # time 
Vabs = Data.v_abs  #xData
P_meas = Data.p_mean #yData

#smoothed p_meas 

p_mean = movavg(P_meas,1000)
p_mean = p_mean.x

v_mean = movavg(Vabs,1000)
v_mean = v_mean.x

#######
#Test
voltage_interpolation = LinearInterpolation(t[:,1], v_mean)

function active_power_load(out,du, u, p, t)
    V_t = 0.48-0.005*sign(t-0.15)+0.005*sign(t-0.35)   # Voltage as a function of time t (you can define your V(t) function)
    V0 = p[1]   # Reference voltage V0
    Npt = p[2]  # Npt parameter
    Nps = p[3]  # Nps parameter
    P0 = p[4]   # P0 parameter
    Tp = p[5]
    xp, P = u
    out[1] = (-u[1] .* (V_t / V0) .^ Npt + P0 .* (V_t / V0) .^ Nps) / Tp .- du[1]
    out[2] = P .- xp.*(V_t/V0).^Npt
end
# Define initial conditions and parameter values
u0 = [0.5,0.03]  # Initial value for x_p
du0 = [0.,0.]
tspan = (0.0, 4.0)  # Time span for simulation
p = [ 0.9, 1., 1., 0.03, 1.]  # Parameter values
differential_vars = [true,false]
f = DAEFunction(active_power_load,syms=[:xp,:P])
dae_prob = DAEProblem(f, du0, u0, tspan,differential_vars = differential_vars,p)
stepsize = 4/40959
sol = DifferentialEquations.solve(dae_prob, IDA();saveat = stepsize)
P_dae = sol[:P]

model = Model(with_optimizer(Ipopt.Optimizer))
@variable(model, V0 >= 0)
@variable(model, Npt)
@variable(model, Nps)
@variable(model, P0)
@variable(model, Tp)


@NLexpression(model, residuals[i=1:length(t)], P_dae[i] - p_mean[i])

@NLobjective(model, Min, sum(residuals[i]^2 for i in 1:length(t)))

optimize!(model)
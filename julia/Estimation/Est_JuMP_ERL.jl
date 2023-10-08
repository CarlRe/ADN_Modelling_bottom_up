using .Data
using JuMP
using Ipopt
using DifferentialEquations
using EasyFit

t = Data.t # time 
Vabs = Data.v_abs  #xData
P_meas = Data.p_mean #yData

#smoothed p_meas 

p_mean = movavg(P_meas,1000)
p_mean = p_mean.x

v_mean = movavg(Vabs,1000)
v_mean = v_mean.x

function ERL!(du, u, p, t)
    a,b,c = p
    du[1] = (1/a)*(-u[1]+ P0*((V/V0))^b - P0*((V/V0)^c))
end

model = Model(with_optimizer(Ipopt.Optimizer))

@variable(model, a)
@variable(model, b)
@variable(model, c) 


@constraint(model, a + b + c == 1.0)
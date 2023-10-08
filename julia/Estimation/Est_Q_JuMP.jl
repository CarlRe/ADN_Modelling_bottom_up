using JuMP
using Ipopt
using Plots
using EasyFit
using .Data 

t = Data.t # time 
Vabs = Data.v_abs  #xData
Q_meas = Data.q_mean #yData

#smoothed p_meas 

q_mean = movavg(Q_meas,1000)
q_mean = q_mean.x

v_mean = movavg(Vabs,1000)
v_mean = v_mean.x

model = Model(with_optimizer(Ipopt.Optimizer));

@variable(model, a)
@variable(model, b)
@variable(model, c) 

@constraint(model, a + b + c == 1.0)


@objective(model, Min, sum(((-0.0172)*(a + b*(x/0.9) +c*(x/0.9)^2) - y)^2 for (x,y) in zip(v_mean,q_mean)))

optimize!(model)

# generate sim Data with fitted values 

function ZIP_reactive(p,V)
    Qq, Iq, Zq = p
    Q0,V0 = [-0.0172, 0.9] 
    Q =  Q0*(Zq*(V/V0).^2 .+Iq*(V/V0) .+Qq)
    return Q
end
param_fit_q = [value(a),value(b),value(c)]
ZIP_Q = ZIP_reactive(param_fit_q,v_mean)

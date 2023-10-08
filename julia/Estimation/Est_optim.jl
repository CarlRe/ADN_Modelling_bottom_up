using Optim
using .Data 

t = Data.t # time 
Vabs = Data.v_abs  #xData
P_meas = Data.p_mean #yData

function objective(params, data_x, data_y)
    a, b, c = params
    V0 = 0.9
    P0 = 0.03
    return sum((P0*(a + b* x/V0 + c * (x/V0)^2 )- y)^2 for (x, y) in zip(data_x, data_y))
end

function constraint(params)
    a, b ,c  = params
    return a + b + c  - 1.0
end

p0 = [1.,-1., 1.]

result = optimize(params -> objective(params, Vabs, P_meas), constraint, p0, NelderMead())

fitted_params = Optim.minimizer(result)

function ZIP_active(p,V)
    Pp, Ip, Zp = p
    P0,V0 = [0.03, 0.9] 
    P =  P0*(Zp*(V/V0).^2 .+Ip*(V/V0) .+Pp)
    return P
end
param_fit = fitted_params
Data_fit = ZIP_active(param_fit,Vabs)
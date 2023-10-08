using .Data
using EasyFit

t = Data.t # time 
Vabs = Data.v_abs  #xData
P_meas = Data.p_mean #yData

#smoothed p_meas 

p_mean = movavg(P_meas,1000)
p_mean = p_mean.x

v_mean = movavg(Vabs,1000)
v_mean = v_mean.x

function ZIP_active(p,V)
    Pp, Ip, Zp = p
    P0,V0 = [0.03, 0.9] 
    P =  P0*(Zp*(V/V0).^2 .+Ip*(V/V0) .+Pp)
    return P
end

p_fit = Float64[]
for i=1:size(v_mean)[1]
    Pp, Ip, Zp = [-3,10,-6]
    P0,V0 = [0.03, 0.9] 
    push!(p_fit,P0*(Zp*(v_mean[i]/V0).^2 .+Ip*(v_mean[i]/V0) .+Pp))
end

param_fit_p = [-3,10,-6]
ZIP_P = ZIP_active(param_fit_p,v_mean)


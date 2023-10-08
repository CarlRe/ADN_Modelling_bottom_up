using .Data 
using LsqFit
using Plots


Vabs = Data.v_abs
P_meas = Data.p_mean

function ZIP_active(p,V)
    P0, V0, Zp, Ip, Pp = p
    P =  P0*(Zp*(V/V0).^2 .+Ip*(V/V0) .+Pp)
    return P
end

function param_constraint(params)
    sum = 0
    for p in params 
        sum += p
    end
    return sum = 0
end

    
    
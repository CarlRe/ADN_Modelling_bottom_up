using PowerDynamics: SlackAlgebraic,WindTurbineGenType4,PVInverterWithFrequencyControl, CSIMinimal, VoltageDependentLoad, ExponentialRecoveryLoad, FourthOrderEq, PQAlgebraic, PiModelLine, StaticLine, Transformer, PowerGrid, write_powergrid, Json
using OrderedCollections: OrderedDict
import PowerDynamics: PVInverterWithFrequencyControl, WindTurbineGenType4
#include(joinpath(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models", raw"LinearERL.jl"))
include(joinpath(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\InductionMotor.jl"))
buses = OrderedDict(
    "bus0" => SlackAlgebraic(U=1),
    "bus1" => ExponentialRecoveryLoad(P0=-0.1, Q0=-0.5, Nps=0.6,Npt = 1,Nqs=1,Nqt=1,
                                      Tp=0.5,Tq=0.5,V0=1.),#Last 20kv
    "bus2" => FourthOrderEq(T_d_dash=6.1, D=0.5, X_d=0.5, 
                             X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.1, 
                             X_q_dash=0.36, P=2, H=1.54, E_f= 1.0),
    "bus3" => CSIMinimal(I_r=0.8),
    "bus4" => ExponentialRecoveryLoad(P0=-0.05, Q0=-0.0025, Nps=0.6, Npt = 3,Nqs=2,Nqt=18,
               Tp=0.01,Tq=0.1,V0=1.), #Last 10kV
    #"bus4" => InductionMotor(P = 0.1, J = 0.1, Tm = 1, B=0.1,ωb=50),
    "bus5" => CSIMinimal(I_r = 0.15), #Solar 10kV
    "bus6" => CSIMinimal(I_r=0.25), #solar 20kv
    "bus7" => FourthOrderEq(T_d_dash=1.1, D=0.5, 
    X_d=0.5, X_q=0.98, Ω=50, X_d_dash=0.185, 
    T_q_dash=0.4, X_q_dash=0.36, P=0.1, H=6.54, E_f= 0.2));

#Leitungsparameter 
R = 0.3;
X = 0.3;
#Leitungslänge als variable 
length1 = 2;
length2 = 0.1;
length3 = 1.9;
length4 = 3;
length5 = 3;
length6 = 3;
length7 = 3;
#Leitungsbeläge
#Y = 1 / (length*R + 1*im*length*X)

branches=OrderedDict(
   
    "branch1" => StaticLine(from= "bus0", to = "bus1", Y = 1 / (length1*R + 1*im*length1*X)),
    "branch2" => StaticLine(from= "bus0", to = "bus2", Y =  1 / (length2*R + 1*im*length2*X)),
    "branch3" => StaticLine(from= "bus0", to = "bus3", Y =  1 / (length3*R + 1*im*length3*X)),
    "branch4" => Transformer(from= "bus0", to = "bus4", y =  1 / (length4*R + 1*im*length4*X), t_ratio = 0.5),
    "branch5" => Transformer(from= "bus0", to = "bus5", y =  1 / (length5*R + 1*im*length5*X),t_ratio = 0.5),
    "branch6" => StaticLine(from = "bus0", to = "bus6",Y =  1 / (length6*R + 1*im*length6*X)),
    "branch7" => Transformer(from= "bus0", to ="bus7", y =  1 / (length7*R + 1*im*length7*X), t_ratio = 0.5));
#

powergrid = PowerGrid(buses, branches)
write_powergrid(powergrid, joinpath(@__DIR__,"grid.json"), Json)


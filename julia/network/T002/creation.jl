using PowerDynamics: SlackAlgebraic, VoltageDependentLoad, ExponentialRecoveryLoad, FourthOrderEq, PQAlgebraic, PiModelLine, StaticLine, Transformer, PowerGrid, write_powergrid, Json
using OrderedCollections: OrderedDict


buses = OrderedDict(
    #"bus1" => SlackAlgebraic(U=1),    #UW Osterburken Trafo 2
   # "bus2" => PQAlgebraic(P=-0.11,Q=-0.05),
    #"bus3" => VoltageDependentLoad(P=0.01,Q=0.01,U=0.1,A=0.01,B=0.01));
    #"bus3" => ExponentialRecoveryLoad(P0=-0.1, Q0=0.01, Nps=0.6,Npt = 3,Nqs=2,Nqt=18,Tp=0.01,Tq=0.1,V0=1.),
    #"bus4" => FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=0.51, H=6.54, E_f= 1),
    #"bus4" => PQAlgebraic(P=-0.0,Q=-0.0),
    #"bus5" => PQAlgebraic(P=-0.3,Q=-0.0));
    "bus1" => SlackAlgebraic(U=1),
    "bus2" => FourthOrderEq(T_d_dash=1.1, D=0.5, X_d=0.5, X_q=0.98, Ω=50, X_d_dash=0.185, T_q_dash=0.4, X_q_dash=0.36, P=0.2, H=6.54, E_f= 1),
    "bus3" => ExponentialRecoveryLoad(P0=-0.1, Q0=0.01, Nps=0.6,Npt = 3,Nqs=2,Nqt=18,Tp=0.01,Tq=0.1,V0=1.),
    "bus4" => ExponentialRecoveryLoad(P0=-0.2, Q0=0.025, Nps=0.6,Npt = 3,Nqs=2,Nqt=18,Tp=0.01,Tq=0.1,V0=1.),
    "bus5"=>  ExponentialRecoveryLoad(P0=-0.2, Q0=0.03, Nps=0.6,Npt = 3,Nqs=2,Nqt=18,Tp=0.01,Tq=0.1,V0=1.));


branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = 2.37403290367663-1im*1.07952208168899),
    "branch2" => StaticLine(from= "bus1", to = "bus3", Y = 2.37403290367663-1im*1.07952208168899),
    "branch3" => StaticLine(from= "bus1", to = "bus4", Y = 4.37403290367663-1im*3.07952208168899),
    "branch4" => StaticLine(from= "bus1", to = "bus5", Y = 8.37403290367663-1im*2.07952208168899));
#

powergrid = PowerGrid(buses, branches)
write_powergrid(powergrid, joinpath(@__DIR__,"gridT002.json"), Json)
using PowerDynamics: SlackAlgebraic, FourthOrderEq, SwingEqLVS, PQAlgebraic, PiModelLine, StaticLine, Transformer, PowerGrid, write_powergrid, Json
using OrderedCollections: OrderedDict


buses = OrderedDict(
    "bus1" => SlackAlgebraic(U=1),   #UW Osterburken Trafo1
    "bus2" => PQAlgebraic(P=-0.2,Q=-0.01), # Osterburken West
    "bus3" => PQAlgebraic(P=-0.35, Q=-0.002), # Osterburken Ost
    "bus4" => PQAlgebraic(P=-0.45,Q=-0.02), #Süden aggregiert
    "bus5" => FourthOrderEq(T_d_dash=4.75, D=2, X_d=1.25, X_q=1.22, Ω=50, X_d_dash=0.232, T_q_dash=1.6, X_q_dash=0.715, P=1.0, H=5.06, E_f= 1));
    



branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = 4.37403290367663-1im*5.07952208168899),
    "branch2" => StaticLine(from = "bus1", to = "bus3", Y = 6.38839373113106-1im*7.41877981679735),
    "branch3" => StaticLine(from = "bus1", to =  "bus4", Y = 0.701451236762956-1im*0.814588533015046),
    "branch4" => StaticLine(from = "bus4", to =  "bus5", Y = 0.701451236762956-1im*0.814588533015046)); #fiktionale Leitung zu Synchrongenerator
   

#

powergrid = PowerGrid(buses, branches)
write_powergrid(powergrid, joinpath(@__DIR__,"gridT001.json"), Json)





using PowerDynamics

@DynamicNode InductionMotor(P,J,Tm,B,ωb) begin
    MassMatrix()
end begin
    @assert P > 0
    @assert J > 0
    @assert Tm > 0 
    @assert B > 0
    @assert ωb > 0
end [[ω, dω],[Te, dTe]] begin
    TI = P/ω

    dω = (Te - TI)/J
    dTe = -(Te + Tm) - (B/ωb)*((ω - ωb)/Tm)
end

export  InductionMotor
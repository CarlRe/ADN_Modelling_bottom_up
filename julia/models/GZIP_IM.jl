


using PowerDynamics

@DynamicNode GZIP_IM(P0,Q0,V0,p1,p2,p3,q1,q2,q3) begin
    MassMatrix()
end begin
    @assert P > 0
end [[ω, dω],[Te, dTe]] begin
    P_ZIP = P0[p1*(abs(u)/V0)^2+p2*(abs(u)/V0)+p3]
    Q_ZIP = Q0[q1*(abs(u)/V0)^2+q2*(abs(u)/V0)+q3]

end
export GZIP_IM
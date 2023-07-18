# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)


using PowerDynamics

@DynamicNode LinearERL(P0, Q0, Nps, Npt, Nqs, Nqt, Tp, Tq, V0) begin
    MassMatrix(m_int = [true, true])
end  begin
    @assert V0 > 0 "Nominal Voltage should be >0"
    @assert Tp > 0 "Load recovery constant should be >0"
    @assert Tq > 0 "Load recovery constant should be >0"
    @assert Nps > 0 "Steady-state load voltage dependence p-axis should be >0"
    @assert Npt > 0 "Transient load voltage dependence p-axis should be >0"
    @assert Nqs > 0 "Steady-state load voltage dependence q-axis should be >0"
    @assert Nqt > 0 "Transient load voltage dependence p-axis should be >0"

end [[x_p, dx_p],[x_q, dx_q]] begin
    Pd = real(u*conj(i))
    Qd = imag(u*conj(i))

    dx_p = (1/Tp)*(-x_p + P0*((Nps-Npt)/P0)*abs(u) )
    dx_q = (1/Tq)*(-x_q + Q0*((Nqs-Nqt)/V0)*abs(u) )

    du = -Pd + x_p + P0*((abs(u)/V0)^Npt) + im*(-Qd + x_q + Q0*((abs(u)/V0)^Nqt))
end

export LinearERL
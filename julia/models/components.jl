

using NetworkDynamics: ODEVertex
import PowerDynamics: dimension, symbolsof, construct_vertex, MassMatrix

@DynamicNode InductionMotor(P,J,Tm,B,ωb) begin
    @assert P > 0
    @assert J > 0
    @assert Tm > 0 
    @assert B > 0
    @assert ωb > 0
end [[ω, dω],[Te, dTe]] begin
    TI = P/ω
    dω = (Te - TI)/J
    dTe = -(Te + Tm) - (B/ωb)*((ω - ωb)/Tm)
    du = u*im*ω
end



@DynamicNode SchmietendorfOriginal(P_m, γ, α, E_f, X) begin
    @assert γ > 0
end [[ω, dω], ] begin
    s = u * conj(i)
    e_q = abs(u)
    de_q = (E_f - e_q - X * imag(s) / e_q) / α
    dω = P_m - γ * ω - real(s)
    du = (de_q / e_q + 1im * ω ) * u
end

@DynamicNode MySwingEq(P, H, D, Ω) begin
MassMatrix(;m_u = true, m_int = [true])
end begin
@assert D > 0 "damping (D) should be >0"
@assert H > 0 "inertia (H) should be >0"
Ω_H = Ω * 2pi / H
end [[ω, dω]] begin
p = real(u * conj(i))
dϕ = ω # dϕ is only a temp variable that Julia should optimize out
du = u * im * dϕ
dω = (P - D*ω - p)*Ω_H
end

@DynamicNode ZIP(P0,Q0,A,B,C,D,E,F) begin
    MassMatrix(m_u = false)
end begin    
    V0 = 1.
end [] begin   
    P  = real(u*conj(i))
    Q  = imag(u*conj(i))
    s = u*conj(i)
    #du = -(P0*(A*(abs(u)/V0)^2+B*(abs(u)/V0)+C)) +P -(1*im*Q0*(D*(abs(u)/V0)^2+E*(abs(u)/V0)+F))+Q
    du = s-(P0*(A*(abs(u)/V0)^2+B*(abs(u)/V0)+C))-(1*im*Q0*(D*(abs(u)/V0)^2+E*(abs(u)/V0)+F))
end

@DynamicNode MotorThird(P,ω0) begin
MassMatrix(m_u = false, m_int = [true,true,true,true])
end begin 
    X_s= 0.15
    X_m = 2.5
    R_r = 0.05
    R_s = 0.2
    X_r = 0.1
    Kₚ = 0.001
    Kᵢ = 0.02
    H=1.
end   [[E_d, dE_d],[E_q,dE_q],[x_ω, dx_ω],[ω,dω],[θ,dθ]] begin
    v_dq = exp(-im*θ)*u
    v_ds = real(v_dq)     
    v_qs = imag(v_dq)
    #PLL Frequenzschätzung
    dx_ω = -v_qs
    σ = ω0 -dx_ω*Kₚ - Kᵢ*x_ω
    
    dθ = σ
    #Schlupf
    #s = σ - ω

    X_t = X_s+(X_r*X_m)/(X_r+X_m)
    X = X_s+X_m
    T0 = (X_s+X_m)/(R_r)

    i_ds = (1/(R_s^2+X_t^2))*(R_s*(v_ds-E_d)+X_t*(v_qs-E_q))
    i_qs = (1/(R_s^2+X_t^2))*(R_s*(v_qs-E_q)-X_t*(v_ds-E_d))

    dE_d = -(1/T0)*(E_d+(X-X_t)*i_qs)-(ω-1)*E_q
    dE_q = -(1/T0)*(E_q-(X-X_t)*i_ds)+(ω-1)*E_d

    P_el = E_d*i_ds + E_q*i_qs 

    dω = (-1/(2*H))*(P-P_el)
    ω = ω0+dω
    i_out = (i_ds+im*i_qs)*exp(im*θ)
    du = i - i_out
end
#=
https://www.ieh.kit.edu/rd_download/industrial.pdf
=#
@DynamicNode MotorThirdKIT(P,ω0) begin
MassMatrix(m_u = false, m_int = [true,true,true,true])
end begin 
    X_s= 0.15
    X_m = 2.5
    R_r = 0.05
    R_s = 0.2
    X_r = 0.1
    Kₚ = 0.001
    Kᵢ = 0.02
    H=1.
end   [[E_d, dE_d],[E_q,dE_q],[x_ω, dx_ω],[ω,dω],[θ,dθ]] begin
    v_dq = exp(-im*θ)*u
    v_ds = real(v_dq)     
    v_qs = imag(v_dq)
    #PLL Frequenzschätzung
    dx_ω = -v_qs
    σ = ω0 -dx_ω*Kₚ - Kᵢ*x_ω
    
    dθ = σ
    #Schlupf
    #s = σ - ω

    X_t = X_s+(X_r*X_m)/(X_r+X_m)
    X = X_s+X_m
    T0 = (X_s+X_m)/(R_r)

    i_ds = (1/(R_s^2+X_t^2))*(R_s*(v_ds-E_d)+X_t*(v_qs-E_q))
    i_qs = (1/(R_s^2+X_t^2))*(R_s*(v_qs-E_q)-X_t*(v_ds-E_d))

    dE_d = -(1/T0)*(E_d+(X-X_t)*i_qs)-(ω0*s)*E_q
    dE_q = -(1/T0)*(E_q-(X-X_t)*i_ds)+(ω0*s)*E_d

    P_el = E_d*i_ds + E_q*i_qs 

    dω = (ω0/(T_A))*(P_el-P_m)

    i_out = (i_ds+im*i_qs)*exp(im*θ)
    du = i - i_out
end

#=
https://www.sciencedirect.com/science/article/pii/S014206150100059X
=#
@DynamicNode MotorThirdSLip(P,ω0) begin
MassMatrix(m_u = false, m_int = [true,true,true,true])
end begin 
    X_s= 0.15
    X_m = 2.5
    R_r = 0.05
    R_s = 0.2
    X_r = 0.1
    Kₚ = 0.001
    Kᵢ = 0.02
    H=1.
end   [[E_d, dE_d],[E_q,dE_q],[x_ω, dx_ω],[ω,dω],[θ,dθ]] begin
    v_dq = exp(-im*θ)*u
    v_ds = real(v_dq)     
    v_qs = imag(v_dq)
    #PLL Frequenzschätzung
    dx_ω = -v_qs
    σ = ω0 -dx_ω*Kₚ - Kᵢ*x_ω
    
    dθ = σ
    #Schlupf
    s = (ω0-ω)/ω0

    X_t = X_s+(X_r*X_m)/(X_r+X_m)
    X = X_s+X_m
    T0 = (X_r+X_m)/(ω0*R_r)

    i_ds = (1/(R_s^2+X_t^2))*(R_s*(v_ds-E_d)+X_t*(v_qs-E_q))
    i_qs = (1/(R_s^2+X_t^2))*(R_s*(v_qs-E_q)-X_t*(v_ds-E_d))

    dE_d = -(1/T0)*(E_d+(X-X_t)*i_qs)-(ω0*s)*E_q
    dE_q = -(1/T0)*(E_q-(X-X_t)*i_ds)+(ω0*s)*E_d

    P_el = E_d*i_ds + E_q*i_qs 

    dω = (ω0/(2*H))*(P_el-P_m)

    i_out = (i_ds+im*i_qs)*exp(im*θ)
    du = i - i_out
end

#=
https://www.sciencedirect.com/science/article/pii/S014206150100059X
=#
@DynamicNode MotorSlipModel1(P, X_s , X_r ,R_r, E, ) begin
MassMatrix(m_u = false, m_int = [true])
end begin 
   H = 1.
end   [[x_p,dx_p]] begin
n=1
X = X_s+X_r
V = abs(u)
E = sqrt((V^2-E^2)*(R_r/X)^2*1/(x_p^2))
P_m = P*(1-x_p)^n
P_el = x_p*(E^2/R_r)
dx_p = (1/2*H)*(P_m-P_el)


Q_el = (X_s + X_r)*((E*x_p)/R_r)^2
du = u*conj(i)-P_el-im*Q_el
end
 

@DynamicNode Motor(P ,ω0, Xs, Xt, T, H ) begin
MassMatrix(m_u = false, m_int = [true,true,true,true,true])
end begin 
    X=0.1
    R=0.01
    Kₚ = 0.001
    Kᵢ = 0.02
end   [[x_d, dx_d],[x_q, dx_q],[x_ω,dx_ω],[ω, dω],[θ,dθ]] begin
    v_dq = exp(-im*θ)*u
    v_d = real(v_dq)     
    v_q = imag(v_dq)
    #PLL Frequenzschätzung
    dx_ω = -v_q
    σ = ω0 -dx_ω*Kₚ - Kᵢ*x_ω
    # σ ist also die Abweichung(abgeschätzt) von der Nennfrequenz resultierend aus der Spannung
    dθ = σ
    #Schlupf
    s = σ - ω

    i_d = (1/(R^2+Xt^2))*(R*(v_d-x_d)+Xt*(v_q-x_q))
    i_q = (1/(R^2+Xt^2))*(R*(v_q-x_q)-Xt*(v_d-x_d))
    P_el = x_d*i_d + x_q*i_q 
    dx_q = -ω0*s*x_d-(1/T)*(x_q-(Xs-Xt)*i_d)
    dx_d = -ω0*s*x_q-(1/T)*(x_d+(Xs-Xt)*i_q)
    dω = (-1/2*H)*((A*ω^2+B*ω+C)*T-P_el) #Änderung des Roterdrehfelds bei Leistungsabweichung 
    du = i -((i_d-v_d*X)-im*(i_q+v_q*X))*exp(im*θ)
    #i = i -((i_d-v_d*X)-im*(i_q+v_q*X))*exp(im*θ)
    #du = u*conj(i)-P_el
end
 #=
SimplifiedSingleCageInductionMachine
 =#
@DynamicNode MotorSingleCageNew(τ_m0 ,ω0) begin
MassMatrix(m_u = false, m_int = [true,true,true])
end begin 
    A = 0.2
    B = 0.0
    C = 0.
    H= 1.
    R_s = 0.1
    R_r = 0.4
    X_m = 0.1
    X_rr = 0.2
    X_p = 0.1
    v_qr = 0.
    v_dr = 0.
end   [[ψ_dr,dψ_dr],[ψ_qr, dψ_qr],[ωr, dωr]] begin
    #v_dq = exp(-im*θ)*u
    v_ds = real(u)     
    v_qs = imag(u)
    #i_dq = exp(-im*θ)*i
    #i_ds = real(i)     
    #i_qs = imag(i)
    
   
    i_qs =
        1 / (R_s^2 + (ω0 * X_p)^2) * (
            (R_s * v_qs - ω0 * X_p * v_ds) -
            (R_s * ω0 * X_m / X_rr * ψ_dr + ω0 * X_p * ω0 * X_m / X_rr * ψ_qr)
        )

    i_ds =
        1 / (R_s^2 + (ω0 * X_p)^2) * (
            (R_s * v_ds + ω0 * X_p * v_qs) -
            (-R_s * ω0 * X_m / X_rr * ψ_qr + ω0 * X_p * ω0 * X_m / X_rr * ψ_dr)
        )
    i_qr = (ψ_qr - X_m * i_qs) / X_rr # derived from 4th row of (4.5.37) in Krause
    i_dr = (ψ_dr - X_m * i_ds) / X_rr # # derived from 5th row of (4.5.37) in Krause
    τ_e = ψ_qr * i_dr - ψ_dr * i_qr 
    dψ_dr =  (v_qr - (ω0 - ωr) * ψ_dr - R_r * i_qr) # (4.5-25) in Krause
    dψ_qr = (v_dr + (ω0 - ωr) * ψ_qr - R_r * i_dr) # (4.5-26) in Krause
    dωr = 1.0 / (2.0 * H) * (τ_e - τ_m0 * (A * ωr^2 + B * ωr + C)) #
    du = i - (i_ds) - (i_qs)
    #ψ_qs = (1/σ)*(v_ds-R_s*i_ds)
    #ψ_ds =  (1/σ)*(R_s*i_qs-v_qs)  
    
   # i_dr = (ψ_ds-L_ss*i_ds)/L_sr
   # i_qr = (ψ_qs-L_ss*i_qs)/L_sr
    #dψ_dr = -(R*i_dr+s*ψ_qr)
    #dψ_qr = -(R*i_qr-s*ψ_dr)

     

    #P_el = ψ_dr*i_qr-ψ_qr*i_dr 
    #P_m = P*(A*ω^2+B*ω+C)
    #dω = (1/(2*H))*(P_el+P_m)

    #i_out = (((i_ds-v_ds*R_s)+im*(i_qs+v_qs*R_s))*exp(im*θ))
    #du = i -i_out
end
 
@DynamicNode ConverterPLL(P,Q,Ω,Kₚ,Kᵢ,K_ω,Kᵥ,Vᵣ) begin
MassMatrix(m_u = false, m_int = [true,true])
end begin

end [[e_Iomega,de_Iomega],[θ,dθ]] begin
u_dq = exp(-im*θ)*u
u_q = imag(u_dq)
e_omega = -u_q
de_Iomega= e_omega
ω = Ω - Kₚ*e_omega - Kᵢ*e_Iomega
dθ = ω
p = -K_ω*(ω-Ω) + P
q = -Kᵥ*(abs(u_dq)-Vᵣ) +Q
s = p + im*q
du = u*conj(i) - s
end

@DynamicNode CSIPLL2(P,Q,Ω,Kₚ,Kᵢ,K_ω,Vᵣ) begin
MassMatrix(m_u = false, m_int = [true,true])
end begin

end [[e_Iomega,de_Iomega],[θ,dθ]] begin
u_dq = exp(-im*θ)*u
u_q = imag(u_dq)
p_node = real(u*conj(i))

e_omega = -u_q
de_Iomega= e_omega
ω = Ω - Kₚ*e_omega - Kᵢ*e_Iomega
dθ = ω
p = -K_ω*(ω-Ω) + P

I_r = p/real(Vᵣ)

du = i - I_r*exp(im*θ)
end




@DynamicNode GridFollowingTecnalia(tau_u,omega_ini,K_pomega,K_iomega,K_omega,K_v,omega_r,V_r,P,Q_r) begin
 MassMatrix(m_u = false, m_int = [true,true,true,true])
 end begin
 @assert tau_u >= 0
 @assert omega_ini >= 0
 @assert K_pomega >= 0
 @assert K_iomega >= 0
 @assert K_omega >= 0
 @assert K_v >= 0
 end [[u_fil_d,du_fil_d],[u_fil_q,du_fil_q],[e_Iomega,de_Iomega],[theta,dtheta]] begin

 u_dq = exp(-im*theta)*u
 u_q = imag(u_dq)
 du_fil_d = 1/tau_u*(-u_fil_d + real(u_dq))
 du_fil_q = 1/tau_u*(-u_fil_q + imag(u_dq))

 u_fil_dq = u_fil_d + im*u_fil_q

 e_omega = -u_q
 de_Iomega= e_omega
 omega = omega_ini - K_pomega*e_omega - K_iomega*e_Iomega
 dtheta = omega

 p = -K_omega*(omega-omega_r) + P
 q = -K_v*(abs(u_fil_dq)-V_r) + Q_r
 i_fil_dq = (p-im*q)/u_fil_dq

 du = i - i_fil_dq*exp(im*theta)
 end


 @DynamicNode GridFormingTecnalia(omega_r,tau_U, tau_I, tau_V, tau_omega,
 tau_P, tau_Q, n_P, n_Q, K_P, K_Q, P, Q, V_r, R_f, X_f) begin
  MassMatrix(m_u = true, m_int = [true,true,true,true,true])
  end begin
  @assert tau_U >= 0
  @assert tau_I >= 0
  @assert tau_P >= 0
  @assert tau_Q >= 0
  @assert n_P >= 0
  @assert n_Q >= 0
  @assert K_P >= 0
  @assert K_Q >= 0
  @assert V_r >= 0
  @assert R_f >= 0
  @assert X_f >= 0
  end [[u_fil_r,du_fil_r],[u_fil_i,du_fil_i],[i_fil_r,di_fil_r],[i_fil_i,
 di_fil_i],[omega, domega]] begin
 
  du_fil_r = 1/tau_U*(-u_fil_r + real(u))
  du_fil_i = 1/tau_U*(-u_fil_i + imag(u))
 
  di_fil_r = 1/tau_I*(-i_fil_r + real(i))
  di_fil_i = 1/tau_I*(-i_fil_i + imag(i))
 
  u_fil = u_fil_r +1im*u_fil_i
  i_fil = i_fil_r +1im*i_fil_i
 
  p = real(u_fil * conj(i_fil))
  dp = du_fil_r*i_fil_r+u_fil_r*di_fil_r+du_fil_i*i_fil_i+du_fil_i*
 di_fil_i
 
  q = imag(u_fil * conj(i_fil))
  dq = -du_fil_r*i_fil_i-u_fil_r*di_fil_i+du_fil_i*i_fil_r+u_fil_i*
 di_fil_r
 
  domega = 1/tau_P*(omega_r-omega) + K_P/tau_P*(P-p) - K_P/n_P*dp
  v = abs(u_fil)
  dv = 1/tau_Q*(V_r-v) + K_Q/tau_Q*(Q-q) - K_Q/n_Q*dq
  du = u*1im*(omega-omega_r)+dv*(u/v)
end


#=
motor_node = Motor(P = -1., ω0 = 50., Xs = 10, Xt = 0.5, T = 1, H = 10 )
gen_node = SchmietendorfOriginal(P_m=2.32, γ=2sqrt(50/(2*5.148)), α=7.4sqrt(50/(2*5.148)), E_f=1, X=0.8979)
display(dimension(motor_node))
symbolsof(motor_node)
display(dimension(gen_node))
symbolsof(gen_node)
=#













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
MassMatrix(m_u = false, m_int = [true,true,true,true,true])
end begin 
    X_s= 0.5
    X_m = 0.5
    R_r = 0.5
    R_s = 0.5
    X_r = 0.5
    T = 1.
    Kₚ = 0.001
    Kᵢ = 0.02
end   [[E_d, dE_d],[E_q,dE_q],[x_ω, dx_ω],[ω,dω],[θ,dθ]] begin
    v_dq = exp(-im*θ)*u
    v_ds = real(v_dq)     
    v_qs = imag(v_dq)
    #PLL Frequenzschätzung
    dx_ω = -v_qs
    σ = ω0 -dx_ω*Kₚ - Kᵢ*x_ω
    
    dθ = σ
    #Schlupf
    s = σ - ω
    X_t = X_s+(X_r*X_m)/(X_r+X_m)
    X = X_s+X_m
    T0 = (X_s+X_m)/(R_r*ω0)
    #i_ds = (1/(R_s^2+X_t^2))*(R_s*(v_ds-E_d)+X_t*(v_qs-E_q))
    #i_qs = (1/(R_s^2+X_t^2))*(R_s*(v_qs-E_q)+X_t*(v_ds-E_d))
    i_ds = (v_qs/X_t)*(ω0/ω)-(E_q/X_t)
    i_qs = -(v_ds/X_t)*(ω0/ω)+(E_d/X_t)
    dE_d = -(1/T0)*(E_d+(X-X_t)*i_qs)+s*E_q
    dE_q = -(1/T0)*(E_q-(X-X_t)*i_ds)-s*E_d
    P_el = E_d*i_ds + E_q*i_qs 
    dω = (ω0/(2*T))*(P-P_el)

    i_out = (i_ds+im*i_qs)*exp(im*θ)
    du = i -i_out
end
 
@DynamicNode Motor(P ,ω0, Xs, Xt, T, H ) begin
MassMatrix(m_u = false, m_int = [true,true,true,true])
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
    i_q = (1/(R^2+Xt^2))*(R*(v_q-x_q)+Xt*(v_d-x_d))
    P_el = x_d*i_d - x_q*i_q 
    dx_q = -ω0*s*x_d-(1/T)*(x_q-(Xs-Xt)*i_d)
    dx_d = -ω0*s*x_q-(1/T)*(x_d+(Xs-Xt)*i_q)
    dω = (1/H)*(P_el+P) #Änderung des Roterdrehfelds bei Leistungsabweichung 
    du = i -((i_d-v_d*X)-im*(i_q+v_q*X))*exp(im*θ)
    #i = i -((i_d-v_d*X)-im*(i_q+v_q*X))*exp(im*θ)
    #du = u*conj(i)-P_el
end

@DynamicNode MotorSingleCage(P ,ω0, H ) begin
MassMatrix(m_u = false, m_int = [true,true,true,true,true])
end begin 
    A = 0.2
    B = 0.0
    C = 0.8
    R_s = 0.013
    R_r = 0.009
    L_ss = 3.867
    L_sr = 3.800
    L_rr = 3.970
    Kₚ = 0.001
    Kᵢ = 0.02
end   [[ψ_dr,dψ_dr],[ψ_qr, dψ_qr],[x_ω,dx_ω],[ω, dω],[θ,dθ]] begin
    v_dq = exp(-im*θ)*u
    v_ds = real(v_dq)     
    v_qs = imag(v_dq)
    i_dq = exp(-im*θ)*i
    i_ds = real(i_dq)     
    i_qs = imag(i_dq)
    #PLL Frequenzschätzung
    dx_ω = -v_qs
    σ = ω0 -dx_ω*Kₚ - Kᵢ*x_ω
    
    dθ = σ
    #Schlupf
    s = σ - ω
    
    ψ_qs = (1/σ)*(v_ds-R_s*i_ds)
    ψ_ds =  (1/σ)*(R_s*i_qs-v_qs)  
    
    i_dr = (ψ_ds-L_ss*i_ds)/L_sr
    i_qr = (ψ_qs-L_ss*i_qs)/L_sr
    dψ_dr = -(R*i_dr+s*ψ_qr)
    dψ_qr = -(R*i_qr-s*ψ_dr)

     

    P_el = ψ_dr*i_qr-ψ_qr*i_dr 
    P_m = P*(A*ω^2+B*ω+C)
    dω = (1/(2*H))*(P_el+P_m)

    i_out = (((i_ds-v_ds*R_s)+im*(i_qs+v_qs*R_s))*exp(im*θ))
    du = i -i_out
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











module Data

using CSV
using MAT
using FileIO
using Statistics
using Plots
using Polynomials
using SpecialFunctions
using LsqFit
using JLD2
include(raw"FrequencyEstimation.jl")
include(raw"import.jl")
using .FrequencyEstimate

#=
ω  = 2*pi*f = 2*pi*T 
T = 1/f = 1*s/50 = 0.02s  
=#

# measurement bestimmt die spezifische Messung
#data = load("data.JLD2","data")
data = Measurement.data
numberOfData = size(data["single_data"]["t"])[1]
measurement = 99
#time  
t = data["single_data"]["t"][measurement]
# timefactor ist 1 oder 2 und wird in FrequencyEstimation bestimmt
timeFactor = checkTimeFactor(t)
t = t.*timeFactor
ω = 2*pi*checkFrequency(data,t,measurement) # frequency wird aus Messung bestimmmt 

# Park Transform (2x3 da dq0)
#K = [cos(θ) cos(θ-((2*pi)/3)) cos(θ+((2*pi)/3)); -sin(θ) -sin(θ-((2*pi)/3)) -sin(θ+((2*pi)/3)) ]
#K_inverse = [cos(θ) -sin(θ);cos(θ-((2*pi)/3)) -sin(θ-((2*pi)/3));cos(θ-((4*pi)/3)) -sin(θ-((4*pi)/3))]


#I_dq,U_dq = zeros(Float64, 2, size(t)[1])
I_dq = Float64[0.,0.]
U_dq = Float64[0.,0.]
I_dq = zeros(Float64, 2, size(t)[1])
U_dq = zeros(Float64, 2, size(t)[1])
U₀ = [data["single_data"]["U1"][measurement]';data["single_data"]["U2"][measurement]';data["single_data"]["U3"][measurement]'] # [V]
U₀ = U₀.*0.001 # [kV]
I₀ = [data["single_data"]["I1"][measurement]';data["single_data"]["I2"][measurement]';data["single_data"]["I3"][measurement]'].*(-1) # [A]
#I₀ = I₀.*0.001 # [kA]
display(ω)
for i=1:(size(t)[1])
    θ = ω*t[i]
    K = [cos(θ) cos(θ-((2*pi)/3)) cos(θ+((2*pi)/3)); -sin(θ) -sin(θ-((2*pi)/3)) -sin(θ+((2*pi)/3)) ]
    U = U₀[:,i]
    U_dq[:,[i]] = K*U
    #push!(U_dq,K*U)
    I = I₀[:,i]
    I_dq[:,[i]] = K*I
    #push!(I_dq,K*I)
end

I_pu = I_dq.*(1/721.7) #I_base =S_base/(U_base*sqrt(3)) = 25MVA/20kV*sqrt(3) =0.7217kA
U_pu = U_dq.*(0.05) #U_base = 20kV

v_d = U_pu[1,:]
v_q = U_pu[2,:]

i_d = I_pu[1,:]
i_q = I_pu[2,:]

# Initialize variables 
v_abs = Float64[]
i_abs = Float64[]
p_active = Float64[]
p_reactive = Float64[]
i_dq = ComplexF64[]
v_dq = ComplexF64[]
power = ComplexF64[]
#Variable declaration
for i=1:(size(t)[1])
    push!(v_abs,sqrt(v_d[i]^2+v_q[i]^2))    # Voltage Magnitude
    push!(v_dq,(v_d[i]+im*v_q[i]))          # Voltage in dq (Complex)
    push!(i_abs,sqrt(i_d[i]^2+i_q[i]^2))    # Current Magnitude
    push!(i_dq,(i_d[i]+im*i_q[i]))          # Current in dq (Complex)
    push!(power,(v_dq[i]*conj(i_dq[i])))    # Power
end





# Smoothing curve
# window = 200
# Calculate mean of Voltage and Power in three stages
function curveSmoothing(t0, value_v, value_s)
window = 200

v_mean = Float64[]
v_mean_secondstage = Float64[]
v_mean_thirdstage = Float64[]
stages_v = [value_v, v_mean, v_mean_secondstage,v_mean_thirdstage]
s_mean=ComplexF64[]
s_mean_secondstage=ComplexF64[]
s_mean_thirdstage=ComplexF64[]
stages_s= [value_s, s_mean, s_mean_secondstage,s_mean_thirdstage]
    for j = 1:3
       window = 100 +j*100
    for i=1:(size(t0)[1])
        if i > ((size(t0)[1])-window)
            push!(stages_v[j+1],stages_v[j+1][((size(t0)[1])-window)])
            push!(stages_s[j+1],stages_s[j+1][((size(t0)[1])-window)])
        else
            push!(stages_v[j+1],mean(stages_v[j][i:(i+window)]))
            push!(stages_s[j+1],mean(stages_s[j][i:(i+window)]))    
        end
    end
    end
    stages = [stages_v[4],stages_s[4]]
end

means = curveSmoothing(t, v_abs, power)
v_mean = real(means[1])
s_mean = means[2]
p_mean = real(s_mean)
q_mean = imag(s_mean)
data_fit = [v_mean real(s_mean) imag(s_mean)]
#fitting
#=
P = Polynomial
xs = range(t[1],t[size(t)[1]],size(t)[1])
ys = v_mean
v_fit1 = fit(P,xs,ys,20)
=#

#step fitting
if abs(v_mean[1]-v_mean[20000])>0.001
    model = model(t, s) = s[1] .+ s[2]*sign.(t.+s[3])
    s0 = [1.,1.,-0.1]
else
    model = model(t, s) = s[1] .+ s[2]*sign.(t.+s[3]).+ s[4]*sign.(t.+s[5])
    s0 = [1.,1.,-0.1,-1.,-0.5]
end
stepfit_v = curve_fit(model, t[:], v_mean[:],s0)
stepfit_p = curve_fit(model, t[:], real(s_mean)[:],s0)
stepfit_q = curve_fit(model, t[:], imag(s_mean)[:],s0)
#test_fit = curve_fit(nonlinear_fit,t[:],v_mean[:])
#Values for simulation
v_start = model(t[1],stepfit_v.param)
v_end = model(t[30000],stepfit_v.param)
v_fit = Float64[]
p_fit = Float64[]

for i=1:size(t)[1]
    push!(v_fit,model(t[i],stepfit_v.param))
    push!(p_fit,model(t[i],stepfit_p.param))
end

# Step Times
function stepTimes(timeseries, values)
    stepTimes = Float64[]
    for i=1:(size(timeseries)[1]-1)
        if ((values[i+1] >  values[i] || values[i+1] <  values[i]))
            push!(stepTimes, t[i])
        end
    end
    return stepTimes
end
# Step times of step fittet voltage
steps = stepTimes(t,v_fit)
v_min = minimum(v_fit[:])
v_max = maximum(v_fit[:])
v_bounds = [v_min, v_max]
p_min = minimum(p_fit[:])
p_max = maximum(p_fit[:])
display("vmin: "*string(v_min)*" ,vmax: "*string(v_max)*" ,pmin: "*string(p_min)*" ,pmax: "*string(p_max))
#Plotting
timespan = (0.,4.)
plot_abs = plot(t,v_abs,xlims= timespan,label = "v_abs")
plot_mean = plot(t,v_mean,xlims= timespan,label = "v_mean")
plot_stepfit_v = plot(t,model(t,stepfit_v.param),xlims= timespan,label = "v_stepfit")
plot_stepfit_p = plot(t,model(t,stepfit_p.param),xlims= timespan,label = "p_stepfit")
plot_v = plot(t,v_abs,xlims= timespan,label = "v_abs");
            plot!(t,v_mean,linecolor=:red,xlims= timespan,label = "v_mean_thirdstage");
            plot!(t,model(t,stepfit_v.param),linecolor=:orange,xlims= timespan,label = "v_stepfit");
plot_p =    plot(t,real(power),xlims= timespan,label = "p");
            plot!(t,real(s_mean),linecolor=:red,xlims= timespan,label = "p_mean");
            plot!(t,model(t,stepfit_p.param),linecolor=:orange,xlims= timespan,label = "p_stepfit")
plot_q =    plot(t,imag(power),xlims= timespan,label = "q");
            plot!(t,imag(s_mean),linecolor=:red,xlims= timespan,label = "q_mean");
            plot!(t,model(t,stepfit_q.param),linecolor=:orange,xlims= timespan,label = "q_stepfit")
plot_pq= plot(real(power),imag(power),label = "pq");
         plot!(real(s_mean),imag(s_mean),label = "pq_mean")
plot(plot_v,plot_p, plot_q;
        layout=(3,1),
        size = (1500, 1000),
        lw=4,
        plot_title = "EMT Daten in dq0",
        xlabel="t[s]")

#savefig("Preprocessing.png")
save("data.jld2", "v_mean", v_mean )
save("data.jld2", "s_mean", s_mean )
save("data.jld2", "steps", steps )
save("data.jld2", "data_fit", data_fit )

 

export v_mean, p_mean, q_mean, v_max, v_min, steps 

end
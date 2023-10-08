
using CSV
using MAT
using FileIO
using Statistics
using Plots
using Polynomials
using SpecialFunctions
using LsqFit
using JLD2
using .FrequencyEstimate
include(raw"FrequencyEstimation.jl")
include(raw"import.jl")


#=
ω  = 2*pi*f = 2*pi*T 
T = 1/f = 1*s/50 = 0.02s  
=#

# measurement bestimmt die spezifische Messung
#data = load("data.JLD2","data")
dataOfInterest = Int64[]
data = Measurement.data
numberOfData = size(data["single_data"]["t"])[1]
for i=1:numberOfData
    if data["single_data"]["wind_W"][i] == 0. && data["single_data"]["pv_W"][i] == 0. 
        push!(dataOfInterest,i)
    end
end


for j in dataOfInterest
for a=1:10
    println("Messreihe"*string(j))
end
measurement = j
#time  
t = data["single_data"]["t"][measurement]
# timefactor ist 1 oder 2 und wird in FrequencyEstimation bestimmt
timeFactor = FrequencyEstimate.checkTimeFactor(t)
t = t.*timeFactor
ω = 2*pi*FrequencyEstimate.checkFrequency(data,t,measurement) # frequency wird aus Messung bestimmmt 

if ω < 300  # statement if ω is not around 50
    ω = 2*ω 
    println("Eingriff")
end

# Park Transform (2x3 da dq0)
#K = [cos(θ) cos(θ-((2*pi)/3)) cos(θ+((2*pi)/3)); -sin(θ) -sin(θ-((2*pi)/3)) -sin(θ+((2*pi)/3)) ]
#K_inverse = [cos(θ) -sin(θ);cos(θ-((2*pi)/3)) -sin(θ-((2*pi)/3));cos(θ-((4*pi)/3)) -sin(θ-((4*pi)/3))]


#I_dq,U_dq = zeros(Float64, 2, size(t)[1])
I_dq = Float64[0.,0.]
U_dq = Float64[0.,0.]
I_dq = zeros(Float64, 2, size(t)[1])
I_dq0 = zeros(Float64, 3, size(t)[1])
U_dq = zeros(Float64, 2, size(t)[1])
U_dq0 = zeros(Float64, 3, size(t)[1])
U₀ = [data["single_data"]["U12"][measurement]';data["single_data"]["U23"][measurement]';data["single_data"]["U31"][measurement]'] # [V]
U₀ = U₀.*0.001 # [kV]
I₀ = [data["single_data"]["I1"][measurement]';data["single_data"]["I2"][measurement]';data["single_data"]["I3"][measurement]'].*(-1) # [A]
#I₀ = I₀.*0.001 # [kA]
display(ω)
for i=1:(size(t)[1])
    θ = ω*t[i]
    #θ = 314.159265*t[i]
    K = sqrt(2/3).*[cos(θ) cos(θ-((2*pi)/3)) cos(θ+((2*pi)/3)); -sin(θ) -sin(θ-((2*pi)/3)) -sin(θ+((2*pi)/3)) ]
    K0 = sqrt(2/3).*[cos(θ) cos(θ-((2*pi)/3)) cos(θ+((2*pi)/3)); -sin(θ) -sin(θ-((2*pi)/3)) -sin(θ+((2*pi)/3)); sqrt(0.5) sqrt(0.5) sqrt(0.5)]
    U = U₀[:,i]
    U_dq0[:,[i]] = K0*U
    U_dq[:,[i]] = K*U
    #push!(U_dq,K*U)
    I = I₀[:,i]
    I_dq0[:,[i]] = K0*I
    I_dq[:,[i]] = K*I
    #push!(I_dq,K*I)
end


#per unit  

I_pu = I_dq0.*(1/721.7) #I_base =S_base/(U_base*sqrt(3)) = 25MVA/20kV*sqrt(3) =0.7217kA
U_pu = U_dq0.*(1/(20)) #U_base = 20kV  

v_d = U_pu[1,:]
v_q = U_pu[2,:]
v_0 = U_pu[3,:]
i_d = I_pu[1,:]
i_q = I_pu[2,:]
i_0 = I_pu[3,:]

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
    push!(p_active,1.5*((v_d[i]*i_d[i]+v_q[i]*i_q[i]))) 
    push!(p_reactive,1.5*((v_q[i]*i_d[i]-v_d[i]*i_q[i])))     
end


# Smoothing curve
# window = 200
# Calculate mean of Voltage and Power in three stages
function curveSmoothing(t0, value_v,value_i ,value_s)
    window = 200

    v_mean = Float64[]
    v_mean_secondstage = Float64[]
    v_mean_thirdstage = Float64[]
    stages_v = [value_v, v_mean, v_mean_secondstage,v_mean_thirdstage]
    i_mean = Float64[]
    i_mean_secondstage = Float64[]
    i_mean_thirdstage = Float64[]
    stages_i = [value_i, i_mean, i_mean_secondstage,i_mean_thirdstage]
    s_mean=ComplexF64[]
    s_mean_secondstage=ComplexF64[]
    s_mean_thirdstage=ComplexF64[]
    stages_s= [value_s, s_mean, s_mean_secondstage,s_mean_thirdstage]
        for j = 1:3
        window = 100 +j*100
            for i=1:(size(t0)[1])
                if i > ((size(t0)[1])-window)
                    push!(stages_v[j+1],stages_v[j+1][((size(t0)[1])-window)])
                    push!(stages_i[j+1],stages_i[j+1][((size(t0)[1])-window)])
                    push!(stages_s[j+1],stages_s[j+1][((size(t0)[1])-window)])
                else
                    push!(stages_v[j+1],mean(stages_v[j][i:(i+window)]))
                    push!(stages_i[j+1],mean(stages_i[j][i:(i+window)]))
                    push!(stages_s[j+1],mean(stages_s[j][i:(i+window)]))    
                end
            end
        end
    stages = [stages_v[4],stages_i[4],stages_s[4]]
    return stages
end

means = curveSmoothing(t, v_abs,i_abs, power)
v_mean = real(means[1])
i_mean = real(means[2])
s_mean = means[3]
p_mean = real(s_mean)
q_mean = imag(s_mean)
data_fit = [v_mean real(s_mean) imag(s_mean)]
#=
#step fitting
#Two case
#Case1: Voltage does one step 
#Case2: Voltage does two steps 
using LsqFit
if abs(v_mean[1]-v_mean[20000])>0.001
    model = model(t, s) = s[1] .+ s[2]*sign.(t.+s[3])
    s0 = [1.,1.,-0.1]
else
    model = model(t, s) = s[1] .+ s[2]*sign.(t.+s[3]).+ s[4]*sign.(t.+s[5])
    s0 = [1.,1.,-0.1,-1.,-0.5]
end
stepfit_v = LsqFit.curve_fit(model, t[:], v_mean[:],s0)
#test_fit = curve_fit(nonlinear_fit,t[:],v_mean[:])
#Values for simulation
v_start = model(t[1],stepfit_v.param)
v_end = model(t[30000],stepfit_v.param)
v_fit = Float64[]


for i=1:size(t)[1]
    push!(v_fit,model(t[i],stepfit_v.param))
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
Δ_v = v_max-v_min
=#
steps = [0.1]
Δ_v = 0.02




##
#Optim
##
using PowerDynamics
using JuMP
using OrderedCollections: OrderedDict
using EasyFit
using Optim
using CSV


parameterNames=[ 
    "Nps",
    "Npt",
    "Nqs", 
    "Nqt",
    "Tp",
    "Tq"]

#dictionary = Dict{String,Any}("Parameter" => parameterNames)
#merge!(dictionary, Dict("test" => [1,1,1,1,1,1]))
function objective_function(parameter)
    
    P_meas = p_mean #yData
    Q_meas = q_mean
    function constraint(parameter)
        sum = 0.
        for i=1:6
            value = max(0.0000001,-parameter[i])
            sum = sum + value
        end
        return sum
    end
    function dotting(smallVector, longVector)
        dottedVector = smallVector
        length_old = length(smallVector)
        length_new = length(longVector)
        missingElements = length_new-length_old 
        dotsPerSpace = Int(round(missingElements/(length_old-1)))
        for i=1:((length_old-1))
            position = i+dotsPerSpace*(i-1)
             mean = (smallVector[position]+smallVector[position+1])/2
             for j=1:dotsPerSpace
                insert!(dottedVector,position+j,mean)
             end
        end
        if length(dottedVector) < 40960
            value = last(dottedVector)
            for i=1:(40960-length(dottedVector))
                push!(dottedVector,value)
            end
            return dottedVector 
        elseif length(dottedVector) > 40960
            for i=1:(length(dottedVector)-40960)
                pop!(dottedVector)
            end
            return dottedVector
        else
            return dottedVector
        end
    end 
    #smoothed p_meas 

    p_mean = movavg(P_meas,1000)
    p_mean = p_mean.x.*(1)
    q_mean = movavg(Q_meas,1000)
    q_mean = q_mean.x.*(1)
    offset = constraint(parameter)
    p = parameter.+offset    
    V₀ = 1. 
    V₁ = 1. - Δ_v
    P_Load = p_mean[1]
    Q_Load = q_mean[1]
    steptime = steps[1]
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U=V₀),
        "bus2" =>  ExponentialRecoveryLoad(P0=P_Load, Q0=Q_Load, Nps=p[1], Npt=p[2], Nqs=p[3], Nqt=p[4], Tp=p[5], Tq=p[6], V0=V₀),
    )

    branches = OrderedDict(
        "branch1" => StaticLine(from = "bus1",to = "bus2",Y=0.5+0.5*im)
    )

    powergrid = PowerGrid(buses,branches)
    parameter_set = [1.;1.;1.;1.;0.1;0.1]
    timespan= (0.0,4.)
    fault = NodeParameterChange(node = "bus1", value = V₁, tspan_fault = (steptime,4.),var=:U)  
    #solution = simulate(fault, powergrid, operationpoint, timespan)
    #λ=0
    function runSimulation()
        try 
            operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
            solution = simulate(fault, powergrid, operationpoint, timespan) 
            λ = 0.000001
        catch 
            λ = 1000.
        end
        return λ
    end
    
    λ = runSimulation()
    if λ < 1 
        operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
        solution = simulate(fault, powergrid, operationpoint, timespan)
    else
        p = [1.;1.;1.;1.;0.1;0.1]
        buses = OrderedDict(
        "bus1" => SlackAlgebraic(U=V₀),
        "bus2" =>  ExponentialRecoveryLoad(P0=P_Load, Q0=Q_Load, Nps=p[1], Npt=p[2], Nqs=p[3], Nqt=p[4], Tp=p[5], Tq=p[6], V0=V₀),
    )

    branches = OrderedDict(
        "branch1" => StaticLine(from = "bus1",to = "bus2",Y=0.5+0.5*im)
    )

    powergrid = PowerGrid(buses,branches)
    parameter_set = [1.;1.;1.;1.;0.1;0.1]
    operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
    timespan= (0.0,4.)
    fault = NodeParameterChange(node = "bus1", value = V₁, tspan_fault = (0.1,4.),var=:U)  
    solution = simulate(fault, powergrid, operationpoint, timespan)
    end
    #p_sim = solution(t,:,:p)[2,:]
    p_sim = solution(solution.dqsol.t,[ "bus2"],:p)[1,:]
    p_sim_new = dotting(p_sim,p_mean)
    p_sim_smooth = movavg(p_sim_new,1516)
    p_sim_smooth = p_sim_smooth.x
    #println(p_sim)
    return sum((p_sim_smooth-p_mean) .^ 2) + λ
end
initial_guess = [1.;1.;1.;1.;0.1;0.1]
result = optimize(objective_function, initial_guess, NelderMead())
optimized_params = Optim.minimizer(result)
for i=1:6
    if optimized_params[i]<0
        optimized_params[i] = optimized_params[i]*(-1)
    end
end

#parameter = runOpt(p,p_mean)
#parameterNames=[ Nps,Npt,Nqs, Nqt,Tp,Tq]
#merge!(dictionary,Dict( "First Evaluation" => optimized_params))
#filePath = "paramFitedAll.csv"
#CSV.write(filePath,dictionary)
using DataFrames
dictionary = DataFrame()
col_name = "Measurement"
dictionary[!, col_name] = [j]
for (idx, element) in enumerate(optimized_params)
    col_name = string(parameterNames[idx])
    dictionary[!, col_name] = [element]
end

filePath = "paramFited.csv"
CSV.write(filePath,dictionary)


end



# Abtastrate = 10kHz -> 10000 Punkte pro Sekunde -> 1000 Punkte pro Zehntel Skeunde -> 100 pro Hundert Sekunde
# Periodendauer bei 50Hz: T=1/f=0.02s -> 2*10^-2 -> 200 Punkte pro Peridode
# Von Nulldurchgang zu Nulldurchgang vergeht eine halbe Periode -> f= 1/T=1/(2*(t₁-t₀))
module FrequencyEstimate

function checkTimeFactor(t0)
    if (t0[12]-t0[11] < 0.00008)
        timeFactor = 2.
    else 
        timeFactor = 1.
    end
    return timeFactor
    end

function checkFrequency(d,t,number)
    #t0 = data["single_data"]["t"][number]'  #Zeivektor
    #t0 = t0.*checkTimeFactor(t0)
    #I0 = data["single_data"]["I1"][number]'  #Phase 1
    t0 = t'  #Zeivektor
    t0 = t0.*checkTimeFactor(t)
    I0 = d["single_data"]["I1"][number]'  #Phase 1
    t₀ = 0. #first zerocrossing 
    t₁ = 1. # second zerocrossing
    f = 0.  #frequency
    n = 0 # used for counting zero crossing
    for i=1:500
        if ((I0[i]*I0[i+1])<0.) && n<1 # Erster Nulldurchgang   (I[i]<2. && I[i]>-2.)
            t₀ = t0[i]
            n = 1
            #println("Erster Nulldruchgang :"*string(t₀))
            #println(string(i))
        elseif ((I0[i]*I0[i+1])<0.) && n<2 && (t0[i]-t₀>0.005) #Zweiter Nulldurchgang
            t₁ = t0[i]
            n = 2
            #println("Zweiter Nulldruchgang: "*string(t₁))
            #println(string(i))
        end
        if n==2
            f  = 1/(2*(t₁-t₀))
            break
        end
    end
    return f
end
end
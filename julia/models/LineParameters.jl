Z₀ = 16 # per unit Impedanz   [Ω] 

#Leitungsparameter 
#NA2XS2Y 1x120 RM/25 12/20 kV  , q = 120mm^2
R = 0.253; #Ω/km  R = 0.2;   0.253 #Ω/km
X = 0.119; #X = 0.35;  0.119
Z = R+1*im*X
Y =1/Z
Y₀=Z₀*Y
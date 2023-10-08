
using Polynomials 
using Plots
#= 
using LsqFit 
# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
model(x, p) = p[1]*exp.(-x.*p[2])

# some example data
# xdata: independent variables
# ydata: dependent variable
xdata = range(0, stop=10, length=20)
ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))
p0 = [0.5, 0.5]

fit = curve_fit(model, xdata, ydata, p0)
=#
xs = range(0, 10, length = 10)
ys = @.exp(-xs)
f = Polynomials.fit(xs, ys) # degree = length(xs) - 1
f2 = Polynomials.fit(xs, ys, 2) # degree = 2
scatter(xs, ys, markerstrokewidth = 0, label = "Data")
plot!(f, extrema(xs)..., label = "Fit")
plot!(f2, extrema(xs)..., label = "Quadratic Fit")
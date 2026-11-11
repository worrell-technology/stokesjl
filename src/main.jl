include("stokesjl.jl")

using .stokesjl
using Plots

x = collect(range(0, 2, length=41))
u = [0.5 <= xi <= 1.0 ? 2.0 : 1.0 for xi in x]
p = plot(x, u)
model = wave1d(u=u, x=x)

convect1d(model, 0.0, .625, .025)
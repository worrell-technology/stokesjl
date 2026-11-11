include("stokesjl.jl")

using .stokesjl

# 1D Convection
x = collect(range(0, 2, length=41))
u = [0.5 <= xi <= 1.0 ? 2.0 : 1.0 for xi in x]
model = wave1d(u=u, x=x)
convect1d(model, 0.0, 1.0, .025)

# Nonlinear Convection
u = [0.5 <= xi <= 1.0 ? 2.0 : 1.0 for xi in x]
model.u = u
convect1d_nl(model, 0.0, 1.0, .025)

# 1D Diffusion
u = [0.5 <= xi <= 1.0 ? 2.0 : 1.0 for xi in x]
model.u = u
σ = .2
Δt = σ * Δx(model)^2 / model.ν 
diffuse1d(model, 0.0, 1.0, Δt)
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
σ = .4
Δt = σ * Δx(model)^2 / model.ν 
diffuse1d(model, 0.0, 1.0, Δt)

# Burgers
f = (x, nu, t) -> -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)^2/(4*nu*(t + 1)))/(4*nu*(t + 1)) - 
                  (-8*t + 2*x - 4*pi)*exp(-(-4*t + x - 2*pi)^2/(4*nu*(t + 1)))/
                  (4*nu*(t + 1)))/(exp(-(-4*t + x - 2*pi)^2/(4*nu*(t + 1))) + 
                  exp(-(-4*t + x)^2/(4*nu*(t + 1)))) + 4
x = collect(range(0, 2 * π, length=101))
model.x = x
model.ν = 0.07
u = [f(xi, model.ν, 0) for xi in x]
model.u = u
dx = Δx(model)
Δt = 0.5 * dx * model.ν
# burgers(model, f, 0.0, 1.0, Δt)
@show minimum(model.u)
@show maximum(model.u)
@show model.ν
@show Δx(model)
@show Δt

# 2D Convection
m = Model(coords=init_coords((0.0, 2.0), 81, (0.0, 2.0), 81))
u = ones(length(m.coords.x), length(m.coords.y))
u[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.u = u
convect2d(m, 0.0, 1.0)

# 2D Nonlinear Convection
#  re-initialize state
v = ones(length(m.coords.x), length(m.coords.y))
v[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.v = v
u = ones(length(m.coords.x), length(m.coords.y))
u[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.u = u
convect2d_nl(m, 0.0, 1.0)

# 2D Diffusion
#  re-initialize state
v = ones(length(m.coords.x), length(m.coords.y))
v[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.v = v
u = ones(length(m.coords.x), length(m.coords.y))
u[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.u = u
diffuse2d(m, 0.0, 1.0)

# 2D Burgers
#  re-initialize state
v = ones(length(m.coords.x), length(m.coords.y))
v[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.v = v
u = ones(length(m.coords.x), length(m.coords.y))
u[round(Int, 0.5/m.coords.dx)+1 : round(Int, 1/m.coords.dx)+1,
  round(Int, 0.5/m.coords.dy)+1 : round(Int, 1/m.coords.dy)+1] .= 2
m.u = u
m.Δt = dt(m, 0.09)
burgers2d(m, 0.0, 1.0)
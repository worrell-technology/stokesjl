module StokesTypes

export Model, wave1d, init_coords
export dt, cfl, n, Δx

Base.@kwdef mutable struct wave1d
    u::Vector{Float64}
    x::Vector{Float64} = collect(range(0, 1, length=length(u)))
    c::Float64 = 1.0
    ν::Float64 = .3
end

n(model::wave1d) = length(model.u)
Δx(model::wave1d) = (model.x[end] - model.x[1]) / (n(model) - 1)

Base.@kwdef struct CoordParams
    mn::Float64 = 0.0
    mx::Float64 = 1.0
    n::Int = 81
end

Base.@kwdef struct Coords2d
    _x_pars::CoordParams
    _y_pars::CoordParams
    x::Vector{Float64} = collect(range(_x_pars.mn, _x_pars.mx, length=_x_pars.n))
    y::Vector{Float64} = collect(range(_y_pars.mn, _y_pars.mx, length=_y_pars.n))
    dx::Float64 = (x[end] - x[1]) / (_x_pars.n - 1)
    dy::Float64 = (y[end] - y[1]) / (_y_pars.n - 1)
end

init_coords(x::Tuple{Float64, Float64}, nx::Int, y::Tuple{Float64, Float64}, ny::Int) = 
    Coords2d(_x_pars=CoordParams(x[1], x[2], nx),
        _y_pars=CoordParams(y[1], y[2], ny))
init_coords(x_pars::CoordParams, y_pars::CoordParams) =
    Coords2d(_x_pars = x_pars, _y_pars = y_pars)
init_coords(x::Vector{Float64}, y::Vector{Float64}) = 
    Coords2d(_x_pars=CoordParams(x[1], x[end], length(x)),
        _y_pars=CoordParams(y[1], y[end], length(y)),
        x=x, y=y)

# Main Model Type
Base.@kwdef mutable struct Model
    coords::Coords2d
    u::Matrix{Float64} = ones(
        length(coords.y),
        length(coords.x)
    )
    cx::Float64 = 1.0
    cy::Float64 = 1.0
    Δt::Float64 = _dt(0.25, cx, cy, coords.dx, coords.dy)
end

cfl(m::Model) = m.Δt * (m.cx/m.coords.dx + m.cy/m.coords.dy) 
dt(m::Model, c::Float64) = c / (m.cx/m.coords.dx + m.cy/m.coords.dy) 
_dt(σ, cx, cy, dx, dy) = σ / (cx/dx + cy/dy) 

end # module StokesTypes
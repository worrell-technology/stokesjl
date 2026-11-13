module stokesjl
include("StokesTypes.jl")

export Model, wave1d, init_coords
export dt, cfl, n, Δx

export convect1d, convect1d_nl, diffuse1d, burgers
export convect2d

using Plots
using .StokesTypes

function convect1d_nl(model::wave1d, start_time, end_time, Δt)
    sz = n(model)
    dx = Δx(model)
    t = start_time
    un = ones(sz)

    anim = @animate while t < end_time
        un = copy(model.u)

        plot(model.x, model.u,
            xlim=(0, 2),
            ylim=(0, 2),
            xlabel="x",
            ylabel="u",
            linewidth=2,
            title="Time = $(round(t, digits=3))",
            legend=false
        )

        for i ∈ 2:sz
            model.u[i] = un[i] - un[i] * (un[i] - un[i-1]) * Δt / dx
        end
        t += Δt
    end
    gif(anim, "results/convect1d_nl.gif", fps=5)
end

function convect1d(model::wave1d, start_time, end_time, Δt)
    sz = n(model)
    dx = Δx(model)
    t = start_time
    un = ones(sz)

    anim = @animate while t <= end_time
        un = copy(model.u)

        plot(model.x, model.u,
            xlim=(0, 2),
            ylim=(0, 2),
            xlabel="x",
            ylabel="u",
            linewidth=2,
            title="Time = $(round(t, digits=3))",
            legend=false
        )

        for i ∈ 2:sz
            model.u[i] = un[i] - model.c * (un[i] - un[i-1]) * Δt / dx
        end

        t += Δt
    end
    gif(anim, "results/convect1d.gif", fps=5)
end

function diffuse1d(model::wave1d, start_time, end_time, Δt)
    sz = n(model)
    dx = Δx(model)
    t = start_time
    un = ones(sz)

    anim = @animate while t <= end_time
        un = copy(model.u)
        
        plot(model.x, model.u,
            xlim=(0, 2),
            ylim=(0, 2),
            xlabel="x",
            ylabel="u",
            linewidth=2,
            title="Time = $(round(t, digits=3))",
            legend=false
        )

        for i ∈ 2:sz-1
            model.u[i] = un[i] + model.ν * 
                (un[i+1] - 2*un[i] + un[i-1]) * Δt / dx^2
        end
        t += Δt
    end every 10
    gif(anim, "results/diffuse1d.gif", fps=5) 
end

function burgers(model::wave1d, analytical, start_time, end_time, Δt)
    sz = n(model)
    dx = Δx(model)
    t = start_time
    un = ones(sz)

    anim = @animate while t <= end_time
        un = copy(model.u)
        u_analytical = [analytical(xi, model.ν, t) for xi in model.x]

        p = plot(model.x, model.u,
                xlim=(0, 2*π),
                ylim=(0, 10),
                xlabel="x",
                ylabel="u",
                linewidth=2,
                title="Time = $(round(t, digits=3))",
                label="Numerical"
        )

        plot!(p, model.x, u_analytical,
            linewidth=2, linestyle=:dash, label="Analytical"
        )

        for i ∈ 2:sz-1
            model.u[i] = un[i] - un[i] * (un[i]-un[i-1]) * Δt / dx +
                model.ν * (un[i+1] - 2 * un[i] + un[i-1]) * Δt / dx^2
        end
        model.u[1] = un[1] - un[1] * (un[1]-un[end-1]) * Δt / dx +
            model.ν * (un[2] - 2 * un[1] + un[end-1]) * Δt / dx^2
        model.u[end] = model.u[1]
        t += Δt
    end
    gif(anim, "results/burgers.gif", fps=5)
end

function convect2d(m::Model, start_time::Float64, end_time::Float64)
    t = start_time
    anim = @animate while t <= end_time
        un = copy(m.u)

        surface(
            m.coords.x,
            m.coords.y,
            m.u,
            color = :viridis,
            xlabel = "X",
            ylabel = "Y",
            zlabel = "U",
            title="Time = $(round(t, digits=3))",
            camera = (45, 45),
            size = (1100, 700),
            dpi=100
        )

        rows, cols = size(un)
        for j in range(2, cols)
            for i in range(2, rows)
                m.u[i, j] = (un[i, j] - (m.cy * (un[i, j] - un[i, j - 1]) * m.Δt / m.coords.dy) - 
                                        (m.cx * (un[i, j] - un[i-1, j]) * m.Δt / m.coords.dx))
                m.u[1,:] .= 1
                m.u[end, :] .= 1
                m.u[:, 1] .= 1
                m.u[:, end] .= 1
            end
        end
        t += m.Δt
    end
    gif(anim, "results/convect2d.gif", fps=5)
end

end # module stokesjl
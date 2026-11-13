module stokesjl
include("StokesTypes.jl")

export wave1d, convect1d, convect1d_nl, diffuse1d, burgers, Δx

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
                label="Numberical"
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
end # module stokesjl
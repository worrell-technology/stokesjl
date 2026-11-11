module stokesjl

export wave1d, convect1d, convect1d_nl, diffuse1d, Δx

using Plots

Base.@kwdef mutable struct wave1d
    u::Vector{Float64}
    x::Vector{Float64} = collect(range(0, 1, length=length(u)))
    c::Float64 = 1.0
    ν::Float64 = .3
end

n(model::wave1d) = length(model.u)
Δx(model::wave1d) = (model.x[end] - model.x[1]) / (n(model) - 1)

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
    end
    gif(anim, "results/diffuse1d.gif", fps=5) 
end

end # module stokesjl

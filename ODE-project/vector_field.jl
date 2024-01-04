using AbstractPlotting
using CairoMakie


array = [(t, y₁, y₂) -> (y₁ * (1.0 + y₂)),
         (t, y₁, y₂) -> (y₂ * (y₁ + 1.0))]

function phase_portrait(𝐹::Array{Function, 1}, P::Number, title::String, xlab::String, ylab::String)
    f(x, y) = Point2f0(𝐹[1](0, x, y), 𝐹[2](0, x, y))
    xs = ys = LinRange(0, P, 3 * ceil(P))
    scene = streamplot(f, xs, ys, arrow_size=0.5, linewidth=1, colormap=:magma)  # arrow_size, colormap
    axis = scene[Axis]
    axis.names.axisnames = (xlab, ylab)
    axis.names.title = title
    axis.names.font = ("Latin Modern Math", "Latin Modern Math")
    axis.ticks.font = ("Latin Modern Math", "Latin Modern Math")
    scene.resolution = (700, 450)
    return scene
end

a = phase_portrait(array, 10, "Phase portrait for mutualism", "Population 1", "Population 2")

save("phase_plot.svg", a)  # pdf, ...

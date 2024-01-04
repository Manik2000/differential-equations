using Plots


"""
    runge_kutta(f::Function, x₀::Real, y₀::Real, x::Real, h::Real=10^(-3))::Array{Float64,1}
Return an array with values of function which is solution to the differential equation
given by function f with intial condition
"""
function runge_kutta(f::Function, x₀::Real, y₀::Real, x::Real, h::Real=10^(-3))::Array{Float64,1}
    n = round(Int64, x / h)
    array = zeros(n)
    for i in 1:n
        k₁ = f(x₀, y₀)
        k₂ = f(x₀ + 1/2*h, y₀ + 1/2*h*k₁)
        k₃ = f(x₀ + 1/2*h, y₀ + 1/2*h*k₂)
        k₄ = f(x₀ + h, y₀ + h*k₃)
        y₀ = y₀ + h / 6 * (k₁ + 2*k₂ + 2*k₃ + k₄)
        x₀ = x₀ + h
        array[i] = y₀
    end
    return array
end

# logistic growth for chosen constants
f(t, y) = 0.6y*(1 - y / 10)
array = runge_kutta(f, 0, 1, 20)
𝐫 = plot(LinRange(0, 20, length(array)), array)
display(𝐫)

# https://en.smath.com/wiki/GetFile.aspx?File=Examples/RK4SystemEquations.pdf
# examples below are taken from this url

array  = [(t, y₁, y₂) -> sin(t) + cos(y₁) + sin(y₂),
        (t, y₁, y₂) -> cos(t) + sin(y₂)]
X = 0
Y = [-1., 1.]



"""
    evaluate(array::Array{Function, 1}, X::Real, Y::Array{Real, 1})::Array{Real, 1}
Return an array with values of function f(t, y₁, y₂) defined in array for given point X
"""
function evaluate(array::Array{Function, 1}, X::Real, Y::Array{T, 1})::Array{T, 1} where T <: Real
    return [i(X, Y[1], Y[2]) for i in array]
end


"""
    runge_kutta(F::Array{Function, 1}, x::Real, X::Real, Y::Array{T, 1}, h::Real=1e-3)::Array{T, 2} where T <: Real
Return an array with values of functions which are solutions to the system of differential equations defined by the array.
X and Y are initial conditions, x is the maximum point.
"""
function runge_kutta(F::Array{Function, 1}, x::Real, X::Real, Y::Array{T, 1}, h::Real=1e-3)::Array{T, 2} where T <: Real
    n = round(Int64, x / h)
    array = zeros((n, 2))
    for i in 1:n
        k₁ = evaluate(F, X, Y)
        k₂ = evaluate(F, X + 1/2 * h, Y .+ 1/2 * h * k₁)
        k₃ = evaluate(F, X + 1/2 * h, Y .+ 1/2 * h * k₂)
        k₄ = evaluate(F, X + h, Y .+ h * k₃)
        Y = @. Y + h / 6 * (k₁ + 2*k₂ + 2*k₃ + k₄)
        X = X + h
        array[i, :] = Y
    end
    return array
end


a = runge_kutta(array, 20, X, Y)

𝑝 = plot(LinRange(0, 5, length(a[:, 1])), a[:, 1])
plot!(LinRange(0, 5, length(a[:, 2])), a[:, 2])
display(𝑝)
# everything works fine


# harmonic oscillator ẍ + x = 0
# can be rewritten as
# ẋ = y, ẏ = -x
# note that y is velocity ⃗v


array2 = [(t, y₁, y₂) -> y₂, (t, y₁, y₂) -> -y₁]
Y2 = [0., 1.]
b = runge_kutta(array2, 30, 0, Y2)

𝑞 = plot(LinRange(0, 5, length(b[:, 1])), b[:, 1])
display(𝑞)
# plot!(LinRange(0, 5, length(b[:, 2])), b[:, 2])
𝑤 = plot(b[:, 1], b[:, 2])  # a kind of phase plane plot
display(𝑤)

# vector field for harmonic oscillator
xs = LinRange(0, 5, length(b[1, :]))

xs = xs[1:150:end]
xs |>println
y₁s = b[:, 1][1:150:end]
y₂s = b[:, 2][1:150:end]

plot(quiver(y₁s, y₂s, quiver=(0.2 * array2[1].(xs, y₁s, y₂s), 0.2 * array2[2].(xs, y₁s, y₂s)), color=:red))

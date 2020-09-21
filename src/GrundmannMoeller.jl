module GrundmannMoeller

using LinearAlgebra
using StaticArrays

export TnScheme
struct TnScheme{N,T}
    name::String
    dim::Int
    weights::Vector{T}
    points::Vector{SVector{N,T}}
    degree::Int
    source::Dict
end

export integrate
"""
    integrate(fun, scheme, vertices::AbstractVector)
    integrate(fun, scheme, vertices::AbstractMatrix)

# Arguments
- `fun`: integrand, should accept an `SVector` as argument
- `scheme`: quadrature scheme
- `vertices`: vertices of the simplex

The vertices need to be passed either as a vector-of-vectors or as a
matrix. In the first case, there need to be `D+1` points with `D`
coordinates each. In the second case, the matrix needs to have size
`D`×`D+1`.
"""
@inbounds function integrate(fun, scheme::TnScheme{N,T},
                             vertices::SMatrix{D,N,U}) where {N,T,D,U}
    @assert N > 0
    @assert N >= D + 1

    ws = scheme.weights
    ps = scheme.points
    @assert length(ws) == length(ps)

    p1 = ps[1]
    x1 = (vertices * p1)::SVector{D}
    X = typeof(x1)
    R = typeof(ws[1] * fun(x1))

    s = zero(R)
    @simd for i in 1:length(ws)
        w = ws[i]
        p = ps[i]
        x = vertices * p
        s += w * fun(x)
    end

    # If `U` is an integer type, then Linearalgebra.det converts to
    # floating-point values; we might want a different type
    vol = R(calc_vol(vertices)) / factorial(D)
    return vol * s
end
@inbounds function integrate(fun, scheme::TnScheme{N,T},
                             vertices::SVector{N,SVector{D,U}}) where {N,T,D,U}
    return integrate(fun, scheme,
                     SMatrix{D,N,U}(vertices[n][d] for d in 1:D, n in 1:N))
end
function integrate(kernel, scheme::TnScheme{N},
                   vertices::AbstractVector) where {N}
    @assert length(vertices) == N
    @assert N > 0
    D = length(vertices[1])
    @assert N >= D + 1
    vertices′ = SVector{N}(map(SVector{D}, vertices))
    return integrate(kernel, scheme, vertices′)
end
function integrate(kernel, scheme::TnScheme{N},
                   vertices::AbstractMatrix) where {N}
    @assert size(vertices, 1) == N
    @assert N > 0
    D = length(vertices[1])
    @assert N >= D + 1
    vertices′ = SMatrix{N,D}(vertices)'
    return integrate(kernel, scheme, vertices)
end

@inbounds function calc_vol(vertices::SMatrix{D,N,T}) where {D,N,T}
    @assert N == D + 1
    X = SMatrix{D,D,T}(vertices[i, j + 1] - vertices[i, 1]
                       for i in 1:D, j in 1:D)
    vol = det(X)
    return vol
end

################################################################################

const source = Dict(:authors => ["A. Grundmann", "H. M. Möller"],
                    :title => "Invariant integration formulas for the n-simplex by combinatorial methods",
                    :journal => "SIAM J. Numer. Anal.", :volume => "15",
                    :year => "1978", :pages => "282-290",
                    :url => "https://doi.org/10.1137/0715019")

export grundmann_moeller
"""
    grundmann_moeller(::Type{T}, ::Val{D}, degree)

# Arguments
- `T`: desired floating point type
- `D`: dimension
- `degree`: desired polynomial degree of accuracy (must be odd)
"""
function grundmann_moeller(::Type{T}, ::Val{D}, degree::Int) where {T,D}
    D::Int
    @assert degree ≥ 0
    @assert isodd(degree)

    N = D + 1

    order = (degree - 1) ÷ 2

    exponents = get_all_exponents(Val(N), order)

    weights = Vector{T}()
    points = Vector{SVector{N,T}}()
    for i in 0:order
        w = T((-1)^i) * big(degree + D - 2 * i)^degree / (big(2)^(2 * order) *
             factorial(big(i)) *
             factorial(big(degree + D - i)))
        for part in exponents[order - i + 1]
            push!(weights, w)
            push!(points, map(p -> T(2 * p + 1) / (degree + D - 2 * i), part))
        end
    end
    weights /= sum(weights)

    name = "Grundmann-Möller(dim=$D, $degree)"
    return TnScheme{N,T}(name, D, weights, points, degree, source)
end

"""
Get all exponent combinations of dimension `dim` and maximum degree
`max_degree`. This method is actually meant for evaluating all
polynomials with these exponents. This problem is similar to the
weak_compositions, e.g.,
<https://stackoverflow.com/a/36748940/353337>. The solution here,
however, only ever adds 1, making it better suited for a possible
extension with actual polynomial evaluation.
"""
function get_all_exponents(::Val{D}, max_degree::Int) where {D}
    # Initialization, level 0
    exponents = [SVector{D,Int}(0 for d in 1:D)]

    all_exponents = Vector{SVector{D,Int}}[]
    push!(all_exponents, exponents)
    for _ in 1:max_degree
        exponents = augment(exponents)
        push!(all_exponents, exponents)
    end

    return all_exponents::Vector{Vector{SVector{D,Int}}}
end

"""
This function takes the values and exponents of a given monomial
level, e.g., [(1,0,0), (0,1,0), (0,0,1)], and augments them by one
level, i.e., [(2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)].
The method works for all dimensions and is based on the observation
that the augmentation can happen by
    1. adding 1 to all exponents from the previous level (i.e.,
       multiplication with x[0]),
    2. adding 1 to the second exponent of all exponent tuples where
       ex[0]==0 (i.e., multiplication with x[1]),
    3. adding 1 to the third exponent of all exponent tuples where
       ex[0]==0, ex[1]=0 (i.e., multiplication with x[2]),
etc. The function call is recursive.
"""
function augment(exponents::Vector{SVector{D,Int}}) where {D}
    if length(exponents) == 0 || D == 0
        return Vector{SVector{D,Int}}()
    end

    idx_leading_zero = [k for k in 1:length(exponents) if exponents[k][1] == 0]
    exponents_with_leading_zero = [popfirst(exponents[k])
                                   for k in idx_leading_zero]

    out1 = augment(exponents_with_leading_zero)
    # increment leading exponent by 1
    out = [SVector(e[1] + 1, popfirst(e)...) for e in exponents]
    append!(out, [SVector(0, e...) for e in out1])

    return out
end

end

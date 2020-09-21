using GrundmannMoeller
using StaticArrays
using Test

const GM = GrundmannMoeller
const SV{N} = SVector{N,Int}

@testset "augment" begin
    @test GM.augment(SV{0}[]) == []
    @test GM.augment(SV{1}[]) == []

    @test GM.augment([SV{1}(0)]) == [SV{1}(1)]
    @test GM.augment([SV{1}(1)]) == [SV{1}(2)]
    @test GM.augment([SV{1}(2)]) == [SV{1}(3)]

    @test GM.augment([SV{2}(0, 0)]) == [SV{2}(1, 0), SV{2}(0, 1)]
    @test GM.augment([SV{2}(1, 0), SV{2}(0, 1)]) ==
          [SV{2}(2, 0), SV{2}(1, 1), SV{2}(0, 2)]
    @test GM.augment([SV{2}(2, 0), SV{2}(1, 1), SV{2}(0, 2)]) ==
          [SV{2}(3, 0), SV{2}(2, 1), SV{2}(1, 2), SV{2}(0, 3)]

    @test GM.augment([SV{3}(0, 0, 0)]) ==
          [SV{3}(1, 0, 0), SV{3}(0, 1, 0), SV{3}(0, 0, 1)]
    @test GM.augment([SV{3}(1, 0, 0), SV{3}(0, 1, 0), SV{3}(0, 0, 1)]) ==
          [SV{3}(2, 0, 0), SV{3}(1, 1, 0), SV{3}(1, 0, 1), SV{3}(0, 2, 0),
           SV{3}(0, 1, 1), SV{3}(0, 0, 2)]
end

@testset "get_all_exponents basic" begin
    @test GM.get_all_exponents(Val(0), 0) == [[SV{0}()]]

    @test GM.get_all_exponents(Val(1), 0) == [[SV{1}(0)]]
    @test GM.get_all_exponents(Val(1), 1) == [[SV{1}(0)], [SV{1}(1)]]
    @test GM.get_all_exponents(Val(1), 2) ==
          [[SV{1}(0)], [SV{1}(1)], [SV{1}(2)]]
    @test GM.get_all_exponents(Val(1), 3) ==
          [[SV{1}(0)], [SV{1}(1)], [SV{1}(2)], [SV{1}(3)]]

    @test GM.get_all_exponents(Val(2), 0) == [[SV{2}(0, 0)]]
    @test GM.get_all_exponents(Val(2), 1) ==
          [[SV{2}(0, 0)], [SV{2}(1, 0), SV{2}(0, 1)]]
    @test GM.get_all_exponents(Val(2), 2) ==
          [[SV{2}(0, 0)], [SV{2}(1, 0), SV{2}(0, 1)],
           [SV{2}(2, 0), SV{2}(1, 1), SV{2}(0, 2)]]
    @test GM.get_all_exponents(Val(2), 3) ==
          [[SV{2}(0, 0)], [SV{2}(1, 0), SV{2}(0, 1)],
           [SV{2}(2, 0), SV{2}(1, 1), SV{2}(0, 2)],
           [SV{2}(3, 0), SV{2}(2, 1), SV{2}(1, 2), SV{2}(0, 3)]]

    @test GM.get_all_exponents(Val(3), 0) == [[SV{3}(0, 0, 0)]]
    @test GM.get_all_exponents(Val(3), 1) ==
          [[SV{3}(0, 0, 0)], [SV{3}(1, 0, 0), SV{3}(0, 1, 0), SV{3}(0, 0, 1)]]
    @test GM.get_all_exponents(Val(3), 2) ==
          [[SV{3}(0, 0, 0)], [SV{3}(1, 0, 0), SV{3}(0, 1, 0), SV{3}(0, 0, 1)],
           [SV{3}(2, 0, 0), SV{3}(1, 1, 0), SV{3}(1, 0, 1), SV{3}(0, 2, 0),
            SV{3}(0, 1, 1), SV{3}(0, 0, 2)]]
end

@testset "get_all_exponents D=$D" for D in 1:5
    exponents = GM.get_all_exponents(Val(D), 0)
    @test exponents == [[SV{D}(0 for d in 1:D)]]
    for n in 1:10
        old_exponents = exponents
        exponents = GM.get_all_exponents(Val(D), n)
        @test length(exponents) == n + 1
        @test exponents[1:(end - 1)] == old_exponents
        @test allunique(exponents[end])
        @test all(x -> sum(x) == n, exponents[end])
        @test all(x -> all(≥(0), x), exponents[end])
        for x in exponents[end - 1]
            for d in 1:D
                y = setindex(x, x[d] + 1, d)
                @test y ∈ exponents[end]
            end
        end
    end
end

@testset "grundmann_moeller basic" begin
    T = Rational{BigInt}

    D = 1
    degree = 1
    scheme = grundmann_moeller(T, Val(D), degree)
    @test scheme.dim == D
    @test scheme.degree == degree
    @test scheme.weights == [1]
    @test scheme.points == [[1//2, 1//2]]

    D = 1
    degree = 3
    scheme = grundmann_moeller(T, Val(D), degree)
    @test scheme.dim == D
    @test scheme.degree == degree
    @test scheme.weights == [2//3, 2//3, -1//3]
    @test scheme.points == [[3//4, 1//4], [1//4, 3//4], [1//2, 1//2]]

    D = 1
    degree = 5
    scheme = grundmann_moeller(T, Val(D), degree)
    @test scheme.dim == D
    @test scheme.degree == degree
    @test scheme.weights == [27//40, 27//40, 27//40, -8//15, -8//15, 1//24]
    @test scheme.points ==
          [[5//6, 1//6], [1//2, 1//2], [1//6, 5//6], [3//4, 1//4], [1//4, 3//4],
           [1//2, 1//2]]

    D = 2
    degree = 1
    scheme = grundmann_moeller(T, Val(D), degree)
    @test scheme.dim == D
    @test scheme.degree == degree
    @test scheme.weights == [1//1]
    @test scheme.points == [[1//3, 1//3, 1//3]]

    D = 2
    degree = 3
    scheme = grundmann_moeller(T, Val(D), degree)
    @test scheme.dim == D
    @test scheme.degree == degree
    @test scheme.weights == [25//48, 25//48, 25//48, -9//16]
    @test scheme.points ==
          [[3//5, 1//5, 1//5], [1//5, 3//5, 1//5], [1//5, 1//5, 3//5],
           [1//3, 1//3, 1//3]]

    D = 3
    degree = 1
    scheme = grundmann_moeller(T, Val(D), degree)
    @test scheme.dim == D
    @test scheme.degree == degree
    @test scheme.weights == [1//1]
    @test scheme.points == [[1//4, 1//4, 1//4, 1//4]]
end

function check_scheme(scheme::TnScheme{2,T}) where {T}
    for p1 in 0:(scheme.degree)
        f(x) = x[1]^p1
        fi = one(T) / (p1 + 1)
        fq = integrate(f, scheme, [(0,), (1,)])
        @test fq == fi
    end
end

function check_scheme(scheme::TnScheme{N,T}) where {N,T}
    D = N - 1
    degree = scheme.degree
    for ps0 in
        CartesianIndex(ntuple(d -> 0, D)):CartesianIndex(ntuple(d -> degree, D))
        ps = SVector(ps0.I)
        if sum(ps) <= degree
            f(xs) = prod(xs[d]^ps[d] for d in 1:D)
            fi = one(T) * prod(factorial(big(p)) for p in ps) /
                 factorial(big(sum(p + 1 for p in ps)))
            vertices = SVector{D,Int}[]
            push!(vertices, ntuple(d -> 0, D))
            for d in 1:D
                push!(vertices, ntuple(==(d), D))
            end
            fq = integrate(f, scheme, vertices)
            @test fq == fi
        end
    end
end

@testset "grundmann_moeller D=$D" for D in 1:5
    T = Rational{BigInt}

    for degree in 1:2:11
        scheme = grundmann_moeller(T, Val(D), degree)
        @test scheme.dim == D
        @test scheme.degree == degree
        check_scheme(scheme)
    end
end

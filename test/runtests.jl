#include("../src/ComplexRationals.jl")
using Test
using ComplexRationals

@testset "ComplexRational Basic Tests" begin
    z1 = ComplexRational(3, 2, 10)
    @test z1 isa ComplexRationals.ComplexRational
    @test z1.a == 3
    @test z1.b == 2
    @test z1.c == 10

    z2 = ComplexRational(3, 2, 1234)
    @test z2 isa Complex{Float64}
    @test isapprox(z2, complex(3 / 1234, 2 / 1234))

    # Normalization
    a = ComplexRational(2, 4, 6)
    @test a == ComplexRational(1, 2, 3)

    # Negative denominator handling
    b = ComplexRational(3, 6, -9)
    @test b == ComplexRational(-1, -2, 3)

    # Multiplying two ComplexRationals
    c = ComplexRational(1, 2, 3)
    d = ComplexRational(3, 4, 5)
    prod = c * d
    expected_prod = ComplexRational(-5, 10, 15) # As you worked out
    @test prod == expected_prod

    # Multiplying with an integer
    e = c * 3
    expected_e = ComplexRational(3, 6, 3)
    @test e == expected_e

    # crationalize on Rational
    r = 3 // 4
    x = crationalize(r, tol=1e-6)
    @test x == ComplexRational(3, 0, 4)

    # crationalize on Complex
    cnum = -1.0 + 0.5im
    y = crationalize(cnum, tol=1e-12)
    @test y == ComplexRational(-2, 1, 2)

    # Test is_negative
    neg1 = ComplexRational(-2, 3, 5)
    @test is_negative(neg1) == true

    neg2 = ComplexRational(0, -3, 5)
    @test is_negative(neg2) == true

    neg3 = ComplexRational(0, 0, 1)
    @test is_negative(neg3) == false
end

@testset "Addition Tests" begin
    x = ComplexRational(1, 2, 3)
    y = ComplexRational(2, -1, 4)

    s = x + y
    # (1/3 + 2i/3) + (2/4 - i/4) = (1/3 + 1/2) + i(2/3 - 1/4)
    # Real: (1/3)*(4/4) + (2/4)*(3/3) = (4/12 + 6/12) = 10/12 = 5/6
    # Imag: (2/3)*(4/4) + (-1/4)*(3/3) = (8/12 - 3/12) = 5/12
    @test s == ComplexRational(10, 5, 12)
end

@testset "Unary and Binary Subtraction Tests" begin
    x = ComplexRational(1, 2, 3)
    y = ComplexRational(2, -1, 4)

    diff = x - y
    # (1/3 + 2i/3) - (2/4 - i/4) = (1/3 - 1/2) + i(2/3 + 1/4)
    # Real: (1/3)*(4/4) - (2/4)*(3/3) = (4/12 - 6/12) = -2/12 = -1/6
    # Imag: (2/3)*(4/4) + (1/4)*(3/3) = (8/12 + 3/12) = 11/12
    @test diff == ComplexRational(-2, 11, 12)

    negx = -x
    @test negx == ComplexRational(-1, -2, 3)
end

@testset "Division Tests" begin
    x = ComplexRational(1, 2, 3)
    y = ComplexRational(2, -1, 4)

    q = x / y
    # messy formula, but just check that (q * y ≈ x)
    @test isapprox(complex(q) * complex(y), complex(x), atol=1e-10)
end

@testset "Identity and Inverse Tests" begin
    z = ComplexRational(3, 4, 5)

    @test one(z) == ComplexRational(1, 0, 1)
    @test zero(z) == ComplexRational(0, 0, 1)

    invz = inv(z)
    # inv((3+4i)/5) = (3 - 4i)/(3^2 + 4^2) = (3 - 4i)/25
    @test invz == ComplexRational(3, -4, 5)

    @test isapprox(complex(z) * complex(invz), 1 + 0im, atol=1e-10)
end

@testset "Exponentiation Tests" begin
    z = ComplexRational(1, 2, 3)

    z2 = z^2
    @test z2 == ComplexRational(-3, 4, 9)

    z3 = z^3
    @test z3 == z * z * z

    @test z^0 == ComplexRational(1, 0, 1)

    println(z)
    @test z^(-1) == inv(z)
end
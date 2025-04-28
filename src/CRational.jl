module CRational

using Preferences, LaTeXStrings


export set_default_digits, get_default_digits, ComplexRational, crationalize, custom_sort_key, sort_mixed, complexrational2str, is_negative

# === Default Coefficient Preferences ===
const DEFAULT_COEFF_PREFS = Dict(
    :DEFAULT_DIGITS => 2,
    :FLOAT_DIGITS => 2,
    :EXP_DIGITS => 2)


function get_default(name::Symbol)
    return @load_preference(String(name), DEFAULT_COEFF_PREFS[name])
end
function set_default(name::Symbol, value)
    @set_preferences!(String(name) => value)
end

"""
    set_default_digits(d::Int)

Sets the default number of digits for displaying coefficients.
"""
function set_default_digits(d::Int)
    set_default(:DEFAULT_DIGITS, d)
end

"""
    get_default_digits() -> Int

Returns the current default number of digits for displaying coefficients.
"""
function get_default_digits()
    return get_default(:DEFAULT_DIGITS)
end

# Define the ComplexRational type with an inner constructor that does the conversion.
struct ComplexRational <: Number
    a::Int
    b::Int
    c::Int  # Denominator, must be nonzero.
    function ComplexRational(a::Int, b::Int, c::Int)
        if c == 0
            throw(ArgumentError("Denominator cannot be zero"))
        end
        # Ensure the denominator is positive.
        if c < 0
            a, b, c = -a, -b, -c
        end
        # Normalize: divide a, b, and c by their greatest common divisor.
        g = gcd(gcd(abs(a), abs(b)), c)
        a, b, c = div(a, g), div(b, g), div(c, g)
        # If the denominator has too many digits, return a Complex{Float64}
        if length(string(c)) > get_default(:DEFAULT_DIGITS)
            # Compute a float complex number from a, b, c.
            return complex(a, b) / c
        end
        new(a, b, c)
    end
end

"""
    is_negative(x::ComplexRational) -> Bool

Assuming `x` is given as (a + i·b)/c with c > 0 (as enforced by its constructor),
this function returns true if the number should be printed with a negative sign.  
The convention here is:
  - If the real part (a) is nonzero, return true if a < 0.
  - Otherwise (if a is zero), return true if the imaginary part (b) is negative.
If both are zero, returns false.
"""
function is_negative(x::ComplexRational)
    if x.a != 0
        return x.a < 0
    else
        return x.b < 0
    end
end

# ------------------------------------------------------------------------------
# Overload conversion to Complex for interactions with Complex numbers.
# ------------------------------------------------------------------------------
import Base: convert, complex
convert(::Type{Complex}, x::ComplexRational) = complex(x.a / x.c, x.b / x.c)
complex(x::ComplexRational) = convert(Complex, x)

# ------------------------------------------------------------------------------
# Basic arithmetic operations. They construct new ComplexRational values and
# rely on the inner constructor to perform the conversion check.
# ------------------------------------------------------------------------------
import Base: +, -, *, /

function +(x::ComplexRational, y::ComplexRational)
    new_a = x.a * y.c + y.a * x.c
    new_b = x.b * y.c + y.b * x.c
    new_c = x.c * y.c
    return ComplexRational(new_a, new_b, new_c)
end

function -(x::ComplexRational, y::ComplexRational)
    new_a = x.a * y.c - y.a * x.c
    new_b = x.b * y.c - y.b * x.c
    new_c = x.c * y.c
    return ComplexRational(new_a, new_b, new_c)
end

function *(x::ComplexRational, y::ComplexRational)
    new_a = x.a * y.a - x.b * y.b
    new_b = x.a * y.b + x.b * y.a
    new_c = x.c * y.c
    return ComplexRational(new_a, new_b, new_c)
end

function /(x::ComplexRational, y::ComplexRational)
    if iszero(y)
        throw(DivideError())
    end
    # For y = (e + i*f)/y.c, we have:
    #   x/y = [(x.a + i*x.b)/x.c] / [(y.a + i*y.b)/y.c]
    #       = (x.a + i*x.b) * (y.c / (y.a + i*y.b)) / x.c
    # Multiply numerator and denominator by the complex conjugate of (y.a + i*y.b):
    #   = ((x.a + i*x.b)*(y.a - i*y.b)*y.c) / (x.c*(y.a^2 + y.b^2))
    e, f, yc = y.a, y.b, y.c
    denom_y = e^2 + f^2
    new_a = x.a * yc * e + x.b * yc * f
    new_b = x.b * yc * e - x.a * yc * f
    new_c = x.c * denom_y
    return ComplexRational(new_a, new_b, new_c)
end

# ------------------------------------------------------------------------------
# iszero, real, and imag functions.
# ------------------------------------------------------------------------------
import Base: iszero, isone, real, imag

iszero(x::ComplexRational) = (x.a == 0 && x.b == 0)
isone(x::ComplexRational) = (x.a == x.c && x.b == 0)

function real(x::ComplexRational)
    return ComplexRational(x.a, 0, x.c)
end

function imag(x::ComplexRational)
    return ComplexRational(x.b, 0, x.c)
end

# For a Rational number, simply convert it to ComplexRational (with zero imaginary part).
function crationalize(x::Rational; tol::Real=10^-12)
    return ComplexRational(numerator(x), 0, denominator(x))
end

# For a Complex number, crationalize the real and imaginary parts separately.
function crationalize(x::Complex; tol::Real=1e-12)
    # Here, Julia's built-in crationalize for Float64 is used for the real parts.
    r = rationalize(real(x), tol=tol)
    i = rationalize(imag(x), tol=tol)
    # Extract numerators and denominators.
    p, q = numerator(r), denominator(r)
    r_im, s = numerator(i), denominator(i)
    # Combine the two rationals to get a common denominator (using q*s).
    a = p * s
    b = r_im * q
    c = q * s
    return ComplexRational(a, b, c)
end

# For a Rational number, simply convert it to ComplexRational (with zero imaginary part).
function crationalize(x::T; tol::Real=10^-12) where {T<:Number}
    return crationalize(complex(x), tol=tol)
end

import Base: +, -

# When adding a ComplexRational and a Complex, try to convert the Complex to a ComplexRational.
+(x::ComplexRational, y::Complex) = x + crationalize(y, tol=1e-12)
+(x::Complex, y::ComplexRational) = crationalize(x, tol=1e-12) + y

# When adding a ComplexRational and a Float64, promote the Float64 to a Complex first.
+(x::ComplexRational, y::Float64) = x + crationalize(Complex(y, 0), tol=1e-12)
+(x::Float64, y::ComplexRational) = crationalize(Complex(x, 0), tol=1e-12) + y

# When adding a ComplexRational and a Rational, convert the Rational to a ComplexRational.
+(x::ComplexRational, y::Rational) = x + crationalize(y, tol=1e-12)
+(x::Rational, y::ComplexRational) = crationalize(x, tol=1e-12) + y

# When adding a ComplexRational and an Int, convert the Int to a Rational first.
+(x::ComplexRational, y::Int) = x + crationalize(Rational(y), tol=1e-12)
+(x::Int, y::ComplexRational) = crationalize(Rational(x), tol=1e-12) + y


# And similarly for subtraction:
-(x::ComplexRational) = ComplexRational(-x.a, -x.b, x.c)

-(x::ComplexRational, y::Complex) = x - crationalize(y, tol=1e-12)
-(x::Complex, y::ComplexRational) = crationalize(x, tol=1e-12) - y

-(x::ComplexRational, y::Float64) = x - crationalize(Complex(y, 0), tol=1e-12)
-(x::Float64, y::ComplexRational) = crationalize(Complex(x, 0), tol=1e-12) - y

-(x::ComplexRational, y::Rational) = x - crationalize(y, tol=1e-12)
-(x::Rational, y::ComplexRational) = crationalize(x, tol=1e-12) - y

-(x::ComplexRational, y::Int) = x - crationalize(Rational(y), tol=1e-12)
-(x::Int, y::ComplexRational) = crationalize(Rational(x), tol=1e-12) - y

import Base: *, /

# --- Multiplication ---

# When multiplying with a Complex, convert the Complex to a ComplexRational.
*(x::ComplexRational, y::Complex) = x * crationalize(y, tol=1e-12)
*(x::Complex, y::ComplexRational) = crationalize(x, tol=1e-12) * y

# When multiplying with a Float64, promote it to a Complex.
*(x::ComplexRational, y::Float64) = x * crationalize(Complex(y, 0), tol=1e-12)
*(x::Float64, y::ComplexRational) = crationalize(Complex(x, 0), tol=1e-12) * y

# When multiplying with a Rational, convert the Rational to a ComplexRational.
*(x::ComplexRational, y::Rational) = x * crationalize(y, tol=1e-12)
*(x::Rational, y::ComplexRational) = crationalize(x, tol=1e-12) * y

# When multiplying with an Int, convert the Int to a Rational first.
*(x::ComplexRational, y::Int) = x * crationalize(Rational(y), tol=1e-12)
*(x::Int, y::ComplexRational) = crationalize(Rational(x), tol=1e-12) * y


# --- Division ---

# When dividing by a Complex, convert the Complex to a ComplexRational.
/(x::ComplexRational, y::Complex) = x / crationalize(y, tol=1e-12)
/(x::Complex, y::ComplexRational) = crationalize(x, tol=1e-12) / y

# When dividing by a Float64, promote it to a Complex.
/(x::ComplexRational, y::Float64) = x / crationalize(Complex(y, 0), tol=1e-12)
/(x::Float64, y::ComplexRational) = crationalize(Complex(x, 0), tol=1e-12) / y

# When dividing by a Rational, convert the Rational to a ComplexRational.
/(x::ComplexRational, y::Rational) = x / crationalize(y, tol=1e-12)
/(x::Rational, y::ComplexRational) = crationalize(x, tol=1e-12) / y

# When dividing by an Int, convert the Int to a Rational first.
/(x::ComplexRational, y::Int) = x / crationalize(Rational(y), tol=1e-12)
/(x::Int, y::ComplexRational) = crationalize(Rational(x), tol=1e-12) / y

import Base: abs, abs2, sqrt
# abs2 for ComplexRational: squared amplitude.
abs2(x::ComplexRational) = ComplexRational(x.a^2 + x.b^2, 0, x.c^2)
sqrt(x::ComplexRational) = crationalize(sqrt(complex(x)))
abs(x::ComplexRational) = crationalize(sqrt(abs2(x)))  # converting to float for sqrt

import Base: adjoint, conj
adjoint(x::ComplexRational) = ComplexRational(x.a, -x.b, x.c)  # conjugate
conj(x::ComplexRational) = ComplexRational(x.a, -x.b, x.c)  # conjugate

import Base: one, zero, inv
# Identity (multiplicative identity)
one(::Type{ComplexRational}) = ComplexRational(1, 0, 1)
one(x::ComplexRational) = ComplexRational(1, 0, 1)

# Zero (additive identity)
zero(::Type{ComplexRational}) = ComplexRational(0, 0, 1)
zero(x::ComplexRational) = ComplexRational(0, 0, 1)

# Inverse (reciprocal)
inv(x::ComplexRational) = ComplexRational(x.c * x.a, -x.c * x.b, x.a^2 + x.b^2)

import Base: ^
function ^(x::ComplexRational, n::Integer)#::ComplexRational
    if n < 0
        return inv(x^(-n))
    end

    result = ComplexRational(1, 0, 1)  # Identity for multiplication
    base = x
    exp = n
    while exp > 0
        if isodd(exp)
            result *= base
        end
        base *= base
        exp ÷= 2
    end
    return result
end

import Base: ==, isless, <, <=, >, >=

# Equality: we use normalized components.
function ==(x::ComplexRational, y::ComplexRational)
    return x.a == y.a && x.b == y.b && x.c == y.c
end

# Define a custom sort key that returns a 3-element tuple.
function custom_sort_key(x)
    # For any real number (including when a Complex or ComplexRational is pure real):
    if iszero(imag(x))
        # Convert to a Complex to extract the real part as a Float.
        # This works even for Rational numbers.
        return (0, real(complex(x)), 0.0)
        # For pure imaginary numbers:
    elseif iszero(real(x))
        # Convert the imaginary part to a Complex and take its real part (its magnitude).
        # We take absolute value to sort by magnitude.
        return (1, real(complex(imag(x))), 0.0)
    else
        # For mixed numbers, use abs(x) as the magnitude.
        # (abs works for built-in Complex; for ComplexRational it is defined to convert via complex)
        return (2, real(complex(abs(x))), 0.0)
    end
end

# A helper function that sorts any vector of numbers (which might include Real, ComplexRational, and Complex).
function sort_mixed(v::AbstractVector)
    return sort(v, by=custom_sort_key)
end

###### printing and display ############################################
function complexrational_plain(x::ComplexRational)
    # Prepare the denominator string: if c==1, it's omitted.
    denom_str = x.c == 1 ? "" : "/" * string(x.c)

    if x.a == 0 && x.b == 0
        return "0"
    elseif x.b == 0
        # Only real part.
        return (x.a < 0 ? "-" * string(abs(x.a)) : string(x.a)) * denom_str
    elseif x.a == 0
        # Only imaginary part.
        return (x.b < 0 ? "-" * string(abs(x.b)) : string(x.b)) * "i" * denom_str
    else
        # Both parts.
        if x.a < 0
            # Global negative: factor out the minus and flip the sign of the imaginary part.
            # For example, if x = (-3 + 2i)/c then we want to display: - (3 - 2i)/c.
            local_im_sign = x.b > 0 ? " - " : " + "
            num_str = "(" * string(abs(x.a)) * local_im_sign * string(abs(x.b)) * "i" * ")"
            return "-" * num_str * denom_str
        else
            local_im_sign = x.b < 0 ? " - " : " + "
            num_str = "(" * string(x.a) * local_im_sign * string(abs(x.b)) * "i" * ")"
            return num_str * denom_str
        end
    end
end

function complexrational_latex(x::ComplexRational)
    if x.a == 0 && x.b == 0
        return "0"
    elseif x.b == 0
        # Only real part.
        if x.a < 0
            return "-" * (x.c == 1 ? string(abs(x.a)) : "\\frac{" * string(abs(x.a)) * "}{" * string(x.c) * "}")
        else
            return x.c == 1 ? string(x.a) : "\\frac{" * string(x.a) * "}{" * string(x.c) * "}"
        end
    elseif x.a == 0
        # Only imaginary part.
        # Now we print the coefficient followed by i.
        if x.b < 0
            return "-" * (x.c == 1 ? string(abs(x.b)) * "i" : "\\frac{" * string(abs(x.b)) * "i}{" * string(x.c) * "}")
        else
            return x.c == 1 ? string(x.b) * "i" : "\\frac{" * string(x.b) * "i}{" * string(x.c) * "}"
        end
    else
        # Both parts.
        if x.a < 0
            # Global negative: factor out the minus sign and flip the sign of the imaginary part.
            local_im_sign = x.b > 0 ? " - " : " + "
            frac_str = x.c == 1 ? "(" * string(abs(x.a)) * local_im_sign * string(abs(x.b)) * "i )" :
                       "\\frac{" * string(abs(x.a)) * local_im_sign * string(abs(x.b)) * "i}{" * string(x.c) * "}"
            return "-" * frac_str
        else
            local_im_sign = x.b < 0 ? " - " : " + "
            return x.c == 1 ? "(" * string(x.a) * local_im_sign * string(abs(x.b)) * "i )" :
                   "\\frac{" * string(x.a) * local_im_sign * string(abs(x.b)) * "i}{" * string(x.c) * "}"
        end
    end
end
function complexrational2str(x::ComplexRational, do_latex::Bool=false)::String
    if do_latex
        return complexrational_latex(x)
    else
        return complexrational_plain(x)
    end
end


import Base: show
function show(io::IO, ::MIME"text/plain", x::ComplexRational)
    print(io, complexrational_plain(x))
end

# For LaTeX-enabled I/O (e.g. in Jupyter notebooks or environments that support it):
function show(io::IO, ::MIME"text/latex", x::ComplexRational)
    print(io, latexstring(complexrational_latex(x)))
end

function show(io::IO, x::ComplexRational)
    show(io::IO, MIME"text/plain"(), x::ComplexRational)
end

end # module CRational
# ComplexRationals.jl

[![Build Status](https://github.com/Ntropic/ComplexRationals.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Ntropic/ComplexRationals.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

**ComplexRationals.jl** defines a `ComplexRational` number type for Julia:  
complex numbers whose real and imaginary parts are **exact integers divided by a shared denominator**.  
This allows for precise arithmetic on small complex numbers without rounding errors,  
while automatically promoting to floating-point `Complex{Float64}` if necessary.

---

## Features

- Exact complex arithmetic with support for:
  - Basic Operations: `+`, `-`, `*`, `/`, `^`, `inv`
  - Comparison operators: `==`, `!=`, `<`, `<=`, `>`, `>=`
  - Conversion from and to Julia's built-in number types
- String representation (plain and LaTeX)

## Installation 
```julia
using Pkg
Pkg.add("ComplexRational")
```
or alternatively
```julia
using Pkg
Pkg.add(url="https://github.com/Ntropic/ComplexRational.jl")
```

## Examples
```julia
using ComplexRational

# Create a ComplexRational number
z = ComplexRational(3, 2, 5)  # (3 + 2i)/5

# Arithmetic
w = ComplexRational(1, -1, 2)
sum = z + w
product = z * w
quotient = z / w

# Exponentiation
squared = z^2
inverse = inv(z)
zero_power = z^0

# Conversion
r = 3//4
cr = crationalize(r)            # Rational -> ComplexRational
cr2 = crationalize(1.5 + 2.0im)  # Complex Float -> ComplexRational

# Access real and imaginary parts
re = real(z)   # (3/5)
im = imag(z)   # (2/5)

# String display
plain_str = string(z)           # "(3 + 2i)/5"
latex_str = string(z, do_latex=true)     # LaTeX code for display
```

## Author 
- [Michael Schilling](https://github.com/Ntropic)
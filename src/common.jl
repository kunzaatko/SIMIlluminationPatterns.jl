# TODO: It should be possible to input missing values in parameters of illumination patterns. <19-11-23> 

@doc raw"""
Super type for `N`-dimensional illumation patterns. Any particular illumination pattern is a subtype of this. See the 
documentation for [`IlluminationPattern` interface](@ref "Unified `IlluminationPattern` Interface") to see how to implement your own illumination pattern subtype.

See also: [`Harmonic2D`](@ref)

# Implementation
When implementing a new subtype `A <: IlluminationPattern{N}` you need to implement the following methods:
+ `(pattern::A)(r::Vararg{Length, N})::Real` representing the intensity of `pattern` at the position `r`

and these optional methods:

TODO

The inner constructor of the subtype `A` should check that `N::Integer` unless your subtype has a specified dimension.
```@example
struct A1{N} <: IlluminationPattern{N}
    function A{N}() where {N}
        @assert N isa Integer
        new{N}()
    end
end

struct A2 <: IlluminationPattern{2}
    function A2()
        # No need for any assertions
        new()
    end
end
```
"""
abstract type IlluminationPattern{N} end
const IP{N} = IlluminationPattern{N}

(ip::IP{N})(::Vararg{Length,N}) where {N} = error("Not Implemented!")

for func in (:(==), :isequal, :isapprox)
    @eval function Base.$func(ip1::A, ip2::B; kwargs...) where {A<:IlluminationPattern,B<:IlluminationPattern}
        nameof(A) === nameof(B) || return false
        fields = fieldnames(A)
        fields === fieldnames(B) || return false

        for f in fields
            isdefined(ip1, f) && isdefined(ip2, f) || return false
            # perform equivalence check to support types that have no defined equality, such
            # as `missing`
            getfield(ip1, f) === getfield(ip2, f) || $func(getfield(ip1, f), getfield(ip2, f); kwargs...) || return false
        end

        return true
    end
end

function Base.hash(ip::IP, h::UInt) where {IP<:IlluminationPattern}
    hashed = hash(IlluminationPattern, h)
    hashed = hash(nameof(IP), hashed)

    for f in fieldnames(IP)
        hashed = hash(getfield(ip, f), hashed)
    end

    return hashed
end

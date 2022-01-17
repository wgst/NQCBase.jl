
export Atoms

"""
    Atoms{S,T<:AbstractFloat}

Basic atomic parameters: element symbols, numbers and masses

Masses are converted to atomic units.
Constructed using either element symbols or masses.

```jldoctest
julia> Atoms(:H)
Atoms{1, Float64}([:H], UInt8[0x01], [1837.4715941070515])

julia> Atoms([:H, :H, :H, :C])
Atoms{4, Float64}([:H, :H, :H, :C], UInt8[0x01, 0x01, 0x01, 0x06], [1837.4715941070515, 1837.4715941070515, 1837.4715941070515, 21894.713607956142])

julia> Atoms([100, 200])
Atoms{2, Float64}([:X, :X], UInt8[0x00, 0x00], [100.0, 200.0])
```
"""
struct Atoms{S,T<:AbstractFloat}
    types::SVector{S,Symbol}
    numbers::SVector{S,UInt8}
    masses::SVector{S,T}
end

function Atoms{T}(atom_types::AbstractVector{Symbol}) where {T}
    S = length(atom_types)
    types = SVector{S,Symbol}(atom_types)
    numbers = SVector{S,UInt8}([element.number for element in elements[atom_types]])
    masses = SVector{S,T}([austrip(element.atomic_mass) for element in elements[atom_types]])
    Atoms{S,T}(types, numbers, masses)
end

function Atoms{T}(masses::AbstractVector) where {T}
    S = length(masses)
    types = SVector{S,Symbol}(fill(:X, S))
    numbers = SVector{S,UInt8}(zeros(S))
    masses = SVector{S,T}(austrip.(masses))
    Atoms{S,T}(types, numbers, masses)
end

Atoms{T}(atom_type) where {T} = Atoms{T}([atom_type])
Atoms(atom_types) = Atoms{Float64}(atom_types)

Base.length(::Atoms{S}) where {S} = S
Base.range(::Atoms{S}) where {S} = range(1; length=S)

Base.IndexStyle(::Atoms) = IndexLinear()
Base.getindex(A::Atoms{S,T}, i::Int) where {S,T} = Atoms{T}(A.types[i])
Base.getindex(A::Atoms{S,T}, i::AbstractRange) where {S,T} = Atoms{T}(A.types[i])

"""
A toy system, representing a phase state (q,p), subject to evolution under a potential V with gradient ∇V. 
The running sum and sum of squares of a given set of observables are recorded throughout the simulation.
"""
mutable struct ToySystem{D,F,dF,BC,O}
    q::SVector{D,Float64}
    p::SVector{D,Float64}

    last_q::SVector{D,Float64}
    
    V::F
    ∇V::dF
    
    boundary_condition!::BC
    
    observables::O
    O_sums::Vector{Float64}
    sq_O_sums::Vector{Float64}
end

ToySystem(q::SVector{D,Float64},p::SVector{D,Float64},V::Function,∇V::Function,bc!::Function,observables) where {D} = ToySystem{D,typeof(V),typeof(∇V),typeof(bc!),typeof(observables)}(q,p,q,V,∇V,bc!,observables,zeros(length(observables)),zeros(length(observables)))
ToySystem(q,p,V,∇V,bc!) = ToySystem(q,p,V,∇V,bc!,[])
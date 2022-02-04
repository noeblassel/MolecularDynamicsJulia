using Unitful

const kᵦ= 1.380649e-23 #Boltzmann's constant (J * K^-1)


#Unitful units of reference

const uₘ=u"u" #mass
const uₗ=u"nm" #length
const uₑ=u"J"  #energy
const uₖ=u"J * K^-1"   #Boltzmann's constant
const uₜ=uₗ*(uₘ/uₑ)^(1//2) #Time
const u_T=uₑ/uₖ #Temperature
const uₚ=uₑ/uₗ^3 #Pressure
const uᵣ=uₗ^(-3) #Density (measured in number of particles per unit space)

lj_params=Dict("argon"=>(m=39.948,σ=0.3405,ϵ=1.66e-21)) #simply add parameters for various species described by the lennard-jones potential


get_reduced_length(species::AbstractString,l::Real)=l/lj_params[species].σ
get_reduced_energy(species::AbstractString,e::Real)=e/lj_params[species].ϵ
get_reduced_time(species::AbstractString,t::Real)=t*sqrt(lj_params[species.ϵ]/(lj_params[species].m*lj_params[species].σ^2))
get_reduced_temperature(species::AbstractString,T::Real)=T*kᵦ/lj_params[species].ϵ
get_reduced_pressure(species::AbstractString,p::Real)=p*lj_params[species].σ^3/lj_params[species].ϵ
get_reduced_density(species::AbstractString,ρ::Real)=ρ*lj_params[species].σ^3

get_reduced_length(species::AbstractString,l::L) where {L<:Unitful.Quantity}=get_reduced_length(species,ustrip(uₗ,l))
get_reduced_energy(species::AbstractString,e::E) where {E<:Unitful.Quantity}=get_reduced_energy(species,ustrip(uₑ,e))
get_reduced_time(species::AbstractString,t::T) where {T<:Unitful.Quantity}=get_reduced_time(species,ustrip(uₜ,t))
get_reduced_temperature(species::AbstractString,T::Θ) where {Θ<:Unitful.Quantity}=get_reduced_temperature(species,ustrip(u_T,T))
get_reduced_pressure(species::AbstractString,p::P) where {P<:Unitful.Quantity}=get_reduced_pressure(species,ustrip(uₚ,p))
get_reduced_density(species::AbstractString,ρ::R) where {R<:Unitful.Quantity}=get_reduced_density(species,ustrip(uᵣ,ρ))

get_physical_length(species::AbstractString,l::Real)=l*lj_params[species].σ*uₗ
get_physical_energy(species::AbstractString,e::Real)=e*lj_params[species].ϵ*uₑ
get_physical_time(species::AbstractString,t::Real)=t*sqrt((lj_params[species].m*lj_params[species].σ^2)/lj_params[species.ϵ])*uₜ
get_physical_temperature(species::AbstractString,T::Real)=T*(lj_params[species].ϵ/kᵦ)*u_T
get_physical_pressure(species::AbstractString,p::Real)=p*(lj_params[species].ϵ/lj_params[species].σ^3)*uₚ
get_physical_density(species::AbstractString,ρ::Real)=(ρ/lj_params[species].σ^3)*uᵣ
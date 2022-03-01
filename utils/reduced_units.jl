const kb= ustrip(Unitful.k) #Boltzmann's constant (J * K^-1)


#Physical units of reference

const uₘ=u"u" #mass
const uₗ=u"nm" #length  
const uₑ=u"eV"  #energy
const uₖ=u"J * K^-1"   #Boltzmann's constant

#Derived
const uₜ=uₗ*(uₘ/uₑ)^(1//2) #Time
const u_T=uₑ/uₖ #Temperature
const uₚ=uₑ/uₗ^3 #Pressure
const uᵣ=uₗ^(-3) #Density (measured in number of particles per unit space)

lj_params=Dict(:Ne=>(m=20.1797,σ=0.2801,ϵ=2.923e-3),:Ar=>(m=39.948,σ=0.341,ϵ=1.03e-2),:Kr=>(m=83.798,σ=0.36274,ϵ=1.401e-2),:Xe=>(m=131.293,σ=0.3949,ϵ=1.949e-2)) #simply add parameters for various species described by the lennard-jones potential


get_reduced_length(species::Symbol,l::Real)=l/lj_params[species].σ
get_reduced_energy(species::Symbol,e::Real)=e/lj_params[species].ϵ
get_reduced_time(species::Symbol,t::Real)=t*sqrt(lj_params[species].ϵ/(lj_params[species].m*lj_params[species].σ^2))
get_reduced_temperature(species::Symbol,T::Real)=T*kb/lj_params[species].ϵ
get_reduced_pressure(species::Symbol,p::Real)=p*lj_params[species].σ^3/lj_params[species].ϵ
get_reduced_density(species::Symbol,ρ::Real)=ρ*lj_params[species].σ^3

get_reduced_length(species::Symbol,l::L) where {L<:Unitful.Quantity}=get_reduced_length(species,ustrip(uₗ,l))
get_reduced_energy(species::Symbol,e::E) where {E<:Unitful.Quantity}=get_reduced_energy(species,ustrip(uₑ,e))
get_reduced_time(species::Symbol,t::T) where {T<:Unitful.Quantity}=get_reduced_time(species,ustrip(uₜ,t))
get_reduced_temperature(species::Symbol,T::Θ) where {Θ<:Unitful.Quantity}=get_reduced_temperature(species,ustrip(u_T,T))
get_reduced_pressure(species::Symbol,p::P) where {P<:Unitful.Quantity}=get_reduced_pressure(species,ustrip(uₚ,p))

function get_reduced_density(species::Symbol,ρ::R) where {R<:Unitful.Quantity}
    if dimension(ρ)==dimension(u"mol*m^-3")
    return get_reduced_density(species,ustrip(uᵣ,ρ*Unitful.Na))
    elseif dimension(ρ)==dimension(u"kg*m^-3")
        return get_reduced_density(species,ustrip(uᵣ,ρ/(lj_params[species].m*uₘ)))
    end
end

get_physical_length(species::Symbol,l::Real)=uconvert(u"nm",l*lj_params[species].σ*uₗ)
get_physical_energy(species::Symbol,e::Real)=uconvert(u"J",e*lj_params[species].ϵ*uₑ)
get_physical_time(species::Symbol,t::Real)=uconvert(u"ns",t*sqrt((lj_params[species].m*lj_params[species].σ^2)/lj_params[species].ϵ)*uₜ)
get_physical_temperature(species::Symbol,T::Real)=uconvert(u"K",T*(lj_params[species].ϵ/kb)*u_T)
get_physical_pressure(species::Symbol,p::Real)=uconvert(u"MPa",p*(lj_params[species].ϵ/lj_params[species].σ^3)*uₚ)

function get_physical_density(species::Symbol,ρ::Real;molar=false)
    if molar
        return uconvert(u"mol * m^-3",(ρ/lj_params[species].σ^3)*uᵣ/Unitful.Na)
    else
        return uconvert(u"kg * m^-3", (ρ/lj_params[species].σ^3)*uᵣ*uₘ*lj_params[species].m)
    end
end

periodic_potential(q::Float64,L::Float64)=sin(2π*q/L)
minus_d_periodic_potential(q::Float64,L::Float64)=-2π*cos(2π*q/L)/L

quadratic_potential(q::Float64,L::Float64=nothing,a::Float64=1.0)=a*q^2/2
minus_d_quadratic_potential(q::Float64,L::Float64=nothing,a::Float64=1.0)=-a*q

double_well_potential(q::Float64,L::Float64=nothing,a::Float64=1.0,b::Float64=2.0)=quadratic_potential(q,a)+b*exp(-q^2/2)
minus_d_double_well_potential(q::Float64,L::Float64=nothing,a::Float64=1.0,b::Float64=2.0)=minus_d_quadratic_potential(q,a)+b*q*exp(-q^2/2)

periodic_potential(q::Float64,L::Float64)=sin(2π*q/L)
minus_d_periodic_potential(q::Float64,L::Float64)=-2π*cos(2π*q/L)/L

quadratic_potential(q::Float64,L::Float64,a::Float64=1.0)=a*q^2/2
minus_d_quadratic_potential(q::Float64,L::Float64,a::Float64=1.0)=-a*q

double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5)=a*q^2/2+b*exp(-q^2/(2c^2))
minus_d_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5)=-a*q+b*q*exp(-q^2/(2c^2))/c^2

tilted_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5,d::Float64=1.0)=a*q^2/2+b*exp(-q^2/(2c^2))+d*q
minus_d_tilted_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5,d::Float64=1.0)=-a*q+b*q*exp(-q^2/(2c^2))/c^2+d

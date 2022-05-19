periodic_potential(q::Float64,L::Float64)=sin(2π*q/L)
grad_periodic_potential(q::Float64,L::Float64)=2π*cos(2π*q/L)/L
d2_periodic_potential(q::Float64,L::Float64)=-4π^2*sin(2π*q/L)/L^2

quadratic_potential(q::Float64,L::Float64,a::Float64=1.0)=a*q^2/2
grad_quadratic_potential(q::Float64,L::Float64,a::Float64=1.0)=a*q
d2_quadratic_potential(q::Float64,L::Float64,a::Float64=1.0)=a

double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5)=a*q^2/2+b*exp(-q^2/(2c^2))
grad_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5)=a*q-b*q*exp(-q^2/(2c^2))/c^2
d2_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5)=a-b*exp(-q^2/(2c^2))/c^2+b*q^2*exp(-q^2/(2c^2))/c^4

tilted_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5,d::Float64=1.0)=a*q^2/2+b*exp(-q^2/(2c^2))+d*q
grad_tilted_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5,d::Float64=1.0)=a*q-b*q*exp(-q^2/(2c^2))/c^2+d
d2_tilted_double_well_potential(q::Float64,L::Float64,a::Float64=1.0,b::Float64=4.0,c::Float64=0.5,d::Float64=1.0)=a-b*exp(-q^2/(2c^2))/c^2+b*q^2*exp(-q^2/(2c^2))/c^4

pathological_potential(q::Float64,L::Float64)=log(1+q^2)-log(1+sin(q^4)^2)
grad_pathological_potential(q::Float64,L::Float64)=-(8q^3*sin(q^4)*cos(q^4))/(1 + sin(q^4)^2)+(2q)/(1 + q^2)

periodic_potential_2D(qx::Float64,qy::Float64,L::Float64=1.0)=cos(2π*qx/L)*sin(2π*qy/L)
minus_grad_periodic_potential_2D(qx::Float64,qy::Float64,L::Float64=1.0)=(-2π/L)*[-sin(2π*qx/L)*sin(2π*qy/L),cos(2π*qy/L)*cos(2π*qx/L)]
d2_periodic_potential_2D(qx::Float64,qy::Float64,L::Float64=1.0)=[-4π^2*cos(2π*qx/L)/L^2 0;0 -4π^2*sin(2π*qy/L)/L^2]




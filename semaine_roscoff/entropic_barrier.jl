using Base.Threads, LinearAlgebra, StaticArrays, Plots

include("ToySystem.jl")
include("ReflectingBoundaryCondition.jl")
include("BAOAB.jl")
#include("GenParRep.jl")

"""
Domain shape: WxH grid of circular states connected by narrow corridors
O--O--O
|  |  |
O--O--O
|  |  |
O--O--O
"""

W=3
H=3

L=2.0
R=3.0
h=0.2

#state_x = i => i*(L+2R)-R-L/2<= q_x <=i*(L+2R)+R+L/2=> i<= (q_x +R +L/2)/(L+2R), i>= (q_x - R- L/2)/(L+2R)=> i∈ [q_x/(L+2R)-1/2,q_x/(L+2R)+1/2] => i =round(q_x/(L+2R))

centers=collect(SVector{2,Float64}[(i*(L+2R),j*(L+2R)) for i=0:W-1,j=0:H-1])
centers=reshape(centers,W*H)

xmax=(W-1)*(L+2R)+R
ymax=(H-1)*(L+2R)+R

urns=[Sphere(center,R) for center in centers]
h_pipes=[OrthoRect(SVector(i*(L+2R)+R-h,j*(L+2R)+0.5h),SVector(i*(L+2R)+R+L+h,j*(L+2R)-0.5h)) for i=0:W-2,j=0:H-1]
h_pipes=reshape(h_pipes,(W-1)*H)

v_pipes=[OrthoRect(SVector(i*(L+2R)-0.5h,j*(L+2R)+R+L+h),SVector(i*(L+2R)+0.5h,j*(L+2R)+R-h)) for i=0:W-1,j=0:H-2]
v_pipes=reshape(v_pipes,W*(H-1))

domain=CompositeShape(vcat(urns,h_pipes,v_pipes))

p0=SVector(0.0,0.0)
q0=rand(centers)

p0=SVector(0.0,0.0)
q0=SVector(0.0,0.0)

obsA(sys)=norm(sys.q-C_A)
obsB(sys)=norm(sys.q-C_B)

V(q)=0.0
∇V(q)=[0.0,0.0]

function get_state(sys)
    state_i=round(Int64,sys.q[1]/(L+2R))
    state_j=round(Int64,sys.q[2]/(L+2R))

    return W*state_i + state_j
end

function bc!(sys)
    if isin(domain,sys.q)
        return nothing
    else
        i=intersection(domain,sys.last_q,sys.q)
        if i!=nothing
            pt,n=i
            q_normal_component=dot(sys.q-pt,n)
            sys.q-=2q_normal_component*n
            p_normal_component=dot(sys.p,n)
            sys.p-=2p_normal_component*n
        end
    end
end

sys=ToySystem(q0,p0,V,∇V,bc!,[])
sim=BAOABIntegrator(0.01,1.0,0.7)

@time for i=1:10000
    simulate!(sys,sim,100)
    println(get_state(sys))
end
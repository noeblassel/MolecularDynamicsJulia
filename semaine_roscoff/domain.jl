using LinearAlgebra, StaticArrays

include("ReflectingBoundaryCondition.jl")

"""
Domain shape: WxH grid of circular states (of radius R) connected by narrow corridors (of width h and length L)
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
h=0.05

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

function get_state(sys,state)
    state_i=round(Int,sys.q[1]/(L+2R))
    state_j=round(Int,sys.q[2]/(L+2R))
    c_x=state_i*(L+2R)
    c_y=state_j*(L+2R)

    if (sys.q[1]-c_x)^2 + (sys.q[2]-c_y)^2 <= R^2
        return W*state_i + state_j
    else
        return state
    end
end

function get_state(sys)
    state_i=round(Int,sys.q[1]/(L+2R))
    state_j=round(Int,sys.q[2]/(L+2R))
    return W*state_i + state_j
end


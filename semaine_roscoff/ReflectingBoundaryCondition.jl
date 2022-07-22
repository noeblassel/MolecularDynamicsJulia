abstract type Shape end

struct Sphere<:Shape
    center::SVector{2,Float64}
    radius::Float64
    sq_radius::Float64
end

Sphere(q,r)=Sphere(q,r,r^2)

isin(shape::Sphere,q)= dot(q-shape.center,q-shape.center) <= shape.sq_radius

struct OrthoRect<:Shape
    tl::SVector{2,Float64}
    tr::SVector{2,Float64}
    bl::SVector{2,Float64}
    br::SVector{2,Float64}
end

OrthoRect(tl,br)=OrthoRect(tl,SVector(br[1],tl[2]),SVector(tl[1],br[2]),br)

isin(shape::OrthoRect,q) = (q[1]<=shape.br[1]) && (q[1]>=shape.tl[1]) && (q[2]>=shape.br[2]) && (q[2]<=shape.tl[2])

function line_intersection(A,B,C,D)
    M=hcat(A-B,D-C)
    (det(M) == 0) && return nothing
    s,t = inv(M)*(D-B)
    if (0<=s<1) && (0<=t<=1)
        return t*A+(1-t)*B
    else
        return nothing
    end
end

function intersection(shape::OrthoRect,q,last_q)
    i=line_intersection(shape.tl,shape.tr,q,last_q)
    if i!=nothing
        return i,SVector(0.0,-1.0)
    end

    i=line_intersection(shape.tr,shape.br,q,last_q)
    if i!=nothing
        return i,SVector(-1.0,0.0)
    end

    i=line_intersection(shape.bl,shape.br,q,last_q)
    if i!=nothing
        return i,SVector(0.0,1.0)
    end
    i=line_intersection(shape.bl,shape.tl,q,last_q)
    if i!=nothing
        return i,SVector(1.0,0.0)
    end
    return nothing
end

struct CompositeShape<:Shape
    subshapes::Vector{Shape}
end

isin(shape::CompositeShape,q)=any(isin(subshape,q) for subshape in shape.subshapes)
"""
Assuming A is inside shape, returns nothing if B is inside the shape, or the first intersection point with the boundary along the ray [A,B] along with the normal vector at that intersection point.
"""
function intersection(shape::Sphere,A,B)
    isin(shape,B) && return nothing
    a=dot(B-A,B-A)
    b=2dot(B-A,A-shape.center)
    c=dot(A-shape.center,A-shape.center)-shape.sq_radius
    Δ=sqrt(b^2-4*a*c)
    s=(-b + Δ)/2a
    i= A +s*(B-A)
    return i,(shape.center-i)/shape.radius
end

function intersection(shape::CompositeShape,A,B)
    isin(shape,B) && return nothing
    last_subshape=first(subshape for subshape in shape.subshapes if isin(subshape,A))
    return intersection(last_subshape,A,B)
end

ref_circle_x=cos.(range(0,2π,100))
ref_circle_y=sin.(range(0,2π,100))

#Plots.plot!(shape::Sphere,color)=plot!(shape.center[1] .+ shape.radius*ref_circle_x, shape.center[2] .+ shape.radius*ref_circle_y,color=color,label="")
#Plots.plot!(shape::OrthoRect,color)=plot!([shape.tl[1],shape.br[1],shape.br[1],shape.tl[1],shape.tl[1]],[shape.tl[2],shape.tl[2],shape.br[2],shape.br[2],shape.tl[2]],color=color,label="")
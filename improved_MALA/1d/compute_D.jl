using Base.Threads

V(q::Float64)=sin(2π*q)
grad_V(q::Float64)=2π*cos(2π*q)

η=1e-3
dq=1e-5
q_range=0.0:dq:1.0
n=length(q_range)
ψ=zero(q_range)


for (i,q) in enumerate(q_range)
    I_threads=[0.0 for i=1:nthreads()]
    @threads for i=1:n-1
        ix=threadid()
        y1=q_range[i]
        y2=q_range[i+1]
        I_threads[ix]+=0.5*(exp(V(q+y1)-η*y1)+exp(V(q+y2)-η*y2))*dq
    end
    ψ[i]=exp(-V(q))*sum(I_threads)

    (i%100 ==0) && (println("$i:$n"))
end

Z=0.0
for (ψ1,ψ2) in zip(ψ,ψ[2:end])
    global Z+=0.5*(ψ1+ψ2)*dq
end

ψ/=Z

D=1.0
for (ψ1,ψ2,q1,q2) in zip(ψ,ψ[2:end],q_range,q_range[2:end])
    global D-=0.5*dq*(ψ1*grad_V(q1)+ψ2*grad_V(q2))/η
end
println(D)
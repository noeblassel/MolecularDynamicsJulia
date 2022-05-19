using LinearAlgebra

potential(q::Vector{Float64})=sin.(2π*q)
grad_potential(q::Vector{Float64})=2π*cos.(2π*q)

function simMALA!(q::Vector{Float64},dt::Float64,n_steps::Int64,rule::String,hist::Matrix{Float64},C::Matrix{Float64})
    M=length(q)
    sd_coords=zero(q)
    msds=zeros(n_steps)
    n_accepted=0
    V=potential(q)
    F=-grad_potential(q)
    σ=sqrt(2dt)
    λ=σ/2

    Vtilde=zero(V)
    Ftilde=zero(F)
    disp=zero(q)

    for i=1:n_steps
        msds[i]=sum(sd_coords .^ 2)/M

        hist[2:end,:]=hist[1:end-1,:]
        hist[1,:]=F
        C.= C+hist .* transpose(hist[1,:])

        G=randn(M)
        qtilde=q+dt*F+σ*G
        disp=qtilde-q
        qtilde.=mod1.(qtilde,(1.0))
        Vtilde.=potential(qtilde)
        Ftilde=-grad_potential(qtilde)
        α=(Vtilde-V)+((G+λ*(F+Ftilde)).^2)/2+(G .^2)/2
        U=rand(M)
        
        if rule=="barker"
            accepted=trues(M)
            @. accepted[α < 0] = ( log(U) < -log( 1 + exp(α)))[α < 0]
            @. accepted[α >= 0] = (log(U) < -α - log(1 + exp(-α)))[α >= 0]
        elseif rule=="metropolis"
            accepted= (log.(U).< -α)
        end

        q[accepted] .= qtilde[accepted]
        V[accepted] .= Vtilde[accepted]
        F[accepted] .= Ftilde[accepted]
        n_accepted+=sum(accepted)

        sd_coords[accepted].+=disp[accepted]
    end
    ts=dt*(0:(n_steps-1))
    Dhat=inv(dot(ts,ts))*dot(ts,msds)/2

    return (Dhat,n_accepted/(n_steps*M))
end



q=zeros(100)
hist=zeros(1000,100)
C=zeros(1000,100)
@time (Dhat,A)=simMALA!(q,0.01,10000,"metropolis",hist,C)

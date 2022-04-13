function update_hist!(bins::Vector{Int64},pt::Float64,lims::Tuple{Float64,Float64})
    N=length(bins)
    (lm,lM)=lims
    
    (pt > lM) && return nothing
    (pt < lm) && return nothing

    bin_ix=Int64(ceil(N*(pt-lm)/(lM-lm)))
    bins[bin_ix]+=1
end


function update_hist2D!(bins::Array{Int64,2},q::Float64,p::Float64,q_lims::Tuple{Float64,Float64},p_lims::Tuple{Float64,Float64})
    (Nq,Np)=size(bins)
    (mq,Mq)=q_lims
    (mp,Mp)=p_lims

    ix_q=Int64(ceil(Nq*(q-mq)/(Mq-mq)))
    ix_p=Int64(ceil(Np*(p-mp)/(Mp-mp)))
    
    ((ix_q<1)||(ix_q>Nq)||(ix_p<1)||(ix_p>Np)) && return nothing

    bins[ix_q,ix_p]+=1
end

function write_hist2D(hist::Array{Int64,2},file::String)
    f=open(file,"w")
    for r in eachrow(hist)
        println(f,join(r,' '))
    end
    close(f)
end

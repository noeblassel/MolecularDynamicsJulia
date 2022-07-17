using Statistics

function asymptotic_var(v)
    data=copy(v)
    avg=mean(v)
    data .-= avg
    N_steps=floor(Int64,log2(length(data)))
    data=data[1:2^N_steps]#crop series to largest possible power of two
    L=1
    N=length(data)
    K=N
    max_var=0.0
    while K>1000
        max_var=max(max_var,varm(data,0.0)*N*inv(K))
        K >>=1
            data=(data[1:2:end]+data[2:2:end])/2
        end
    return max_var
end

run(`./scp_sv_results.sh`)

etas=collect(vcat(collect(0.008:0.008:0.072),collect(0.08:0.08:0.8)))
normalizing_csts=Dict("SINUSOIDAL"=>1/2,"CONSTANT"=>2/π,"LINEAR"=>-4/π^2)

methods=["SINUSOIDAL","CONSTANT","LINEAR"]
for method in methods
    output_file="norton_$method.dat"
    f_output=open(output_file,"w")
    norm=normalizing_csts[method]
    println(f_output,"eta norton_response est_normalized_response N_samps  AV_response, AV_normalized_response")
    println(method)
    for eta in etas 
        println("\t",eta)
        f=open("norton_forcing_$(method)_$(eta).out")
        C=reinterpret(Float64,read(f))
        N=length(C)
        m=mean(C)
        σ2=asymptotic_var(C)
        println(f_output,join([eta,m,inv(m*norm),N,σ2,σ2/(norm^2*m^4)]," ")) #get other variance by delta method (see comment below)
        close(f)
        rm("norton_forcing_$(method)_$(eta).out")
    end
    close(f_output)
end

#φ(m)=1/(norm*m) => φ'(m)=-1/(norm * m^2)
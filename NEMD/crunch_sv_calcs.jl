using Statistics

function asymptotic_var(v::Vector{Float64})
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

#run(`./scp_sv_results.sh`)
etas=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

output_file="sv_sinusoidal.dat"
f_output=open(output_file,"w")
println(f_output,"eta Im(response) Re(response) N_samps AV_im_response AV_re_response")


for eta in etas 
    println(eta)
    f=open("fourier_response_SINUSOIDAL_$(eta).out")
    C=read(f)
    C=reinterpret(ComplexF64,C)
    N=length(C)
    m=mean(C)
    v_im=asymptotic_var(imag(C))
    v_re=asymptotic_var(real(C))
    println(f_output,join([eta,imag(m),real(m),N,v_im,v_re]," "))
    close(f)
    rm("fourier_response_SINUSOIDAL_$(eta).out")
end

close(f_output)
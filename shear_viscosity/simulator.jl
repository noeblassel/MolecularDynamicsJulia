using Molly

mutable struct NortonShearViscosityTest{rngType}
    v::Float64
    dt::Float64
    β::Float64
    γ::Float64
    rseed::Uint32
    rng::rngType

end

function NortonShearViscosityTest(;v::Float64,dt::Float64,T::Float64,γ::Float64,rseed::UInt32=floor(UInt32,time()))
    return NortonShearViscosityTest{typeof(rng)}(v,inv(T),γ,rseed,rng)
end

function simulate_norton!(sys::System{D},sim::NortonShearViscositySimulator,n_steps::Int64;parallel::Bool=true) where {D}
    α= exp(-sim.γ*sim.dt)*ones(length(sys))
    σ=zeros(α)
    @. σ=sqrt((1-α^2)/β)
    neighbors=find_neighbors(sys,sys.neighbor_finder;parallel=parallel)
    accels=accelerations(sys,neighbors;parallel=parallel)
    G=zeros(sys.velocities)

    ### views into state vector
    accels_array=reinterpret(reshape,Float64,accels)
    coords_array=reinterpret(reshape,Float64,sys.coords)
    velocities_array=reinterpret(reshape,Float64,sys.velocities)

    y_accels=view(reshape,)


    Ly=sys.box_size[2]
    F_q=

    for step_n=1:n_steps
        run_loggers!(sys,neighbors,step_n)
        #B step
        #A step
        #B step
        #O step
    end
end
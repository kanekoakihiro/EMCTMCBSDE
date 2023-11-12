using Random, Distributions, LinearAlgebra
@doc """
Arguments:
    - BSDE : struct defined in ./equation.jl
    - basis : basis functions used for regression
    - M : number of samples
    - N : number of temporal steps
"""
function LSMC1D(BSDE, basis::Function, M::Int, N::Int)
    Δt = BSDE.T/N; tj = range(0.0, step =Δt, length = N+1)
    ΔWj = zeros(M, N+1); Wj = zeros(M, N+1)
    @inbounds for j = 1:M
        ΔWj[j, :] .= [0; sqrt(Δt) .* randn(N)]
        Wj[j, :] .= cumsum(ΔWj[j, :])
    end

    Sj = zeros(M, N+1)
    Sj[:,1] .= BSDE.X0
    @inbounds for j = 2:N+1
        Sj[:, j] = Sj[:,j-1] .+ BSDE.drift.(Sj[:,j-1]).*Δt .+ BSDE.diffusion.(Sj[:,j-1]).*ΔWj[:,j]
    end

    Zj = zeros(M,N+1); Yj = zeros(M,N+1)
    Yj[1:end,end] = BSDE.terminal.(Sj[1:end,end])
    reg_coeff(X,y)=X\y
    @inbounds for n in N:-1:1
        A = reduce(hcat,basis.(Sj[1:end,n]))'
        b = Yj[1:end,n+1] .* ΔWj[1:end,n+1]
        α = reg_coeff(A,b)
        Zj[1:end,n] = A * α ./ Δt

        c = Yj[1:end,n+1] .+ Δt .* BSDE.driver.(tj[n], Sj[1:end,n], Yj[1:end,n+1], Zj[1:end,n])
        β = reg_coeff(A,c)
        Yj[1:end,n] = A * β
    end;

    return mean(Yj[1:end,2] .+ Δt .* BSDE.driver.(tj[2], Sj[1:end,1], Yj[1:end,2], Zj[1:end,2]))
end

function LSMC(BSDE, basis::Function, M::Int, N::Int)
    Δt = BSDE.T/N; tj = range(0.0, step =Δt, length = N+1)
    d = length(BSDE.X0)
    ΔWj = zeros(M, N+1, d); Wj = zeros(M, N+1, d)

    ΔWj[:, 2:end, :] .= sqrt(Δt) .* randn(M, N, d)
    Wj .= cumsum(ΔWj, dims=2)
    
    Sj = zeros(M, N+1, d)
    Sj[:,1,:] .= reshape(BSDE.X0, (1,d))
    @inbounds for j = 2:N+1
        @views Sj[:,j,:] = (
            Sj[:,j-1,:] .+ reduce(hcat, BSDE.drift.(eachslice(Sj[:,j-1,:],dims=1)) .*Δt)'
            .+ reduce(hcat, BSDE.diffusion.(eachslice(Sj[:,j-1,:],dims=1)) .* eachslice(ΔWj[:,j,:],dims=1))'
        )
    end

    Zj = zeros(M,N+1,d); Yj = zeros(M,N+1)
    @views Yj[:,end] = BSDE.terminal.(eachslice(Sj[:,end,:], dims=1))
    reg_coeff(X,y)=X\y
    @inbounds for n in N:-1:1
        @views A = reduce(hcat,basis.(eachslice(Sj[:,n,:],dims=1)))'
        @inbounds for p in 1:d
            @views b = Yj[:,n+1] .* ΔWj[:,n+1,p]
            α = reg_coeff(A,b)
            Zj[:,n,p] = A * α ./ Δt
        end
        @views c = Yj[:,n+1] .+ Δt .* BSDE.driver.(tj[n], eachslice(Sj[:,n,:],dims=1), Yj[:,n+1], eachslice(Zj[:,n,:],dims=1))
        β = reg_coeff(A,c)
        Yj[:,n] = A * β
    end;
    return @views mean(Yj[1:end,2] .+ Δt .* BSDE.driver.(tj[2], eachslice(Sj[:,1,:],dims=1), Yj[:,2], eachslice(Zj[:,2,:],dims=1)))
end
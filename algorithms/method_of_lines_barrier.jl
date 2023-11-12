using DifferentialEquations, ExponentialUtilities

function MethodOfLinesBarrier(bsde::AbstractBSDE, grid::Any, Nₜ::Int, scheme::OrdinaryDiffEqAlgorithm, EXPINT::Bool)
    attrs = generateDiffAttributes(grid, bsde.drift, bsde.diffusion, bsde.driver, bsde.terminal1)
    typeof(grid) == Array{AbstractGrid, 1} ? states = eachslice(attrs.grid,dims=1) : states = attrs.grid

    Ξ = bsde.domain_cond.(states)
    Ξᶜ = .~(Ξ)
    Ξ = Diagonal(float.(Ξ))
    Ξᶜ = Diagonal(float.(Ξᶜ))

    # Ξ = zeros(attrs.N)
    # Ξ[bsde.domain_cond.(states)] .= 1.0
    # Ξᶜ = 1.0 .- Ξ
    # Ξ = Diagonal(Ξ)
    # Ξᶜ = Diagonal(Ξᶜ)

    compensat_func(t::Float64) = Ξᶜ*bsde.terminal2.(bsde.T-t, states)
    Nonlinear(v, p, t) = Ξ*attrs.Nonlinear(v, p, bsde.T-t)
    
    dt = (bsde.T-0.0) / Nₜ
    sol = zeros(attrs.N, Nₜ+1)
    sol[:,1] = Ξ*attrs.Terminal + Ξᶜ*bsde.terminal2.(bsde.T, states)
    for tₙ in 1:Nₜ
        δ = compensat_func(dt*tₙ) - compensat_func(dt*(tₙ-1))
        iszero(δ) ? compensat = δ : compensat = expv(dt*tₙ, Ξ*attrs.Q, δ; m=100)
        if EXPINT
            Linear = DiffEqArrayOperator(Ξ*attrs.Q);
            prob = SplitODEProblem(Linear, Nonlinear, sol[:,tₙ], (0.0,dt));
        else
            _Nonlinear(v, p, t) = Nonlinear(v, p, t) + Ξ*attrs.Q * v
            prob = ODEProblem(_Nonlinear, sol[:,tₙ], (0.0,dt))
        end
        expint = solve(prob, scheme, adaptive = false, dt = dt).u[2]
        sol[:,tₙ+1] = expint + compensat
    end
    return [sol, attrs.grid, attrs.Q]
end
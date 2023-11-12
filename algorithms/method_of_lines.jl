using DifferentialEquations
function MethodOfLines(bsde::AbstractBSDE, grid::Any, Nₜ::Int, scheme::OrdinaryDiffEqAlgorithm, EXPINT::Bool)
    attrs = generateDiffAttributes(grid, bsde.drift, bsde.diffusion, bsde.driver, bsde.terminal)
    dt = (bsde.T-0.0) / Nₜ
    if EXPINT
        Linear = DiffEqArrayOperator(attrs.Q);
        prob = SplitODEProblem(Linear, attrs.Nonlinear, attrs.Terminal, (0.0, bsde.T));
    else
        _Nonlinear(v, p, t) = attrs.Nonlinear(v, p, t) + attrs.Q * v
        prob = ODEProblem(_Nonlinear, attrs.Terminal, (0.0, bsde.T))
    end
    expint = solve(prob, scheme, adaptive = false, dt = dt)
    sol = reduce(hcat,expint.u)
    return [sol, attrs.grid, attrs.Q]
end
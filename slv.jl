module slv

export SABRSLV, HestonSLV, HestonSABRSLV, RootQuadraticSLV, RootQuadraticSLV, TanHypSLV, HypHypSLV, SLV, αHyperSV, FourOverTwoSV, SteinSteinSV

@doc raw"""
```
dS_t = \omega(S_t,v_t)dt + m(v_t)\Gamma(S_t)dW_t^{(1)},
dv_t = \mu(v_t)dt + \sigma(v_t)dW_t^{(2)},
\langle W^{(1)}, W^{(2)}\rangle_t = \rho t
 with \rho\in(-1,1).
```
"""
mutable struct SLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function SLV(S0, v0, ρ, ω, m, Γ, μ, σ)
        new(S0, v0, ρ, ω, m, Γ, μ, σ)
    end
end

struct SABRSLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function SABRSLV(S0::Float64, v0::Float64, α::Float64, β::Float64, ρ::Float64, r::Float64)
        #### coefficients
        ω(s::Float64,v::Float64) = 0.0
        Γ(s::Float64) = s^β
        m(v::Float64) = v
        μ(v::Float64) = 0.0
        σ(v::Float64) = α*v
        new(S0, v0, ρ, ω, m, Γ, μ, σ)
    end
end

# struct JacobiSLV
#     ## initial values
#     S0::Float64
#     v0::Float64
#     ## corralation coefficient
#     ρ::Float64
#     ## coefficients of the price process
#     ω::Function
#     m::Function
#     Γ::Function
#     ## coefficients of the volatility process
#     μ::Function
#     σ::Function
#     function HestonSLV(S0::Float64, v0::Float64, ρ::Float64, r::Float64, η::Float64, α::Float64, θ₀::Float64)
#         Q(v) = (v-vₘᵢₙ)*(v-vₘₐₓ)/(sqrt(vₘₐₓ)-sqrt(vₘᵢₙ))^2
#         #### coefficients
#         ω(s::Float64,v::Float64) = (r-v/2)
#         Γ(s::Float64) = 
#         sv = HestonSV(v0, η, θ₀, α)
#         # m(v::Float64) = sqrt(abs(v))
#         # μ(v::Float64) = η*(θ₀-abs(v))
#         # σ(v::Float64) = α*sqrt(abs(v))
#         new(S0, v0, ρ, ω, sv.m, Γ, sv.μ, sv.σ)
#     end
# end

struct HestonSLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function HestonSLV(S0::Float64, v0::Float64, ρ::Float64, r::Float64, η::Float64, α::Float64, θ₀::Float64)
        #### coefficients
        ω(s::Float64,v::Float64) = r*s
        Γ(s::Float64) = s
        sv = HestonSV(v0, η, θ₀, α)
        # m(v::Float64) = sqrt(abs(v))
        # μ(v::Float64) = η*(θ₀-abs(v))
        # σ(v::Float64) = α*sqrt(abs(v))
        new(S0, v0, ρ, ω, sv.m, Γ, sv.μ, sv.σ)
    end
end

struct HestonSABRSLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function HestonSABRSLV(S0::Float64, v0::Float64, ρ::Float64, r::Float64, β::Float64, η::Float64, α::Float64, θ₀::Float64)
        #### coefficients
        ω(s::Float64,v::Float64) = r*s
        Γ(s::Float64) = s^β

        sv = HestonSV(v0, η, θ₀, α)
        # m(v::Float64) = sqrt(abs(v))
        # μ(v::Float64) = η*(θ₀-abs(v))
        # σ(v::Float64) = α*sqrt(abs(v))
        new(S0, v0, ρ, ω, sv.m, Γ, sv.μ, sv.σ)
    end
end


struct RootQuadraticSLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function RootQuadraticSLV(S0::Float64, v0::Float64, ρ::Float64, r::Float64, a::Float64, b::Float64, c::Float64, β::Float64, sv)
        #### coefficients
        ω(s::Float64,v::Float64) = r*s
        Γ(s::Float64) = sqrt(a*s^2+b*s+c)
        new(S0, v0, ρ, ω, sv.m, Γ, sv.μ, sv.σ)
    end
end


struct TanHypSLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function TanHypSLV(S0::Float64, v0::Float64, ρ::Float64, r::Float64, β::Float64, sv)
        #### coefficients
        ω(s::Float64,v::Float64) = r*s
        Γ(s::Float64) = tanh(β*s)
        new(S0, v0, ρ, ω, sv.m, Γ, sv.μ, sv.σ)
    end
end

struct HypHypSLV
    ## initial values
    S0::Float64
    v0::Float64
    ## corralation coefficient
    ρ::Float64
    ## coefficients of the price process
    ω::Function
    m::Function
    Γ::Function
    ## coefficients of the volatility process
    μ::Function
    σ::Function
    function HypHypSLV(S0::Float64, v0::Float64, ρ::Float64, b::Float64, β::Float64, σ₀::Float64, η::Float64, α::Float64)
        #### coefficients
        ω(s::Float64,v::Float64) = b*s
        Γ(s::Float64) = σ₀*((1-β+β^2)*s + (β-1)*(sqrt(s^2+β^2*((1-s)^2))-β))/β
        m(v::Float64) = v + sqrt(v^2+1)
        μ(v::Float64) = -η*v
        σ(v::Float64) = α*sqrt(2*η)
        new(S0, v0, ρ, ω, m, Γ, μ, σ)
    end
end

struct SV
    v0::Float64
    m::Function
    μ::Function
    σ::Function
    function SV(v0::Float64, m::Function, μ::Function, σ::Function)
        new(v0, m, μ, σ)
    end
end


struct HestonSV
    v0::Float64
    m::Function
    μ::Function
    σ::Function
    function HestonSV(v0::Float64, η::Float64, θ₀::Float64, α::Float64)
        m(v) = sqrt(v)
        μ(v) = η*(θ₀-v)
        σ(v) = α*sqrt(v)
        new(v0, m, μ, σ)
    end
end

struct FourOverTwoSV
    v0::Float64
    m::Function
    μ::Function
    σ::Function
    function FourOverTwoSV(v0::Float64, a::Float64, b::Float64, η::Float64, θ₀::Float64, σᵥ::Float64, α::Float64)
        m(v) = a*sqrt(v) + b/sqrt(v)
        μ(v) = η*(θ₀-v)
        σ(v) = σᵥ*sqrt(v)
        new(v0, m, μ, σ)
    end
end

struct SteinSteinSV
    v0::Float64
    m::Function
    μ::Function
    σ::Function
    function SteinSteinSV(v0::Float64, η::Float64, θ₀::Float64, α::Float64)
        m(v) = v
        μ(v) = η*(θ₀-v)
        σ(v) = α
        new(v0, m, μ, σ)
    end
end

struct αHyperSV
    v0::Float64
    m::Function
    μ::Function
    σ::Function
    function αHyperSV(v0::Float64, η::Float64, θ::Float64, aᵥ::Float64, α::Float64)
        m(v) = exp(v)
        μ(v) = η-θ*exp(aᵥ*v)
        σ(v) = α
        new(v0, m, μ, σ)
    end
end

end
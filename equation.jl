module equation
export AbstractBSDE, BSDE, BSDEwithBddStoppingTime

abstract type AbstractBSDE end

@doc raw"""
```
Description:
    Represents Markov BSDEs such as
        X_t = X_0 + \int_0^t\mu(X_s)ds + \int_0^t\sigma(X_s)dW_s
        Y_t = g(X_T) + \int_t^Tf(s,X_s,Y_s,Z_s)ds - \int_t^TZ_sdW_s
    for t∈[0,r].

Attributes:
    - T : terminal time
    - X0 : initial value X_0 of X_t
    - drift : drift coefficient \mu(x)
    - diffusion : diffusion coefficient \sigma(x)
    - driver : driver f(t,x,y,z)
    - terminal : terminal condition g(x) of Y_t
```
"""
mutable struct BSDE <: AbstractBSDE
    T
    X0
    drift
    diffusion
    driver
    terminal
end

mutable struct BSDEwithBddStoppingTime <: AbstractBSDE
    T
    domain_cond
    X0
    drift
    diffusion
    driver
    terminal1
    terminal2
end
end
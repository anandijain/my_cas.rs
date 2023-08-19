using ModelingToolkit, Graphs, DiffEqBase, Test, UnPack

# https://github.com/SciML/ModelingToolkit.jl/blob/master/test/structural_transformation/index_reduction.jl
# seems like ss handles ode_order_loweringq, but i don't see where it gets called, so it might not 


# Define some variables
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t) xˍt(t) yˍt(t) xˍˍt(t) yˍˍt(t)
D = Differential(t)

eqs2 = [D(D(x)) ~ T * x,
    D(D(y)) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs2, t, [x, y, T], [L, g], name = :pendulum)
lowered_sys = ModelingToolkit.ode_order_lowering(pendulum2)

lowered_eqs = [
    D(xˍt) ~ T * x,
    D(yˍt) ~ T * y - g,
    D(x) ~ xˍt,
    D(y) ~ yˍt,
    0 ~ x^2 + y^2 - L^2]
@test ODESystem(lowered_eqs, t, [xˍt, yˍt, x, y, T], [L, g], name = :pendulum) ==
      lowered_sys
@test isequal(equations(lowered_sys), lowered_eqs)


using DifferentialEquations, Symbolics, ModelingToolkit, Sundials, Plots,  DataFrames, CSV

# Define the ODE system for D(x) ~ -p*x
function simple_system!(du, u, p, t)
    du[1] = -p[1] * u[1]
end

# Initial conditions and parameters
u0 = [1.0]
p = [2.0]
tspan = (0.0, 10.0)

# Solve the system
prob = ODEProblem(simple_system!, u0, tspan, p)
sol = solve(prob, CVODE_BDF(), reltol=1e-8, abstol=1e-8, adaptive=false, saveat=LinRange(0, 10, 100))

# Plot the result
# plot(sol, xlabel="Time", ylabel="x", title="Simple ODE system simulation")
df = DataFrame(sol)
CSV.write("simple.csv", df)
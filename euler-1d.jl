using Plots
using Printf
using Profile
using BenchmarkTools
using StaticArrays

include("fluxes.jl")
include("RiemannSolver.jl")


struct EOS
    gamma::Float64

    function EOS(_gamma)
        new(_gamma)
    end
end


function getp(eos::EOS, rho::Float64, e::Float64)::Float64 
    (eos.gamma-1)*rho*e
end

function gete(eos::EOS, rho::Float64, p::Float64)::Float64
    p/(eos.gamma-1)/rho
end

function getc(eos::EOS, rho::Float64, p::Float64)::Float64
    sqrt(eos.gamma*p/rho)
end

function calcdt(eos::EOS, dx::Float64, CFL::Float64, U::Vector{SVector{3, Float64}})
    to_primitive = x::SVector{3, Float64} -> SA[x[1], x[2]/x[1], getp(eos, x[1], x[3]/x[1]-x[2]*x[2]/x[1]/x[1]/2.0)]
    local_dt = x::SVector{3, Float64} -> CFL*dx/(abs(x[2]) + getc(eos, x[1], x[3]))

    dt_opt = U[1] |> to_primitive |> local_dt

    for vect in U
        dt_new = vect |> to_primitive |> local_dt
        if dt_new < dt_opt
            dt_opt = dt_new
        end
    end
    dt_opt
end



function calc_hllc_flux(eos::EOS, Ul::SVector{3, Float64}, Ur::SVector{3, Float64})::SVector{3, Float64}
    rhol = Ul[1]; ul = Ul[2]/Ul[1]; El = Ul[3]/Ul[1]
    el = El - 0.5*ul*ul; pl = getp(eos, rhol, el); cl = getc(eos, rhol, pl); hl = El + pl/rhol
    Fl = SA[rhol*ul, pl + rhol*ul*ul, ul*(pl + rhol*El)]
    rhor = Ur[1]; ur = Ur[2]/Ur[1]; Er = Ur[3]/Ur[1]
    er = Er - 0.5*ur*ur; pr = getp(eos, rhor, er); cr = getc(eos, rhor, pr); hr = Er + pr/rhor    
    Fr = SA[rhor*ur, pr + rhor*ur*ur, ur*(pr + rhor*Er)] 
    
    p_pvrs = .5*(pl + pr) - 0.5*(ur - ul) * 0.5*(rhol + rhor) * 0.5*(cl + cr)
    p_star = max(0., p_pvrs)

    rho_av = sqrt(rhol*rhor)
    u_av = (ul*sqrt(rhol) + ur*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhor))
    h_av = (hl*sqrt(rhol) + hr*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhor))
    E_av = (El*sqrt(rhol) + Er*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhor))
    p_av = (h_av - E_av) * rho_av;
	c_av = getc(eos, rho_av, p_av);
    
    sl = u_av - c_av; sr = u_av + c_av
    s_star = (pr - pl + rhol*ul*(sl - ul) - rhor*ur*(sr - ur))/(rhol*(sl - ul) - rhor*(sr - ur))
    D = SA[0.0, 1.0, s_star]
    Ul_star = (sl*Ul - Fl + p_star*D) / (sl - s_star)
    Ur_star = (sr*Ur - Fr + p_star*D) / (sr - s_star)
    Fl_star = Fl + sl*(Ul_star - Ul)
    Fr_star = Fr + sr*(Ur_star - Ur)
 
    Fhllc = zeros(SVector{3})
    if 0. <= sl  
        Fhllc = Fl
    elseif sl <= 0.0 <= s_star
        Fhllc = Fl_star
    elseif s_star <= 0.0 <= sr
        Fhllc = Fr_star
    else
        Fhllc = Fr
    end       
end


function main()
        
    println("fsLA: simple laser hydrodynamic code v0.1, (c) Vadim V. Shepelev, Ph.D.")

    eos = EOS(1.4)
    # i.c.
    N = 1000
    x_min = 0.0 
    x_max = 1.0
    x_0 = 0.3
    x = collect(range(x_min, x_max, N+1))
    U = Vector[]

    (t_min, t_max) = 0.0, 0.2
    global t = t_min
    global dx = (x_max - x_min)/N

    global CFL = 0.9

    rho_l = 1.0
    u_l = 0.75
    p_l = 1.0

    rho_r = 0.125
    u_r = 0.0
    p_r = 0.1

    U = [x[i] < x_0 ? 
        SA[rho_l, rho_l*u_l, rho_l*(gete(eos, rho_l, p_l) + u_l*u_l/2.0)] : 
        SA[rho_r, rho_r*u_r, rho_r*(gete(eos, rho_r, p_r) + u_r*u_r/2.0)] for i in 1:N ]
    pushfirst!(U, U[begin])
    pushfirst!(U, U[begin])
    push!(U, U[end])
    push!(U, U[end])
    U_new = []
    copy!(U_new, U)
    F = [zeros(SVector{3}) for _ = 1:N+2+2]

    println()
    println("size(x) = ", size(x))
    println("size(U) = ", size(U))
    println()

    println("Starting simulation...")

    global counter = 1
    while t < t_max
        tau = calcdt(eos, dx, CFL, U)
        if counter <= 5 
            tau *= 0.2
        end
 
        # @printf("iter=%d t=%.3g dt=%.3g CFL=%g\n", counter, t, tau, CFL)
        @inbounds for i = 2:N+2+1
            F[i] = calc_hllc_flux(eos, U[i-1], U[i])
        end
        @inbounds for i = 1+2:N+2
            U_new[i] = U[i] - tau/dx*(F[i+1] - F[i])        
        end
        copy!(U, U_new)

        t += tau
        counter += 1
    end


    println("Exact solution:")
    y = sample_riemann(x .- 0.3, 0.2)   # Analytical solution, vector of HydroStatusSimple structures
    rho_ex = [elem.rho for elem in y]
end

@time main()

# Profile.clear()
# @profile main()
# Profile.print()
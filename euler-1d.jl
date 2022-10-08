# include("fluxes.jl")

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


function calc_hllc_flux(eos::EOS, Ul::Vector{Float64}, Ur::Vector{Float64})::Vector{Float64}
    (rhol, ul, El) = (Ul[1], Ul[2]/Ul[1], Ul[3]/Ul[1])
    el = El - 0.5*ul*ul
    pl = getp(eos, rhol, el)
    cl = getc(eos, rhol, pl) 
    hl = El + pl/rhol
    Fl = [rhol*ul, pl + rhol*ul*ul, ul*(pl + rhol*El)]
    (rhor, ur, Er) = (Ur[1], Ur[2]/Ur[1], Ur[3]/Ur[1])
    er = Er - 0.5*ur*ur
    pr = getp(eos, rhor, er)
    cr = getc(eos, rhor, pr) 
    hr = Er + pr/rhor    
    Fr = [rhor*ur, pr + rhor*ur*ur, ur*(pr + rhor*Er)] 
    
    p_pvrs = .5*(pl + pr) - 0.5*(ur - ul) * 0.5*(rhol + rhor) * 0.5*(cl + cr)
    p_star = max(0., p_pvrs)

    rho_av = sqrt(rhol*rhor)
    u_av = (ul*sqrt(rhol) + ur*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhol))
    h_av = (hl*sqrt(rhol) + hr*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhol))
    E_av = (El*sqrt(rhol) + Er*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhol))
    p_av = (pl*sqrt(rhol) + pr*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhol))
    c_av = (cl*sqrt(rhol) + cr*sqrt(rhor)) / (sqrt(rhol)+sqrt(rhol))
    
    (sl, sr) = (u_av - c_av, u_av + c_av)
    s_star = (pr - pl + rhol*ul*(sl - ul) - rhor*ur*(sr - ur))/(rhol*(sl - ul) - rhor*(sr - ur))
    D = [0.0, 1.0, s_star]
    Ul_star = (sl*Ul - Fl + p_star*D) / (sl - s_star)
    Ur_star = (sr*Ur - Fr + p_star*D) / (sr - s_star)
    Fl_star = Fl + sl*(Ul_star - Ul)
    Fr_star = Fr + sr*(Ur_star - Ur)
 
    Fhllc = zeros(3)
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

println("fsLA: simple laser hydrodynamic code v0.1, (c) Vadim V. Shepelev, Ph.D.")

eos = EOS(1.4
)
# i.c.
N = 5
x_min = 0 
x_max = 1
x_0 = 0.3
x = collect(range(x_min, x_max, N+1))
U = collect(zeros(3, N+2+2))

dt = 0.01
dx = (x_max - x_min)/N

rho_l = 1.0
u_l = 0.75
p_l = 1.0

rho_r = 1.0
u_r = 0.0
p_r = 0.1

#=for i in 1:N
    if x[i] < x_0 
        U[:, i] = [rho_l, u_l, p_l]
    else 
        U[:, i] = [rho_r, u_r, p_r]
    end
=#
U = [x[i] < x_0 ? [rho_l, u_l, p_l] : [rho_r, u_r, p_r]  for i in 1:N ]
pushfirst!(U, U[begin])
pushfirst!(U, U[begin])
push!(U, U[end])
push!(U, U[end])
U_new = []
copy!(U_new, U)
F = [zeros(3) for _ = 1:N+2+2]

println("x = ", x)
println()
println()

println("U = ",  U)
println()
println()

println("Starting calc_hllc_flux()...")

for i = 2:N+1
    F[i] = calc_hllc_flux(eos, U[i-1], U[i])
end

println("F = ",  F, " ", typeof(F), " ", size(F))
println()
println()

for i = 1+2:N+2
    U_new[i] = U[i] - dt/dx*(F[i+1] - F[i])
end


println("U_new = ",  U_new)
println()
println()

copy!(U, U_new)

println("U = ",  U)
println()
println()
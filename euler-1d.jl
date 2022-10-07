# include("fluxes.jl")


function hllc(nx,gamma,uL,uR,f,fL,fR)
    gm = gamma-1.0
    Ds = Array{Float64}(undef,3)
    Ds[1], Ds[2] = 0.0, 1.0
    
    for i = 1:nx+1
        # left state
        rhLL = uL[i,1]
        uuLL = uL[i,2]/rhLL
        eeLL = uL[i,3]/rhLL
        ppLL = gm*(eeLL*rhLL - 0.5*rhLL*(uuLL*uuLL))
        aaLL = sqrt(abs(gamma*ppLL/rhLL))
    
        # right state
        rhRR = uR[i,1]
        uuRR = uR[i,2]/rhRR
        eeRR = uR[i,3]/rhRR
        ppRR = gm*(eeRR*rhRR - 0.5*rhRR*(uuRR*uuRR))
        aaRR = sqrt(abs(gamma*ppRR/rhRR))
    
        # compute SL and Sr
        SL = min(uuLL,uuRR) - max(aaLL,aaRR)
        SR = max(uuLL,uuRR) + max(aaLL,aaRR)
    
        # compute compound speed
        SP = (ppRR - ppLL + rhLL*uuLL*(SL-uuLL) - rhRR*uuRR*(SR-uuRR))/
             (rhLL*(SL-uuLL) - rhRR*(SR-uuRR)) #never get zero
    
        # compute compound pressure
        PLR = 0.5*(ppLL + ppRR + rhLL*(SL-uuLL)*(SP-uuLL)+
                   rhRR*(SR-uuRR)*(SP-uuRR))
    
        # compute D
        Ds[3] = SP
    
        if (SL >= 0.0)
            for m = 1:3
                f[i,m] = fL[i,m]
            end
        elseif (SR <= 0.0)
            for m =1:3
                f[i,m] = fR[i,m]
            end
        elseif ((SP >=0.0) & (SL <= 0.0))
            for m = 1:3
                f[i,m] = (SP*(SL*uL[i,m]-fL[i,m]) + SL*PLR*Ds[m])/(SL-SP)
            end
        elseif ((SP <= 0.0) & (SR >= 0.0))
            for m = 1:3
                f[i,m] = (SP*(SR*uR[i,m]-fR[i,m]) + SR*PLR*Ds[m])/(SR-SP)
            end
        end
    end
    end

println("Simple laser hydrodynamic code in Julia v0.1, (c) Vadim V. Shepelev, Ph.D.")

# i.c.
N = 100
x_min = 0 
x_max = 1
x_0 = 0.3
x = collect(range(x_min, x_max, N+1))
U = collect(zeros(3, N))

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

println(x, U)


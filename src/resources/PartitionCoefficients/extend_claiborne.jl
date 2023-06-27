using StatGeochem, LsqFit, Plots

ree3 = ["Pr","Pr","Nd","Sm","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",] # 3+ rare earth elements
r = [0.106,0.101,0.100,0.096,0.094,0.092,0.091,0.089,0.088,0.087,0.086,0.085,] # Ionic Radii
rp = 0.085:0.001:0.107 # Radius range to plot over

# Equation we're fitting to (from Blundy and Wood, 1994):
# lnD0 + a * (r0/2*(x-r0)^2 + 1/3*(x-r0)^3)
# where a = 4π E Na / RT
@. blundy_wood(x,param) = param[1] + param[2] * (param[3]/2*(x-param[3])^2 + 1/3*(x-param[3])^3)

h = plot(framestyle=:box, xlabel="Ionic radius", ylabel="log10 kD")

T = 500:10:1000
kd_La = zeros(length(T))
kd_Pr = zeros(length(T))
rD = r[3:end]
kD = zeros(length(rD))
for j in eachindex(T)

    for i in 1:10
        kD[i] = log10(claiborne_zircon_kd(ree3[2+i], T[j]))
    end
    plot!(h, rD, kD, seriestype=:scatter, label="")

    # Fit to Blundy and Wood curve
    param = [maximum(kD), -10000, 0.095] # initial guess
    fobj = curve_fit(blundy_wood, rD, kD, param) # Fit

    plot!(h, r, blundy_wood(r,fobj.param), label="") # Plot

    kd_La[j] = 10.0^blundy_wood(r[1],fobj.param)
    kd_Pr[j] = 10.0^blundy_wood(r[2],fobj.param)
end
display(h)

@. f(T, param) = param[1]*exp(param[2]/(T+273.15))
param = [nanmean(kd_La), 5000] # Initial guess

## --- La vs temp

h1 = plot(T, kd_La, xlabel="Temperature (C)", ylabel="La kD", seriestype=:scatter, framestyle=:box, label="")
fobj = curve_fit(f, T, kd_La, param)
plot!(h1, T, f(T, fobj.param), label="")
display(h1)

@info "La zrn/melt kd = $(fobj.param[1]) * exp($(fobj.param[2])/T)"

## -- Fit kd Pr vs temp

h1 = plot(T, kd_Pr, xlabel="Temperature (C)", ylabel="Pr kD", seriestype=:scatter, framestyle=:box, label="")
fobj = curve_fit(f, T, kd_Pr, param)
plot!(h1, T, f(T, fobj.param), label="")
display(h1)

@info "Pr zrn/melt kd = $(fobj.param[1]) * exp($(fobj.param[2])/T)"

# Test creating and normalizing Composition objects

x = NCKFMASHTOtrace((1.0:length(fieldnames(NCKFMASHTOtrace)))...,)
@test x isa NCKFMASHTOtrace{Float64}
@test x == NCKFMASHTOtrace((1:40.)...,)
xn = normalize(x) 
@test sum(e->xn[e], majorelements(xn)) ≈ 99.86110228500358 
xn = normalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.83028850953377
@test propertynames(xn)[1:10] == majorelements(xn) == majorelements(NCKFMASHTOtrace)
@test propertynames(xn)[11:40] == traceelements(xn) == traceelements(NCKFMASHTOtrace)
@test xn isa NCKFMASHTOtrace{Float64}
@test 0.5*xn + 0.5*xn == xn

xl = NCKFMASHTOlogtrace(x)
@test xl isa NCKFMASHTOlogtrace{Float64}
@test xl == NCKFMASHTOlogtrace((1:10.)..., log.(11:40.)...,)
xn = normalize(xl) 
@test xn isa NCKFMASHTOlogtrace{Float64}
@test sum(e->xn[e], majorelements(xn)) ≈ 99.86110228500358 
xn = normalize(xl, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.83028850953377
@test propertynames(xn)[1:10] == majorelements(xn) == majorelements(NCKFMASHTOlogtrace)
@test propertynames(xn)[11:40] == traceelements(xn) == traceelements(NCKFMASHTOlogtrace)
@test xn isa NCKFMASHTOlogtrace{Float64}
@test 0.5*xn + 0.5*xn == xn

x = NCKFMASHTOCrtrace((1.0:length(fieldnames(NCKFMASHTOCrtrace)))...,)
@test x isa NCKFMASHTOCrtrace{Float64}
@test x == NCKFMASHTOCrtrace((1:40.)...,)
xn = normalize(x) 
@test sum(e->xn[e], majorelements(xn)) ≈ 99.8858879401411
xn = normalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:10]) ≈ 99.86309677278784
@test propertynames(xn)[1:11] == majorelements(xn) == majorelements(NCKFMASHTOCrtrace)
@test propertynames(xn)[12:40] == traceelements(xn) == traceelements(NCKFMASHTOCrtrace)
@test xn isa NCKFMASHTOCrtrace{Float64}
@test 0.5*xn + 0.5*xn == xn

xl = NCKFMASHTOCrlogtrace(x)
@test xl isa NCKFMASHTOCrlogtrace{Float64}
@test xl == NCKFMASHTOCrlogtrace((1:11.)..., log.(12:40.)...,)
xn = normalize(xl) 
@test sum(e->xn[e], majorelements(xn)) ≈ 99.8858879401411
xn = normalize(xl, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:10]) ≈ 99.86309677278784
@test propertynames(xn)[1:11] == majorelements(xn) == majorelements(NCKFMASHTOCrlogtrace)
@test propertynames(xn)[12:40] == traceelements(xn) == traceelements(NCKFMASHTOCrlogtrace)
@test xn isa NCKFMASHTOCrlogtrace{Float64}
@test 0.5*xn + 0.5*xn == xn

# Composition arrays
ca = CompositionArray{NCKFMASHTOCrtrace{Float64}}(undef, 100)
@test ca isa CompositionArray{NCKFMASHTOCrtrace{Float64}}
@test ca[1] isa NCKFMASHTOCrtrace{Float64}
@test ca[1:10] isa CompositionArray{NCKFMASHTOCrtrace{Float64}}
@test length(ca) === 100
@test eachindex(ca) === Base.OneTo(100)
@test ca.SiO2 isa Vector{Float64}
@test length(ca.SiO2) === 100
@test eachindex(ca.SiO2) === Base.OneTo(100)
@test ca.SiO2 === ca[:SiO2]
@test ca.La isa Vector{Float64}
@test length(ca.La) === 100
@test eachindex(ca.La) === Base.OneTo(100)
@test ca.La === ca[:La]

# Randomize and normalize composition arrays
StatGeochem.rand!(ca)
@test nanmean(ca.SiO2) ≈ 10 atol = 2
@test sum(e->ca[1][e], majorelements(ca)) ≈ 99.99 atol=0.02
renormalize!(ca; anhydrous=true)
@test sum(e->ca[1][e], filter(x->!(x===:H2O), majorelements(ca)))  ≈ 99.99 atol=0.02


ca = CompositionArray{NCKFMASHTOlogtrace{Float64}}(undef, 100)
StatGeochem.rand!(ca)
@test nanmean(ca.SiO2) ≈ 10 atol = 2
@test sum(e->ca[1][e], majorelements(ca)) ≈ 99.99 atol=0.02
renormalize!(ca; anhydrous=true)
@test sum(e->ca[1][e], filter(x->!(x===:H2O), majorelements(ca)))  ≈ 99.99 atol=0.02

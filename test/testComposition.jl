# Test creating and normalizing Composition objects

x = NCKFMASHTOtrace((1.0:length(fieldnames(NCKFMASHTOtrace)))...,)
@test x isa NCKFMASHTOtrace{Float64}
@test x == NCKFMASHTOtrace((1:40.)...,)
xn = normalize(x) 
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.86110228500358 
xn = normalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.83028850953377
@test propertynames(xn)[1:10] == StatGeochem.majors(xn) == StatGeochem.majors(NCKFMASHTOtrace)
@test propertynames(xn)[11:40] == StatGeochem.traces(xn) == StatGeochem.traces(NCKFMASHTOtrace)
@test xn isa NCKFMASHTOtrace{Float64}

xl = NCKFMASHTOlogtrace(x)
@test xl isa NCKFMASHTOlogtrace{Float64}
@test xl == NCKFMASHTOlogtrace((1:10.)..., log.(11:40.)...,)
xn = normalize(xl) 
@test xn isa NCKFMASHTOlogtrace{Float64}
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.86110228500358 
xn = normalize(xl, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.83028850953377
@test propertynames(xn)[1:10] == StatGeochem.majors(xn) == StatGeochem.majors(NCKFMASHTOlogtrace)
@test propertynames(xn)[11:40] == StatGeochem.traces(xn) == StatGeochem.traces(NCKFMASHTOlogtrace)
@test xn isa NCKFMASHTOlogtrace{Float64}

x = NCKFMASHTOCrtrace((1.0:length(fieldnames(NCKFMASHTOCrtrace)))...,)
@test x isa NCKFMASHTOCrtrace{Float64}
@test x == NCKFMASHTOCrtrace((1:40.)...,)
xn = normalize(x) 
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.8858879401411
xn = normalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:10]) ≈ 99.86309677278784
@test propertynames(xn)[1:11] == StatGeochem.majors(xn) == StatGeochem.majors(NCKFMASHTOCrtrace)
@test propertynames(xn)[12:40] == StatGeochem.traces(xn) == StatGeochem.traces(NCKFMASHTOCrtrace)
@test xn isa NCKFMASHTOCrtrace{Float64}

xl = NCKFMASHTOCrlogtrace(x)
@test xl isa NCKFMASHTOCrlogtrace{Float64}
@test xl == NCKFMASHTOCrlogtrace((1:11.)..., log.(12:40.)...,)
xn = normalize(xl) 
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.8858879401411
xn = normalize(xl, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:10]) ≈ 99.86309677278784
@test propertynames(xn)[1:11] == StatGeochem.majors(xn) == StatGeochem.majors(NCKFMASHTOCrlogtrace)
@test propertynames(xn)[12:40] == StatGeochem.traces(xn) == StatGeochem.traces(NCKFMASHTOCrlogtrace)
@test xn isa NCKFMASHTOCrlogtrace{Float64}

# Composition arrays
ca = CompositionArray{NCKFMASHTOCrtrace{Float64}}(undef, 100)
@test ca isa CompositionArray{NCKFMASHTOCrtrace{Float64}}
@test length(ca) == 100
@test eachindex(ca) == 1:100
@test ca[1] isa NCKFMASHTOCrtrace{Float64}
@test ca[1:10] isa CompositionArray{NCKFMASHTOCrtrace{Float64}}
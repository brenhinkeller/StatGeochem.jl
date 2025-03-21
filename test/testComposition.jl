# Test creating and normalizing Composition objects

x = NCKFMASHTOtrace((1.0:length(fieldnames(NCKFMASHTOtrace)))...,)
@test x isa NCKFMASHTOtrace{Float64}
@test x == NCKFMASHTOtrace((1:40.)...,)
xn = normalize(x) 
@test xn isa NCKFMASHTOtrace{Float64}
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.86110228500358 
@test propertynames(xn)[1:10] == StatGeochem.majors(xn)

xl = NCKFMASHTOlogtrace(x)
@test xl isa NCKFMASHTOlogtrace{Float64}
@test xl == NCKFMASHTOlogtrace((1:10.)..., log.(11:40.)...,)
xn = normalize(xl) 
@test xn isa NCKFMASHTOlogtrace{Float64}
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.86110228500358 
@test propertynames(xn)[1:10] == StatGeochem.majors(xn)

x = NCKFMASHTOCrtrace((1.0:length(fieldnames(NCKFMASHTOCrtrace)))...,)
@test x isa NCKFMASHTOCrtrace{Float64}
@test x == NCKFMASHTOCrtrace((1:40.)...,)
xn = normalize(x) 
@test xn isa NCKFMASHTOCrtrace{Float64}
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.8858879401411
@test propertynames(xn)[1:11] == StatGeochem.majors(xn)

xl = NCKFMASHTOCrlogtrace(x)
@test xl isa NCKFMASHTOCrlogtrace{Float64}
@test xl == NCKFMASHTOCrlogtrace((1:11.)..., log.(12:40.)...,)
xn = normalize(xl) 
@test xn isa NCKFMASHTOCrlogtrace{Float64}
@test sum(e->xn[e], StatGeochem.majors(xn)) ≈ 99.8858879401411
@test propertynames(xn)[1:11] == StatGeochem.majors(xn)
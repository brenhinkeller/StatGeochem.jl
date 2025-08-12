## --- Test creating and normalizing Composition objects

x = NCKFMASHTOtrace((1.0:length(fieldnames(NCKFMASHTOtrace)))...,)
@test x isa NCKFMASHTOtrace{Float64}
@test x === NCKFMASHTOtrace{Float64}(1:40)
@test haskey(x, "SiO2") && haskey(x, :SiO2)
@test x["SiO2"] ≈ x[:SiO2] ≈ 1
@test ntuple(x) === ntuple(Float64, 40)
@test NamedTuple(x) === NamedTuple{fieldnames(NCKFMASHTOtrace)}(1:40.)
@test majorelementvalues(x) === ((1:10.)...,)
@test traceelementvalues(x) === ((11:40.)...,)
xn = renormalize(x) 
@test sum(majorelementvalues(xn)) ≈ 99.86110228500358 
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.83028850953377
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.83028850953377
@test propertynames(xn)[1:10] == majorelements(xn) == majorelements(NCKFMASHTOtrace)
@test propertynames(xn)[11:40] == traceelements(xn) == traceelements(NCKFMASHTOtrace)
@test xn isa NCKFMASHTOtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASHTOlogtrace(x)
@test x isa NCKFMASHTOlogtrace{Float64}
@test x === NCKFMASHTOlogtrace{Float64}([1:10; log.(11:40)])
@test haskey(x, "SiO2") && haskey(x, :SiO2)
@test x["SiO2"] ≈ x[:SiO2] ≈ 1
@test ntuple(x) === ((1:10.)..., log.(11:40)...,)
@test NamedTuple(x) === NamedTuple{fieldnames(NCKFMASHTOlogtrace)}([1:10; log.(11:40)])
@test majorelementvalues(x) === ((1:10.)...,)
@test traceelementvalues(x) === (log.(11:40.)...,)
xn = renormalize(x) 
@test xn isa NCKFMASHTOlogtrace{Float64}
@test sum(majorelementvalues(xn)) ≈ 99.86110228500358 
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.83028850953377
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.83028850953377
@test propertynames(xn)[1:10] == majorelements(xn) == majorelements(NCKFMASHTOlogtrace)
@test propertynames(xn)[11:40] == traceelements(xn) == traceelements(NCKFMASHTOlogtrace)
@test xn isa NCKFMASHTOlogtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASHTOCrtrace((1.0:length(fieldnames(NCKFMASHTOCrtrace)))...,)
@test x isa NCKFMASHTOCrtrace{Float64}
@test x === NCKFMASHTOCrtrace{Float64}(1:40)
@test ntuple(x) === ntuple(Float64, 40)
@test majorelementvalues(x) === ((1:11.)...,)
@test traceelementvalues(x) === ((12:40.)...,)
xn = renormalize(x) 
@test sum(majorelementvalues(xn)) ≈ 99.8858879401411
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.86309677278784
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:10]) ≈ 99.86309677278784
@test propertynames(xn)[1:11] == majorelements(xn) == majorelements(NCKFMASHTOCrtrace)
@test propertynames(xn)[12:40] == traceelements(xn) == traceelements(NCKFMASHTOCrtrace)
@test xn isa NCKFMASHTOCrtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASHTOCrlogtrace(x)
@test x isa NCKFMASHTOCrlogtrace{Float64}
@test x === NCKFMASHTOCrlogtrace{Float64}([1:11; log.(12:40)])
@test ntuple(x) === ((1:11.)..., log.(12:40)...,)
@test majorelementvalues(x) === ((1:11.)...,)
@test traceelementvalues(x) === (log.(12:40.)...,)
xn = renormalize(x) 
@test sum(majorelementvalues(xn)) ≈ 99.8858879401411
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.86309677278784
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:10]) ≈ 99.86309677278784
@test propertynames(xn)[1:11] == majorelements(xn) == majorelements(NCKFMASHTOCrlogtrace)
@test propertynames(xn)[12:40] == traceelements(xn) == traceelements(NCKFMASHTOCrlogtrace)
@test xn isa NCKFMASHTOCrlogtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASTtrace((1.0:length(fieldnames(NCKFMASTtrace)))...,)
@test x isa NCKFMASTtrace{Float64}
@test x === NCKFMASTtrace{Float64}(1:38)
@test ntuple(x) === ntuple(Float64, 38)
@test majorelementvalues(x) === ((1:8.)...,)
@test traceelementvalues(x) === ((9:38.)...,)
xn = renormalize(x) 
@test sum(majorelementvalues(xn)) ≈ 99.8045494240446
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.8045494240446
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:8]) ≈ 99.8045494240446
@test propertynames(xn)[1:8] == majorelements(xn) == majorelements(NCKFMASTtrace)
@test propertynames(xn)[9:38] == traceelements(xn) == traceelements(NCKFMASTtrace)
@test xn isa NCKFMASTtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASTlogtrace(x)
@test x isa NCKFMASTlogtrace{Float64}
@test x === NCKFMASTlogtrace{Float64}([1:8; log.(9:38)])
@test ntuple(x) === ((1:8.)..., log.(9:38)...,)
@test majorelementvalues(x) === ((1:8.)...,)
@test traceelementvalues(x) === (log.(9:38.)...,)
xn = renormalize(x) 
@test xn isa NCKFMASTlogtrace{Float64}
@test sum(majorelementvalues(xn)) ≈ 99.8045494240446
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.8045494240446
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:8]) ≈ 99.8045494240446
@test propertynames(xn)[1:8] == majorelements(xn) == majorelements(NCKFMASTlogtrace)
@test propertynames(xn)[9:38] == traceelements(xn) == traceelements(NCKFMASTlogtrace)
@test xn isa NCKFMASTlogtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASTCrtrace((1.0:length(fieldnames(NCKFMASTCrtrace)))...,)
@test x isa NCKFMASTCrtrace{Float64}
@test x === NCKFMASTCrtrace{Float64}(1:38)
@test ntuple(x) === ntuple(Float64, 38)
xn = renormalize(x) 
@test sum(majorelementvalues(xn)) ≈ 99.8455721816923
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.8455721816923
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.8455721816923
@test propertynames(xn)[1:9] == majorelements(xn) == majorelements(NCKFMASTCrtrace)
@test propertynames(xn)[10:38] == traceelements(xn) == traceelements(NCKFMASTCrtrace)
@test xn isa NCKFMASTCrtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

x = NCKFMASTCrlogtrace(x)
@test x isa NCKFMASTCrlogtrace{Float64}
@test x === NCKFMASTCrlogtrace{Float64}([1:9; log.(10:38)])
@test ntuple(x) === ((1:9.)..., log.(10:38)...,)
xn = renormalize(x) 
@test sum(majorelementvalues(xn)) ≈ 99.8455721816923
xn = renormalize(dehydrate(x)) 
@test sum(majorelementvalues(xn)) ≈ 99.8455721816923
xn = renormalize(x, anhydrous=true) 
@test sum(e->xn[e], propertynames(xn)[1:9]) ≈ 99.8455721816923
@test propertynames(xn)[1:9] == majorelements(xn) == majorelements(NCKFMASTCrlogtrace)
@test propertynames(xn)[10:38] == traceelements(xn) == traceelements(NCKFMASTCrlogtrace)
@test xn isa NCKFMASTCrlogtrace{Float64}
@test 0.5*xn + 0.5*xn == xn
@test !isnan(xn)
@test isnan(xn * NaN)

# Additional composition types
x = NCKFMASHTOCO2trace((1.0:length(fieldnames(NCKFMASHTOCO2trace)))...,)
@test x isa NCKFMASHTOCO2trace{Float64}
@test x === NCKFMASHTOCO2trace{Float64}(1:41)
@test ntuple(x) === ntuple(Float64, 41)
@test majorelementvalues(x) === ((1:11.)...,)
@test traceelementvalues(x) === ((12:41.)...,)

x = NCKFMASHTOCO2logtrace(x)
@test x isa NCKFMASHTOCO2logtrace{Float64}
@test x === NCKFMASHTOCO2logtrace{Float64}([1:11; log.(12:41)])
@test ntuple(x) === ((1:11.)..., log.(12:41)...,)
@test majorelementvalues(x) === ((1:11.)...,)
@test traceelementvalues(x) === (log.(12:41.)...,)

# Test alternate constructors
elements =       [ "SiO2", "Al2O3",  "FeO",  "MgO",  "CaO", "Na2O",  "K2O",  "H2O",  "CO2",]
concentrations = [50.0956, 15.3224, 8.5103, 9.2520, 9.6912, 2.5472, 0.8588, 2.0000, 0.6000,]
data = NamedTuple{(Symbol.(elements)...,)}((concentrations...,))
c1 = NCKFMASHTOtrace{Float64}(data)
@test c1 isa NCKFMASHTOtrace{Float64}
@test c1["SiO2"] ≈ c1[:SiO2] ≈ 50.0956
data = Dict(elements .=> concentrations)
c2 = NCKFMASHTOtrace{Float64}(data)
@test c2 isa NCKFMASHTOtrace{Float64}
@test c2["SiO2"] ≈ c2[:SiO2] ≈ 50.0956
c3 = NCKFMASHTOtrace{Float64}(concentrations, elements)
@test c3 isa NCKFMASHTOtrace{Float64}
@test c3["SiO2"] ≈ c3[:SiO2] ≈ 50.0956
c4 = NCKFMASHTOtrace{Float64}(concentrations, elements)
@test c4 isa NCKFMASHTOtrace{Float64}
@test c4["SiO2"] ≈ c4[:SiO2] ≈ 50.0956
@test c1 === c2 === c3 === c4
@test StatGeochem.normconst(c1) ≈ 0.9827750000000002

## --- Composition distributions
μ = renormalize(NCKFMASHTOCrtrace{Float64}(40:-1:1))
Σ = [Float64(i==j) for i in 1:40, j in 1:40]
d = CompositionNormal(μ, Σ)
@test d isa CompositionNormal{Float64, NCKFMASHTOCrtrace{Float64}}
@test d isa StatGeochem.CompositionDistribution{NCKFMASHTOCrtrace{Float64}}
@test rand(d) isa NCKFMASHTOCrtrace{Float64}
@test rand(d,10) isa Vector{NCKFMASHTOCrtrace{Float64}}
@test mean(d) === μ
@test var(d) === NCKFMASHTOCrtrace{Float64}(diag(Σ))
@test std(d) === NCKFMASHTOCrtrace{Float64}(sqrt.(diag(Σ)))
@test cov(d) == Σ
@test pdf(d, μ) ≈ pdf(d, mean(d)) ≈ 1.0874333119089363e-16
@test logpdf(d, μ) ≈ logpdf(d, mean(d)) ≈ -36.75754132818691
@test pdf(d, μ) ≈ pdf(d, logtrace(mean(d))) ≈ 1.0874333119089363e-16
@test logpdf(d, μ) ≈ logpdf(d, logtrace(mean(d))) ≈ -36.75754132818691

d = CompositionNormal(logtrace(μ), Σ)
@test d isa CompositionNormal{Float64, NCKFMASHTOCrlogtrace{Float64}}
@test d isa StatGeochem.CompositionDistribution{NCKFMASHTOCrlogtrace{Float64}}
@test rand(d) isa NCKFMASHTOCrlogtrace{Float64}
@test rand(d,10) isa Vector{NCKFMASHTOCrlogtrace{Float64}}
@test mean(d) === logtrace(μ)
@test var(d) === NCKFMASHTOCrlogtrace{Float64}(diag(Σ))
@test std(d) === NCKFMASHTOCrlogtrace{Float64}(sqrt.(diag(Σ)))
@test cov(d) == Σ
@test pdf(d, μ) ≈ pdf(d, mean(d)) ≈ 1.0874333119089363e-16
@test logpdf(d, μ) ≈ logpdf(d, mean(d)) ≈ -36.75754132818691
@test pdf(d, μ) ≈ pdf(d, lineartrace(mean(d))) ≈ 1.0874333119089363e-16
@test logpdf(d, μ) ≈ logpdf(d, lineartrace(mean(d))) ≈ -36.75754132818691

## --- Composition arrays
ca = CompositionArray{NCKFMASHTOCrtrace{Float64}}(undef, 99)
@test ca isa CompositionVector{NCKFMASHTOCrtrace{Float64}}
@test ca[1] isa NCKFMASHTOCrtrace{Float64}
@test ca[1:10] isa CompositionVector{NCKFMASHTOCrtrace{Float64}}
@test view(ca, 1:10) isa CompositionVector{NCKFMASHTOCrtrace{Float64}}
@test copy(ca[1:10]) isa CompositionVector{NCKFMASHTOCrtrace{Float64}}
@test length(ca) === 99
@test eachindex(ca) === Base.OneTo(99)
@test ca.SiO2 isa Vector{Float64}
@test length(ca.SiO2) === 99
@test eachindex(ca.SiO2) === Base.OneTo(99)
@test ca.SiO2 === ca[:SiO2]
@test ca.La isa Vector{Float64}
@test length(ca.La) === 99
@test eachindex(ca.La) === Base.OneTo(99)
@test ca.La === ca[:La]

# Randomize and renormalize CompositionArrays
StatGeochem.rand!(ca)
@test nanmean(ca.SiO2) ≈ 10 atol = 3
@test sum(e->ca[1][e], majorelements(ca)) ≈ 99.99 atol=0.02
car = renormalize(ca; anhydrous=true)
@test sum(e->car[1][e], filter(x->!(x===:H2O), majorelements(car)))  ≈ 99.99 atol=0.02
renormalize!(ca; anhydrous=true)
@test sum(e->ca[1][e], filter(x->!(x===:H2O), majorelements(ca)))  ≈ 99.99 atol=0.02
cad = dehydrate(ca)
@test sum(e->cad[1][e], majorelements(cad))  ≈ 99.99 atol=0.02
can = naninf(cad)
@test sum(e->can[1][e], majorelements(can))  ≈ 99.99 atol=0.02

# Other conversions and transformations of CompositionArrays
calog = logtrace(ca)
@test calog isa CompositionVector{NCKFMASHTOCrlogtrace{Float64}}
calin = lineartrace(calog)
@test calin ≈ ca
catup = TupleDataset(ca)
@test catup isa NamedTuple
@test CompositionArray{NCKFMASHTOCrtrace{Float64}}(catup) ≈ ca
cadict = DictDataset(ca)
@test cadict isa Dict
@test CompositionArray{NCKFMASHTOCrtrace{Float64}}(cadict) ≈ ca
cmat = unelementify(ca, floatout=true)
@test cmat isa Matrix{Float64}

cam = partiallymix!(copy(ca), 1)
@test cam[50] ≈ 0.5*ca[1] + 0.5*ca[99]
@test isnormalized(cam[50], anhydrous=true)
@test !isnormalized(cam[50], anhydrous=false)
cam = partiallymix!(copy(ca), 0.5)
@test cam[50] ≈ 0.5ca[50] + 0.5(0.5*ca[1] + 0.5*ca[99])
@test isnormalized(cam[50], anhydrous=true)

ca = CompositionArray{NCKFMASHTOlogtrace{Float64}}(undef, 99)
StatGeochem.rand!(ca)
@test nanmean(ca.SiO2) ≈ 10 atol = 2
@test sum(e->ca[1][e], majorelements(ca)) ≈ 99.99 atol=0.02
renormalize!(ca; anhydrous=true)
@test sum(e->ca[1][e], filter(x->!(x===:H2O), majorelements(ca)))  ≈ 99.99 atol=0.02

cam = partiallymix!(copy(ca), 1)
@test cam[50] ≈ 0.5*ca[1] + 0.5*ca[99]
@test isnormalized(cam[50], anhydrous=true)
@test !isnormalized(cam[50], anhydrous=false)
cam = partiallymix!(copy(ca), 0.5)
@test cam[50] ≈ 0.5ca[50] + 0.5(0.5*ca[1] + 0.5*ca[99])
@test isnormalized(cam[50], anhydrous=true)

# Statistical properties of composition arrays
μ = renormalize(NCKFMASTtrace{Float64}(38:-1:1))
Σ = [Float64(i==j) for i in 1:38, j in 1:38]
d = CompositionNormal(μ, Σ)
ca = CompositionArray{NCKFMASTtrace{Float64}}(undef, 10000)
StatGeochem.rand!(ca, d)
@test nanmean(ca) ≈ μ rtol=0.1
@test nanvar(ca) ≈ var(d) rtol=0.3
@test nanstd(ca) ≈ std(d) rtol=0.1
@test nansem(ca) ≈ std(d)/sqrt(length(ca)) rtol=0.1
@test nancov(ca) ≈ cov(d) rtol=0.3
@test nancovem(ca) ≈ cov(d)/length(ca) rtol=0.3

# Test implicit conversion (filling array with different type)
@test compositiondist isa CompositionNormal{Float64, NCKFMASHTOlogtrace{Float64}}
calog = CompositionArray(rand(compositiondist, 10000))
@test calog isa CompositionVector{NCKFMASHTOlogtrace{Float64}}
μ = nanmean(calog)
@test μ.SiO2 ≈ 60 atol=5
@test μ.Lu ≈ -1.15 atol=0.5

# Test conversion to different eltype
ca = CompositionArray{NCKFMASHTOtrace}(calog)
@test ca isa CompositionVector{NCKFMASHTOtrace{Float64}}
@test nanmean(ca) ≈ nanmean(lineartrace(calog))

# Test implicit conversion of comparable Compositions
StatGeochem.rand!(ca, compositiondist)
μ = nanmean(ca)
@test μ.SiO2 ≈ 60 atol=5
@test μ.Lu ≈ 0.315 atol=0.5
## -- Changepoint.jl tests

    A = [randn(100).-2; randn(100).+2]
    nsteps = 10000
    burnin = 4000

    dist = changepoint(A, nsteps; np=1)
    @test isapprox(nanmean(dist[burnin:end]), 101, atol=4)
    dist = changepoint(A, nsteps; npmin=1, npmax=5)[burnin:end,:]
    @test isapprox(nanmean(dist[dist.>0]), 101, atol=4)

    dist = changepoint(A, ones(200), nsteps; np=1)
    @test isapprox(nanmean(dist[burnin:end]), 101, atol=4)
    dist = changepoint(A, ones(200), nsteps; npmin=1, npmax=5)[burnin:end,:]
    @test isapprox(nanmean(dist[dist.>0]), 101, atol=4)


## --

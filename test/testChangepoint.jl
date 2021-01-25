## -- Changepoint.jl tests

    A = [randn(100).-2; randn(100).+2]
    nsteps = 10000
    burnin = 1000

    dist = changepoint(A, 10000; np=1)
    @test isapprox(nanmean(dist[burnin:end]), 101, atol=2)
    dist = changepoint(A, 10000; npmin=1, npmax=5)[burnin:end,:]
    @test isapprox(nanmean(dist[dist.>0])), 101, atol=2)

    dist = changepoint(A, ones(200), 10000; np=1)
    @test isapprox(nanmean(dist[burnin:end]), 101, atol=2)
    dist = changepoint(A, ones(200), 10000; npmin=1, npmax=5)[burnin:end,:]
    @test isapprox(nanmean(dist[dist.>0])), 101, atol=2)


## --

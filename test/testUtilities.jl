## --- Geochronology.jl

    @test eHf(0.2818792, 0.001009289, 1424.878) ≈ -0.9639408985107067

    (c, m, el, eu) = bin_bsr_eHf(1:4000, fill(0.2818792, 4000), fill(0.001009289,4000), 1:4000, 0, 4000, 4, fill(10,4000), fill(0.001,4000), fill(1e-6,4000), fill(10,4000), 500)
    @test c == 500:1000:3500
    isapprox(m, [-21.316, 0.724, 23.187, 45.97], atol=5)
    isapprox(el, [2.254, 2.216, 2.28, 2.319], atol=0.5)
    isapprox(eu, [2.245, 2.216, 2.273, 2.348], atol=0.5)

## --- GIS.jl

    # Calculate slope on a sphere
    A = [8 7 10 -1 -1; 5 2 2 -1 -1; 8 5 3 1 6; 8 10 3 9 5; -1 8 4 10 4]
    @test aveslope(A*1000, 1:4, 1:4, 1, Int64, minmatval=0, km_per_lat=111.1) ==
        [0 20 0 0 0; 19 23 8 0 0; 18 15 11 22 15; 4 0 15 35 16; 0 14 15 0 0]
    @test aveslope(A*1000, 1:4, 1:4, 1, Float16, minmatval=0, km_per_lat=111.1) ==
        Float16[0.0 19.72 0.0 0.0 0.0; 18.75 22.72 8.29 0.0 0.0; 18.45 15.375 11.2 22.45 14.63; 3.824 0.0 15.06 35.4 15.914; 0.0 13.59 15.39 0.0 0.0]
    @test maxslope(A*1000, 1:4, 1:4, 1, minmatval=0, km_per_lat=111.1) ==
        [0 45 72 0 0; 27 16 32 0 6; 38 36 22 14 45; 18 23 22 41 36; 0 32 38 45 0]
    @test maxslope(A*1000, 1:4, 1:4, 1, Float16, minmatval=0, km_per_lat=111.1) ==
        Float16[0.0 45.0 72.0 0.0 0.0; 27.02 15.914 31.5 0.0 6.367; 38.22 36.0 22.3 13.52 45.06; 18.05 22.56 22.3 40.5 36.1; 0.0 31.86 38.22 44.6 0.0]


    # Test ESRI arc/info ASCII grid import function
    grid = """
    ncols         4
    nrows         6
    xllcorner     0.0
    yllcorner     0.0
    cellsize      50.0
    NODATA_value  -9999
    -9999 -9999 5 2
    -9999 20 100 36
    3 8 35 10
    32 42 50 6
    88 75 27 9
    13 5 1 -9999
    """

    f = open("grid.asc","w")
    print(f, grid)
    close(f)

    (data, metadata) = importAAIGrid("grid.asc", Int64, undefval=-999)
    @test eltype(data) == Int64
    @test data == [-9999 -9999 5 2; -9999 20 100 36; 3 8 35 10; 32 42 50 6; 88 75 27 9; 13 5 1 -9999]
    @test metadata["nodata"] == -9999



    # Random lat-lon generation
    @test isa(randlatlon(), Tuple{Float64,Float64})
    lat,lon = randlatlon(100)
    @test maximum(lat) <= 90
    @test minimum(lat) >= -90
    @test maximum(lon) <= 180
    @test minimum(lon) >= -180

    lat,lon = randlatlon(100,land=true)
    @test maximum(lat) <= 90
    @test minimum(lat) >= -90
    @test maximum(lon) <= 180
    @test minimum(lon) >= -180
    @test all(find_land(lat, lon))
    @test length(lat) == length(lon) == 100

    # Calculate arc-degree distance
    @test isapprox(haversine(1, 0, 0, 0), 1)
    @test isapprox(haversine(0, 1, 0, 0), 1)
    @test isapprox(haversine(0, 0, 1, 0), 1)
    @test isapprox(haversine(0, 0, 0, 1), 1)
    @test isapprox(haversine(0, 0, 0, 0), 0)

    @test isapprox(haversine(90, 2, 0, 0), 90)
    @test isapprox(haversine(0, 0, 90, 2), 90)

    # Centroid of a set of lats and lons on a sphere
    lats, lons = [-1, 1, 0, 0.], [0, 0, -1, 1.]
    @test all(centroid(lats, lons) .≈ (0.0, 0.0))
    lats, lons = [43.69852352,43.69852385,43.69944918,43.69945593,], [-116.0948702,-116.0936334,-116.0936182,-116.0948765,];
    @test all(centroid(lats, lons) .≈ (43.698988121696146, -116.0942495750004))

    # Maximum arc-degree distance between a list of points on a sphere
    lats, lons = [-1, 1, 0, 0.], [0, 0, -1, 1.]
    @test all(dist_uncert(lats, lons) .≈ (0.0, 0.0, 1.0))
    lats, lons = [0, 0, 0, 0], [0, 30, 23, 90]
    @test all(dist_uncert(lats, lons) .≈ (0.0, 34.15788270532762, 45.0))


## --- Etc.jl

    # Test digitize_plotmarkers
    img = load("assets/xyscatter.png")
    C = eltype(img)
    (x,dx,y,dy) = digitize_plotmarkers(img, C(0,0.604,0.976,1), (0,10), (0,10))
    @test isapprox(x, 1:10, atol = 0.1)
    @test isapprox(y, 1:10, atol = 0.1)


    # Test digitize_plotline
    img = load("assets/xysin.png")
    C = eltype(img)
    (x,y) = digitize_plotline(img, C(0,0.604,0.976,1), (0,2pi), (-1.1,1.1))
    @test isapprox(sin.(x), y, atol=0.1)


## ---

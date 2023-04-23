## --- tc1/tc1.jl

   @test find_tc1_crust(33,-100) == 40.0
   @test find_tc1_crust([33,-30,70],[-100,20,-130]) == [40.0, 38.0, 46.0]

   @test find_tc1_lith(33,-100) == 132.0
   @test find_tc1_lith([33,-30,70],[-100,20,-130]) == [132.0, 123.0, 167.0]

   @test find_tc1_age(33,-100) == (1400.0, 1100.0, 1700.0)
   @test find_tc1_age([33,-30,70],[-100,20,-130]) ==
      ([1400.0, 1400.0, 2100.0], [1100.0, 1100.0, 1700.0], [1700.0, 1700.0, 2500.0])

## --- PartitionCoefficients/PartitionCoefficients.jl

   @test Float16.(claiborne_zircon_kd.(["Hf","Th","U","Y","Nb","Nd","Sm","Tb","Eu","Dy","Gd","Ho","Er","Yb","Tm","Lu"], 800)) ==
      Float16[1.303e3, 6.46, 36.97, 33.47, 0.2556, 0.02368, 0.4363, 10.695, 0.8784, 19.56, 3.293, 40.34, 62.88, 101.75, 94.94, 126.4]

## --- Elevation.jl

   lat = [43.7022,-26.2041,-19.5723,-34.9285,46.4908]
   lon = [-72.2896,28.0473,65.7550,138.6007,9.8355]
   @test find_geolcont(lat, lon) == [3, 1, 7, 5, 2]
   @test find_geolcont(43.702245, -72.0929) == fill(3)

   @test find_geolprov(lat, lon) == [10, 31, 0, 10, 10]
   @test find_geolprov(43.702245, -72.0929) == fill(10)

   @test find_land(lat, lon) == Bool[1, 1, 0, 1, 1]
   @test find_land(43.702245, -72.0929) == fill(true)

   A = (1:200)*(1:200)'
   @test find_etopoelev(A, -90:-89, -180:-179) == [1,3721]
   @test find_srtm15plus(A, -90:0.1:-89.5, -180:0.1:-179.5) == [1, 625, 2401, 5329, 9409, 14641]
   @test find_seafloorage(A, 80.738:-0.1:80, 0:0.1:0.7) == [1, 80, 266, 570, 975, 1472, 2090, 2667]


## ---

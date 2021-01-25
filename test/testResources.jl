## --- tc1/tc1.jl

   @test find_tc1_crust(33,-100) == 40.0
   @test find_tc1_crust([33,-30,70],[-100,20,-130]) == [40.0, 38.0, 46.0]

   @test find_tc1_lith(33,-100) == 132.0
   @test find_tc1_lith([33,-30,70],[-100,20,-130]) == [132.0, 123.0, 167.0]

   @test find_tc1_age(33,-100) == (1400.0, 1100.0, 1700.0)
   @test find_tc1_age([33,-30,70],[-100,20,-130]) ==
      ([1400.0, 1400.0, 2100.0], [1100.0, 1100.0, 1700.0], [1700.0, 1700.0, 2500.0])

## --- PartitionCoefficients/PartitionCoefficients.jl

   @test claiborne_zircon_kd.(["Sm","Yb"], 800) ≈ [0.4362887941802957, 101.77562766576213]


## --- Colormaps.jl

   using Colors: Color

   cmap = resize_colormap(viridis, 10)
   @test length(cmap) == 10
   @test isa(cmap, Array{<:Color,1})

   matrix = rand(10,10)
   # Specifiying limits
   img1 = imsc(matrix, viridis, 0, 1)
   @test isa(img1, Array{<:Color,2})
   img2 = imsci(matrix, viridis, 0, 1)
   @test isa(img2, AbstractArray{<:Color,2})
   @test all(img1 .== img2)

   # Auto-ranging
   img1 = imsc(matrix, viridis)
   @test isa(img1, Array{<:Color,2})
   img2 = imsci(matrix, viridis)
   @test isa(img2, AbstractArray{<:Color,2})
   @test all(img1 .== img2)


## ---

## --- Get Crust 1.0

	# Have to download crust 1 before testing it!
	@test get_crust1() == 0

	# Note: Available layers in Crust 1:
	# 1) water
	# 2) ice
	# 3) upper sediments   (VP, VS, rho not defined in all cells)
	# 4) middle sediments  "
	# 5) lower sediments   "
	# 6) upper crystalline crust
	# 7) middle crystalline crust
	# 8) lower crystalline crust

	lats = [43.70, 39.2508, NaN]
	lons = [-72.29, -106.2925, NaN]

## -- Test seismic data

	# Test single lat-lon-depth point from crust1
	@test [x[1] for x in find_crust1_seismic(50,50,8)] ≈ [7.1, 4.05, 3000.0]
	@test [x[1] for x in find_crust1_seismic(50,NaN,8)] ≈ [NaN, NaN, NaN] nans=true

	# Test fetching lat-lon array from crust1
	@test all(isapprox.(find_crust1_seismic(lats,lons,6), ([6.3, 6.0, NaN], [3.63, 3.5, NaN], [2790.0, 2720.0, NaN]), nans=true))
	@test all(isapprox.(find_crust1_seismic(lats,lons,7), ([6.6, 6.6, NaN], [3.80, 3.8, NaN], [2860.0, 2860.0, NaN]), nans=true))
	@test all(isapprox.(find_crust1_seismic(lats,lons,8), ([7.0, 7.2, NaN], [3.99, 4.1, NaN], [2950.0, 3030.0, NaN]), nans=true))
	@test all(isapprox.(find_crust1_seismic(lats,lons,:upper_crust), ([6.3, 6.0, NaN], [3.63, 3.5, NaN], [2790.0, 2720.0, NaN]), nans=true))
	@test all(isapprox.(find_crust1_seismic(lats,lons,:middle_crust), ([6.6, 6.6, NaN], [3.80, 3.8, NaN], [2860.0, 2860.0, NaN]), nans=true))
	@test all(isapprox.(find_crust1_seismic(lats,lons,:lower_crust), ([7.0, 7.2, NaN], [3.99, 4.1, NaN], [2950.0, 3030.0, NaN]), nans=true))

## --  Test layer thickness data

	@test find_crust1_thickness(lats,lons,6) ≈ [15.4, 17.73, NaN] nans=true
	@test find_crust1_thickness(lats,lons,7) ≈ [13.57, 16.21, NaN] nans=true
	@test find_crust1_thickness(lats,lons,8) ≈ [7.70, 16.72, NaN] nans=true


## -- Test cumulative thickness data

	@test find_crust1_base(lats,lons,6) ≈ [-14.99, -14.91, NaN] nans=true
	@test find_crust1_base(lats,lons,7) ≈ [-28.56, -31.12, NaN] nans=true
	@test find_crust1_base(lats,lons,8) ≈ [-36.26, -47.84, NaN] nans=true

## -- Test complete data

	@test all(isapprox.(find_crust1_layer(lats,lons,6), ([6.3, 6.0, NaN], [3.63, 3.5, NaN], [2790., 2720., NaN], [15.4, 17.73, NaN]), nans=true))
	@test all(isapprox.(find_crust1_layer(lats,lons,7), ([6.6, 6.6, NaN], [3.8, 3.8, NaN], [2860., 2860., NaN], [13.57, 16.21, NaN]), nans=true))
	@test all(isapprox.(find_crust1_layer(lats,lons,8), ([7.0, 7.2, NaN], [3.99, 4.1, NaN], [2950., 3030., NaN], [7.70, 16.72, NaN]), nans=true))

## ---

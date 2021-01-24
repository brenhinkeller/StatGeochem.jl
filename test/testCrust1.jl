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

	const lats = [43.70, 39.2508]
	const lons = [-72.29, -106.2925]

## -- Test seismic data

	# Test single lat-lon-depth point from crust1
	@test [x[1] for x in find_crust1_seismic(50,50,8)] == [7.1, 4.05, 3000.0]

	# Test fetching lat-lon array from crust1
	@test find_crust1_seismic(lats,lons,6) == ([6.3, 6.0], [3.63, 3.5], [2790.0, 2720.0])
	@test find_crust1_seismic(lats,lons,7) == ([6.6, 6.6], [3.80, 3.8], [2860.0, 2860.0])
	@test find_crust1_seismic(lats,lons,8) == ([7.0, 7.2], [3.99, 4.1], [2950.0, 3030.0])

## --  Test layer thickness data

	@test find_crust1_thickness(lats,lons,6) ≈ [15.4, 17.73]
	@test find_crust1_thickness(lats,lons,7) ≈ [13.57, 16.21]
	@test find_crust1_thickness(lats,lons,8) ≈ [7.70, 16.72]


## -- Test cumulative thickness data

	@test find_crust1_base(lats,lons,6) ≈ [-14.99, -14.91]
	@test find_crust1_base(lats,lons,7) ≈ [-28.56, -31.12]
	@test find_crust1_base(lats,lons,8) ≈ [-36.26, -47.84]

## -- Test complete data

	@test find_crust1_layer(lats,lons,6) ≈ ([6.3, 6.0], [3.63, 3.5], [2.79, 2.72], [15.4, 17.73])
	@test find_crust1_layer(lats,lons,7) ≈ ([6.6, 6.6], [3.8, 3.8], [2.86, 2.86], [13.57, 16.21])
	@test find_crust1_layer(lats,lons,8) ≈ ([7.0, 7.2], [3.99, 4.1], [2.95, 3.03], [7.70, 16.72])

## ---

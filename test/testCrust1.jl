## -- Test getting seismic data

	# Test single lat-lon-depth point from crust1
	@test [x[1] for x in find_crust1_seismic(50,50,8)] == [7.1, 4.05, 3000.0]

	# Test fetching lat-lon array from crust1
	@test find_crust1_seismic([50,51],[50,51],8) == ([7.1, 7.1], [4.05, 4.05], [3000.0, 3000.0])

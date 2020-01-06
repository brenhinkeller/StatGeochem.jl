## Test getting seismic data (bug fix: previously did not work with arrays)

using StatGeochem

# Test single value 
out = find_crust1_seismic(50,50,8)
out = [x[1] for x in out]
if out != [7.1, 4.05, 3000.0]
	throw(AssertionError("Incorrect value returned"))
end 

out = find_crust1_seismic([50,51],[50,51],8)
if out != ([7.1, 7.1], [4.05, 4.05], [3000.0, 3000.0])
	throw(AssertionError("Incorrect value returned"))
end 
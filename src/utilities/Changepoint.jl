function update_changepoint_model!(m, σ, d, boundaries, np)
	@inbounds for i=1:np+1
		r = boundaries[i]:boundaries[i+1]
		for col = 1:size(d,2)
			m[r,col] .= nanmean(view(d,r,col))
			σ[r,col] .= nanstd(view(d,r,col))
		end
	end
end

function update_changepoint_mu!(m, d, boundaries, np)
	@inbounds for i=1:np+1
		r = boundaries[i]:boundaries[i+1]
		for col = 1:size(d,2)
			m[r,col] .= nanmean(view(d,r,col))
		end
	end
end

# # Not currently used
# function update_changepoint_sigma!(σ, d, boundaries, np)
# 	@inbounds for i=1:np+1
# 		r = boundaries[i]:boundaries[i+1]
# 		for col = 1:size(d,2)
# 			σ[r,col] .= nanstd(view(d,r,col))
# 		end
# 	end
# end


"""
```julia
changepoint(data, [sigma], nsteps; np, npmin, npmax)
```
Given an ordered array of `data` points, optionally with uncertainties `sigma`,
use a Markov chain Monte Carlo approach based on that of Gallagher et al., 2010
(10.1016/j.epsl.2011.09.015) to estimate the position (by index) and optionally
number of changepoints that best explain the `data`. Will return the results for
`nsteps` steps of the Markov chain.

Optionally, you may also specify as keyword arguments either `np` or
`npmin` and `npmax` to constrain the number of allowed changepoints.

Currently prints results to the terminal (stdout), one line per step of the Markov
chain, but a more efficent output format is probably desirable in a future version
of this function.

### Examples
```julia
julia> A = [randn(100).-2; randn(100).+2];

julia> dist = changepoint(A, 10000; np=1);

julia> dist[9000:end]
1001-element Vector{Int64}:
 101
 101
 101
 101
   ⋮
 101
 101
 101
```
"""
function changepoint(data::AbstractArray, nsteps::Integer; np::Integer=0, npmin::Integer=0, npmax::Integer=min(size(data,1) ÷ 2, 11))

    MOVE = 0.70
    BIRTH = 0.15
    DEATH = 0.15

	DEBUG = false
	FORMATTED = true

	T = float(eltype(data))
	nrows = size(data,1)
	ncolumns = size(data,2)
	m = similar(data, T)
	σ = similar(data, T) #Array{T}(undef, ncolumns)
	σₚ = similar(data, T) #Array{T}(undef, ncolumns)

	# Number of possible changepoint locations
	K = nrows-1

	# Parse provided options
	if 0 < np <= K
		# If np is specified, use that
		npmin = npmax = np
	else
		# Otherwise, ensure all provided values are plausible and go with that
		npmax > K && (npmax = K)
		npmin < 0 && (npmin = 0)
		np = min(max(npmin, 2), npmax)
	end

	# Allocate output array of changepoints
	result = fill(0, nsteps, npmax)

	# Create and fill initial boundary point array
	boundaries = Array{Int}(undef, K+2)
	boundaries[1] = 1
	boundaries[np+2] = nrows
	boundaries[2:np+1] .= rand(2:nrows-1, np)
	boundariesₚ = similar(boundaries)
	boundary_sigma = nrows/np

	np = count_unique!(view(boundaries,1:np+2)) - 2

	# Calculate initial proposal and log likelihood
	update_changepoint_model!(m, σ, data, boundaries, np)
	ll = normpdf_ll(m, σ, data)

	# The actual loop
	@inbounds for i = 1:nsteps

		# Randomly choose a type of modification to the model
		r = rand()
		u = rand()

		# Update the model with the chosen modification
		if r < MOVE && np>0
			# Move a changepoint
			copyto!(boundariesₚ,1,boundaries,1,np+2)
			# Pick which changepoint to move
			pick = rand(2:np+1)
			# Move the changepoint
			boundary_adj = randn()*boundary_sigma
			boundary_prop = boundariesₚ[pick] + round(Int, boundary_adj)

			# Treat ends of array as periodic boundary conditions
			boundariesₚ[pick] = mod(boundary_prop - 1, nrows-1) + 1

			# Check if this has caused any redundancies
			npₚ = count_unique!(view(boundariesₚ,1:np+2)) - 2

			# Update the model
			update_changepoint_mu!(m, data, boundariesₚ, npₚ)

			# Calculate log likelihood for proposal
			if (1 < boundary_prop < nrows)
				llₚ = normpdf_ll(m, σ, data)
			else
				llₚ = -Inf
			end

			DEBUG && println("Move: llₚ-ll = $llₚ - $ll")
			if log(u) < llₚ-ll
				DEBUG && println("Accepted!")
				ll = llₚ
				boundary_sigma = abs(boundary_adj)*2.9
				# println("sigma: $boundary_sigma")
				copyto!(boundaries,1,boundariesₚ,1,np+2)
				# for n=1:np
				# 	print("$(boundariesₚ[n+1]),")
 				# end
				# FORMATTED && print("\n")
			end

		elseif r < MOVE+BIRTH
			# Add a changepoint
			if np < npmax
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Propose a new changepoint
				boundariesₚ[np+3] = rand(2:nrows-1)
				npₚ = count_unique!(view(boundariesₚ,1:np+3)) - 2

				# Update the model
				# update_changepoint_model!(m, σ, data, boundariesₚ, npₚ)
				update_changepoint_mu!(m, data, boundariesₚ, npₚ)

				# Calculate log likelihood for proposal
				lqz = sum(1 ./ (2*σ.*σ))
				llₚ = normpdf_ll(m, σ, data)
				DEBUG && println("Birth: -lqz+llₚ-ll = $(-lqz) + $llₚ - $ll")
				if log(u) < llₚ-lqz-ll
					DEBUG && println("Accepted!")
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					# for n=1:np
					# 	print("$(boundariesₚ[n+1]),")
					# end
					# FORMATTED && print("\n")
				end
			end
		elseif r < MOVE+BIRTH+DEATH
			# Delete a changepoint
			if np > npmin
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Pick which changepoint to delete
				pick = rand(2:np+1)
				boundariesₚ[pick]=boundariesₚ[pick+1]
				npₚ = count_unique!(view(boundariesₚ,1:np+2)) - 2

				# Update the model
				# update_changepoint_model!(m, σ, data, boundariesₚ, npₚ)
				update_changepoint_mu!(m, data, boundariesₚ, npₚ)

				# Calculate log likelihood for proposal
				llₚ = normpdf_ll(m, σ, data)

				lqz = sum(1 ./ (2*σ.*σ))
				DEBUG && println("Death: lqz+llₚ-ll = $lqz + $llₚ - $ll")
				if log(u) < llₚ+lqz-ll
					DEBUG && println("Accepted!")
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					# for n=1:np
					# 	print("$(boundariesₚ[n+1]),")
					# end
					# FORMATTED && print("\n")
				end
			end
		end

		# Record results
		result[i, 1:np] .= boundaries[2:(np+1)]
	end
	return result
end
function changepoint(data::AbstractArray, sigma::AbstractArray, nsteps::Integer; np::Integer=0, npmin::Integer=0, npmax::Integer=min(size(data,1) ÷ 2, 11))

    MOVE = 0.70
    BIRTH = 0.15
    DEATH = 0.15

	DEBUG = false
	FORMATTED = true

	T = float(eltype(data))
	nrows = size(data,1)
	ncolumns = size(data,2)
	m = similar(data, T)

	# Number of possible changepoint locations
	K = nrows-1

	# Parse provided options
	if 0 < np <= K
		# If np is specified, use that
		npmin = npmax = np
	else
		# Otherwise, ensure all provided values are plausible and go with that
		npmax > K && (npmax = K)
		npmin < 0 && (npmin = 0)
		np = min(max(npmin, 2), npmax)
	end

	# Allocate output array of changepoints
	result = fill(0, nsteps, npmax)

	# Create and fill initial boundary point array
	boundaries = Array{Int}(undef, K+2)
	boundaries[1] = 1
	boundaries[np+2] = nrows
	boundaries[2:np+1] .= rand(2:nrows-1, np)
	boundariesₚ = similar(boundaries)
	boundary_sigma = nrows/np

	np = count_unique!(view(boundaries,1:np+2)) - 2

	# Calculate initial proposal and log likelihood
	update_changepoint_mu!(m, data, boundaries, np)
	ll = normpdf_ll(m, sigma, data)

	# The actual loop
	@inbounds for i = 1:nsteps

		# Randomly choose a type of modification to the model
		r = rand()
		u = rand()

		# Update the model with the chosen modification
		if r < MOVE && np>0
			# Move a changepoint
			copyto!(boundariesₚ,1,boundaries,1,np+2)
			# Pick which changepoint to move
			pick = rand(2:np+1)
			# Move the changepoint
			boundary_adj = randn()*boundary_sigma
			boundary_prop = boundariesₚ[pick] + round(Int, boundary_adj)

			# Treat ends of array as periodic boundary conditions
			boundariesₚ[pick] = mod(boundary_prop - 1, nrows-1) + 1

			# Check if this has caused any redundancies
			npₚ = count_unique!(view(boundariesₚ,1:np+2)) - 2

			# Update the model
			update_changepoint_mu!(m, data, boundariesₚ, npₚ)

			# Calculate log likelihood for proposal
			if (1 < boundary_prop < nrows)
				llₚ = normpdf_ll(m, sigma, data)
			else
				llₚ = -Inf
			end

			DEBUG && println("Move: llₚ-ll = $llₚ - $ll")
			if log(u) < llₚ-ll
				DEBUG && println("Accepted!")
				ll = llₚ
				boundary_sigma = abs(boundary_adj)*2.9
				# println("sigma: $boundary_sigma")
				copyto!(boundaries,1,boundariesₚ,1,np+2)
				# for n=1:np
				# 	print("$(boundariesₚ[n+1]),")
				# end
				# FORMATTED && print("\n")
			end
		elseif r < MOVE+BIRTH
			# Add a changepoint
			if np < npmax
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Propose a new changepoint
				boundariesₚ[np+3] = rand(2:nrows-1)
				npₚ = count_unique!(view(boundariesₚ,1:np+3)) - 2

				# Update the model
				update_changepoint_mu!(m, data, boundariesₚ, npₚ)

				# Calculate log likelihood for proposal
				lqz = sum(1 ./ (2*sigma.*sigma))
				llₚ = normpdf_ll(m, sigma, data)
				DEBUG && println("Birth: -lqz+llₚ-ll = $(-lqz) + $llₚ - $ll")
				if log(u) < llₚ-lqz-ll
					DEBUG && println("Accepted!")
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					# for n=1:np
					# 	print("$(boundariesₚ[n+1]),")
					# end
					# FORMATTED && print("\n")
				end
			end
		elseif r < MOVE+BIRTH+DEATH
			# Delete a changepoint
			if np > npmin
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Pick which changepoint to delete
				pick = rand(2:np+1)
				boundariesₚ[pick]=boundariesₚ[pick+1]
				npₚ = count_unique!(view(boundariesₚ,1:np+2)) - 2

				# Update the model
				update_changepoint_mu!(m, data, boundariesₚ, npₚ)

				# Calculate log likelihood for proposal
				llₚ = normpdf_ll(m, sigma, data)

				lqz = sum(1 ./ (2*sigma.*sigma))
				DEBUG && println("Death: lqz+llₚ-ll = $lqz + $llₚ - $ll")
				if log(u) < llₚ+lqz-ll
					DEBUG && println("Accepted!")
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					# for n=1:np
					# 	print("$(boundariesₚ[n+1]),")
					# end
					# FORMATTED && print("\n")
				end
			end
		end

		# Record results
		result[i, 1:np] .= boundaries[2:(np+1)]
	end
	return result
end
export changepoint

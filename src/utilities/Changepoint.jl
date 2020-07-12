function update_changepoint_model!(m, σ, d, boundaries, np)
	@inbounds for i=1:np+1
		r = boundaries[i]:boundaries[i+1]
		for col = 1:size(m,2)
			m[r,col] .= nanmean(view(d,r,col))
			σ[r,col] .= nanstd(view(d,r,col))
		end
	end
end


function changepoint(data::AbstractArray, nsims::Integer, npoints::Integer; npmin::Integer=0, npmax::Integer=0)

    MOVE = 0.30
    UPDATE = 0.00
    SIGMA = 0.20
    BIRTH = 0.25
    DEATH = 0.25

	DEBUG = false
	FORMATTED = false

	T = float(eltype(data))
	model = similar(data, T)
	nrows = size(data,1)
	ncolumns = size(data,2)
	σ = Array{mtype}(T, ncolumns) # TODO: Decide whether to actually do something with these sigmas
	σₚ = Array{mtype}(T, ncolumns)

	# Number of possible changepoint locations
	K = nrows-1
	if !(0 < npmax < K)
		npmax = K
	end

	# Number of initial changepoints
	np = npₚ = 2

	# Fill initial boundary point array
	boundaries = Array{Int}(undef, K+2)
	boundaries[1] = 1
	boundaries[np+2] = nrows
	boundaries[2:np+1] .= rand(2:nrows-1, np)
	boundariesₚ = similar(boundaries)
	boundary_sigma = nrows/np

	np = count_unique!(view(boundaries,np+2)) - 2
	copyto!(boundariesₚ,1,boundaries,1,np+2)


	# The actual loop
	@inbounds for i = 1:nsims

		# Randomly choose a type of modification to the model
		r = rand()
		u = rand()

		if r < MOVE && np>0
			# Move a changepoint
			copyto!(boundariesₚ,1,boundaries,1,np+2)
			# Pick which changepoint to move
			pick = rand(2:np+1)
			# Move the changepoint
			boundariesₚ[pick] .+= round(Int, randn()*boundary_sigma)
			# Treat ends of array as periodic boundary conditions
			boundariesₚ[pick] = mod(boundariesₚ[pick] - 1, nrows-1) + 1
			# Check if this has caused any redundancies
			npₚ = count_unique!(view(boundariesₚ,np+2)) - 2

			# Update the model
			update_changepoint_model!(m, σ, d, boundariesₚ, npₚ)

			# Calculate log likelihood for proposal
			llₚ = normpdf_ll(m, σ, d)

			DEBUG && print("Move: llₚ-ll = %g - %g\n",llₚ,ll)
			if log(u) < llₚ-ll
				DEBUG && print("Accepted!\n")
				ll = llₚ
				boundary_sigma = abs(boundariesₚ[pick]-boundaries[pick])*2.9
				# printf("%g\n",boundary_sigma)
				copyto!(boundaries,1,boundariesₚ,1,np+2)
				for n=1:np
					printf("%u,",boundariesₚ[n+1])
				end
				FORMATTED && print("\n")
				# printf("\n%u\n",np)
			end

		elseif r < MOVE+SIGMA
			# Adjust σ
			copyto!(σₚ,σ)

			# Choose a new sigma
			j = rand(1:columns)
			σₚ[j] = rand()
			sll = log(σₚ[j]/σ[j]) * rows

			# Calculate log likelihood for proposal
			llₚ = normpdf_ll(m, σₚ, d)

			DEBUG && println("Sigma $(σₚ[j]): llₚ-ll = $llₚ - $ll")
			if log(u) < llₚ+sll-ll
				# If accepted
				DEBUG && println("Accepted!")
				ll = llₚ
				copyto!(σ,σₚ)
			end

		elseif r < MOVE+SIGMA+BIRTH
			# Add a changepoint
			if np < npmax
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Propose a new changepoint
				boundariesₚ[np+3] = rand(2:nrows-1)
				npₚ = count_unique!(view(boundariesₚ,np+3)) - 2

				# Update the model
				update_changepoint_model!(m, σ, d, boundariesₚ, npₚ)

				# Calculate log likelihood for proposal
				lqz=lqxz(d, m, &boundariesₚ[pick], rows, columns); #TODO: do something with / make replacement for lqxz
				llₚ = normpdf_ll(m, σ, d)
				DEBUG && print("Birth: -lqz+llₚ-ll = %g + %g - %g\n",-lqz,llₚ,ll)
				if log(u) < llₚ-lqz-ll
					DEBUG && print("Accepted!\n")
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					for n=1:np
						printf("%u,",boundariesₚ[n+1])
					end
					FORMATTED && print("\n")
				end
			end
		elseif r < MOVE+SIGMA+BIRTH+DEATH
			# Delete a changepoint
			if np > npmin
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Pick which changepoint to delete
				pick = rand(2:np+1)
				boundariesₚ[pick]=boundariesₚ[pick+1]
				npₚ = count_unique!(view(boundariesₚ,np+2) - 2

				# Update the model
				update_changepoint_model!(m, σ, d, boundariesₚ, npₚ)

				# Calculate log likelihood for proposal
				llₚ = normpdf_ll(m, σ, d)

				lqz=lqxz(d, m, &boundaries[pick], rows, columns);
				DEBUG && println("Death: lqz+llₚ-ll = $lqz + $llₚ - $ll")
				if log(u) < llₚ+lqz-ll
					DEBUG && println("Accepted!")
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					for n=1:np
						printf("%u,",boundariesₚ[n+1])
					end
					FORMATTED && print("\n")
				end
			end

		else
			# Update the model
			npₚ = count_unique!(view(boundariesₚ,np+2)) - 2
			update_changepoint_model!(m, σ, d, boundariesₚ, npₚ)

			# Calculate log likelihood for proposal
			llₚ = normpdf_ll(m, σ, d)
			if log(u) < llₚ-ll
				ll = llₚ
			end
		end
	end
end

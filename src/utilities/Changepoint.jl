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
	σ = Array{mtype}(T, ncolumns)
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
			# Move the changepoint between its boundaries
			nextValueUniform = pcg32_random_r(&rng) / ( RAND_MAX_U32 / ((boundariesₚ[pick+1]-1)-(boundariesₚ[pick-1]+1)+1) ) + (boundariesₚ[pick-1]+1);
			nextValueGaussian = (uint32_t)(pcg_gaussian_ziggurat(&rng, lastDifference) + boundaries[pick]);
			if (abs((int)nextValueUniform-(int)boundariesₚ[pick]) > abs((int)nextValueGaussian-(int)boundariesₚ[pick]) && nextValueGaussian>boundariesₚ[pick-1] && nextValueGaussian<boundariesₚ[pick+1] && nextValueGaussian != boundariesₚ[pick]){
				# printf("GaussianCloser: %u %u %u ", nextValueUniform, boundariesₚ[pick], nextValueGaussian);
				boundariesₚ[pick] = nextValueGaussian;
			else
				# printf("UniformCloser: %u %u %u ", nextValueUniform, boundariesₚ[pick], nextValueGaussian);
				boundariesₚ[pick] = nextValueUniform;
			end
			# Update the model
			update_model(boundariesₚ, np, rows, columns, &rng, d, σ, m);

			# Calculate log likelihood for proposal
			llₚ=log_likelihood(d, m, σ, rows, columns);

			DEBUG && print("Move: llₚ-ll = %g - %g\n",llₚ,ll);}
			if log(u) < llₚ-ll
				DEBUG && print("Accepted!\n");}
				ll = llₚ
				lastDifference = fabs(((double)boundariesₚ[pick]-(double)boundaries[pick])*2.9);
				# printf("%g\n",lastDifference);
				copyto!(boundaries,1,boundariesₚ,1,np+2)
				for n=1:np
					printf("%u,",boundariesₚ[n+1]);
				end
				FORMATTED && print("\n")
				# printf("\n%u\n",np);
			end

		elseif r < MOVE+SIGMA
			# Adjust σ
			copyto!(σₚ,σ)

			# Choose a new sigma
			j=rand(1:columns)
			σₚ[j]=rand()
			sll=log(pow((σₚ[j]/σ[j]),rows));

			# Calculate log likelihood for proposal
			llₚ=log_likelihood(d, m, σₚ, rows, columns);

			DEBUG && println("Sigma $(σₚ[j]): llₚ-ll = $llₚ - $ll")
			if log(u) < llₚ+sll-ll
				# If accepted
				DEBUG && println("Accepted!")
				ll = llₚ
				copyto!(σ,σₚ)
				# for (int n=0; n<np+2; n++){
				# 	printf("%u,",boundariesₚ[n]);
				# }
				# printf("\n%u\n",np);
				# printf("sigma: ");
				# for (int n=0; n<columns; n++){
				# 	printf("%g,",σ[n]);
				# }
				# printf("\n");
			end

		elseif r < MOVE+SIGMA+BIRTH
			# Add a changepoint
			if np < npmax
				copyto!(boundariesₚ,1,boundaries,1,np+2)

				# Pick which changepoint to add right of
				pick=rand(1:np+1)
				boundariesₚ[np+3] = rand(boundariesₚ[pick]:boundariesₚ[pick+1])
				npₚ = count_unique!(view(boundariesₚ,np+3)) - 2

				# Update the model
				update_model(boundariesₚ, npₚ, rows, columns, &rng, d, σ, m);

				# Calculate log likelihood for proposal
				lqz=lqxz(d, m, &boundariesₚ[pick], rows, columns);
				llₚ=log_likelihood(d, m, σ, rows, columns);
				DEBUG && print("Birth: -lqz+llₚ-ll = %g + %g - %g\n",-lqz,llₚ,ll)
				if log(u) < llₚ-lqz-ll
					DEBUG && print("Accepted!\n");}
					ll = llₚ
					np = npₚ
					copyto!(boundaries,1,boundariesₚ,1,np+2)
					for n=1:np
						printf("%u,",boundariesₚ[n+1])
					end
					FORMATTED && print("\n")
				# printf("\n%u\n",np);
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
				update_model(boundariesₚ, np, rows, columns, &rng, d, σ, m);

				# Calculate log likelihood for proposal
				llₚ=log_likelihood(d, m, σ, rows, columns);
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
				# printf("\n%u\n",np);
				end
			end

		else
			# Update the model
			update_model(boundariesₚ, np, rows, columns, &rng, d, σ, m)

			# Calculate log likelihood for proposal
			llₚ=log_likelihood(d, m, σ, rows, columns)
			if log(u) < llₚ-ll
				# printf("Update: llₚ-ll = %g - %g\n",llₚ,ll);
				# printf("Accepted!\n");
				ll = llₚ
			end
		end
		# printf("%u\n",i);
	end
end

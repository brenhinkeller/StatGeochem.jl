function changepoint(data::AbstractArray, nsims::Integer, npoints::Integer; npmin::Integer=0, npmax::Integer=0)

    MOVE = 0.30
    UPDATE = 0.00
    SIGMA = 0.20
    BIRTH = 0.25
    DEATH = 0.25

	DEBUG = false

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
	boundarypoints = Array{Int}(undef, K+2)
	boundarypoints[1] = 1
	boundarypoints[np+2] = nrows
	boundarypoints[2:np+1] .= rand(2:nrows-1, np)
	boundarypointsₚ = copy(boundarypoints)



	np = count_unique!(view(boundarypoints,np+2)) - 2
	copyto!(boundarypoints,1,boundarypointsₚ,1,np+2)
	copyto!(boundarypoints,1,boundarypointsₚ,1,np+2)


	# The actual loop
	for (i=0;i<nsims;i++){

		# Randomly choose a type of modification to the model
		r = rand()
		u = rand()

		if r < MOVE && np>0
			# Move a changepoint

			copyto!(boundarypoints,1,boundarypointsₚ,1,np+2)
			# Pick which changepoint to move
			pick = rand(2:np+1)
			# Move the changepoint between its boundaries
			nextValueUniform = pcg32_random_r(&rng) / ( RAND_MAX_U32 / ((boundarypointsₚ[pick+1]-1)-(boundarypointsₚ[pick-1]+1)+1) ) + (boundarypointsₚ[pick-1]+1);
			nextValueGaussian = (uint32_t)(pcg_gaussian_ziggurat(&rng, lastDifference) + boundarypoints[pick]);
			if (abs((int)nextValueUniform-(int)boundarypointsₚ[pick]) > abs((int)nextValueGaussian-(int)boundarypointsₚ[pick]) && nextValueGaussian>boundarypointsₚ[pick-1] && nextValueGaussian<boundarypointsₚ[pick+1] && nextValueGaussian != boundarypointsₚ[pick]){
				# printf("GaussianCloser: %u %u %u ", nextValueUniform, boundarypointsₚ[pick], nextValueGaussian);
				boundarypointsₚ[pick] = nextValueGaussian;
			else
				# printf("UniformCloser: %u %u %u ", nextValueUniform, boundarypointsₚ[pick], nextValueGaussian);
				boundarypointsₚ[pick] = nextValueUniform;
			end
			# Update the model
			update_model(boundarypointsₚ, np, rows, columns, &rng, d, s, m);

			# Calculate log likelihood for proposal
			llP=log_likelihood(d, m, s, rows, columns);

			DEBUG && print("Move: llP-ll = %g - %g\n",llP,ll);}
			if (u<exp(llP-ll)){
				DEBUG && print("Accepted!\n");}
				ll=llP;
				lastDifference = fabs(((double)boundarypointsₚ[pick]-(double)boundarypoints[pick])*2.9);
				# printf("%g\n",lastDifference);
				copyArrayUint(boundarypointsₚ,np+2,boundarypoints);
				for n=1:np
					printf("%u,",boundarypointsₚ[n+1]);
				end
				FORMATTED && print("\n")
				# printf("\n%u\n",np);
			end

		elseif r < MOVE+SIGMA
			# Adjust σ
			copyArray(s,columns,sP);

			# Choose a new sigma
			j=pcg32_random_r(&rng)/(RAND_MAX_U32/columns);
			sP[j]=pcg32_random_r(&rng)/4294967295.0;
			sll=log(pow((sP[j]/s[j]),rows));

			# Calculate log likelihood for proposal
			llP=log_likelihood(d, m, sP, rows, columns);

			DEBUG && println("Sigma $(sP[j]): llP-ll = $llP - $ll")
			if (u<exp(sll+llP-ll)){ # If accepted
				DEBUG && println("Accepted!")
				ll=llP;
				copyArray(sP,columns,s);
				# for (int n=0; n<np+2; n++){
				# 	printf("%u,",boundarypointsₚ[n]);
				# }
				# printf("\n%u\n",np);
				# printf("sigma: ");
				# for (int n=0; n<columns; n++){
				# 	printf("%g,",s[n]);
				# }
				# printf("\n");
			end

		elseif r < MOVE+SIGMA+BIRTH
			# Add a changepoint
			if np < npmax
				copyto!(boundarypoints,1,boundarypointsₚ,1,np+2)

				# Pick which changepoint to add right of
				pick=pcg32_random_r(&rng)/(RAND_MAX_U32/(np+1));
				boundarypointsₚ[np+2] = pcg32_random_r(&rng) / ( RAND_MAX_U32 / ((boundarypointsₚ[pick+1]-1)-(boundarypointsₚ[pick]+1)+1) ) + (boundarypointsₚ[pick]+1);
				npₚ = count_unique!(view(boundarypointsₚ,np+3)) - 2

				# Update the model
				update_model(boundarypointsₚ, npₚ, rows, columns, &rng, d, s, m);

				# Calculate log likelihood for proposal
				lqz=lqxz(d, m, &boundarypointsₚ[pick], rows, columns);
				llP=log_likelihood(d, m, s, rows, columns);
				DEBUG && print("Birth: -lqz+llP-ll = %g + %g - %g\n",-lqz,llP,ll);}
				if u < exp(-lqz+llP-ll)
					DEBUG && print("Accepted!\n");}
					ll=llP;
					np=npₚ;
					copyArrayUint(boundarypointsₚ,np+2,boundarypoints);
					for n=1:np
						printf("%u,",boundarypointsₚ[n+1]);
					end
					FORMATTED && print("\n")
				# printf("\n%u\n",np);
				end
			end
		elseif r < MOVE+SIGMA+BIRTH+DEATH
			# Delete a changepoint
			if np > npmin
				copyto!(boundarypoints,1,boundarypointsₚ,1,np+2)

				# Pick which changepoint to delete
				pick = rand(2:np+1)
				boundarypointsₚ[pick]=boundarypointsₚ[pick+1];
				npₚ = count_unique!(view(boundarypointsₚ,np+2) - 2

				# Update the model
				update_model(boundarypointsₚ, np, rows, columns, &rng, d, s, m);

				# Calculate log likelihood for proposal
				llP=log_likelihood(d, m, s, rows, columns);
				lqz=lqxz(d, m, &boundarypoints[pick], rows, columns);
				DEBUG && println("Death: lqz+llP-ll = $lqz + $llP - $ll")
				if (u<exp(lqz+llP-ll))
					DEBUG && println("Accepted!")
					ll=llP;
					np=npₚ;
					copyArrayUint(boundarypointsₚ,np+2,boundarypoints);
					for n=1:np
						printf("%u,",boundarypointsₚ[n+1]);
					end
					FORMATTED && print("\n")
				# printf("\n%u\n",np);
				end
			end

		else
			# Update the model
			update_model(boundarypointsₚ, np, rows, columns, &rng, d, s, m);

			# Calculate log likelihood for proposal
			llP=log_likelihood(d, m, s, rows, columns);
			if u < exp(llP-ll)
				# printf("Update: llP-ll = %g - %g\n",llP,ll);
				# printf("Accepted!\n");
				ll=llP;

			end
		end
		# printf("%u\n",i);
	end
end

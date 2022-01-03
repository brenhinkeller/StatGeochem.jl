var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = StatGeochem","category":"page"},{"location":"#StatGeochem","page":"Home","title":"StatGeochem","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for StatGeochem.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [StatGeochem]","category":"page"},{"location":"#StatGeochem.Ayers_tsphene-NTuple{10, Any}","page":"Home","title":"StatGeochem.Ayers_tsphene","text":"TC = Ayers_tsphene(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)\n\nCalculate sphene saturation temperature in degrees Celsius Following the sphene saturation calibration of Ayers et al., 2018 (doi: 10.1130/abs/2018AM-320568)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Ayers_tspheneTiO2-NTuple{11, Any}","page":"Home","title":"StatGeochem.Ayers_tspheneTiO2","text":"TiO2Sat = Ayers_tspheneTiO2(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)\n\nCalculate sphene saturation TiO2 concentration (in wt. %) for a given temperature (in C) following the sphene saturation calibration of Ayers et al., 2018 (doi: 10.1130/abs/2018AM-320568)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Boehnke_tzirc-NTuple{11, Any}","page":"Home","title":"StatGeochem.Boehnke_tzirc","text":"T = Boehnke_tzirc(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, Zr)\n\nCalculate zircon saturation temperature in degrees Celsius Following the zircon saturation calibration of Boehnke, Watson, et al., 2013 (doi: 10.1016/j.chemgeo.2013.05.028)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Boehnke_tzircM-NTuple{10, Number}","page":"Home","title":"StatGeochem.Boehnke_tzircM","text":"M = Boehnke_tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)\n\nCalculate zircon saturation M-value based on major element concentrations Following the zircon saturation calibration of Boehnke, Watson, et al., 2013 (doi: 10.1016/j.chemgeo.2013.05.028)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Boehnke_tzircZr-NTuple{11, Any}","page":"Home","title":"StatGeochem.Boehnke_tzircZr","text":"ZrSat = Boehnke_tzircZr(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)\n\nCalculate zircon saturation Zr concentration for a given temperature (in C) Following the zircon saturation calibration of Boehnke, Watson, et al., 2013 (doi: 10.1016/j.chemgeo.2013.05.028)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Harrison_tapatite-Union{Tuple{T}, Tuple{T, T}} where T<:Number","page":"Home","title":"StatGeochem.Harrison_tapatite","text":"TC = Harrison_tapatite(SiO2, P2O5)\n\nCalculate apatite saturation temperature in degrees Celcius following the apatite saturation model of Harrison and Watson 1984 (doi: 10.1016/0016-7037(84)90403-4)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Harrison_tapatiteP-Tuple","page":"Home","title":"StatGeochem.Harrison_tapatiteP","text":"As Harrison_tapatiteP2O5, but returns saturation phosphorus concentration in PPM P\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Harrison_tapatiteP2O5-Union{Tuple{T}, NTuple{6, T}} where T<:Number","page":"Home","title":"StatGeochem.Harrison_tapatiteP2O5","text":"P2O5 = Harrison_tapatiteP2O5(SiO2, Al2O3, CaO, Na2O, K2O, T)\n\nCalculate P2O5 concentration (in wt.%) required for apatite saturation at a given T (in C) following the apatite saturation model of Harrison and Watson 1984 (doi: 10.1016/0016-7037(84)90403-4) with the correction of Bea et al. 1992 (doi: 10.1016/0024-4937(92)90033-U) where applicable\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.LREEmolwt-NTuple{6, Any}","page":"Home","title":"StatGeochem.LREEmolwt","text":"LREEmolwt(La, Ce, Pr, Nd, Sm, Gd)\n\nReturns the average molecular weight of the LREE considered in the REEt value from the monazite saturation model of Montel 1993 (doi: 10.1016/0009-2541(93)90250-M)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.LREEt-Union{Tuple{T}, NTuple{6, T}} where T<:Number","page":"Home","title":"StatGeochem.LREEt","text":"LREEt(La, Ce, Pr, Nd, Sm, Gd)\n\nReturns the sum of the LREE concentrations divided by their respective molar masses. If REE are input in parts per million by weight (ppmw), the result is in units of moles per megagram. This is equivalent to the REEt value from the monazite saturation model of Montel 1993 (doi: 10.1016/0009-2541(93)90250-M)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Montel_tmonazite-NTuple{16, Any}","page":"Home","title":"StatGeochem.Montel_tmonazite","text":"TC = Montel_tmonazite(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O, La, Ce, Pr, Nd, Sm, Gd)\n\nCalculate monazite saturation temperature in degrees Celcius following the monazite saturation model of Montel 1993 (doi: 10.1016/0009-2541(93)90250-M)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Montel_tmonaziteREE-NTuple{11, Any}","page":"Home","title":"StatGeochem.Montel_tmonaziteREE","text":"REEt = Montel_tmonaziteREE(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O, T)\n\nCalculate monazite saturation REEt value (in [ppm/mol.wt.]) for a given temperature (in C) following the monazite saturation model of Montel 1993 (doi: 10.1016/0009-2541(93)90250-M), where:\n\nD = (Na + K + Li + 2Ca) / Al * 1/(Al + Si)) # all as molar cation fractions (not at. %!) ln(REEt) = 9.50 + 2.34D + 0.3879√H2O - 13318/T # H2O as wt.% REEt = Σ REEᵢ(ppm) / at. weight (g/mol)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Rusiecka_tmonaziteREE-Tuple{Any, Any}","page":"Home","title":"StatGeochem.Rusiecka_tmonaziteREE","text":"LREEt = Rusiecka_tmonaziteREE(P_ppm, TC)\n\nCalculate the LREEt (mol/Megagram) value required for monazite saturation at a temperature of TC degrees celcius and P ppmw phosphorous present, following the solubility model of Rusiecka & Baker, 2019 (doi: 10.2138/am-2019-6931)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Rusiecka_txenotimeY-Tuple{Any, Any}","page":"Home","title":"StatGeochem.Rusiecka_txenotimeY","text":"LREEt = Rusiecka_txenotimeY(P_ppm, TC)\n\nCalculate the Y (ppmw) concentration required for xenotime saturation at a temperature of TC degrees celcius and P ppmw phosphorous present, following the solubility model of Rusiecka & Baker, 2019 (doi: 10.2138/am-2019-6931)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Tollari_tapatite-NTuple{9, Any}","page":"Home","title":"StatGeochem.Tollari_tapatite","text":"TC = Tollari_tapatite(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5)\n\nCalculate apatite saturation temperature in degrees Celcius following the apatite saturation model of Tollari et al. 2006 (doi: 10.1016/j.gca.2005.11.024)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.Tollari_tapatiteP2O5-Union{Tuple{T}, Tuple{T, T, T}} where T<:Number","page":"Home","title":"StatGeochem.Tollari_tapatiteP2O5","text":"P2O5 = Tollari_tapatiteP2O5(SiO2, CaO, T)\n\nCalculate P2O5 concentration (in wt.%) required for apatite saturation at a given T (in C) following the apatite saturation model of Tollari et al. 2006 (doi: 10.1016/j.gca.2005.11.024)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.bin_bsr-Tuple{Function, AbstractVector, AbstractVector, Any, Any, Integer}","page":"Home","title":"StatGeochem.bin_bsr","text":"bin_bsr([f!::Function], x::Vector, y::VecOrMat, xmin, xmax, nbins, [w];\n    \tx_sigma = zeros(size(x)),\n    \ty_sigma = zeros(size(y)),\n    \tnresamplings = 1000,\n    \tsem = :sigma,\n    \tp = 0.2\n)\n\nReturns the bincenters c, means or medians m, and uncertainties of the mean or median for a variable y binned by independent variable x into nbins equal bins between xmin and xmax, after nresamplings boostrap resamplings with acceptance probability p.\n\nIf a 2-d array (matrix) of y values is provided, each column will be treated as a separate variable, means and uncertainties will be returned column-wise.\n\nOptional keyword arguments and defaults:\n\nx_sigma = zeros(size(x))\n\nA vector representing the uncertainty (standard deviation) of each x value\n\ny_sigma = zeros(size(y))\n\nA vector representing the uncertainty (standard deviation) of each y value\n\nnresamplings = 1000\n\nThe number of resamplings to conduct\n\nsem = :sigma\n\nFormat of the uncertainty estimate of the distribution of the mean. If :sigma is chosen, a tuple of three vectors (c, m, e) will be returned, where e is the standard error of the mean. If :CI or :pctile is chosen, a tuple of four vectors (c, m, el, eu) will be returned, where el and eu are the lower and upper bounds of the 95% confidence interval.\n\np = 0.2\n\nResampling probabilities, either as a scalar or a vector of the same length as x\n\nExamples:\n\n(c,m,e) = bin_bsr(nanbinmedian!, x, y, 0, 4000, 40, x_sigma=0.05x, p=probability, sem=:sigma)\n\n(c,m,el,eu) = bin_bsr(nanbinmean!, x, y, 0, 4000, 40, x_sigma=0.05x, p=probability, sem=:pctile)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.bin_bsr_ratios-Tuple{Function, AbstractVector, AbstractVector, AbstractVector, Any, Any, Integer}","page":"Home","title":"StatGeochem.bin_bsr_ratios","text":"(c, m, el, eu) = bin_bsr_ratios([f!::Function], x::Vector, num::Vector, denom::Vector, xmin, xmax, nbins, [w];\n    \tx_sigma = zeros(size(x)),\n    \tnum_sigma = zeros(size(num)),\n    \tdenom_sigma = zeros(size(denom)),\n    \tnresamplings = 1000,\n    \tp::Union{Number,Vector} = 0.2\n)\n\nReturns the bincenters c, means m, as well as upper (el) and lower (eu) 95% CIs of the mean for a ratio num/den binned by x into nbins equal bins between xmin and xmax, after nresamplings boostrap resamplings with acceptance probability p.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.bincounts-Tuple{AbstractArray, Number, Number, Integer}","page":"Home","title":"StatGeochem.bincounts","text":"(bincenters, N) = bincounts(x::AbstractArray, xmin::Number, xmax::Number, nbins::Integer)\n\nTally the number of samples that fall into each of nbins equally spaced x bins between xmin and xmax, aligned with bin edges as xmin:(xmax-xmin)/nbins:xmax\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.binmeans-Tuple{AbstractArray, AbstractArray, Number, Number, Integer}","page":"Home","title":"StatGeochem.binmeans","text":"(c,m,e) = binmeans(x, y, xmin, xmax, nbins, [weight]; resamplingratio::Number=1)\n\nThe means (ignoring NaNs) of y values binned by x, into each of nbins equally spaced x bins between xmin and xmax, returning bincenters, means, and standard errors of the mean.\n\nTo calculate binned medians only (without uncertainties), see nanmean\n\nExamples\n\n(c,m,e) = binmeans(x, y, 0, 4000, 40)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.binmedians-Tuple{AbstractArray, AbstractArray, Number, Number, Integer}","page":"Home","title":"StatGeochem.binmedians","text":"binmedians(x::AbstractArray, y::AbstractArray, min::Number, max::Number, nbins::Integer;\n    \tresamplingratio::Number=1\n)\n\nThe medians (ignoring NaNs) of y values binned by x, into each of nbins equally spaced x bins between xmin and xmax, returning bincenters, medians, and equivalent standard errors of the mean (1.4828 * median abolute deviation)\n\nTo calculate binned medians only (without uncertainties), see nanmedian\n\nExamples\n\n(c,m,e) = binmedians(x, y, 0, 4000, 40)\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.bsr!-Tuple{Function, Array, Vector{Int64}, AbstractArray, Number, Number}","page":"Home","title":"StatGeochem.bsr!","text":"bsr!([f::Function], resampled::Array, index::Vector{Int}, data, sigma, p;\n    \trng::AbstractRNG=MersenneTwister()\n)\n\nFill resampled with data boostrap resampled from a (sample-per-row / element-per-column) dataset data with uncertainties sigma and resampling probabilities p, optionally using random numbers generated by f where f is a function of the form f(rng, data[i], sigma[i])\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.bsresample","page":"Home","title":"StatGeochem.bsresample","text":"resampled = bsresample(dataset::Dict, nrows, [elements], [p];\n    \t kernel = gaussian,\n    \t rng = MersenneTwister()\n)\n\nBootstrap resample a dictionary-based dataset with uncertainties stored either in dataset[\"err\"] or dataset[\"[variable]_sigma\"]\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.bsresample-2","page":"Home","title":"StatGeochem.bsresample","text":"resampled = bsresample(data::AbstractArray, sigma, nrows, [p];\n    \t kernel = gaussian,\n    \t rng = MersenneTwister(),\n    \t return_index = false\n)\n\nBootstrap resample a (sample-per-row / element-per-column) array of data with uncertainties sigma and resampling probabilities p\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.changepoint-Tuple{AbstractArray, Integer}","page":"Home","title":"StatGeochem.changepoint","text":"changepoint(data, [sigma], nsteps; np, npmin, npmax)\n\nGiven an ordered array of data points, optionally with uncertainties sigma, use a Markov chain Monte Carlo approach based on that of Gallagher et al., 2010 (10.1016/j.epsl.2011.09.015) to estimate the position (by index) and optionally number of changepoints that best explain the data. Will return the results for nsteps steps of the Markov chain.\n\nOptionally, you may also specify as keyword arguments either np or npmin and npmax to constrain the number of allowed changepoints.\n\nCurrently prints results to the terminal (stdout), one line per step of the Markov chain, but a more efficent output format is probably desirable in a future version of this function.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.concatenatedatasets-Tuple{AbstractDict, AbstractDict}","page":"Home","title":"StatGeochem.concatenatedatasets","text":"concatenatedatasets(d1::AbstractDict, d2::AbstractDict)\n\nVertically concatenate two Dict-based datasets, variable-by-variable\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.delim_string_function-Tuple{Function, AbstractString, Char, Type}","page":"Home","title":"StatGeochem.delim_string_function","text":"delim_string_function(f, str, delim, T;\n    \tmerge::Bool=false,\n\nParse a delimited string str with delimiter delim into substrings that will then be operated upon by function f. The results of f will be returned in an array with eltype T.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.delim_string_parse","page":"Home","title":"StatGeochem.delim_string_parse","text":"delim_string_parse(str, delim, T;\n    \tmerge::Bool=false,\n    \tundefval=NaN)\n\nParse a delimited string str with delimiter delim into values of type T and return the answers as an array with eltype T\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.delim_string_parse!","page":"Home","title":"StatGeochem.delim_string_parse!","text":"delim_string_parse!(result, str, delim, [T];\n    \toffset::Integer=0,\n    \tmerge::Bool=false,\n    \tundefval=NaN)\n\nParse a delimited string str with delimiter delim into values of type T and return the answers in a pre-allocated result array provided as input.\n\nIf T is not specified, it the eltype of the result array will be used by default.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.digitize_plotline-NTuple{4, Any}","page":"Home","title":"StatGeochem.digitize_plotline","text":"(x,y) = digitize_plotline(img, line_color, xlims, ylims; atol=0.16)\n\nCalculate approximate x and y positions for a colored line in an image\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.digitize_plotmarkers-NTuple{4, Any}","page":"Home","title":"StatGeochem.digitize_plotmarkers","text":"(x,dx,y,dy) = digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)\n\nCalculate approximate x and y positions and uncertainties for colored markers in an image\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.elementify","page":"Home","title":"StatGeochem.elementify","text":"elementify(data::Array, elements::Array=data[1,:];\n    \timportas=:Dict,\n    \tstandardize::Bool=true,\n    \tfloattype=Float64,\n    \tskipstart::Integer=1,\n    \tskipnameless::Bool=true\n)\n\nConvert a flat array data into a dictionary (importas=:Dict) or named tuple (importas=:Tuple) with each column as a variable. Tuples are substantially more efficient, so should be favored where possible.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.eustar-NTuple{4, Number}","page":"Home","title":"StatGeochem.eustar","text":"eustar(Nd::Number, Sm::Number, Gd::Number, Tb::Number)\n\nCalculate expected europium concentration, Eu*, based on abundance of adjacent rare earths. Full four-element log-linear interpolation, using ionic radii\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.eustar-Tuple{Number, Number}","page":"Home","title":"StatGeochem.eustar","text":"eustar(Sm::Number, Gd::Number)\n\nCalculate expected europium concentration, Eu*, based on abundance of adjacent rare earths. Simple geometric mean interpolation from Sm and Gd alone\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.exportdataset-Tuple{Union{Dict, NamedTuple}, AbstractString, AbstractChar}","page":"Home","title":"StatGeochem.exportdataset","text":"exportdataset(dataset, [elements], filepath, delim;\n    \tfloatout::Bool=false,\n    \tfindnumeric::Bool=false,\n    \tskipnan::Bool=true,\n    \tdigits::Integer,\n    \tsigdigits::Integer\n    \trows=:\n)\n\nConvert a dict or named tuple of vectors into a 2-D array with variables as columns Export a dataset (in the form of either a Dict or a NamedTuple), optionally specifying which elements to export, as a delimited ASCII text file with the name specified by filepath and delimiter delim.\n\nPossible keyword arguments include:\n\n\tdigits\n\tsigdigits\n\nSpecify a number of absolute or significant digits to which to round the printed output. Default is no rounding.\n\n\tskipnan\n\nLeave NaNs as empty cells in the delimited output file. Boolean; true by default.\n\n\tfloatout\n\nForce all output to be represented as a floating-point number, or else NaN. Boolean; false by default.\n\n\tfindnumeric\n\nExport only numeric columns. Boolean; false by default.\n\n\trows\n\nspecify which rows of the dataset to export. Default : exports all rows.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.feoconversion","page":"Home","title":"StatGeochem.feoconversion","text":"feoconversion(FeO::Number=NaN, Fe2O3::Number=NaN, FeOT::Number=NaN, Fe2O3T::Number=NaN)\n\nCompiles data from FeO, Fe2O3, FeOT, and Fe2O3T into  a single FeOT value.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.find_tc1_age-Tuple{Number, Number}","page":"Home","title":"StatGeochem.find_tc1_age","text":"find_tc1_age(lat::Number,lon::Number)\n\nReturn a tuple (age, age_min, age_max) containing the nominal, upper, and lower tc1 age bounds for the 1x1 arc degree grid cell containing lat and lon\n\nfind_tc1_age(lat::AbstractArray,lon::AbstractArray)\n\nReturn a tuple (age, age_min, age_max) where age, age_min, and age_max are arrays containing the nominal, upper and lower tc1 age bounds for each location pair lat[i], lon[i]\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.find_tc1_crust-Tuple{Number, Number}","page":"Home","title":"StatGeochem.find_tc1_crust","text":"find_tc1_crust(lat::Number,lon::Number)\n\nFind the depth to the 550C isotherm for the 1x1 arc degree grid cell containing lat and lon\n\nfind_tc1_crust(lat::AbstractArray,lon::AbstractArray)\n\nFor each pair of latitudes and longitudes given by lat and lon, find the depth to the 550C isotherm for the 1x1 arc degree grid cell containing lat[i] and lon[i]\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.find_tc1_lith-Tuple{Number, Number}","page":"Home","title":"StatGeochem.find_tc1_lith","text":"find_tc1_lith(lat::Number,lon::Number)\n\nFind the depth to the 1300C isotherm for the 1x1 arc degree grid cell containing lat and lon\n\nfind_tc1_lith(lat::AbstractArray,lon::AbstractArray)\n\nFor each pair of latitudes and longitudes given by lat and lon, find the depth to the 1300C isotherm for the 1x1 arc degree grid cell containing lat[i] and lon[i]\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.floatify","page":"Home","title":"StatGeochem.floatify","text":"floatify(x, T::Type=Float64)\n\nConvert x to a floating-point number (default Float64) by any means necessary\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.importdataset-Tuple{AbstractString, AbstractChar}","page":"Home","title":"StatGeochem.importdataset","text":"function importdataset(filepath, delim;\n    \timportas=:Dict,\n    \tstandardize::Bool=true,\n    \tfloattype=Float64,\n    \tskipstart::Integer=0,\n    \tskipnameless::Bool=true,\n    \tmindefinedcolumns::Integer=0\n)\n\nImport a delimited file specified by filepath with delimiter delim as a dataset in the form of either a Dict or a NamedTuple.\n\nPossible keyword arguments include:\n\n\timportas\n\nSpecify the format of the imported dataset. Options include :Dict and :Tuple\n\n\tstandardize\n\nConvert columns to uniform type wherever possible. Boolean; true by default.\n\n\tfloattype\n\nPreferred floating-point type for numerical data. Float64 by default.\n\n\tskipstart\n\nIgnore this many rows at the start of the input file (useful if input file has a header or other text before the column names). 0 by default.\n\n\tskipnameless\n\nSkip columns with no column name. Boolean; true by default\n\n\tmindefinedcolumns\n\nSkip rows with fewer than this number of delimiters. 0 by default.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.invweight-Tuple{AbstractArray, AbstractArray, AbstractArray}","page":"Home","title":"StatGeochem.invweight","text":"k = invweight(lat::AbstractArray, lon::AbstractArray, age::AbstractArray;\n    \tlp::Number=2,\n    \tspatialscale=1.8,\n    \tagescale=38.0\n)\n\nFind the inverse weights k (proportional to spatiotemporal sample density) for a set of geological samples with specified latitude (lat), logitude (lon), and age (of crystallization, deposition, etc.).\n\nThe default spatialscale and agescale are taken from Keller and Schoene 2012. However, alternative scalings can be supplied. If an array is supplied for either spatialscale, agescale, or both, a 3-d matrix of k values will be returned, with dimensions length(spatialscale)length(agescale)nrows.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.invweight-Tuple{AbstractArray, Number}","page":"Home","title":"StatGeochem.invweight","text":"k = invweight(nums::AbstractArray, scale::Number; lp=2)\n\nFind the inverse weights for a single array nums for a given scale, and exponent lp (default lp = 2).\n\nReturns an array k where k[i] is the \"inverse weight\" for element i of the input array.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.invweight_age-Tuple{AbstractArray}","page":"Home","title":"StatGeochem.invweight_age","text":"k = invweight_age(age::AbstractArray; lp::Number=2, agescale::Number=38.0)\n\nFind the inverse weights k (proportional to temporal sample density) for a set of geological samples with specified age (of crystallization, deposition, etc.).\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.invweight_location-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"StatGeochem.invweight_location","text":"k = invweight_location(lat::AbstractArray, lon::AbstractArray;\n    \tlp::Number=2,\n    \tspatialscale::Number=1.8\n)\n\nFind the inverse weights k (proportional to spatial sample density) for a set of geological samples with specified latitude (lat), and logitude (lon).\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.melts_clean_modes-Tuple{String}","page":"Home","title":"StatGeochem.melts_clean_modes","text":"melts_clean_modes(scratchdir::String; index=1)\n\nRead and parse / clean-up modal phase proportions from specified MELTS run directory Returns an elementified dictionary\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.melts_configure","page":"Home","title":"StatGeochem.melts_configure","text":"melts_configure(meltspath::String, scratchdir::String, composition::Array{Float64},\n    \telements::Array,\n    \tT_range::Array=[1400, 600],\n    \tP_range::Array=[10000,10000];)\n\nConfigure and run a MELTS simulation using alphaMELTS. Optional keyword arguments and defaults:\n\nbatchstring::String=\"1 sc.melts 10 1 3 1 liquid 1 1.0 0 10 0 4 0 \"\n\ndT = -10\n\ndP = 0\n\nindex = 1\n\nversion = \"pMELTS\"\n\nmode = \"isobaric\"\n\nfo2path = \"FMQ\" Oxygen fugacity buffer to follow, e.g., FMQ or NNO+1\n\nfractionatesolids::Bool = false Fractionate all solids\n\nsuppress::Array{String} = [] Supress individual phases (specify as strings in array, i.e. [\"leucite\"])\n\nverbose::Bool = true Print verbose MELTS output to terminal (else, write it to melts.log)\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.melts_query-Tuple{String}","page":"Home","title":"StatGeochem.melts_query","text":"melts_query_modes(scratchdir::String; index=1)\n\nRead all phase proportions from Phase_main_tbl.txt in specified MELTS run directory Returns an elementified dictionary\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.melts_query_liquid-Tuple{String}","page":"Home","title":"StatGeochem.melts_query_liquid","text":"melts_query_liquid(scratchdir::String; index=1)\n\nRead liquid composition from Liquid_comp_tbl.txt in specified MELTS run directory Returns an elementified dictionary\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.melts_query_modes-Tuple{String}","page":"Home","title":"StatGeochem.melts_query_modes","text":"melts_query_modes(scratchdir::String; index=1)\n\nRead modal phase proportions from Phase_mass_tbl.txt in specified MELTS run Returns an elementified dictionary\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.melts_query_solid-Tuple{String}","page":"Home","title":"StatGeochem.melts_query_solid","text":"melts_query_solid(scratchdir::String; index=1)\n\nRead solid composition from Solid_comp_tbl.txt in specified MELTS run directory Returns an elementified dictionary\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.nonnumeric-Tuple{Any}","page":"Home","title":"StatGeochem.nonnumeric","text":"nonnumeric(x)\n\nReturn true for if x is not missing but cannot be parsed as a number\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.oxideconversion!-Tuple{Dict}","page":"Home","title":"StatGeochem.oxideconversion!","text":"dataset = oxideconversion!(dataset::Dict; unitratio::Number=10000)\n\nConvert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3, ...) in place.\n\nIf metals are as PPM, set unitratio=10000 (default); if metals are as wt%, set unitratio = 1\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.oxideconversion-Tuple{Dict}","page":"Home","title":"StatGeochem.oxideconversion","text":"dataset = oxideconversion(dataset::Dict; unitratio::Number=10000)\n\nConvert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3, ...).\n\nIf metals are as PPM, set unitratio=10000 (default); if metals are as wt%, set unitratio = 1\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.parsedlm","page":"Home","title":"StatGeochem.parsedlm","text":"parsedlm(str::AbstractString, delimiter::Char, T::Type=Float64; rowdelimiter::Char='\\n')\n\nParse a string delimited by both row and column into a single (2-D) matrix. Default column delimiter is newline. Similar to readdlm, but operating on a string instead of a file.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.perplex_configure_geotherm","page":"Home","title":"StatGeochem.perplex_configure_geotherm","text":"perplex_configure_geotherm(perplexdir::String, scratchdir::String, composition::Array{<:Number},\n    \telements::String=[\"SIO2\",\"TIO2\",\"AL2O3\",\"FEO\",\"MGO\",\"CAO\",\"NA2O\",\"K2O\",\"H2O\"],\n    \tP_range::Array{<:Number}=[280,28000], T_surf::Number=273.15, geotherm::Number=0.1;\n    \tdataset::String=\"hp02ver.dat\",\n    \tindex::Integer=1,\n    \tnpoints::Integer=100,\n    \tsolution_phases::String=\"O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n\",\n    \texcludes::String=\"ts\\nparg\\ngl\\nged\\nfanth\\ng\\n\",\n    \tmode_basis::String=\"vol\",  #[\"vol\", \"wt\", \"mol\"]\n    \tcomposition_basis::String=\"wt\",  #[\"vol\", \"wt\", \"mol\"]\n    \tfluid_eos::Integer=5)\n\nSet up a PerpleX calculation for a single bulk composition along a specified geothermal gradient and pressure (depth) range. P specified in bar and T_surf in Kelvin, with geothermal gradient in units of Kelvin/bar\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.perplex_configure_isobar","page":"Home","title":"StatGeochem.perplex_configure_isobar","text":"perplex_configure_isobar(perplexdir::String, scratchdir::String, composition::Array{<:Number},\n    \telements::String=[\"SIO2\",\"TIO2\",\"AL2O3\",\"FEO\",\"MGO\",\"CAO\",\"NA2O\",\"K2O\",\"H2O\"]\n    \tP::Number=10000, T::Array{<:Number}=[500+273.15, 1500+273.15];\n    \tdataset::String=\"hp11ver.dat\",\n    \tindex::Integer=1,\n    \tnpoints::Integer=100,\n    \tsolution_phases::String=\"O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n\",\n    \texcludes::String=\"ts\\nparg\\ngl\\nged\\nfanth\\ng\\n\",\n    \tmode_basis::String=\"vol\",  #[\"vol\", \"wt\", \"mol\"]\n    \tcomposition_basis::String=\"wt\",  #[\"vol\", \"wt\", \"mol\"]\n    \tfluid_eos::Integer=5)\n\nSet up a PerpleX calculation for a single bulk composition along a specified isobaric temperature gradient. P specified in bar and T_range in Kelvin\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.perplex_configure_pseudosection","page":"Home","title":"StatGeochem.perplex_configure_pseudosection","text":"perplex_configure_pseudosection(perplexdir::String, scratchdir::String, composition::Array{<:Number},\n    \telements::Array{String}=[\"SIO2\",\"TIO2\",\"AL2O3\",\"FEO\",\"MGO\",\"CAO\",\"NA2O\",\"K2O\",\"H2O\"],\n    \tP::Array{<:Number}=[280, 28000], T::Array{<:Number}=[273.15, 1500+273.15];\n    \tdataset::String=\"hp11ver.dat\",\n    \tindex::Integer=1,\n    \txnodes::Integer=42,\n    \tynodes::Integer=42,\n    \tsolution_phases::String=\"O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n\",\n    \texcludes::String=\"ts\\nparg\\ngl\\nged\\nfanth\\ng\\n\",\n    \tmode_basis::String=\"vol\", #[\"vol\", \"wt\", \"mol\"]\n    \tcomposition_basis::String=\"wt\", #[\"wt\", \"mol\"]\n    \tfluid_eos::Number=5)\n\nSet up a PerpleX calculation for a single bulk composition across an entire 2d P-T space. P specified in bar and T in Kelvin.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.perplex_query_modes-Tuple{String, String, Array{<:Number}, Array{<:Number}}","page":"Home","title":"StatGeochem.perplex_query_modes","text":"perplex_query_modes(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};\n    \tindex::Integer=1, npoints::Integer=200, include_fluid=\"y\")\n\nQuery modal mineralogy (mass proportions) along a specified P-T path using a pre-computed pseudosection. Results are returned as a dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_modes-Tuple{String, String}","page":"Home","title":"StatGeochem.perplex_query_modes","text":"perplex_query_modes(perplexdir::String, scratchdir::String;\n    \tdof::Integer=1, index::Integer=1, include_fluid=\"y\")\n\nQuery modal mineralogy (mass proportions) along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d grid / pseudosection (dof=2). Results are returned as a dictionary.\n\nCurrently returns vol %\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_phase-Tuple{String, String, String, Array{<:Number}, Array{<:Number}}","page":"Home","title":"StatGeochem.perplex_query_phase","text":"perplex_query_phase(perplexdir::String, scratchdir::String, phase::String, P::Array{<:Number}, T::Array{<:Number};\n    \tindex::Integer=1, npoints::Integer=200, include_fluid=\"y\", clean_units::Bool=true)\n\nQuery all perplex-calculated properties for a specified phase (e.g. \"Melt(G)\") along a specified P-T path using a pre-computed pseudosection. Results are returned as a dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_phase-Tuple{String, String, String}","page":"Home","title":"StatGeochem.perplex_query_phase","text":"perplex_query_phase(perplexdir::String, scratchdir::String, phase::String;\n    \tdof::Integer=1, index::Integer=1, include_fluid=\"y\", clean_units::Bool=true)\n\nQuery all perplex-calculated properties for a specified phase (e.g. \"Melt(G)\") along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d grid / pseudosection (dof=2). Results are returned as a dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_point-Tuple{String, String, Number, Number}","page":"Home","title":"StatGeochem.perplex_query_point","text":"perplex_query_point(perplexdir::String, scratchdir::String, P::Number, T::Number; index::Integer=1)\n\nQuery perplex results at a single P,T point in a pseudosection. Results are returned as a string.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_point-Tuple{String, String, Number}","page":"Home","title":"StatGeochem.perplex_query_point","text":"perplex_query_point(perplexdir::String, scratchdir::String, indvar::Number; index::Integer=1)\n\nQuery perplex results at a single temperature on an isobar or single pressure on a geotherm. Results are returned as a string.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_seismic-Tuple{String, String, Array{<:Number}, Array{<:Number}}","page":"Home","title":"StatGeochem.perplex_query_seismic","text":"perplex_query_seismic(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};\n    \tindex::Integer=1, npoints::Integer=200, include_fluid=\"n\")\n\nQuery perplex seismic results along a specified P-T path using a pre-computed pseudosection. Results are returned as a dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_seismic-Tuple{String, String}","page":"Home","title":"StatGeochem.perplex_query_seismic","text":"perplex_query_seismic(perplexdir::String, scratchdir::String;\n    \tdof::Integer=1, index::Integer=1, include_fluid=\"n\")\n\nQuery perplex seismic results along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d grid / pseudosection (dof=2). Results are returned as a dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_system-Tuple{String, String, Array{<:Number}, Array{<:Number}}","page":"Home","title":"StatGeochem.perplex_query_system","text":"function perplex_query_system(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};\n    \tindex::Integer=1, npoints::Integer=200, include_fluid=\"y\",clean_units::Bool=true)\n\nQuery all perplex-calculated properties for the system (with or without fluid) along a specified P-T path using a pre-computed pseudosection. Results are returned as a dictionary. Set include_fluid=\"n\" to return solid+melt only.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.perplex_query_system-Tuple{String, String}","page":"Home","title":"StatGeochem.perplex_query_system","text":"perplex_query_system(perplexdir::String, scratchdir::String;\n    \tindex::Integer=1, include_fluid=\"y\", clean_units::Bool=true)\n\n?\n\nQuery all perplex-calculated properties for the system (with or without fluid) along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d grid / pseudosection (dof=2). Results are returned as a dictionary. Set include_fluid=\"n\" to return solid+melt only.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.plausiblynumeric-Tuple{Any}","page":"Home","title":"StatGeochem.plausiblynumeric","text":"plausiblynumeric(x)\n\nReturn true if x can be parsed as a number, else false\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.randsample","page":"Home","title":"StatGeochem.randsample","text":"randsample(dataset::Dict, nrows, [elements], [p])\n\nBootstrap resample (without uncertainty) a dataset dict to length nrows. Optionally provide weights p either as a vector (one-weight-per-sample) or scalar.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.randsample-2","page":"Home","title":"StatGeochem.randsample","text":"randsample(data::Array, nrows, [p])\n\nBootstrap resample (without uncertainty) a data array to length nrows. Optionally provide weights p either as a vector (one-weight-per-sample) or scalar.\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.renormalize!","page":"Home","title":"StatGeochem.renormalize!","text":"renormalize!(dataset, [elements]; total=1.0)\n\nNormalize in-place a (i.e., compositional) dataset defined by a Dict or NamedTuple of one-dimensional numerical arrays, such that all the elements (i.e., variables – by default all keys in the datset) sum to a given total (by default, 1.0).\n\nNote that the arrays representing each element or variable are assumed to be of uniform length\n\n\n\n\n\n","category":"function"},{"location":"#StatGeochem.renormalize!-Tuple{AbstractArray}","page":"Home","title":"StatGeochem.renormalize!","text":"renormalize!(A::AbstractArray; dim, total=1.0)\n\nNormalize an array A in place such that it sums to total. Optionally may specify a dimension dim along which to normalize.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.system-Tuple{AbstractString}","page":"Home","title":"StatGeochem.system","text":"system(cmdstr::AbstractString)\n\nDirect access to the command line through C's system function – without stripping/sanitizing special characters,  in contrast to Julia's safer run() function. This allows pipelining, etc. in shell commands. Returns 0 on success.\n\n\n\n\n\n","category":"method"},{"location":"#StatGeochem.unelementify","page":"Home","title":"StatGeochem.unelementify","text":"unelementify(dataset, elements;\n    \tfloatout::Bool=false,\n    \tfloattype=Float64,\n    \tfindnumeric::Bool=false,\n    \tskipnan::Bool=false,\n    \trows=:\n)\n\nConvert a dict or named tuple of vectors into a 2-D array with variables as columns\n\n\n\n\n\n","category":"function"}]
}

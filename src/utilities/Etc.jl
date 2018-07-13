## --- Two-sided linear exponential functions

    # Two-sided linear exponential distribution joined by an atan sigmoid.
    function doubleLinearExponential(x,p)
        # If to a normal-esque PDF, parameters p roughly correspond to:
        # p[1] = pre-exponential (normaliation constant)
        # p[2] = mean (central moment)
        # p[3] = standard deviation
        # p[4] = sharpness
        # p[5] = skew
        xs = (x-p[2])./p[3]; # X scaled by mean and variance
        v = 1/2-atan.(xs)/pi; # Sigmoid (positive on LHS)
        f = p[1] .* exp.((p[4].^2).*(p[5].^2).*xs.*v - (p[4].^2)./(p[5].^2).*xs.*(1-v));
        return f
    end
    export doubleLinearExponential

    # Log of two-sided linear exponential distribution joined by an atan sigmoid.
    function doubleLinearExponentialLL(x,p)
        # If to a normal-esque PDF, parameters p roughly correspond to:
        # p[1] = pre-exponential (normaliation constant)
        # p[2] = mean (central moment)
        # p[3] = standard deviation
        # p[4] = sharpness
        # p[5] = skew
        xs = (x-p[2,:])./p[3,:]; # X scaled by mean and variance
        v = 1/2-atan.(xs)/pi; # Sigmoid (positive on LHS)
        f = log.(p[1,:]) + (p[4,:].^2).*(p[5,:].^2).*xs.*v - (p[4,:].^2)./(p[5,:].^2).*xs.*(1-v);
        return f
    end
    export doubleLinearExponentialLL

## --- Some useful distributions

    UniformDistribution = [1.0, 1.0];
    export UniformDistribution

    TriangularDistribution = [2.0,1.0,0.0];
    export TriangularDistribution

    TruncatedNormalDistribution = [0.241971, 0.251742, 0.261481, 0.271153, 0.280724, 0.29016, 0.299423, 0.308478, 0.317288, 0.325817, 0.334031, 0.341892, 0.349368, 0.356425, 0.363032, 0.369157, 0.374774, 0.379856, 0.384378, 0.38832, 0.391662, 0.394389, 0.396487, 0.397946, 0.398759, 0.398922, 0.398434, 0.397297, 0.395518, 0.393104, 0.390067, 0.386423, 0.382188, 0.377383, 0.372031, 0.366156, 0.359787, 0.352951, 0.345681, 0.338008, 0.329966, 0.32159, 0.312916, 0.303978, 0.294815, 0.285461, 0.275953, 0.266327, 0.256617, 0.246858, 0.237083, 0.227324, 0.21761, 0.207972, 0.198437, 0.18903, 0.179775, 0.170694, 0.161808, 0.153134, 0.144689, 0.136486, 0.128538, 0.120856, 0.113448, 0.10632, 0.0994771, 0.092923, 0.0866592, 0.0806857, 0.0750015, 0.069604, 0.0644895, 0.0596534, 0.05509, 0.0507926, 0.0467541, 0.0429665, 0.0394214, 0.0361097, 0.0330223, 0.0301496, 0.0274819, 0.0250094, 0.0227222, 0.0206105, 0.0186646, 0.0168748, 0.0152318, 0.0137263, 0.0123494, 0.0110926, 0.00994734, 0.00890582, 0.00796034, 0.00710363, 0.00632878, 0.00562925, 0.00499887, 0.00443185];
    export TruncatedNormalDistribution

    HalfNormalDistribution = reverse([0.00308455799258221,0.00344513787810736,0.00384359593870402,0.00428337694406554,0.00476817640292968,0.00530195190886020,0.00588893424076672,0.00653363811239984,0.00724087245603857,0.00801575011663268,0.00886369682387602,0.00979045930117787,0.0108021123623887,0.0119050648395517,0.0131060641780270,0.0144121995292184,0.0158309031659599,0.0173699500415578,0.0190374553106744,0.0208418696288452,0.0227919720475949,0.0248968603240084,0.0271659384673712,0.0296089013512526,0.0322357162272980,0.0350566009871371,0.0380819990313005,0.0413225506189564,0.0447890605896858,0.0484924623684547,0.0524437781874190,0.0566540754832023,0.0611344194557710,0.0658958218049086,0.0709491856924629,0.0763052470128359,0.0819745120904443,0.0879671919608544,0.0942931334317431,0.100961747160455,0.107981933026376,0.115362003118266,0.123109602698694,0.131231629549353,0.139734152141830,0.148622327117986,0.157900316601788,0.167571205899929,0.177636922181184,0.188098154753774,0.198954277585497,0.210203274732500,0.221841669358911,0.233864457040601,0.246265044051699,0.259035191331784,0.272164964824556,0.285642692865021,0.299454931271490,0.313586436771017,0.328020149351987,0.342737184095615,0.357716832989081,0.372936577167131,0.388372109966426,0.403997371108118,0.419784592249449,0.435704354065101,0.451725654934249,0.467815991220326,0.483941449038287,0.500066807309317,0.516155651806475,0.532170499797510,0.548072934794057,0.563823750820605,0.579383105522966,0.594710681345605,0.609765853921034,0.624507866733523,0.638896011044705,0.652889810001052,0.666449205783599,0.679534748609516,0.692107786353846,0.704130653528599,0.715566858335938,0.726381266502836,0.736540280606647,0.746012013614657,0.754766455385986,0.762775630921048,0.770013749192028,0.776457341447094,0.782085387950912,0.786879432203880,0.790823681771635,0.793905094954024,0.796113452627904,0.797441414709954,0.797884560802865]);
    export HalfNormalDistribution

    EllisDistribution = reverse([0.0704953225358788,0.0729834190955530,0.0758052035283703,0.0789774877136435,0.0825170835306852,0.0864408028588081,0.0907654575773249,0.0955078595655482,0.100684820702791,0.106313153211279,0.112445604872168,0.119269623079656,0.127004210008982,0.135868367835387,0.146081098734114,0.157861404880401,0.171428288449491,0.187000751616623,0.204797796557039,0.225038425445979,0.247941640458685,0.273726443770397,0.302611837556356,0.334816823991803,0.370560405251978,0.410061583512122,0.453537559207078,0.501069896764001,0.552512350023633,0.607696867685180,0.666455398447850,0.728619891010849,0.794022294073386,0.862494556334665,0.933868626493896,1.00797645325028,1.08464998530304,1.16372117135136,1.24517058719607,1.33534668147077,1.44922848391669,1.60240934801522,1.81048262724771,2.08904167509554,2.45367984504008,2.91457155069222,3.45095399773449,4.03114563097533,4.68723738141821,5.64309998743355,7.15849979749155]);
    export EllisDistribution

    MeltsZirconDistribution = reverse([0.000250810042630289,0.209249326535399,0.418247843028168,0.627246359520935,0.836373869090645,1.04374357155004,1.24322922919823,1.43183209158837,1.61129844550338,1.77688486888092,1.91661079082960,2.02685847965723,2.09922776245012,2.14214951000388,2.16049597963371,2.15970102681196,2.14446347280944,2.11932089388414,2.08857186945567,2.04960073596920,2.01041078899474,1.96673182187778,1.92474422194058,1.87915814526374,1.83509405646190,1.79113931931073,1.74610988807671,1.70247028643145,1.66040197025845,1.61746239633646,1.57521575240718,1.53439000675662,1.49515411445647,1.45616189914587,1.41776712500444,1.37956742674343,1.34262162411021,1.30694333629830,1.27238055574041,1.23881668637306,1.20536150622499,1.17185978304267,1.13930634959523,1.10770181433460,1.07702403102922,1.04755288935075,1.01891846438717,0.991307020548356,0.964597559566676,0.937935447447592,0.912547254896304,0.888432981912822,0.864397309803580,0.841090957191169,0.818607007784804,0.796806321622935,0.775811959280348,0.755412283398949,0.735655522179200,0.716614802353539,0.698283859993860,0.680538233670485,0.663358513619274,0.646723719544709,0.630552743575445,0.615079934324004,0.599677772774101,0.585285935569870,0.570974188360921,0.557289930703817,0.544044974206469,0.531081757193001,0.518528309854683,0.506357798936355,0.494356584164449,0.482751168373747,0.471302277398808,0.460045602813051,0.449055447325410,0.438327785548115,0.427997270253091,0.417817958928498,0.407979196700429,0.398511598892507,0.389449671179037,0.380749356382449,0.372391900054805,0.364322271587003,0.356549119619268,0.349118200705518,0.341882681922588,0.334822491050461,0.327796643262725,0.320830884764955,0.313845271864816,0.307017155215140,0.300244336445030,0.293473757055415,0.286703177665800,0.279932598276178]);
    export MeltsZirconDistribution

    MeltsVolcanicZirconDistribution = reverse([0.000126311537071135,0.138678719074772,0.277231126612473,0.415780749218191,0.554244033653465,0.692478823982408,0.828041603755384,0.958465872031531,1.08363596175037,1.20336837717604,1.31086441368344,1.39918689258733,1.46297882113932,1.50876500985997,1.53936076666456,1.55766152364649,1.56721340549175,1.57055368507083,1.57017608482977,1.56375900729345,1.55673917794689,1.54528204105819,1.53446510415179,1.51888561033322,1.50413627708598,1.48878532213817,1.47093626533172,1.45367531631901,1.43722401954538,1.41884953441590,1.40020284908090,1.38217173159262,1.36491496000665,1.34781339259003,1.32918673600548,1.31064062341720,1.29264559362098,1.27506639070296,1.25789091287708,1.24124868413795,1.22415364593791,1.20619847793269,1.18871940548722,1.17166291592868,1.15507107892548,1.13905814837756,1.12329852032616,1.10796490309432,1.09296239013487,1.07797063027672,1.06271754003706,1.04720311991123,1.03171712350443,1.01658666090717,1.00179914992514,0.987288063694791,0.973203816081034,0.959396520209570,0.945904001861541,0.932723654268405,0.919949056056459,0.907472933652423,0.895233778569715,0.883320760619869,0.871740160686158,0.860440380738532,0.849189937637448,0.836826698333265,0.824476857182168,0.812345940043929,0.800432420550727,0.788729129991484,0.777279410925614,0.766113420680074,0.755111910752756,0.744481086182191,0.734045338222697,0.723854066251978,0.713944076764060,0.704174797556201,0.694697443892042,0.685378323765953,0.676278206907816,0.667385068616828,0.658707640374667,0.650208560617337,0.641977872085481,0.633875778399037,0.625836373705725,0.617898448430796,0.610049502743681,0.602380243936629,0.594892672217814,0.587607589519472,0.580506613961880,0.573451694478375,0.566417381408414,0.559389964999820,0.552362549285121,0.545335133570427]);
    export MeltsVolcanicZirconDistribution


## --- Parse a delimited string

    function parseDelimString!(parsed::Array, str::String, delim::Char, parseType::Type; offset=0)

        # Ignore initial delimiter
        last_delim_pos=0;
        if str[1]==delim
            last_delim_pos=1;
        end

        # Cycle through string parsing text betweeen delims
        delim_pos=0;
        n = offset;
        for i=1:length(str)
            if str[i] == delim
                delim_pos = i;
                    if delim_pos>last_delim_pos
                    n += 1;
                    if delim_pos>last_delim_pos+1
                        parsed[n] = parse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                    end
                    last_delim_pos = delim_pos;
                end
            end
        end

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1;
            parsed[n] = parse(parseType, str[(last_delim_pos+1):length(str)])
        end

        # Return the number of parsed values
        return n-offset
    end
    export parseDelimString

## --- Read Arc ASCII grid files

    function parseAAIGrid(fname, parseType)
        f = open(fname);

            metadata = Dict();
            metadata["ncols"] = parse(match(r"  *(.*?)$", readline(f))[1])
            metadata["nrows"] = parse(match(r"  *(.*?)$", readline(f))[1])
            metadata["xll_corner"] = parse(match(r"  *(.*?)$", readline(f))[1])
            metadata["yll_corner"] = parse(match(r"  *(.*?)$", readline(f))[1])
            metadata["cellsize"] = parse(match(r"  *(.*?)$", readline(f))[1])
            metadata["nodata"] = parse(match(r"  *(.*?)$", readline(f))[1])

            nrows = metadata["nrows"];
            ncols = metadata["ncols"];

            data = Array{Int16}(ncols,nrows)
            for i=1:nrows
                l = readline(f);
                parseDelimString!(data, l, ' ', Int16, offset=(i-1)*ncols);
            end

        close(f)

        return (data', metadata)
    end
    export parseAAIGrid

## --- Calculate slope from a DEM

    function quickMaxSlopeEarth(matrix, x_lon_cntr, y_lat_cntr, cellsize; minmatval=-500)
        # Returns slope in units/kilometer given a latitude-longitude grid of z-values

        # Allocate output array
        slope = Array{UInt16}(size(matrix))

        # Average size of a degree on Earth
        km_per_lat = 111.1;

        # Fill in the center first
        distNS = 2*cellsize*km_per_lat;
        for i=2:size(matrix,1)-1
            # Distance between grid cell centers
            km_per_lon = cos(y_lat_cntr[i]*pi/180) * km_per_lat;
            distEW = 2*cellsize*km_per_lon;
            distDiag = sqrt(distNS^2+distEW^2);

            for j = 2:size(matrix,2)-1
                # Gradients, in matrix units per km
                if matrix[i,j]<minmatval
                    slope[i,j] = 0;
                else
                    if matrix[i+1,j]<minmatval || matrix[i-1,j]<minmatval
                        NS = 0;
                    else
                        NS = abs(matrix[i+1,j]-matrix[i-1,j])/distNS;
                    end
                    if matrix[i,j+1]<minmatval || matrix[i,j-1]<minmatval
                        EW = 0;
                    else
                        EW = abs(matrix[i,j+1]-matrix[i,j-1])/distEW
                    end
                    if matrix[i+1,j-1]<minmatval || matrix[i-1,j+1]<minmatval
                        NESW = 0;
                    else
                        NESW = abs(matrix[i+1,j-1]-matrix[i-1,j+1])/distDiag
                    end
                    if matrix[i+1,j+1]<minmatval || matrix[i-1,j-1]<minmatval
                        NWSE = 0;
                    else
                        NWSE = abs(matrix[i+1,j+1]-matrix[i-1,j-1])/distDiag
                    end

                    # Record the steepest slope
                    slope[i,j] = round(UInt16, min(max(max(NS,EW), max(NESW,NWSE)),0xffff));
                end
            end

            # Fill in edges too
            distEW = cellsize*km_per_lon;
            distDiag = sqrt((distNS/2)^2+distEW^2);

            # Left edge
            if matrix[i+1,1]<minmatval || matrix[i-1,1]<minmatval
                NS = 0;
            else
                NS = abs(matrix[i+1,1]-matrix[i-1,1])/distNS
            end
            if matrix[i,2]<minmatval || matrix[i,1]<minmatval
                EW = 0;
            else
                EW = abs(matrix[i,2]-matrix[i,1])/distEW
            end
            if matrix[i+1,1]<minmatval || matrix[i-1,2]<minmatval
                NESW = 0;
            else
                NESW = abs(matrix[i+1,1]-matrix[i-1,2])/distDiag
            end
            if matrix[i+1,2]<minmatval || matrix[i-1,1]<minmatval
                NWSE = 0;
            else
                NWSE = abs(matrix[i+1,2]-matrix[i-1,1])/distDiag
            end
            slope[i,1] = round(UInt16, min(max(max(NS,EW), max(NESW,NWSE)),0xffff));

            # Right edge
            if matrix[i+1,end]<minmatval || matrix[i-1,end]<minmatval
                NS = 0;
            else
                NS = abs(matrix[i+1,end]-matrix[i-1,end])/distNS
            end
            if matrix[i,end]<minmatval || matrix[i,end-1]<minmatval
                EW = 0;
            else
                EW = abs(matrix[i,end]-matrix[i,end-1])/distEW
            end
            if matrix[i+1,end-1]<minmatval || matrix[i-1,end]<minmatval
                NEWS = 0;
            else
                NESW = abs(matrix[i+1,end-1]-matrix[i-1,end])/distDiag
            end
            if matrix[i+1,end]<minmatval || matrix[i-1,end-1]<minmatval
                NWSE = 0;
            else
                NWSE = abs(matrix[i+1,end]-matrix[i-1,end-1])/distDiag
            end
            slope[i,end] = round(UInt16, min(max(max(NS,EW), max(NESW,NWSE)),0xffff));
        end

        # Fill in the top and bottom row
        distNS = cellsize*km_per_lat;

        # Top row
        km_per_lon = cos(y_lat_cntr[1]*pi/180) * km_per_lat;
        distEW = 2*cellsize*km_per_lon;
        distDiag = sqrt(distNS^2+(distEW/2)^2)
        for j = 2:size(matrix,2)-1
            # Gradients, in meters per km
            if matrix[2,j]<minmatval || matrix[1,j]<minmatval
                NS = 0;
            else
                NS = abs(matrix[2,j]-matrix[1,j])/distNS
            end
            if matrix[1,j+1]<minmatval || matrix[1,j-1]<minmatval
                EW = 0;
            else
                EW = abs(matrix[1,j+1]-matrix[1,j-1])/distEW
            end
            if matrix[2,j-1]<minmatval || matrix[1,j]<minmatval
                NESW = 0;
            else
                NESW = abs(matrix[2,j-1]-matrix[1,j])/distDiag
            end
            if matrix[2,j+1]<minmatval || matrix[1,j]<minmatval
                NWSE = 0;
            else
                NWSE = abs(matrix[2,j+1]-matrix[1,j])/distDiag
            end
            slope[1,j] = round(UInt16, min(max(max(NS,EW), max(NESW,NWSE)),0xffff));
        end
        slope[1,1] = 0;
        slope[1,end] = 0;

        # Bottom row
        km_per_lon = cos(y_lat_cntr[end]*pi/180) * km_per_lat;
        distEW = 2*cellsize*km_per_lon;
        distDiag = sqrt(distNS^2+(distEW/2)^2)
        for j = 2:size(matrix,2)-1
            # Gradients, in meters per Km
            if matrix[end-1,j]<minmatval || matrix[end,j]<minmatval
                NS = 0;
            else
                NS = abs(matrix[end-1,j]-matrix[end,j])/distNS
            end
            if matrix[end,j+1]<minmatval || matrix[end,j-1]<minmatval
                EW = 0;
            else
                EW = abs(matrix[end,j+1]-matrix[end,j-1])/distEW
            end
            if matrix[end-1,j-1]<minmatval || matrix[end,j]<minmatval
                NESW = 0;
            else
                NESW = abs(matrix[end-1,j-1]-matrix[end,j])/distDiag
            end
            if matrix[end-1,j+1]<minmatval || matrix[end,j]<minmatval
                NWSE = 0;
            else
                NWSE = abs(matrix[end-1,j+1]-matrix[end,j])/distDiag
            end
            slope[end,j] = round(UInt16, min(max(max(NS,EW), max(NESW,NWSE)),0xffff));
        end
        slope[end,1] = 0;
        slope[end,end] = 0;

        return slope
    end
    export quickMaxSlopeEarth

    function quickAveSlopeEarth(matrix, x_lon_cntr, y_lat_cntr, cellsize; minmatval=-500, maxmatval=9000)
        # Returns slope in units/kilometer given a latitude-longitude grid of z-values

        # Allocate intermediate and output arrays
        distance = Array{Float64}(8);
        local_slopes = Array{Float64}(8);
        slope = Array{UInt16}(size(matrix))

        # Index offsets to cycle through:
        #         [N,NE,E,SE,S,SW,W,NW]
        ioffset = [-1,-1,0,1,1,1,0,-1];
        joffset = [0,1,1,1,0,-1,-1,-1];
        #
        # i.e. Layout:
        # 8 1 2
        # 7 x 3
        # 6 5 4

        # Average size of a degree on Earth
        km_per_lat = 111.1;

        # Distance between grid cell centers
        # N, S
        distance[[1,5]] = cellsize*km_per_lat;

        # Fill in the center first
        for i=2:size(matrix,1)-1
            # Distance between grid cell centers
            km_per_lon = cos(y_lat_cntr[i]*pi/180) * km_per_lat;
            distance[[3,7]] = cellsize*km_per_lon; #E, W
            distance[[2,4,6,8]] = sqrt(distance[1]^2+distance[3]^2);  # Diagonals

            # Center
            for j = 2:size(matrix,2)-1
                # Gradients, in matrix z-units per km
                here = matrix[i,j];
                if here<minmatval || here>maxmatval
                    slope[i,j] = 0;
                else
                    for k=1:8
                        there = matrix[i+ioffset[k],j+joffset[k]];
                        if there<minmatval || there>maxmatval
                            local_slopes[k] = 0;
                        else
                            local_slopes[k] = abs(there-here)/distance[k];
                        end
                    end
                    # Record the average slope
                    slope[i,j] = round(UInt16, min(mean(local_slopes),0xffff));
                end
            end

            # Left edge
            here = matrix[i,1];
            if here<minmatval || here>maxmatval
                slope[i,1] = 0;
            else
                for k=1:5
                    there = matrix[i+ioffset[k],1+joffset[k]];
                    if there<minmatval || there>maxmatval
                        local_slopes[k] = 0;
                    else
                        local_slopes[k] = abs(there-here)/distance[k];
                    end
                end
                slope[i,1] = round(UInt16, min(mean(local_slopes[1:5]),0xffff));
            end

            # Right edge
            here = matrix[i,end];
            if here<minmatval || here>maxmatval
                slope[i,end] = 0;
            else
                for k=[5,6,7,8,1]
                    there = matrix[i+ioffset[k],end+joffset[k]];
                    if there<minmatval || there>maxmatval
                        local_slopes[k] = 0;
                    else
                        local_slopes[k] = abs(there-here)/distance[k];
                    end
                end
                slope[i,end] = round(UInt16, min(mean(local_slopes[[5,6,7,8,1]]),0xffff));
            end
        end

        # Top row
        km_per_lon = cos(y_lat_cntr[1]*pi/180) * km_per_lat;
        distance[[3,7]] = cellsize*km_per_lon; #E, W
        distance[[2,4,6,8]] = sqrt(distance[1]^2+distance[3]^2);  # Diagonals
        for j = 2:size(matrix,2)-1
            # Gradients, in matrix units per km
            here = matrix[1,j];
            if here<minmatval || here>maxmatval
                slope[1,j] = 0;
            else
                for k=3:7
                    there = matrix[1+ioffset[k],j+joffset[k]];
                    if there<minmatval || there>maxmatval
                        local_slopes[k] = 0;
                    else
                        local_slopes[k] = abs(there-here)/distance[k];
                    end
                end
                slope[1,j] = round(UInt16, min(mean(local_slopes[3:7]),0xffff));
            end
        end
        slope[1,1] = 0;
        slope[1,end] = 0;

        # Bottom row
        km_per_lon = cos(y_lat_cntr[end]*pi/180) * km_per_lat;
        distance[[3,7]] = cellsize*km_per_lon; #E, W
        distance[[2,4,6,8]] = sqrt(distance[1]^2+distance[3]^2);  # Diagonals
        for j = 2:size(matrix,2)-1
            # Gradients, in matrix units per km
            here = matrix[end,j];
            if here<minmatval || here>maxmatval
                slope[end,j] = 0;
            else
                for k=[7,8,1,2,3]
                    there = matrix[end+ioffset[k],j+joffset[k]];
                    if there<minmatval || there>maxmatval
                        local_slopes[k] = 0;
                    else
                        local_slopes[k] = abs(there-here)/distance[k];
                    end
                end
                slope[end,j] = round(UInt16, min(mean(local_slopes[[7,8,1,2,3]]),0xffff));
            end
        end
        slope[end,1] = 0;
        slope[end,end] = 0;

        return slope
    end
    export quickAveSlopeEarth

## --- End of File

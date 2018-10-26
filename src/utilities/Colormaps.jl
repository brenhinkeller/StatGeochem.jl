## --- Matplotlib colormaps

    # Viridis
    viridis = RGB{N0f8}.(
        [0.267004, 0.26851, 0.269944, 0.271305, 0.272594, 0.273809, 0.274952, 0.276022, 0.277018, 0.277941, 0.278791, 0.279566, 0.280267, 0.280894, 0.281446, 0.281924, 0.282327, 0.282656, 0.28291, 0.283091, 0.283197, 0.283229, 0.283187, 0.283072, 0.282884, 0.282623, 0.28229, 0.281887, 0.281412, 0.280868, 0.280255, 0.279574, 0.278826, 0.278012, 0.277134, 0.276194, 0.275191, 0.274128, 0.273006, 0.271828, 0.270595, 0.269308, 0.267968, 0.26658, 0.265145, 0.263663, 0.262138, 0.260571, 0.258965, 0.257322, 0.255645, 0.253935, 0.252194, 0.250425, 0.248629, 0.246811, 0.244972, 0.243113, 0.241237, 0.239346, 0.237441, 0.235526, 0.233603, 0.231674, 0.229739, 0.227802, 0.225863, 0.223925, 0.221989, 0.220057, 0.21813, 0.21621, 0.214298, 0.212395, 0.210503, 0.208623, 0.206756, 0.204903, 0.203063, 0.201239, 0.19943, 0.197636, 0.19586, 0.1941, 0.192357, 0.190631, 0.188923, 0.187231, 0.185556, 0.183898, 0.182256, 0.180629, 0.179019, 0.177423, 0.175841, 0.174274, 0.172719, 0.171176, 0.169646,0.168126, 0.166617, 0.165117, 0.163625, 0.162142, 0.160665, 0.159194, 0.157729, 0.15627, 0.154815, 0.153364,0.151918, 0.150476, 0.149039, 0.147607, 0.14618, 0.144759, 0.143343, 0.141935, 0.140536, 0.139147, 0.13777, 0.136408, 0.135066, 0.133743, 0.132444, 0.131172, 0.129933, 0.128729, 0.127568, 0.126453, 0.125394, 0.124395,0.123463, 0.122606, 0.121831, 0.121148, 0.120565, 0.120092, 0.119738, 0.119512, 0.119423, 0.119483, 0.119699, 0.120081, 0.120638, 0.12138, 0.122312, 0.123444, 0.12478, 0.126326, 0.128087, 0.130067, 0.132268, 0.134692,0.137339, 0.14021, 0.143303, 0.146616, 0.150148, 0.153894, 0.157851, 0.162016, 0.166383, 0.170948, 0.175707,0.180653, 0.185783, 0.19109, 0.196571, 0.202219, 0.20803, 0.214, 0.220124, 0.226397, 0.232815, 0.239374, 0.24607, 0.252899, 0.259857, 0.266941, 0.274149, 0.281477, 0.288921, 0.296479, 0.304148, 0.311925, 0.319809, 0.327796, 0.335885, 0.344074, 0.35236, 0.360741, 0.369214, 0.377779, 0.386433, 0.395174, 0.404001, 0.412913, 0.421908, 0.430983, 0.440137, 0.449368, 0.458674, 0.468053, 0.477504, 0.487026, 0.496615, 0.506271, 0.515992, 0.525776, 0.535621, 0.545524, 0.555484, 0.565498, 0.575563, 0.585678, 0.595839, 0.606045, 0.616293, 0.626579, 0.636902, 0.647257, 0.657642, 0.668054, 0.678489, 0.688944, 0.699415, 0.709898, 0.720391, 0.730889, 0.741388, 0.751884, 0.762373, 0.772852, 0.783315, 0.79376, 0.804182, 0.814576, 0.82494, 0.83527, 0.845561, 0.85581, 0.866013, 0.876168, 0.886271, 0.89632, 0.906311, 0.916242, 0.926106, 0.935904, 0.945636, 0.9553, 0.964894, 0.974417, 0.983868, 0.993248],
        [0.004874, 0.009605, 0.014625, 0.019942, 0.025563, 0.031497, 0.037752, 0.044167, 0.050344, 0.056324, 0.062145, 0.067836, 0.073417, 0.078907, 0.08432, 0.089666, 0.094955, 0.100196, 0.105393, 0.110553, 0.11568, 0.120777, 0.125848, 0.130895, 0.13592, 0.140926, 0.145912, 0.150881, 0.155834, 0.160771, 0.165693, 0.170599, 0.17549,0.180367, 0.185228, 0.190074, 0.194905, 0.199721, 0.20452, 0.209303, 0.214069, 0.218818, 0.223549, 0.228262,0.232956, 0.237631, 0.242286, 0.246922, 0.251537, 0.25613, 0.260703, 0.265254, 0.269783, 0.27429, 0.278775, 0.283237, 0.287675, 0.292092, 0.296485, 0.300855, 0.305202, 0.309527, 0.313828, 0.318106, 0.322361, 0.326594,0.330805, 0.334994, 0.339161, 0.343307, 0.347432, 0.351535, 0.355619, 0.359683, 0.363727, 0.367752, 0.371758, 0.375746, 0.379716, 0.38367, 0.387607, 0.391528, 0.395433, 0.399323, 0.403199, 0.407061, 0.41091, 0.414746,0.41857, 0.422383, 0.426184, 0.429975, 0.433756, 0.437527, 0.44129, 0.445044, 0.448791, 0.45253, 0.456262, 0.459988, 0.463708, 0.467423, 0.471133, 0.474838, 0.47854, 0.482237, 0.485932, 0.489624, 0.493313, 0.497, 0.500685, 0.504369, 0.508051, 0.511733, 0.515413, 0.519093, 0.522773, 0.526453, 0.530132, 0.533812, 0.537492, 0.541173, 0.544853, 0.548535, 0.552216, 0.555899, 0.559582, 0.563265, 0.566949, 0.570633, 0.574318, 0.578002, 0.581687, 0.585371, 0.589055, 0.592739, 0.596422, 0.600104, 0.603785, 0.607464, 0.611141, 0.614817, 0.61849, 0.622161, 0.625828, 0.629492, 0.633153, 0.636809, 0.640461, 0.644107, 0.647749, 0.651384, 0.655014, 0.658636, 0.662252, 0.665859, 0.669459, 0.67305, 0.676631, 0.680203, 0.683765, 0.687316, 0.690856, 0.694384, 0.6979, 0.701402, 0.704891, 0.708366, 0.711827, 0.715272, 0.718701, 0.722114, 0.725509, 0.728888, 0.732247, 0.735588, 0.73891, 0.742211, 0.745492, 0.748751, 0.751988, 0.755203, 0.758394, 0.761561, 0.764704, 0.767822, 0.770914, 0.77398, 0.777018, 0.780029, 0.783011, 0.785964, 0.788888, 0.791781, 0.794644, 0.797475, 0.800275, 0.803041, 0.805774, 0.808473, 0.811138, 0.813768, 0.816363, 0.818921, 0.821444, 0.823929, 0.826376, 0.828786, 0.831158, 0.833491, 0.835785, 0.838039, 0.840254, 0.84243, 0.844566, 0.846661, 0.848717, 0.850733, 0.852709, 0.854645, 0.856542, 0.8584, 0.860219, 0.861999, 0.863742, 0.865448, 0.867117, 0.868751, 0.87035, 0.871916, 0.873449, 0.874951, 0.876424, 0.877868, 0.879285, 0.880678, 0.882046, 0.883393, 0.88472, 0.886029, 0.887322, 0.888601, 0.889868, 0.891125, 0.892374, 0.893616, 0.894855, 0.896091, 0.89733, 0.89857, 0.899815, 0.901065, 0.902323, 0.90359, 0.904867, 0.906157],
        [0.329415, 0.335427, 0.341379, 0.347269, 0.353093, 0.358853, 0.364543, 0.370164, 0.375715, 0.381191, 0.386592, 0.391917, 0.397163, 0.402329, 0.407414, 0.412415, 0.417331, 0.42216, 0.426902, 0.431554, 0.436115, 0.440584, 0.44496, 0.449241, 0.453427, 0.457517, 0.46151, 0.465405, 0.469201, 0.472899, 0.476498, 0.479997, 0.483397, 0.486697, 0.489898, 0.493001, 0.496005, 0.498911, 0.501721, 0.504434, 0.507052, 0.509577, 0.512008, 0.514349, 0.516599, 0.518762, 0.520837, 0.522828, 0.524736, 0.526563, 0.528312, 0.529983, 0.531579, 0.533103, 0.534556, 0.535941, 0.53726, 0.538516, 0.539709, 0.540844, 0.541921, 0.542944, 0.543914, 0.544834, 0.545706, 0.546532, 0.547314, 0.548053, 0.548752, 0.549413, 0.550038, 0.550627, 0.551184, 0.55171, 0.552206, 0.552675, 0.553117, 0.553533, 0.553925, 0.554294, 0.554642, 0.554969, 0.555276, 0.555565, 0.555836, 0.556089, 0.556326, 0.556547, 0.556753, 0.556944, 0.55712, 0.557282, 0.55743, 0.557565, 0.557685, 0.557792, 0.557885, 0.557965, 0.55803, 0.558082, 0.558119, 0.558141, 0.558148, 0.55814, 0.558115, 0.558073, 0.558013, 0.557936, 0.55784, 0.557724, 0.557587, 0.55743, 0.55725, 0.557049, 0.556823, 0.556572, 0.556295, 0.555991, 0.555659, 0.555298, 0.554906,0.554483, 0.554029, 0.553541, 0.553018, 0.552459, 0.551864, 0.551229, 0.550556, 0.549841, 0.549086, 0.548287, 0.547445, 0.546557, 0.545623, 0.544641, 0.543611, 0.54253, 0.5414, 0.540218, 0.538982, 0.537692, 0.536347, 0.534946, 0.533488, 0.531973, 0.530398, 0.528763, 0.527068, 0.525311, 0.523491, 0.521608, 0.519661, 0.517649,0.515571, 0.513427, 0.511215, 0.508936, 0.506589, 0.504172, 0.501686, 0.499129, 0.496502, 0.493803, 0.491033, 0.488189, 0.485273, 0.482284, 0.479221, 0.476084, 0.472873, 0.469588, 0.466226, 0.462789, 0.459277, 0.455688, 0.452024, 0.448284, 0.444467, 0.440573, 0.436601, 0.432552, 0.428426, 0.424223, 0.419943, 0.415586, 0.411152, 0.40664, 0.402049, 0.397381, 0.392636, 0.387814, 0.382914, 0.377939, 0.372886, 0.367757, 0.362552, 0.357269, 0.35191, 0.346476, 0.340967, 0.335384, 0.329727, 0.323998, 0.318195, 0.312321, 0.306377, 0.300362, 0.294279, 0.288127, 0.281908, 0.275626, 0.269281, 0.262877, 0.256415, 0.249897, 0.243329, 0.236712, 0.230052, 0.223353, 0.21662, 0.209861, 0.203082, 0.196293, 0.189503, 0.182725, 0.175971, 0.169257, 0.162603, 0.156029, 0.149561, 0.143228, 0.137064, 0.131109, 0.125405, 0.120005, 0.114965, 0.110347, 0.106217, 0.102646, 0.099702, 0.097452, 0.095953, 0.09525, 0.095374, 0.096335, 0.098125, 0.100717, 0.104071, 0.108131, 0.112838, 0.118128, 0.123941, 0.130215, 0.136897, 0.143936]
        )
    export viridis

    # Plasma
    plasma = RGB{N0f8}.(
        [0.050383, 0.063536, 0.075353, 0.086222, 0.096379, 0.10598, 0.115124, 0.123903, 0.132381, 0.140603, 0.148607, 0.156421, 0.16407, 0.171574, 0.17895, 0.186213, 0.193374, 0.200445, 0.207435, 0.21435, 0.221197, 0.227983, 0.234715, 0.241396, 0.248032, 0.254627, 0.261183, 0.267703, 0.274191, 0.280648, 0.287076, 0.293478, 0.299855,0.30621, 0.312543, 0.318856, 0.32515, 0.331426, 0.337683, 0.343925, 0.35015, 0.356359, 0.362553, 0.368733, 0.374897, 0.381047, 0.387183, 0.393304, 0.399411, 0.405503, 0.41158, 0.417642, 0.423689, 0.429719, 0.435734, 0.441732, 0.447714, 0.453677, 0.459623, 0.46555, 0.471457, 0.477344, 0.48321, 0.489055, 0.494877, 0.500678, 0.506454, 0.512206, 0.517933, 0.523633, 0.529306, 0.534952, 0.54057, 0.546157, 0.551715, 0.557243, 0.562738, 0.568201, 0.573632, 0.579029, 0.584391, 0.589719, 0.595011, 0.600266, 0.605485, 0.610667, 0.615812, 0.620919, 0.625987, 0.631017, 0.636008, 0.640959, 0.645872, 0.650746, 0.65558, 0.660374, 0.665129, 0.669845, 0.674522, 0.67916, 0.683758, 0.688318, 0.69284, 0.697324, 0.701769, 0.706178, 0.710549, 0.714883, 0.719181, 0.723444, 0.72767, 0.731862, 0.736019, 0.740143, 0.744232, 0.748289, 0.752312, 0.756304, 0.760264, 0.764193, 0.76809, 0.771958, 0.775796, 0.779604, 0.783383, 0.787133, 0.790855, 0.794549, 0.798216, 0.801855, 0.805467, 0.809052, 0.812612, 0.816144, 0.819651, 0.823132, 0.826588, 0.830018, 0.833422, 0.836801, 0.840155, 0.843484, 0.846788, 0.850066, 0.853319, 0.856547, 0.85975, 0.862927, 0.866078, 0.869203, 0.872303, 0.875376, 0.878423, 0.881443, 0.884436, 0.887402, 0.89034, 0.89325, 0.896131, 0.898984, 0.901807, 0.904601, 0.907365, 0.910098, 0.9128, 0.915471, 0.918109, 0.920714, 0.923287, 0.925825, 0.928329, 0.930798, 0.933232, 0.93563, 0.93799, 0.940313, 0.942598, 0.944844, 0.947051, 0.949217, 0.951344, 0.953428, 0.95547, 0.957469, 0.959424, 0.961336, 0.963203, 0.965024, 0.966798, 0.968526, 0.970205, 0.971835, 0.973416, 0.974947, 0.976428, 0.977856, 0.979233, 0.980556, 0.981826, 0.983041, 0.984199, 0.985301, 0.986345, 0.987332, 0.98826, 0.989128, 0.989935, 0.990681, 0.991365, 0.991985, 0.992541, 0.993032, 0.993456, 0.993814, 0.994103, 0.994324, 0.994474, 0.994553, 0.994561, 0.994495, 0.994355, 0.994141, 0.993851, 0.993482, 0.993033, 0.992505, 0.991897, 0.991209, 0.990439, 0.989587, 0.988648, 0.987621, 0.986509, 0.985314, 0.984031, 0.982653, 0.98119, 0.979644, 0.977995, 0.976265, 0.974443, 0.97253, 0.970533, 0.968443, 0.966271, 0.964021, 0.961681, 0.959276, 0.956808, 0.954287, 0.951726, 0.949151, 0.946602, 0.944152, 0.941896, 0.940015],
        [0.029803, 0.028426, 0.027206, 0.026125, 0.025165, 0.024309, 0.023556, 0.022878, 0.022258, 0.021687, 0.021154, 0.020651, 0.020171, 0.019706, 0.019252, 0.018803, 0.018354, 0.017902, 0.017442, 0.016973, 0.016497, 0.016007, 0.015502, 0.014979, 0.014439, 0.013882, 0.013308, 0.012716, 0.012109, 0.011488, 0.010855, 0.010213, 0.009561, 0.008902, 0.008239, 0.007576, 0.006915, 0.006261, 0.005618, 0.004991, 0.004382, 0.003798, 0.003243, 0.002724, 0.002245, 0.001814, 0.001434, 0.001114, 0.000859, 0.000678, 0.000577, 0.000564, 0.000646, 0.000831, 0.001127, 0.00154, 0.00208, 0.002755, 0.003574, 0.004545, 0.005678, 0.00698, 0.00846, 0.010127, 0.01199, 0.014055, 0.016333, 0.018833, 0.021563, 0.024532, 0.027747, 0.031217, 0.03495, 0.038954, 0.043136, 0.047331, 0.051545, 0.055778, 0.060028, 0.064296, 0.068579, 0.072878, 0.07719, 0.081516, 0.085854, 0.090204, 0.094564, 0.098934, 0.103312, 0.107699, 0.112092, 0.116492, 0.120898, 0.125309, 0.129725, 0.134144, 0.138566, 0.142992, 0.147419, 0.151848, 0.156278, 0.160709, 0.165141, 0.169573, 0.174005, 0.178437, 0.182868, 0.187299, 0.191729, 0.196158, 0.200586, 0.205013, 0.209439, 0.213864, 0.218288, 0.222711, 0.227133, 0.231555, 0.235976, 0.240396, 0.244817, 0.249237, 0.253658, 0.258078, 0.2625, 0.266922, 0.271345, 0.27577, 0.280197, 0.284626, 0.289057, 0.293491, 0.297928, 0.302368, 0.306812, 0.311261, 0.315714, 0.320172, 0.324635, 0.329105, 0.33358, 0.338062, 0.342551, 0.347048, 0.351553, 0.356066, 0.360588, 0.365119, 0.36966, 0.374212, 0.378774, 0.383347, 0.387932, 0.392529, 0.397139, 0.401762, 0.406398, 0.411048, 0.415712, 0.420392, 0.425087, 0.429797, 0.434524, 0.439268, 0.444029, 0.448807, 0.453603, 0.458417, 0.463251, 0.468103, 0.472975, 0.477867, 0.48278, 0.487712, 0.492667, 0.497642, 0.502639, 0.507658, 0.512699, 0.517763, 0.52285, 0.52796, 0.533093, 0.53825, 0.543431, 0.548636, 0.553865, 0.559118, 0.564396, 0.5697, 0.575028, 0.580382, 0.585761, 0.591165, 0.596595, 0.602051, 0.607532, 0.613039, 0.618572, 0.624131, 0.629718, 0.63533, 0.640969, 0.646633, 0.652325, 0.658043, 0.663787, 0.669558, 0.675355, 0.681179, 0.68703, 0.692907, 0.69881, 0.704741, 0.710698, 0.716681, 0.722691, 0.728728, 0.734791, 0.74088, 0.746995, 0.753137, 0.759304, 0.765499, 0.77172, 0.777967, 0.784239, 0.790537, 0.796859, 0.803205, 0.809579, 0.815978, 0.822401, 0.828846, 0.835315, 0.841812, 0.848329, 0.854866, 0.861432, 0.868016, 0.874622, 0.88125, 0.887896, 0.894564, 0.901249, 0.90795, 0.914672, 0.921407, 0.928152, 0.934908, 0.941671, 0.948435, 0.95519, 0.961916, 0.96859, 0.975158],
        [0.527975, 0.533124, 0.538007, 0.542658, 0.547103, 0.551368, 0.555468, 0.559423, 0.56325, 0.566959, 0.570562, 0.574065, 0.577478, 0.580806, 0.584054, 0.587228, 0.59033, 0.593364, 0.596333, 0.599239, 0.602083, 0.604867, 0.607592, 0.610259, 0.612868, 0.615419, 0.617911, 0.620346, 0.622722, 0.625038, 0.627295, 0.62949, 0.631624, 0.633694, 0.6357, 0.63764, 0.639512, 0.641316, 0.643049, 0.64471, 0.646298, 0.64781, 0.649245, 0.650601, 0.651876, 0.653068, 0.654177, 0.655199, 0.656133, 0.656977, 0.65773, 0.65839, 0.658956, 0.659425, 0.659797, 0.660069, 0.66024, 0.66031, 0.660277, 0.660139, 0.659897, 0.659549, 0.659095, 0.658534, 0.657865, 0.657088, 0.656202, 0.655209, 0.654109, 0.652901, 0.651586, 0.650165, 0.64864, 0.64701, 0.645277, 0.643443, 0.641509, 0.639477, 0.637349, 0.635126, 0.632812, 0.630408, 0.627917, 0.625342, 0.622686, 0.619951, 0.61714, 0.614257, 0.611305, 0.608287, 0.605205, 0.602065, 0.598867, 0.595617, 0.592317, 0.588971, 0.585582, 0.582154, 0.578688, 0.575189, 0.57166, 0.568103, 0.564522, 0.560919, 0.557296, 0.553657, 0.550004, 0.546338, 0.542663, 0.538981, 0.535293, 0.531601, 0.527908, 0.524216, 0.520524, 0.516834, 0.513149, 0.509468, 0.505794, 0.502126, 0.498465, 0.494813, 0.491171, 0.487539, 0.483918, 0.480307, 0.476706, 0.473117, 0.469538, 0.465971, 0.462415, 0.45887, 0.455338, 0.451816, 0.448306, 0.444806, 0.441316, 0.437836, 0.434366, 0.430905, 0.427455, 0.424013, 0.420579, 0.417153, 0.413734, 0.410322, 0.406917, 0.403519, 0.400126, 0.396738, 0.393355, 0.389976, 0.3866, 0.383229, 0.37986, 0.376494, 0.37313, 0.369768, 0.366407, 0.363047, 0.359688, 0.356329, 0.35297, 0.34961, 0.346251, 0.34289, 0.339529, 0.336166, 0.332801, 0.329435, 0.326067, 0.322697, 0.319325, 0.315952, 0.312575, 0.309197, 0.305816, 0.302433, 0.299049, 0.295662, 0.292275, 0.288883, 0.28549, 0.282096, 0.278701, 0.275305, 0.271909, 0.268513, 0.265118, 0.261721, 0.258325, 0.254931, 0.25154, 0.248151, 0.244767, 0.241387, 0.238013, 0.234646, 0.231287, 0.227937, 0.224595, 0.221265, 0.217948, 0.214648, 0.211364, 0.2081, 0.204859, 0.201642, 0.198453, 0.195295, 0.19217, 0.189084, 0.186041, 0.183043, 0.180097, 0.177208, 0.174381, 0.171622, 0.168938, 0.166335, 0.163821, 0.161404, 0.159092, 0.156891, 0.154808, 0.152855, 0.151042, 0.149377, 0.14787, 0.146529, 0.145357, 0.144363, 0.143557, 0.142945, 0.142528, 0.142303, 0.142279, 0.142453, 0.142808, 0.143351, 0.144061, 0.144923, 0.145919, 0.147014, 0.14818, 0.14937, 0.15052, 0.151566, 0.152409, 0.152921, 0.152925, 0.152178, 0.150328, 0.146861, 0.140956, 0.131326]
        )
    export plasma

    # Magma
    magma = RGB{N0f8}.(
        [0.001462, 0.002258, 0.003279, 0.004512, 0.00595, 0.007588, 0.009426, 0.011465, 0.013708, 0.016156, 0.018815, 0.021692, 0.024792, 0.028123, 0.031696, 0.03552, 0.039608, 0.04383, 0.048062, 0.05232, 0.056615, 0.060949, 0.06533, 0.069764, 0.074257, 0.078815, 0.083446, 0.088155, 0.092949, 0.097833, 0.102815, 0.107899, 0.113094, 0.118405, 0.123833, 0.12938, 0.135053, 0.140858, 0.146785, 0.152839, 0.159018, 0.165308, 0.171713, 0.178212, 0.184801, 0.19146, 0.198177, 0.204935, 0.211718, 0.218512, 0.225302, 0.232077, 0.238826, 0.245543, 0.25222, 0.258857, 0.265447, 0.271994, 0.278493, 0.284951, 0.291366, 0.29774, 0.304081, 0.310382, 0.316654, 0.322899, 0.329114, 0.335308, 0.341482, 0.347636, 0.353773, 0.359898, 0.366012, 0.372116, 0.378211, 0.384299, 0.390384, 0.396467, 0.402548, 0.408629, 0.414709, 0.420791, 0.426877, 0.432967, 0.439062, 0.445163, 0.451271, 0.457386,0.463508, 0.46964, 0.47578, 0.481929, 0.488088, 0.494258, 0.500438, 0.506629, 0.512831, 0.519045, 0.52527, 0.531507, 0.537755, 0.544015, 0.550287, 0.556571, 0.562866, 0.569172, 0.57549, 0.581819, 0.588158, 0.594508, 0.600868, 0.607238, 0.613617, 0.620005, 0.626401, 0.632805, 0.639216, 0.645633, 0.652056, 0.658483, 0.664915, 0.671349, 0.677786, 0.684224, 0.690661, 0.697098, 0.703532, 0.709962, 0.716387, 0.722805, 0.729216, 0.735616,0.742004, 0.748378, 0.754737, 0.761077, 0.767398, 0.773695, 0.779968, 0.786212, 0.792427, 0.798608, 0.804752, 0.810855, 0.816914, 0.822926, 0.828886, 0.834791, 0.840636, 0.846416, 0.852126, 0.857763, 0.86332, 0.868793, 0.874176, 0.879464, 0.884651, 0.889731, 0.8947, 0.899552, 0.904281, 0.908884, 0.913354, 0.917689, 0.921884,0.925937, 0.929845, 0.933606, 0.937221, 0.940687, 0.944006, 0.94718, 0.95021, 0.953099, 0.955849, 0.958464, 0.960949, 0.96331, 0.965549, 0.967671, 0.96968, 0.971582, 0.973381, 0.975082, 0.97669, 0.97821, 0.979645, 0.981, 0.982279, 0.983485, 0.984622, 0.985693, 0.9867, 0.987646, 0.988533, 0.989363, 0.990138, 0.990871, 0.991558, 0.992196, 0.992785, 0.993326, 0.993834, 0.994309, 0.994738, 0.995122, 0.99548, 0.99581, 0.996096, 0.996341, 0.99658, 0.996775, 0.996925, 0.997077, 0.997186, 0.997254, 0.997325, 0.997351, 0.997351, 0.997341, 0.997285, 0.997228, 0.997138, 0.997019, 0.996898, 0.996727, 0.996571, 0.996369, 0.996162, 0.995932, 0.99568, 0.995424, 0.995131, 0.994851, 0.994524, 0.994222, 0.993866, 0.993545, 0.99317, 0.992831, 0.99244, 0.992089, 0.991688,0.991332, 0.99093, 0.99057, 0.990175, 0.989815, 0.989434, 0.989077, 0.988717, 0.988367, 0.988033, 0.987691, 0.987387, 0.987053],
        [0.000466, 0.001295, 0.002305, 0.00349, 0.004843, 0.006356, 0.008022, 0.009828, 0.011771, 0.01384, 0.016026,0.01832, 0.020715, 0.023201, 0.025765, 0.028397, 0.03109, 0.03383, 0.036607, 0.039407, 0.04216, 0.044794, 0.047318, 0.049726, 0.052017, 0.054184, 0.056225, 0.058133, 0.059904, 0.061531, 0.06301, 0.064335, 0.065492, 0.066479, 0.067295, 0.067935, 0.068391, 0.068654, 0.068738, 0.068637, 0.068354, 0.067911, 0.067305, 0.066576, 0.065732, 0.064818, 0.063862, 0.062907, 0.061992, 0.061158, 0.060445, 0.059889, 0.059517, 0.059352, 0.059415, 0.059706, 0.060237, 0.060994, 0.061978, 0.063168, 0.064553, 0.066117, 0.067835, 0.069702, 0.07169, 0.073782, 0.075972, 0.078236, 0.080564, 0.082946, 0.085373, 0.087831, 0.090314, 0.092816, 0.095332, 0.097855, 0.100379,0.102902, 0.10542, 0.10793, 0.110431, 0.11292, 0.115395, 0.117855, 0.120298, 0.122724, 0.125132, 0.127522, 0.129893, 0.132245, 0.134577, 0.136891, 0.139186, 0.141462, 0.143719, 0.145958, 0.148179, 0.150383, 0.152569, 0.154739, 0.156894, 0.159033, 0.161158, 0.163269, 0.165368, 0.167454, 0.16953, 0.171596, 0.173652, 0.175701, 0.177743, 0.179779, 0.181811, 0.18384, 0.185867, 0.187893, 0.189921, 0.191952, 0.193986, 0.196027, 0.198075, 0.200133, 0.202203, 0.204286, 0.206384, 0.208501, 0.210638, 0.212797, 0.214982, 0.217194, 0.219437, 0.221713,0.224025, 0.226377, 0.228772, 0.231214, 0.233705, 0.236249, 0.238851, 0.241514, 0.244242, 0.24704, 0.249911,0.252861, 0.255895, 0.259016, 0.262229, 0.26554, 0.268953, 0.272473, 0.276106, 0.279857, 0.283729, 0.287728,0.291859, 0.296125, 0.30053, 0.305079, 0.309773, 0.314616, 0.31961, 0.324755, 0.330052, 0.3355, 0.341098, 0.346844, 0.352734, 0.358764, 0.364929, 0.371224, 0.377643, 0.384178, 0.39082, 0.397563, 0.4044, 0.411324, 0.418323, 0.42539, 0.432519, 0.439703, 0.446936, 0.45421, 0.46152, 0.468861, 0.476226, 0.483612, 0.491014, 0.498428, 0.505851, 0.51328, 0.520713, 0.528148, 0.535582, 0.543015, 0.550446, 0.557873, 0.565296, 0.572706, 0.580107, 0.587502, 0.594891, 0.602275, 0.609644, 0.616999, 0.62435, 0.631696, 0.639027, 0.646344, 0.653659, 0.660969, 0.668256, 0.675541, 0.682828, 0.690088, 0.697349, 0.704611, 0.711848, 0.719089, 0.726324, 0.733545, 0.740772, 0.747981, 0.75519, 0.762398, 0.769591, 0.776795, 0.783977, 0.791167, 0.798348, 0.805527, 0.812706, 0.819875, 0.827052, 0.834213, 0.841387, 0.84854, 0.855711, 0.862859, 0.870024, 0.877168, 0.88433, 0.89147, 0.898627, 0.905763, 0.912915, 0.920049, 0.927196, 0.934329, 0.94147, 0.948604, 0.955742, 0.962878, 0.970012, 0.977154, 0.984288, 0.991438],
        [0.013866, 0.018331, 0.023708, 0.029965, 0.03713, 0.044973, 0.052844, 0.06075, 0.068667, 0.076603, 0.084584,0.09261, 0.100676, 0.108787, 0.116965, 0.125209, 0.133515, 0.141886, 0.150327, 0.158841, 0.167446, 0.176129,0.184892, 0.193735, 0.20266, 0.211667, 0.220755, 0.229922, 0.239164, 0.248477, 0.257854, 0.267289, 0.276784,0.286321, 0.295879, 0.305443, 0.315, 0.324538, 0.334011, 0.343404, 0.352688, 0.361816, 0.370771, 0.379497, 0.387973, 0.396152, 0.404009, 0.411514, 0.418647, 0.425392, 0.431742, 0.437695, 0.443256, 0.448436, 0.453248, 0.45771, 0.46184, 0.46566, 0.46919, 0.472451, 0.475462, 0.478243, 0.480812, 0.483186, 0.48538, 0.487408, 0.489287, 0.491024, 0.492631, 0.494121, 0.495501, 0.496778, 0.49796, 0.499053, 0.500067, 0.501002, 0.501864, 0.502658, 0.503386, 0.504052, 0.504662, 0.505215, 0.505714, 0.50616, 0.506555, 0.506901, 0.507198, 0.507448, 0.507652, 0.507809, 0.507921, 0.507989, 0.508011, 0.507988, 0.50792, 0.507806, 0.507648, 0.507443, 0.507192, 0.506895, 0.506551, 0.506159, 0.505719, 0.50523, 0.504692, 0.504105, 0.503466, 0.502777, 0.502035, 0.501241, 0.500394, 0.499492, 0.498536, 0.497524, 0.496456, 0.495332, 0.49415, 0.49291, 0.491611, 0.490253, 0.488836, 0.487358, 0.485819, 0.484219, 0.482558, 0.480835, 0.479049, 0.477201, 0.47529, 0.473316, 0.471279, 0.46918, 0.467018, 0.464794, 0.462509, 0.460162, 0.457755, 0.455289, 0.452765, 0.450184, 0.447543, 0.444848, 0.442102, 0.439305, 0.436461, 0.433573, 0.430644, 0.427671, 0.424666, 0.421631, 0.418573, 0.415496, 0.412403, 0.409303, 0.406205, 0.403118, 0.400047, 0.397002, 0.393995, 0.391037, 0.388137, 0.385308, 0.382563, 0.379915, 0.377376, 0.374959, 0.372677, 0.370541, 0.368567, 0.366762, 0.365136, 0.363701, 0.362468, 0.361438, 0.360619, 0.360014, 0.35963, 0.359469, 0.359529, 0.35981, 0.360311, 0.36103, 0.361965, 0.363111, 0.364466, 0.366025, 0.367783, 0.369734, 0.371874, 0.374198, 0.376698, 0.379371, 0.38221, 0.38521, 0.388365, 0.391671, 0.395122, 0.398714, 0.402441, 0.406299, 0.410283, 0.41439, 0.418613, 0.42295, 0.427397, 0.431951, 0.436607, 0.441361, 0.446213, 0.45116, 0.456192, 0.461314, 0.466526, 0.471811, 0.477182, 0.482635, 0.488154, 0.493755, 0.499428, 0.505167, 0.510983, 0.516859, 0.522806, 0.528821, 0.534892, 0.541039, 0.547233, 0.553499, 0.55982, 0.566202, 0.572645, 0.57914, 0.585701, 0.592307, 0.598983, 0.605696, 0.612482, 0.619299, 0.626189, 0.633109, 0.640099, 0.647116, 0.654202, 0.661309, 0.668481, 0.675675, 0.682926, 0.690198, 0.697519, 0.704863, 0.712242, 0.719649, 0.727077, 0.734536, 0.742002, 0.749504]
        )
    export magma

    # Inferno
    inferno = RGB{N0f8}.(
        [0.001462, 0.002267, 0.003299, 0.004547, 0.006006, 0.007676, 0.009561, 0.011663, 0.013995, 0.016561, 0.019373, 0.022447, 0.025793, 0.029432, 0.033385, 0.037668, 0.042253, 0.046915, 0.051644, 0.056449, 0.06134, 0.066331, 0.071429, 0.076637, 0.081962, 0.087411, 0.09299, 0.098702, 0.104551, 0.110536, 0.116656, 0.122908, 0.129285, 0.135778, 0.142378, 0.149073, 0.15585, 0.162689, 0.169575, 0.176493, 0.183429, 0.190367, 0.197297, 0.204209, 0.211095, 0.217949, 0.224763, 0.231538, 0.238273, 0.244967, 0.25162, 0.258234, 0.26481, 0.271347, 0.27785,0.284321, 0.290763, 0.297178, 0.303568, 0.309935, 0.316282, 0.32261, 0.328921, 0.335217, 0.3415, 0.347771, 0.354032, 0.360284, 0.366529, 0.372768, 0.379001, 0.385228, 0.391453, 0.397674, 0.403894, 0.410113, 0.416331, 0.422549, 0.428768, 0.434987, 0.441207, 0.447428, 0.453651, 0.459875, 0.4661, 0.472328, 0.478558, 0.484789, 0.491022, 0.497257, 0.503493, 0.50973, 0.515967, 0.522206, 0.528444, 0.534683, 0.54092, 0.547157, 0.553392, 0.559624, 0.565854, 0.572081, 0.578304, 0.584521, 0.590734, 0.59694, 0.603139, 0.60933, 0.615513, 0.621685, 0.627847, 0.633998, 0.640135, 0.64626, 0.652369, 0.658463, 0.66454, 0.670599, 0.676638, 0.682656, 0.688653, 0.694627, 0.700576, 0.7065, 0.712396, 0.718264, 0.724103, 0.729909, 0.735683, 0.741423, 0.747127, 0.752794, 0.758422, 0.76401, 0.769556, 0.775059, 0.780517, 0.785929, 0.791293, 0.796607, 0.801871, 0.807082, 0.812239, 0.817341, 0.822386, 0.827372, 0.832299, 0.837165, 0.841969, 0.846709, 0.851384, 0.855992, 0.860533, 0.865006, 0.869409, 0.873741, 0.878001, 0.882188, 0.886302, 0.890341, 0.894305, 0.898192, 0.902003, 0.905735, 0.90939, 0.912966, 0.916462, 0.919879, 0.923215, 0.92647, 0.929644, 0.932737, 0.935747, 0.938675, 0.941521, 0.944285, 0.946965, 0.949562, 0.952075, 0.954506, 0.956852, 0.959114, 0.961293, 0.963387, 0.965397, 0.967322, 0.969163, 0.970919, 0.97259, 0.974176, 0.975677, 0.977092, 0.978422, 0.979666, 0.980824, 0.981895, 0.982881, 0.983779, 0.984591, 0.985315, 0.985952, 0.986502, 0.986964, 0.987337, 0.987622, 0.987819, 0.987926, 0.987945, 0.987874, 0.987714, 0.987464, 0.987124, 0.986694, 0.986175, 0.985566, 0.984865, 0.984075, 0.983196, 0.982228, 0.981173, 0.980032, 0.978806, 0.977497, 0.976108, 0.974638, 0.973088, 0.971468, 0.969783, 0.968041, 0.966243, 0.964394, 0.962517, 0.960626, 0.95872, 0.956834, 0.954997, 0.953215, 0.951546, 0.950018, 0.948683, 0.947594, 0.946809, 0.946392, 0.946403, 0.946903, 0.947937, 0.949545, 0.95174, 0.954529, 0.957896, 0.961812, 0.966249, 0.971162, 0.976511, 0.982257, 0.988362],
        [0.000466, 0.00127, 0.002249, 0.003392, 0.004692, 0.006136, 0.007713, 0.009417, 0.011225, 0.013136, 0.015133, 0.017199, 0.019331, 0.021503, 0.023702, 0.025921, 0.028139, 0.030324, 0.032474, 0.034569, 0.03659, 0.038504, 0.040294, 0.041905, 0.043328, 0.044556, 0.045583, 0.046402, 0.047008, 0.047399, 0.047574, 0.047536, 0.047293, 0.046856, 0.046242, 0.045468, 0.044559, 0.043554, 0.042489, 0.041402, 0.040329, 0.039309, 0.0384, 0.037632, 0.03703, 0.036615, 0.036405, 0.036405, 0.036621, 0.037055, 0.037705, 0.038571, 0.039647, 0.040922, 0.042353, 0.043933, 0.045644, 0.04747, 0.049396, 0.051407, 0.05349, 0.055634, 0.057827, 0.06006, 0.062325, 0.064616, 0.066925, 0.069247, 0.071579, 0.073915, 0.076253, 0.078591, 0.080927, 0.083257, 0.08558, 0.087896, 0.090203, 0.092501, 0.09479, 0.097069, 0.099338, 0.101597, 0.103848, 0.106089, 0.108322, 0.110547, 0.112764, 0.114974, 0.117179, 0.119379, 0.121575, 0.123769, 0.12596, 0.12815, 0.130341, 0.132534, 0.134729, 0.136929, 0.139134, 0.141346, 0.143567, 0.145797, 0.148039, 0.150294, 0.152563, 0.154848, 0.157151, 0.159474, 0.161817, 0.164184, 0.166575, 0.168992, 0.171438, 0.173914, 0.176421, 0.178962, 0.181539, 0.184153, 0.186807, 0.189501, 0.192239,0.195021, 0.197851, 0.200728, 0.203656, 0.206636, 0.20967, 0.212759, 0.215906, 0.219112, 0.222378, 0.225706,0.229097, 0.232554, 0.236077, 0.239667, 0.243327, 0.247056, 0.250856, 0.254728, 0.258674, 0.262692, 0.266786, 0.270954, 0.275197, 0.279517, 0.283913, 0.288385, 0.292933, 0.297559, 0.30226, 0.307038, 0.311892, 0.316822, 0.321827, 0.326906, 0.33206, 0.337287, 0.342586, 0.347957, 0.353399, 0.358911, 0.364492, 0.37014, 0.375856,0.381636, 0.387481, 0.393389, 0.399359, 0.405389, 0.411479, 0.417627, 0.423831, 0.430091, 0.436405, 0.442772, 0.449191, 0.45566, 0.462178, 0.468744, 0.475356, 0.482014, 0.488716, 0.495462, 0.502249, 0.509078, 0.515946, 0.522853, 0.529798, 0.53678, 0.543798, 0.55085, 0.557937, 0.565057, 0.572209, 0.579392, 0.586606, 0.593849,0.601122, 0.608422, 0.61575, 0.623105, 0.630485, 0.63789, 0.64532, 0.652773, 0.66025, 0.667748, 0.675267, 0.682807, 0.690366, 0.697944, 0.70554, 0.713153, 0.720782, 0.728427, 0.736087, 0.743758, 0.751442, 0.759135, 0.766837, 0.774545, 0.782258, 0.789974, 0.797692, 0.805409, 0.813122, 0.820825, 0.828515, 0.836191, 0.843848, 0.851476, 0.859069, 0.866624, 0.874129, 0.881569, 0.888942, 0.896226, 0.903409, 0.910473, 0.917399, 0.924168, 0.930761, 0.937159, 0.943348, 0.949318, 0.955063, 0.960587, 0.965896, 0.971003, 0.975924, 0.980678, 0.985282,0.989753, 0.994109, 0.998364],
        [0.013866, 0.01857, 0.024239, 0.030909, 0.038558, 0.046836, 0.055143, 0.06346, 0.071862, 0.080282, 0.088767,0.097327, 0.10593, 0.114621, 0.123397, 0.132232, 0.141141, 0.150164, 0.159254, 0.168414, 0.177642, 0.186962,0.196354, 0.205799, 0.215289, 0.224813, 0.234358, 0.243904, 0.25343, 0.262912, 0.272321, 0.281624, 0.290788,0.299776, 0.308553, 0.317085, 0.325338, 0.333277, 0.340874, 0.348111, 0.354971, 0.361447, 0.367535, 0.373238, 0.378563, 0.383522, 0.388129, 0.3924, 0.396353, 0.400007, 0.403378, 0.406485, 0.409345, 0.411976, 0.414392,0.416608, 0.418637, 0.420491, 0.422182, 0.423721, 0.425116, 0.426377, 0.427511, 0.428524, 0.429425, 0.430217, 0.430906, 0.431497, 0.431994, 0.4324, 0.432719, 0.432955, 0.433109, 0.433183, 0.433179, 0.433098, 0.432943,0.432714, 0.432412, 0.432039, 0.431594, 0.43108, 0.430498, 0.429846, 0.429125, 0.428334, 0.427475, 0.426548,0.425552, 0.424488, 0.423356, 0.422156, 0.420887, 0.419549, 0.418142, 0.416667, 0.415123, 0.413511, 0.411829, 0.410078, 0.408258, 0.406369, 0.404411, 0.402385, 0.40029, 0.398125, 0.395891, 0.393589, 0.391219, 0.388781, 0.386276, 0.383704, 0.381065, 0.378359, 0.375586, 0.372748, 0.369846, 0.366879, 0.363849, 0.360757, 0.357603, 0.354388, 0.351113, 0.347777, 0.344383, 0.340931, 0.337424, 0.333861, 0.330245, 0.326576, 0.322856, 0.319085, 0.315266, 0.311399, 0.307485, 0.303526, 0.299523, 0.295477, 0.29139, 0.287264, 0.283099, 0.278898, 0.274661, 0.27039, 0.266085, 0.26175, 0.257383, 0.252988, 0.248564, 0.244113, 0.239636, 0.235133, 0.230606, 0.226055, 0.221482, 0.216886, 0.212268, 0.207628, 0.202968, 0.198286, 0.193584, 0.18886, 0.184116, 0.17935, 0.174563, 0.169755, 0.164924, 0.16007, 0.155193, 0.150292, 0.145367, 0.140417, 0.13544, 0.130438, 0.125409, 0.120354,0.115272, 0.110164, 0.105031, 0.099874, 0.094695, 0.089499, 0.084289, 0.079073, 0.073859, 0.068659, 0.063488, 0.058367, 0.053324, 0.048392, 0.043618, 0.03905, 0.034931, 0.031409, 0.028508, 0.02625, 0.024661, 0.02377, 0.023606, 0.024202, 0.025592, 0.027814, 0.030908, 0.034916, 0.039886, 0.045581, 0.05175, 0.058329, 0.065257, 0.072489, 0.07999, 0.087731, 0.095694, 0.103863, 0.112229, 0.120785, 0.129527, 0.138453, 0.147565, 0.156863, 0.166353, 0.176037, 0.185923, 0.196018, 0.206332, 0.216877, 0.227658, 0.238686, 0.249972, 0.261534, 0.273391,0.285546, 0.29801, 0.31082, 0.323974, 0.337475, 0.351369, 0.365627, 0.380271, 0.395289, 0.410665, 0.426373, 0.442367, 0.458592, 0.47497, 0.491426, 0.50786, 0.524203, 0.540361, 0.556275, 0.571925, 0.587206, 0.602154, 0.61676, 0.631017, 0.644924]
        )
    export inferno

## --- Other colormaps

    # Fire colormap
    fire = RGB{N0f8}.(
        vcat(fill(1,120),linsp(0.992,0.05,136)), # r
        vcat(linsp(0.9,0,120),fill(0,136)), # g
        vcat(linsp(0.9,0,120),fill(0,136)) #b
        )
    export fire

    # Water colormap
    water = RGB{N0f8}.(
        vcat(linsp(0.9,0,136),fill(0,120)), # r
        vcat(linsp(0.9,0,136),fill(0,120)), # g
        vcat(fill(1,136),linsp(0.992,0.05,120)) #b
        )
    export water

## --- Resize and interpolate colormaps

    # Linearly interpolate cmap at positions xq
    function linterp_colormap(x,cmap,xq)
        # Extract red, green, and blue vectors
        cmap_r = cmap .|> c -> c.r
        cmap_g = cmap .|> c -> c.g
        cmap_b = cmap .|> c -> c.b
        # Interpolate
        r_interp = linterp1(x,cmap_r,xq)
        g_interp = linterp1(x,cmap_g,xq)
        b_interp = linterp1(x,cmap_b,xq)
        # Convert back to a color
        return RGB.(r_interp,g_interp,b_interp)
    end
    export linterp_colormap

    function resize_colormap(cmap,n)
        cNum = length(cmap)
        return linterp_colormap(1:cNum,cmap,linsp(1,cNum,n))
    end
    export resize_colormap

## --- Map colormaps to images

    function imsc(matrix::Array,colormap::Array,cmin::Number=0,cmax::Number=0)
        Nc = length(colormap) - 1
        if cmin>=cmax
            cmin = nanminimum(matrix)
            cmax = nanmaximum(matrix)
        end
        crange = cmax - cmin
        return  matrix .|> x -> colormap[isnan(x) ? 1 : floor(UInt, min(max((x-cmin)/crange*Nc,0), Nc))+1]
    end
    export imsc

    function imsc_log10f(matrix::Array,from::Number,colormap::Array,cmin::Number=0,cmax::Number=0)
        Nc = length(colormap) - 1
        if cmin>=cmax
            cmin = log10f(nanminimum(matrix),from)
            cmax = log10f(nanmaximum(matrix),from)
        end
        crange = cmax - cmin
        return  matrix .|> x -> colormap[isnan(x) ? 1 : floor(UInt, min(max((log10f(x,from)-cmin)/crange*Nc,0), Nc))+1]
    end
    export imsc_log10f

## -- End of File

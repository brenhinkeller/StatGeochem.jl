@test find_litho1_property(40:45, fill(-100,6), :ice, :thickness) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@test find_litho1_property(40:45, fill(-100,6), :water, :thickness) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@test find_litho1_property(40:45, fill(-100,6), :upper_sediments, :vp) ≈ [2500.0, 2500.0, 2500.0, 2500.0, 2500.0, 2500.0]
@test find_litho1_property(40:45, fill(-100,6), :upper_sediments, :vs) ≈ [1070.0, 1070.0, 1070.0, 1070.0, 1070.0, 1070.0]
@test find_litho1_property(40:45, fill(-100,6), :upper_sediments, :rho) ≈ [2110.0, 2110.0, 2110.0, 2110.0, 2110.0, 2110.0]
@test find_litho1_property(40:45, fill(-100,6), :upper_sediments, :bottom) ≈ [-0.25, -0.352, -0.301, -0.208, -0.063, 0.045]
@test find_litho1_property(40:45, fill(-100,6), :upper_sediments, :thickness) ≈ [0.5, 0.5, 0.5, 0.5, 0.5, 0.608]

@test find_litho1_property(40:45, fill(-100,6), :middle_sediments, :vp) ≈ [4000.0, 4000.0, 4000.0, 4000.0, 4000.0, 4000.0]
@test find_litho1_property(40:45, fill(-100,6), :middle_sediments, :vs) ≈ [2130.0, 2130.0, 2130.0, 2130.0, 2130.0, 2130.0]
@test find_litho1_property(40:45, fill(-100,6), :middle_sediments, :rho) ≈ [2370.0, 2370.0, 2370.0, 2370.0, 2370.0, 2370.0]
@test find_litho1_property(40:45, fill(-100,6), :middle_sediments, :bottom) ≈ [0.662, 0.382, 0.453, 0.297, 0.235, 0.37]
@test find_litho1_property(40:45, fill(-100,6), :middle_sediments, :thickness) ≈ [0.912, 0.734, 0.754, 0.505, 0.298, 0.325]

@test find_litho1_property(40:45, fill(-100,6), :lower_sediments, :thickness) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@test find_litho1_property(40:45, fill(-100,6), :upper_crust, :vp) ≈ [6252.38, 6225.32, 6094.56, 6233.38, 6180.0, 6219.24]
@test find_litho1_property(40:45, fill(-100,6), :upper_crust, :vs) ≈ [3626.89, 3610.68, 3534.85, 3615.36, 3584.4, 3607.79]
@test find_litho1_property(40:45, fill(-100,6), :upper_crust, :rho) ≈ [2791.52, 2780.64, 2722.24, 2784.24, 2760.4, 2776.46]
@test find_litho1_property(40:45, fill(-100,6), :upper_crust, :bottom) ≈ [16.877, 16.043, 14.313, 16.055, 15.497, 15.499]
@test find_litho1_property(40:45, fill(-100,6), :upper_crust, :thickness) ≈ [16.215, 15.661, 13.86, 15.758, 15.262, 15.129]

@test find_litho1_property(40:45, fill(-100,6), :middle_crust, :vp) ≈ [6641.63, 6640.34, 6500.87, 6648.94, 6592.0, 6600.29]
@test find_litho1_property(40:45, fill(-100,6), :middle_crust, :vs) ≈ [3821.52, 3818.19, 3738.0, 3823.14, 3790.4, 3798.31]
@test find_litho1_property(40:45, fill(-100,6), :middle_crust, :rho) ≈ [2888.84, 2884.4, 2823.81, 2888.13, 2863.4, 2871.72]
@test find_litho1_property(40:45, fill(-100,6), :middle_crust, :bottom) ≈ [34.838, 33.669, 29.909, 33.783, 32.659, 32.24]
@test find_litho1_property(40:45, fill(-100,6), :middle_crust, :thickness) ≈ [17.961, 17.626, 15.596, 17.728, 17.162, 16.741]

@test find_litho1_property(40:45, fill(-100,6), :lower_crust, :vp) ≈ [7147.59, 7159.11, 7008.75, 7168.39, 7107.0, 7100.2]
@test find_litho1_property(40:45, fill(-100,6), :lower_crust, :vs) ≈ [4074.5, 4077.58, 3991.94, 4082.86, 4047.9, 4048.27]
@test find_litho1_property(40:45, fill(-100,6), :lower_crust, :rho) ≈ [3028.95, 3029.65, 2966.02, 3033.58, 3007.6, 3009.81]
@test find_litho1_property(40:45, fill(-100,6), :lower_crust, :bottom) ≈ [51.119, 49.339, 43.779, 49.55, 47.917, 47.426]
@test find_litho1_property(40:45, fill(-100,6), :lower_crust, :thickness) ≈ [16.281, 15.67, 13.87, 15.767, 15.258, 15.186]

@test find_litho1_property(40:45, fill(-100,6), :sclm, :vp) ≈ [8171.22, 8173.6, 8194.44, 8242.08, 8251.74, 8220.29]
@test find_litho1_property(40:45, fill(-100,6), :sclm, :vs) ≈ [4655.96, 4657.32, 4669.19, 4696.34, 4701.85, 4683.93]
@test find_litho1_property(40:45, fill(-100,6), :sclm, :rho) ≈ [3300.0, 3300.0, 3300.0, 3300.0, 3300.0, 3300.0]
@test find_litho1_property(40:45, fill(-100,6), :sclm, :bottom) ≈ [182.204, 234.025, 241.833, 246.228, 243.389, 251.825]
@test find_litho1_property(40:45, fill(-100,6), :sclm, :thickness) ≈ [131.085, 184.686, 198.054, 196.678, 195.472, 204.399]

@test find_litho1_property(40:45, fill(-100,6), :asthenosphere, :vp) ≈ [8007.79, 8010.13, 8030.55, 8077.24, 8086.7, 8055.88]
@test find_litho1_property(40:45, fill(-100,6), :asthenosphere, :vs) ≈ [4352.06, 4353.33, 4364.43, 4389.8, 4394.94, 4378.2]
@test find_litho1_property(40:45, fill(-100,6), :asthenosphere, :rho) ≈ [3300.0, 3300.0, 3300.0, 3300.0, 3300.0, 3300.0]
@test find_litho1_property(40:45, fill(-100,6), 10, :vp) ≈ [8007.79, 8010.13, 8030.55, 8077.24, 8086.7, 8055.88]

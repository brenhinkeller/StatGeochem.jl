## -- Common elemental masses in g/mol (AMU/atom)

    const molarmass = Dict{String,Float64}("Pd" => 106.421,"Fl" => 289.0,"Nb" => 92.906372,"C" => 12.011,"P" => 30.9737619985,"Ag" => 107.86822,"Gd" => 157.253,"Si" => 28.085,"Ru" => 101.072,"At" => 210.0,"Uus" => 294.0,"Sb" => 121.7601,"Cs" => 132.905451966,"Cn" => 285.0,"Uut" => 285.0,"Be" => 9.01218315,"Ac" => 227.0,"Cf" => 251.0,"Bh" => 270.0,"Sr" => 87.621,"Ga" => 69.7231,"Ta" => 180.947882,"Te" => 127.603,"Np" => 237.0,"Lr" => 262.0,"Pu" => 244.0,"U" => 238.028913,"Kr" => 83.7982,"Y" => 88.905842,"Sg" => 271.0,"Ca" => 40.0784,"Au" => 196.9665695,"K" => 39.09831,"Rn" => 222.0,"Ra" => 226.0,"Ce" => 140.1161,"Uuo" => 294.0,"V" => 50.94151,"Fr" => 223.0,"Mo" => 95.951,"Pr" => 140.907662,"Th" => 232.03774,"Br" => 79.904,"Zn" => 65.382,"He" => 4.0026022,"Sc" => 44.9559085,"H" => 1.008,"Al" => 26.98153857,"S" => 32.06,"Uup" => 289.0,"Ar" => 39.9481,"Ge" => 72.6308,"Er" => 167.2593,"Fe" => 55.8452,"Mg" => 24.305,"F" => 18.9984031636,"La" => 138.905477,"Rf" => 267.0,"W" => 183.841,"Li" => 6.94,"Dy" => 162.5001,"O" => 15.999,"B" => 10.81,"Bi" => 208.980401,"Mn" => 54.9380443,"Re" => 186.2071,"Db" => 270.0,"Hf" => 178.492,"Cm" => 247.0,"Cl" => 35.45,"In" => 114.8181,"Ds" => 281.0,"Rb" => 85.46783,"Po" => 209.0,"Lv" => 293.0,"Am" => 243.0,"Pa" => 231.035882,"Se" => 78.9718,"Ba" => 137.3277,"Nd" => 144.2423,"Pm" => 145.0,"Rh" => 102.905502,"Ti" => 47.8671,"Tb" => 158.925352,"Hs" => 277.0,"Zr" => 91.2242,"Sm" => 150.362,"Cr" => 51.99616,"Cu" => 63.5463,"Fm" => 257.0,"Tc" => 97.0,"Tl" => 204.38,"Ne" => 20.17976,"Hg" => 200.5923,"Mt" => 276.0,"N" => 14.007,"Es" => 252.0,"Yb" => 173.0455,"Lu" => 174.96681,"Eu" => 151.9641,"Na" => 22.989769282,"No" => 259.0,"Os" => 190.233,"Ni" => 58.69344,"Ho" => 164.930332,"Co" => 58.9331944,"Md" => 258.0,"Ir" => 192.2173,"Pt" => 195.0849,"Tm" => 168.934222,"As" => 74.9215956,"Sn" => 118.7107,"Xe" => 131.2936,"I" => 126.904473,"Cd" => 112.4144,"Pb" => 207.21,"Rg" => 282.0,"Bk" => 247.0)
    export molarmass

    const molarmasspercation = Dict{String,Float64}("MgO" => 40.304,"SO2" => 64.058,"TiO2" => 79.8651,"CaO" => 56.077400000000004,"P2O5" => 70.9712619985,"NiO" => 74.69244,"CoO" => 74.9321944,"SiO2" => 60.083,"Al2O3" => 50.980038570000005,"CO2" => 44.009,"Cr2O3" => 75.99466000000001,"FeO" => 71.8442,"K2O" => 47.097809999999996,"Fe2O3" => 79.8437,"H2O" => 9.0075,"MnO" => 70.9370443,"Na2O" => 30.989269282000002)
    export molarmasspercation

    """
    julia```
    ionicradius::NamedTuple
    ```
    A named tuple containing the ionic radii, in picometers, for the naturally ocurring 
    elements in their geologically common redox states. Where multiple redox states may 
    be expected for a single element in nature, they are suffixed to the atomic symbol, 
    e.g. `Fe2` vs `Fe3`.

    Radii are those reported as "crystal" ionic radii by Shannon [1],  as tabulated at
    https://en.wikipedia.org/wiki/Ionic_radius#Tables

    [1]  R. D. Shannon (1976). "Revised effective ionic radii and systematic studies of interatomic 
    distances in halides and chalcogenides". Acta Crystallogr A. 32 (5): 751-767.
    doi:10.1107/S0567739476001551
    """
    const ionicradius = (;
        H = -4.,
        Li = 90.,
        Be = 59.,
        B = 41.,
        C = 30.,
        N = 27.,
        O = 126.,
        F = 119.,
        Na = 116.,
        Mg = 86.,
        Al = 67.5,
        Si = 54.,
        P = 52.,
        S = 170.,
        Cl = 167.,
        K = 152.,
        Ca = 114.,
        Sc = 88.5,
        Ti = 74.5,
        V3 = 78.,
        V5 = 68.,
        Cr3 = 75.5,
        Cr6 = 58.,
        Mn2 = 81.,
        Mn4 = 67.,
        Fe2 = 75.,
        Fe3 = 69.,
        Co = 79.,
        Ni = 83.,
        Cu1 = 91.,
        Cu2 = 87.,
        Zn = 88.,
        Ga = 76.,
        Ge = 67.,
        As = 60.,
        Se = 184.,
        Br = 182.,
        Rb = 166.,
        Sr = 132.,
        Y = 104.,
        Zr = 86.,
        Nb = 78.,
        Mo4 = 79.,
        Mo6 = 73.,
        Ru = 82.,
        Rh = 80.5,
        Pd = 100.,
        Ag = 129.,
        Cd = 109.,
        In = 94.,
        Sn = 83.,
        Sb = 90.,
        Te = 207.,
        I = 206.,
        Cs = 167.,
        Ba = 149.,
        La = 117.2,
        Ce3 = 115.,
        Ce4 = 101.,
        Pr = 113.,
        Nd = 112.3,
        Sm = 109.8,
        Eu2 = 131,
        Eu3 = 108.7,
        Gd = 107.8,
        Tb = 106.3,
        Dy = 105.2,
        Ho = 104.1,
        Er = 103.,
        Tm = 102.,
        Yb = 100.8,
        Lu = 100.1,
        Hf = 85.,
        Ta = 78.,
        W = 74.,
        Re4 = 77.,
        Re7 = 67.,
        Os4 = 77.,
        Os8 = 53.,
        Ir = 76.5,
        Pt = 76.5,
        Au = 151.,
        Hg = 116.,
        Tl = 164.,
        Pb = 133.,
        Bi = 117.,
        Th = 108.,
        U4 = 103.,
        U6 = 87.,
    )
    export ionicradius

    """
    julia```
    ioniccharge::NamedTuple
    ```
    A named tuple containing the ionic charges corresponding to the ionic radii in `ionicradius`.
    """
    const ioniccharge = (;
        H = +1,
        Li = +1,
        Be = +2,
        B = +3,
        C = +4,
        N = +5,
        O = -2,
        F = -1,
        Na = +1,
        Mg = +2,
        Al = +3,
        Si = +4,
        P = +5,
        S = -2,
        Cl = -1,
        K = +1,
        Ca = +2,
        Sc = +3,
        Ti = +4,
        V3 = +3,
        V5 = +5,
        Cr3 = +3,
        Cr6 = +6,
        Mn2 = +2,
        Mn4 = +4,
        Fe2 = +2,
        Fe3 = +3,
        Co = +2,
        Ni = +2,
        Cu1 = +1,
        Cu2 = +2,
        Zn = +2,
        Ga = +3,
        Ge = +4,
        As = +5,
        Se = -2,
        Br = -1,
        Rb = +1,
        Sr = +2,
        Y = +3,
        Zr = +4,
        Nb = +5,
        Mo4 = +4,
        Mo6 = +6,
        Ru = +3,
        Rh = +3,
        Pd = +2,
        Ag = +1,
        Cd = +2,
        In = +3,
        Sn = +4,
        Sb = +3,
        Te = -2,
        I = -1,
        Cs = +1,
        Ba = +2,
        La = +3,
        Ce3 = +3,
        Ce4 = +4,
        Pr = +3,
        Nd = +3,
        Sm = +3,
        Eu2 = +2,
        Eu3 = +3,
        Gd = +3,
        Tb = +3,
        Dy = +3,
        Ho = +3,
        Er = +3,
        Tm = +3,
        Yb = +3,
        Lu = +3,
        Hf = +4,
        Ta = +5,
        W = +6,
        Re4 = +4,
        Re7 = +7,
        Os4 = +4,
        Os8 = +8,
        Ir = +4,
        Pt = +4,
        Au = +1,
        Hg = +2,
        Tl = +1,
        Pb = +2,
        Bi = +3,
        Th = +4,
        U4 = +4,
        U6 = +6,
    )
    export ioniccharge

    """
    julia```
    atomicnumber::NamedTuple
    ```
    A named tuple containing the atomic number of each element.
    """
    const atomicnumber = (;
        H = 1,
        He = 2,
        Li = 3,
        Be = 4,
        B = 5,
        C = 6,
        N = 7,
        O = 8,
        F = 9,
        Ne = 10,
        Na = 11,
        Mg = 12,
        Al = 13,
        Si = 14,
        P = 15,
        S = 16,
        Cl = 17,
        Ar = 18,
        K = 19,
        Ca = 20,
        Sc = 21,
        Ti = 22,
        V = 23,
        Cr = 24,
        Mn = 25,
        Fe = 26,
        Co = 27,
        Ni = 28,
        Cu = 29,
        Zn = 30,
        Ga = 31,
        Ge = 32,
        As = 33,
        Se = 34,
        Br = 35,
        Kr = 36,
        Rb = 37,
        Sr = 38,
        Y = 39,
        Zr = 40,
        Nb = 41,
        Mo = 42,
        Tc = 43,
        Ru = 44,
        Rh = 45,
        Pd = 46,
        Ag = 47,
        Cd = 48,
        In = 49,
        Sn = 50,
        Sb = 51,
        Te = 52,
        I = 53,
        Xe = 54,
        Cs = 55,
        Ba = 56,
        La = 57,
        Ce = 58,
        Pr = 59,
        Nd = 60,
        Pm = 61,
        Sm = 62,
        Eu = 63,
        Gd = 64,
        Tb = 65,
        Dy = 66,
        Ho = 67,
        Er = 68,
        Tm = 69,
        Yb = 70,
        Lu = 71,
        Hf = 72,
        Ta = 73,
        W = 74,
        Re = 75,
        Os = 76,
        Ir = 77,
        Pt = 78,
        Au = 79,
        Hg = 80,
        Tl = 81,
        Pb = 82,
        Bi = 83,
        Po = 84,
        At = 85,
        Rn = 86,
        Fr = 87,
        Ra = 88,
        Ac = 89,
        Th = 90,
        Pa = 91,
        U = 92,
        Np = 93,
        Pu = 94,
        Am = 95,
        Cm = 96,
        Bk = 97,
        Cf = 98,
        Es = 99,
        Fm = 100,
        Md = 101,
        No = 102,
        Lr = 103,
        Rf = 104,
        Db = 105,
        Sg = 106,
        Bh = 107,
        Hs = 108,
        Mt = 109,
        Ds = 110,
        Rg = 111,
        Cn = 112,
        Nh = 113,
        Fl = 114,
        Mc = 115,
        Lv = 116,
        Ts = 117,
        Og = 118,
    )
    export atomicnumber
    
## ---

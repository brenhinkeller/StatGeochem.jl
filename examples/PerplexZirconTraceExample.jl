## --- Load the StatGeochem package which has the resampling functions we'll want

    using StatGeochem
    using Plots

## --- Download and unzip Keller and Schoene (2012) dataset

    if ~isfile("ign.h5") # Unless it already exists
        download("https://storage.googleapis.com/statgeochem/ign.h5.gz","./ign.h5.gz")
        download("https://storage.googleapis.com/statgeochem/err2srel.csv","./err2srel.csv")
        run(`gunzip -f ign.h5.gz`) # Unzip file
    end

    # Read HDF5 file
    using HDF5
    ign = TupleDataset(h5read("ign.h5","vars"))

## --- Clean up and calculate average granitoid
    ign.FeOT .= feoconversion.(ign.FeO, ign.Fe2O3, ign.FeOT, ign.Fe2O3)
    ign = (;ign..., P=fill(NaN, size(ign.P2O5)), Mn=fill(NaN, size(ign.P2O5)), H2O=fill(4.0, size(ign.P2O5)))
    metalconversion!(ign)
    oxideconversion!(ign)

    igncomp = CompositionArray{NCKFMASHTOtrace{Float64}}((;ign..., FeO=ign.FeOT, O2=0.18/4*ign.FeOT))
    t = (62 .< ign.SiO2 .< 74) .& (ign.Age .< 541)
    avegranite = nanmean(igncomp[t])

## --- Prepare for Perple_X run

    scratchdir = "./"
    dataset="hp633ver.dat"
    melt_model = "melt(HGPH)"
    fsp_model = "Fsp(C1)"
    plag_model = "Pl(I1,HP)"
    solution_phases = melt_model*"\n"*fsp_model*"\n"*plag_model*"\n"*"Sp(HGP)\nGt(HGP)\nO(HGP)\nOpx(HGP)\nCpx(HGP)\nCrd(HGP)\nBi(HGP)\nMica(W)\nEp(HP)\ncAmph(G)\nIlm(WPH)\nChl(W)\n"

    excludes ="ged\nfanth\ngl\n" # Use solutions instead of endmembers
    excludes *= "zrc\n" # Also exclude zircon since we'll calculate it ourselves

## --- Run PerpleX

    scratchdir = "./" # Where to store perplex files and output 
    index = 1 
    perplex_configure_isobar(scratchdir, avegranite, 4000., (500+273.15, 1000+273.15);
        dataset,
        solution_phases,
        excludes,
        index,
        npoints = 256,
        nonlinear_subdivision = true,
    )

    results = perplextrace_query(scratchdir, avegranite; 
        index,
        npoints=256,
    )

## --- End of File
## -- Perplex interface: 1. Configuration

"""
```julia
perplex_configure_geotherm(scratchdir::String, composition, [elements],
    \tP_range=(280,28000), T_surf::Number=273.15, geotherm::Number=0.1;
    \tdataset::String="hp633ver.dat",
    \tindex::Integer=1,
    \tnpoints::Integer=100,
    \tsolution_phases::String="melt(HGPH)\\nPl(I1,HP)\\nFsp(C1)\\nSp(HGP)\\nGt(HGP)\\nO(HGP)\\nOpx(HGP)\\nCpx(HGP)\\nCrd(HGP)\\nBi(HGP)\\nMica(W)\\nEp(HP)\\ncAmph(G)\\nIlm(WPH)\\nChl(W)\\n",
    \texcludes::String="ged\\nfanth\\ng\\n",
    \tmode_basis::String="wt",  #["vol", "wt", "mol"]
    \tcomposition_basis::String="wt",  #["vol", "wt", "mol"]
    \tfluid_eos::Integer=5)
```

Set up a PerpleX calculation for a single bulk composition along a specified
geothermal gradient and pressure (depth) range. P specified in bar and T_surf
in Kelvin, with geothermal gradient in units of Kelvin/bar
"""
function perplex_configure_geotherm(scratchdir, composition::AbstractComposition, args...; dataset="hp633ver.dat", kwargs...)
    data = majorelementvalues(composition)
    elements = String.(majorelements(composition))
    perplex_dataset_uppercase(dataset) && (elements = uppercase.(elements))
    return perplex_configure_geotherm(scratchdir, data, elements, args...; dataset, kwargs...)
end
function perplex_configure_geotherm(scratchdir::String, composition::Collection{Number}, elements::Collection{<:AbstractString},
        P_range::NTuple{2,Number}=(280,28000), T_surf::Number=273.15, geotherm::Number=0.1;
        dataset::String="hp633ver.dat",
        index::Integer=1,
        npoints::Integer=100,
        solution_phases::String="melt(HGPH)\nPl(I1,HP)\nFsp(C1)\nSp(HGP)\nGt(HGP)\nO(HGP)\nOpx(HGP)\nCpx(HGP)\nCrd(HGP)\nBi(HGP)\nMica(W)\nEp(HP)\ncAmph(G)\nIlm(WPH)\nChl(W)\n",
        excludes::String="ged\nfanth\ngl\n",
        mode_basis::String="wt",
        composition_basis::String="wt",
        fluid_eos::Integer=5
    )
    @assert eachindex(composition) == eachindex(elements) "`composition` and `elements` must match"
    build = joinpath(Perple_X_jll.PATH[], "build")# path to PerpleX build
    vertex = joinpath(Perple_X_jll.PATH[], "vertex")# path to PerpleX vertex

    #Configure working directory
    prefix = joinpath(scratchdir, "out$(index)/")
    system("rm -rf $prefix; mkdir -p $prefix")

    # Place required data files
    artifact_path = artifact"perplex-datafiles"
    system("cp $(joinpath(artifact_path, "perplex-datafiles/$dataset")) $prefix") 
    system("cp $(joinpath(artifact_path,"perplex-datafiles/perplex_option.dat")) $prefix")
    system("cp $(joinpath(artifact_path,"perplex-datafiles/solution_model.dat")) $prefix")

    # Edit perplex_option.dat to specify number of nodes at which to solve
    system("sed -e \"s/1d_path .*|/1d_path                   $npoints $npoints |/\" -i.backup $(prefix)perplex_option.dat")

    # Edit perplex_option.dat to output all seismic properties
    #println("editing perplex options ")
    system("sed -e \"s/seismic_output .*|/seismic_output                   all |/\" -i.backup $(prefix)perplex_option.dat")

    # Specify whether we want volume or weight percentages
    system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

    # Create build batch file.
    fp = open(prefix*"build.bat", "w")

    # Name, components, and basic options. P-T conditions.
    elementstring = ""
    for i in eachindex(composition)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(first(P_range))\n$(last(P_range))\ny\n") # v7.1.6/7.1.8
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n5\n3\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(first(P_range))\n$(last(P_range))\ny\n") # v6.8.1

    # Whole-rock composition
    for i ∈ eachindex(composition)
        if !isnan(composition[i])
            write(fp,"$(composition[i]) ")
        end
    end

    # Solution model
    if length(excludes) > 0
        write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\n$fluid_eos\nGeothermal") 
    else
        write(fp,"\nn\nn\ny\nsolution_model.dat\n$(solution_phases)\n$fluid_eos\nGeothermal") 
    end
    close(fp)

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # build PerpleX problem definition
    system("export $ldpathadjust; cd $prefix; $build < build.bat > build.log")

    # Run PerpleX vertex calculations
    result = system("export $ldpathadjust; cd $prefix; echo $index | $vertex > vertex.log")
    return result
end
export perplex_configure_geotherm

"""
```julia
perplex_configure_isobar(scratchdir::String, composition, [elements],
    \tP::Number=10000, T_range::NTuple{2,Number}=(500+273.15, 1500+273.15);
    \tdataset::String="hp633ver.dat",
    \tindex::Integer=1,
    \tnpoints::Integer=100,
    \tsolution_phases::String="melt(HGPH)\\nPl(I1,HP)\\nFsp(C1)\\nSp(HGP)\\nGt(HGP)\\nO(HGP)\\nOpx(HGP)\\nCpx(HGP)\\nCrd(HGP)\\nBi(HGP)\\nMica(W)\\nEp(HP)\\ncAmph(G)\\nIlm(WPH)\\nChl(W)\\n",
    \texcludes::String="ged\\nfanth\\ng\\n",
    \tmode_basis::String="wt",  #["vol", "wt", "mol"]
    \tcomposition_basis::String="wt",  #["wt", "mol"]
    \tnonlinear_subdivision::Bool=false,
    \tfluid_eos::Integer=5)
```

Set up a PerpleX calculation for a single bulk composition along a specified
isobaric temperature gradient. P specified in bar and T_range in Kelvin
"""
function perplex_configure_isobar(scratchdir, composition::AbstractComposition, args...; dataset="hp633ver.dat", kwargs...)
    data = majorelementvalues(composition)
    elements = String.(majorelements(composition))
    perplex_dataset_uppercase(dataset) && (elements = uppercase.(elements))
    return perplex_configure_isobar(scratchdir, data, elements, args...; dataset, kwargs...)
end
function perplex_configure_isobar(scratchdir::String, composition::Collection{Number}, elements::Collection{<:AbstractString},
        P::Number=10000, T_range::NTuple{2,Number}=(500+273.15, 1500+273.15);
        dataset::String="hp633ver.dat",
        index::Integer=1,
        npoints::Integer=100,
        solution_phases::String="melt(HGPH)\nPl(I1,HP)\nFsp(C1)\nSp(HGP)\nGt(HGP)\nO(HGP)\nOpx(HGP)\nCpx(HGP)\nCrd(HGP)\nBi(HGP)\nMica(W)\nEp(HP)\ncAmph(G)\nIlm(WPH)\nChl(W)\n",
        excludes::String="ged\nfanth\ngl\n",
        mode_basis::String="wt",
        composition_basis::String="wt",
        nonlinear_subdivision::Bool=false,
        fluid_eos::Integer=5
    )
    @assert eachindex(composition) == eachindex(elements) "`composition` and `elements` must match"
    build = joinpath(Perple_X_jll.PATH[], "build")# path to PerpleX build
    vertex = joinpath(Perple_X_jll.PATH[], "vertex")# path to PerpleX vertex

    #Configure working directory
    prefix = joinpath(scratchdir, "out$(index)/")
    system("rm -rf $prefix; mkdir -p $prefix")

    # Place required data files
    artifact_path = artifact"perplex-datafiles"
    system("cp $(joinpath(artifact_path, "perplex-datafiles/$dataset")) $prefix") 
    system("cp $(joinpath(artifact_path,"perplex-datafiles/perplex_option.dat")) $prefix")
    system("cp $(joinpath(artifact_path,"perplex-datafiles/solution_model.dat")) $prefix")

    # Edit perplex_option.dat to specify number of nodes at which to solve
    system("sed -e \"s/1d_path .*|/1d_path                   $npoints $npoints |/\" -i.backup $(prefix)perplex_option.dat")

    # Specify whether we want volume or weight percentages
    system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

    # Turn on nonlinear subdivision and change resolution
    if nonlinear_subdivision
        system("sed -e \"s/non_linear_switch .*|/non_linear_switch              T |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s:initial_resolution .*|:initial_resolution        1/2 1/4 |:\" -i.backup $(prefix)perplex_option.dat")
    end

    # Create build batch file
    # Options based on Perplex v7.1.6
    fp = open(prefix*"build.bat", "w")

    # Name, components, and basic options. P-T conditions.
    elementstring = ""
    for i in eachindex(composition)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\nn\nn\n2\n$(first(T_range))\n$(last(T_range))\n$P\ny\n") # v7.1.8
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n$fluid_eos\nn\nn\n2\n$(first(T_range))\n$(last(T_range))\n$P\ny\n") # v7.1.6
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n$fluid_eos\n3\nn\nn\n2\n$(first(T_range))\n$(last(T_range))\n$P\ny\n") # v6.8.1

    # Whole-rock composition
    for i ∈ eachindex(composition)
        if !isnan(composition[i])
            write(fp,"$(composition[i]) ")
        end
    end

    # Solution model
    write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\n$fluid_eos\nIsobaric\n")
    close(fp)

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"
    # build PerpleX problem definition
    system("export $ldpathadjust; cd $prefix; $build < build.bat > build.log")

    # Run PerpleX vertex calculations
    result = system("export $ldpathadjust; cd $prefix; printf \"$index\n\" | $vertex > vertex.log")
    return result
end
export perplex_configure_isobar

"""
```julia
perplex_configure_path(scratchdir::String, composition::Collection{Number}, PTdir::String="", PTfilename::String="",
    \telements::String=("SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O"),
    \tT_range::NTuple{2,Number}=(500+273.15, 1050+273.15);
    \tdataset::String="hp633ver.dat",
    \tindex::Integer=1,
    \tsolution_phases::String="melt(HGPH)\\nPl(I1,HP)\\nFsp(C1)\\nSp(HGP)\\nGt(HGP)\\nO(HGP)\\nOpx(HGP)\\nCpx(HGP)\\nCrd(HGP)\\nBi(HGP)\\nMica(W)\\nEp(HP)\\ncAmph(G)\\nIlm(WPH)\\nChl(W)\\n",
    \texcludes::String="ged\\nfanth\\ng\\n",
    \tmode_basis::String="wt",  #["vol", "wt", "mol"]
    \tcomposition_basis::String="wt",  #["wt", "mol"]
    \tnonlinear_subdivision::Bool=false,
    \tfluid_eos::Integer=5,
    \tfractionate::Integer=0)
```

Set up a PerpleX calculation for a single bulk composition along a specified
pressure–temperature path with T as the independent variable. 

P specified in bar and T_range in Kelvin
"""
function perplex_configure_path(scratchdir::String, composition::Collection{Number}, PTdir::String="", PTfilename = "",
    elements::Collection{String}=("SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O"),
    T_range::NTuple{2,Number}=(500+273.15, 1050+273.15);
    dataset::String="hp633ver.dat",
    index::Integer=1,
    solution_phases::String="melt(HGPH)\nPl(I1,HP)\nFsp(C1)\nSp(HGP)\nGt(HGP)\nO(HGP)\nOpx(HGP)\nCpx(HGP)\nCrd(HGP)\nBi(HGP)\nMica(W)\nEp(HP)\ncAmph(G)\nIlm(WPH)\nChl(W)\n",
    excludes::String="ged\nfanth\ngl\n",
    mode_basis::String="wt",  #["vol", "wt", "mol"]
    composition_basis::String="wt",  #["wt", "mol"]
    nonlinear_subdivision::Bool=false,
    fluid_eos::Integer=5,
    fractionate::Integer=0,
    )
    @assert eachindex(composition) == eachindex(elements) "`composition` and `elements` must match"
    build = joinpath(Perple_X_jll.PATH[], "build")# path to PerpleX build
    vertex = joinpath(Perple_X_jll.PATH[], "vertex")# path to PerpleX vertex

    # Configure working directory
    prefix = joinpath(scratchdir, "out$(index)/")
    system("rm -rf $prefix; mkdir -p $prefix")

    # Place required data files
    artifact_path = artifact"perplex-datafiles"
    system("cp $(joinpath(artifact_path, "perplex-datafiles/$dataset")) $prefix") 
    system("cp $(joinpath(artifact_path,"perplex-datafiles/perplex_option.dat")) $prefix")
    system("cp $(joinpath(artifact_path,"perplex-datafiles/solution_model.dat")) $prefix")
    
    # Specify whether we want volume or weight percentages
    system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
    
    # Turn on nonlinear subdivision and change resolution
    if nonlinear_subdivision
        system("sed -e \"s/non_linear_switch .*|/non_linear_switch              T |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s:initial_resolution .*|:initial_resolution        1/2 1/4 |:\" -i.backup $(prefix)perplex_option.dat")
    end

    # Create default P–T.dat path if one is not provided
    # TODO: change this to be an interpolated PT path?
    if PTdir == ""
        # Input parameters
        P_range = (2000, 6000, 10000, 14000, 18000) #bar
        T_range = (550+273.15, 1050+273.15) #K
        T_int = 10 #Interval for T 

        T = T_range[1]:T_int:T_range[2]
        P = zeros(length(T))

        for i in 1:length(T)
            if i == length(T)
                P[i] = P_range[end]
            else
                P[i] = P_range[floor(Int64, (i/length(T)) * length(P_range) + 1)]
            end 
        end

        # Save P–T path as .dat file
        # Apparently you need to have it as T and then P despite what Perplex tells you
        PTfilename = "P–T.dat"
        PTfile = joinpath(prefix, PTfilename)
        open(PTfile, "w") do file
            for i in zip(T, P)
                write(file, "$(i[1])\t$(i[2])\n")
            end
        end
    else 
        system("cp $(PTdir) $prefix")
    end

    # Create build batch file
    fp = open(prefix*"build.bat", "w")

    # Name, components, and basic options. P-T conditions.
    elementstring = ""
    for i in eachindex(composition)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\ny\n$PTfilename\ny\n") #7.1.8
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n5\ny\n$PTfilename\ny\n") #7.1.6
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n$fluid_eos\ny\n$PTfilename\n2\ny\n") #6.8.7

    # Whole-rock composition
    for i ∈ eachindex(composition)
        if !isnan(composition[i])
            write(fp,"$(composition[i]) ")
        end
    end

    # Solution model
    write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\n5\nP-T Path") #7.1.8
    # write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\nP-T Path") #7.1.6

    close(fp)

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"
  
    # build PerpleX problem definition
    system("export $ldpathadjust; cd $prefix; $build < build.bat > build.log")

    # Run PerpleX vertex calculations
    result = system("export $ldpathadjust; cd $prefix; printf \"$index\n$fractionate\n\" | $vertex > vertex.log")

    return result
end
export perplex_configure_path

"""
```julia
perplex_configure_pseudosection(scratchdir::String, composition, [elements::Collection{String}],
    \tP::NTuple{2,Number}=(280, 28000), T::NTuple{2,Number}=(273.15, 1500+273.15);
    \tdataset::String="hp633ver.dat",
    \tindex::Integer=1,
    \txnodes::Integer=42,
    \tynodes::Integer=42,
    \tsolution_phases::String="melt(HGPH)\\nPl(I1,HP)\\nFsp(C1)\\nSp(HGP)\\nGt(HGP)\\nO(HGP)\\nOpx(HGP)\\nCpx(HGP)\\nCrd(HGP)\\nBi(HGP)\\nMica(W)\\nEp(HP)\\ncAmph(G)\\nIlm(WPH)\\nChl(W)\\n",
    \texcludes::String="ged\\nfanth\\ng\\n",
    \tmode_basis::String="wt", #["vol", "wt", "mol"]
    \tcomposition_basis::String="wt", #["wt", "mol"]
    \tfluid_eos::Number=5)
```

Set up a PerpleX calculation for a single bulk composition across an entire
2d P-T space. P specified in bar and T in Kelvin.
"""
function perplex_configure_pseudosection(scratchdir, composition::AbstractComposition, args...; dataset="hp633ver.dat", kwargs...)
    data = majorelementvalues(composition)
    elements = String.(majorelements(composition))
    perplex_dataset_uppercase(dataset) && (elements = uppercase.(elements))
    return perplex_configure_pseudosection(scratchdir, data, elements, args...; dataset, kwargs...)
end
function perplex_configure_pseudosection(scratchdir::String, composition::Collection{Number}, elements::Collection{<:AbstractString},
        P::NTuple{2,Number}=(280, 28000), T::NTuple{2,Number}=(273.15, 1500+273.15);
        dataset::String="hp633ver.dat",
        index::Integer=1,
        xnodes::Integer=42,
        ynodes::Integer=42,
        solution_phases::String="melt(HGPH)\nPl(I1,HP)\nFsp(C1)\nSp(HGP)\nGt(HGP)\nO(HGP)\nOpx(HGP)\nCpx(HGP)\nCrd(HGP)\nBi(HGP)\nMica(W)\nEp(HP)\ncAmph(G)\nIlm(WPH)\nChl(W)\n",
        excludes::String="ged\nfanth\ngl\n",
        mode_basis::String="wt",
        composition_basis::String="wt",
        fluid_eos::Number=5
    )        
    @assert eachindex(composition) == eachindex(elements) "`composition` and `elements` must match"
    build = joinpath(Perple_X_jll.PATH[], "build")# path to PerpleX build
    vertex = joinpath(Perple_X_jll.PATH[], "vertex")# path to PerpleX vertex

    #Configure working directory
    prefix = joinpath(scratchdir, "out$(index)/")
    system("rm -rf $prefix; mkdir -p $prefix")

    # Place required data files
    artifact_path = artifact"perplex-datafiles"
    system("cp $(joinpath(artifact_path, "perplex-datafiles/$dataset")) $prefix") 
    system("cp $(joinpath(artifact_path,"perplex-datafiles/perplex_option.dat")) $prefix")
    system("cp $(joinpath(artifact_path,"perplex-datafiles/solution_model.dat")) $prefix")

    # Edit data files to specify number of nodes at which to solve
    system("sed -e \"s/x_nodes .*|/x_nodes                   $xnodes $xnodes |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/y_nodes .*|/y_nodes                   $ynodes $ynodes |/\" -i.backup $(prefix)perplex_option.dat")

    # Specify whether we want volume or weight percentages
    system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

    # Create build batch file
    # Options based on Perplex v7.1.6
    fp = open(prefix*"build.bat", "w")

    # Name, components, and basic options. P-T conditions.
    elementstring = ""
    for i in eachindex(composition)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n2\nn\nn\nn\n$elementstring\nn\n2\n$(first(T))\n$(last(T))\n$(first(P))\n$(last(P))\ny\n") # v6.8.7

    # Whole-rock composition
    for i ∈ eachindex(composition)
        if !isnan(composition[i])
            write(fp,"$(composition[i]) ")
        end
    end

    # Solution models
    write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\n$fluid_eos\nPseudosection")
    close(fp)

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # build PerpleX problem definition
    system("export $ldpathadjust; cd $prefix; $build < build.bat > build.log")

    # Run PerpleX vertex calculations
    result = system("export $ldpathadjust; cd $prefix; echo $index | $vertex > vertex.log")
    return result
end
export perplex_configure_pseudosection

## -- Perplex interface: 2. queries

"""
```julia
perplex_query_point(scratchdir::String, indvar::Number; index::Integer=1)
```

Query perplex results at a single temperature on an isobar or single pressure
on a geotherm. Results are returned as a string.
"""
function perplex_query_point(scratchdir::String, indvar::Number; index::Integer=1)
    werami = joinpath(Perple_X_jll.PATH[], "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

    # Sanitize T inputs to avoid PerpleX escape sequence
    if indvar == 999
        indvar = 999.001
    end

    # Create werami batch file
    # Options based on Perplex v7.1.6
    fp = open(prefix*"werami.bat", "w")
    write(fp,"$index\n1\n$indvar\n999\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.txt")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Read results and return them if possible
    data = ""
    try
        # Read entire output file as a string
        fp = open("$(prefix)$(index)_1.txt", "r")
        data = read(fp, String)
        close(fp)
    catch
        # Return empty string if file doesn't exist
        @warn "$(prefix)$(index)_1.txt could not be parsed, perplex may not have run"

    end
    return data
end
"""
```julia
perplex_query_point(scratchdir::String, P::Number, T::Number; index::Integer=1)
```

Query perplex results at a single P,T point in a pseudosection.
Results are returned as a string.
"""
function perplex_query_point(scratchdir::String, P::Number, T::Number; index::Integer=1)
    werami = joinpath(Perple_X_jll.PATH[], "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

    # Sanitize T inputs to avoid PerpleX escape sequence
    if P == 99
        P = 99.001
    end
    if T == 99
        T = 99.001
    end

    # Create werami batch file
    # Options based on Perplex v6.7.2
    fp = open(prefix*"werami.bat", "w")
    write(fp,"$index\n1\n$T\n$P\n99\n99\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.txt")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Read results and return them if possible
    data = ""
    try
        # Read entire output file as a string
        fp = open("$(prefix)$(index)_1.txt", "r")
        data = read(fp, String)
        close(fp)
    catch
        # Return empty string if file doesn't exist
        @warn "$(prefix)$(index)_1.txt could not be parsed, perplex may not have run"

    end
    return data
end
export perplex_query_point

"""
```julia
perplex_query_seismic(scratchdir::String;
    \tindex::Integer=1, 
    \tdof::Integer=1, 
    \tinclude_fluid::String="n", 
    \tnpoints::Integer=0
    \tmanual_grid::Bool=npoints>0, 
    \timportas=:Dict,
)
```

Query perplex seismic results along a previously configured 1-d path (dof=1,
isobar or geotherm) or 2-d grid / pseudosection (dof=2).
Results are returned as a dictionary.
"""
function perplex_query_seismic(scratchdir::String;
        index::Integer=1, 
        dof::Integer=1, 
        include_fluid::String="n", 
        npoints::Integer=0,
        manual_grid::Bool=npoints>0, 
        importas=:Dict,
    )
    # Query a pre-defined path (isobar or geotherm)
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    if manual_grid
        # Edit perplex_option.dat to specify a manual grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   F |/\" -i.backup $(prefix)perplex_option.dat")
        if dof == 1
            write(fp,"$index\n3\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\nn\n$npoints\n0\n") # v7.1.8+ 1d path
            # write(fp,"$index\n3\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\n0\n") # v6.7.8 1d path
        elseif dof == 2
            # v7.1.8+ 2d grid
            write(fp,"$index\n2\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\nn\n$npoints\n0\n")
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
    else 
        # Edit perplex_option.dat to specify an automatic grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   T |/\" -i.backup $(prefix)perplex_option.dat")
        if dof == 1
            write(fp,"$index\n3\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\nn\n1\n0\n") # v7.1.6 1d path
        elseif dof==2 
            # v6.7.8 2d grid
            write(fp,"$index\n2\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\nn\n1\n0\n")
        else 
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
    end
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    data = nothing
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary
        data = elementify(data, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        data = importas==:Dict ? Dict() : ()
    end
    return data
end
"""
```julia
perplex_query_seismic(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1,
    \tnpoints::Integer=200, 
    \tinclude_fluid="n", 
    \timportas=:Dict,
)
```
```julia
perplex_query_seismic(scratchdir::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, 
    \tinclude_fluid="n", 
    \timportas=:Dict,
)
```
Query perplex seismic results along a specified P-T path using a pre-computed
pseudosection. Results are returned as a dictionary.
"""
function perplex_query_seismic(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
        index::Integer=1,
        npoints::Integer=200, 
        include_fluid="n", 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8 pseudosection
    write(fp,"$index\n3\nn\n$(first(T))\n$(first(P))\n$(last(T))\n$(last(P))\n$npoints\n2\nn
            $include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    data = nothing
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary
        data = elementify(data, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        data = importas==:Dict ? Dict() : ()
    end
    return data
end
function perplex_query_seismic(scratchdir::String, P::AbstractArray, T::AbstractArray;
        index::Integer=1, 
        include_fluid="n", 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Write TP data to file
    fp = open(prefix*"TP.tsv", "w")
    for i in eachindex(T,P)
        write(fp,"$(T[i])\t$(P[i])\n")
    end
    close(fp)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8 pseudosection
    write(fp,"$index\n4\n2\nTP.tsv\n1\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    data = nothing
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary
        data = elementify(data, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        data = importas==:Dict ? Dict() : ()
    end
    return data
end
export perplex_query_seismic


"""
```julia
perplex_query_phase(scratchdir::String, phase::String;
    \tindex::Integer=1, 
    \tdof::Integer=1,
    \tinclude_fluid="y", 
    \tclean_units::Bool=true, 
    \tnpoints::Integer=0,
    \tmanual_grid::Bool=npoints>0, 
    \timportas=:Dict,
)
```

Query all perplex-calculated properties for a specified phase (e.g. "Melt(G)")
along a previously configured 1-d path (dof=1, isobar, geotherm, or P–T path) or 2-d
grid / pseudosection (dof=2). Results are returned as a dictionary.
"""
function perplex_query_phase(scratchdir::String, phase::String;
        index::Integer=1, 
        dof::Integer=1,
        include_fluid="y", 
        clean_units::Bool=true, 
        npoints::Integer=0,
        manual_grid::Bool=npoints>0, 
        importas=:Dict,
    )
    # Query a pre-defined path (isobar or geotherm)
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")

    if manual_grid
        # Edit perplex_option.dat to specify a manual grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   F |/\" -i.backup $(prefix)perplex_option.dat")

        #v7.1.8+ (same for dof=1 & 2??)
        write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\nn\n$npoints\n5\n0\n") #5 is for the case of immiscible phases

    else 
        # Edit perplex_option.dat to specify an automatic grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   T |/\" -i.backup $(prefix)perplex_option.dat")

        #v7.1.8+ (same for dof=1 & 2??)
        write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\nn\n1\n5\n0\n") #5 is for the case of immiscible phases
    end

    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        elements = data[1,:]

        # Renormalize weight percentages
        t = contains.(elements,"wt%")
        total_weight = nansum(Float64.(data[2:end,t]),dim=2)
        hasdata = count(x->x>0, total_weight) 
        if (hasdata > 0) && !(50 < nansum(total_weight)/hasdata < 150)
            @warn "Perple_X may be reporting incorrect or unnormalized phase compositions"
        end
        data[2:end,t] .*= 100 ./ total_weight

        # Clean up element names
        if clean_units
            elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
            elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
        end

        # Convert to a dictionary
        result = elementify(data, elements, skipstart=1, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
"""
```julia
perplex_query_phase(scratchdir::String, phase::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1, 
    \tnpoints::Integer=200, 
    \tinclude_fluid="y", 
    \tclean_units::Bool=true, 
    \timportas=:Dict,
)
```
```julia
perplex_query_phase(scratchdir::String, phase::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, 
    \tinclude_fluid="y", 
    \tclean_units::Bool=true, 
    \timportas=:Dict,
)
```

Query all perplex-calculated properties for a specified phase (e.g. "Melt(G)")
along a specified P-T path using a pre-computed pseudosection. Results are
returned as a dictionary.
"""
function perplex_query_phase(scratchdir::String, phase::String, P::NTuple{2,Number}, T::NTuple{2,Number};
        index::Integer=1, 
        npoints::Integer=200, 
        include_fluid="y", 
        clean_units::Bool=true, 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    
    # If a named phase (e.g. feldspar) has multiple immiscible phases, average them (5)
    write(fp,"$index\n3\nn\n$(first(T))\n$(first(P))\n$(last(T))\n$(last(P))\n$npoints\n36\n2\n$phase\n$include_fluid\n5\n0\n") # v6.7.8, v7.1.9+ pseudosection
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        elements = data[1,:]

        # Renormalize weight percentages
        t = contains.(elements,"wt%")
        total_weight = nansum(Float64.(data[2:end,t]),dim=2)
        hasdata = count(x->x>0, total_weight) 
        if (hasdata > 0) && !(50 < nansum(total_weight)/hasdata < 150)
            @warn "Perple_X may be reporting incorrect or unnormalized phase compositions"
        end
        data[2:end,t] .*= 100 ./ total_weight

        # Clean up element names
        if clean_units
            elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
            elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
        end

        # Convert to a dictionary
        result = elementify(data, elements, skipstart=1, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
function perplex_query_phase(scratchdir::String, phase::String, P::AbstractArray, T::AbstractArray;
        index::Integer=1, 
        include_fluid="y", 
        clean_units::Bool=true, 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Write TP data to file
    fp = open(prefix*"TP.tsv", "w")
    for i in eachindex(T,P)
        write(fp,"$(T[i])\t$(P[i])\n")
    end
    close(fp)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")

    # If a named phase (e.g. feldspar) has multiple immiscible phases, average them (5).
    write(fp,"$index\n4\n2\nTP.tsv\n1\n36\n2\n$phase\n$include_fluid\n5\n0\n") # v6.7.8, 7.1.9+
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        elements = data[1,:]

        # Renormalize weight percentages
        t = contains.(elements,"wt%")
        total_weight = nansum(Float64.(data[2:end,t]),dim=2)
        hasdata = count(x->x>0, total_weight) 
        if (hasdata > 0) && !(50 < nansum(total_weight)/hasdata < 150)
            @warn "Perple_X may be reporting incorrect or unnormalized phase compositions"
        end
        data[2:end,t] .*= 100 ./ total_weight

        # Clean up element names
        if clean_units
            elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
            elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
        end

        # Convert to a dictionary
        result = elementify(data, elements, skipstart=1, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
export perplex_query_phase


"""
```julia
perplex_query_modes(scratchdir::String;
    \tindex::Integer=1, 
    \tdof::Integer=1, 
    \tinclude_fluid="y", 
    \tnpoints::Integer=0,
    \tmanual_grid::Bool=npoints>0,
    \timportas=:Dict,
)
```

Query modal mineralogy (mass proportions) along a previously configured 1-d
path (dof=1, isobar, geotherm, or P–T path) or 2-d grid / pseudosection (dof=2).
Results are returned as a dictionary.

Currently returns wt% 
"""
function perplex_query_modes(scratchdir::String;
        index::Integer=1, 
        dof::Integer=1, 
        include_fluid="y", 
        npoints::Integer=0,
        manual_grid::Bool=npoints>0,
        importas=:Dict,
    )
    # Query a pre-defined path (isobar or geotherm)
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    if manual_grid
         # Edit perplex_option.dat to specify a manual grid size
         system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   F |/\" -i.backup $(prefix)perplex_option.dat")
        if dof == 1 
            write(fp,"$index\n3\n38\n3\n$include_fluid\n37\n0\nn\n$npoints\n0\n") # v7.1.8+ 1d path
            # write(fp,"$index\n3\n38\n3\n$include_fluid\n37\n0\n0\n1\n0\n") # v7.1.6 1d path
            # write(fp,"$index\n3\n38\n3\nn\n37\n0\n0\n") # v6.7.8 1d path
        elseif dof == 2
            # v6.7.8 2d grid 
            #TODO: check this with a pseudosection example
            write(fp,"$index\n2\n25\nn\n$include_fluid\nn\n1\n0\n")
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
    else
         # Edit perplex_option.dat to specify an automatic grid size
         system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   T |/\" -i.backup $(prefix)perplex_option.dat")

        if dof == 1 
            write(fp,"$index\n3\n38\n3\n$include_fluid\n37\n0\n0\n1\n0\n") # v7.1.6 1d path
        elseif dof == 2
            # v6.7.8 2d grid
            write(fp,"$index\n2\n25\nn\n$include_fluid\nn\n1\n0\n")
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
    end
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.phm*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.phm")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.phm")
    # Replace "Missing data" with just "Missing" 
    file_content = read("$(prefix)$(index)_1.phm", String)
    modified_content = replace(file_content, "Missing data" => replace("Missing data", "Missing data" => "Missing"))
    write("$(prefix)$(index)_1.phm", modified_content)

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    
    if dof == 1 
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.phm", skipstart=8)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.phm could not be parsed, perplex may not have run"
        end
        # Convert to a dictionary.
        table = elementify(data, importas=importas)
        # Create results dictionary
        phase_names = unique(table["Name"])

        if haskey(table, "node#") #PT path
            nodes = unique(table["node#"])

            # Create result dictionary
            result = Dict{String, Vector{Float64}}(i => zeros(length(nodes)) for i in phase_names)
            result["P(bar)"] = zeros(length(nodes))

            # Loop through table
            for n in nodes
                # Index table 
                n_idx = table["node#"] .== n
                # Index phase name and weight(kg) 
                name = table["Name"][n_idx]
                kg = table["phase,kg"][n_idx]
                #  Calculate wt% and add to results dictionary
                for i in zip(name, kg)
                    result[i[1]][floor(Int64, n)] = (i[2]/nansum(kg)) * 100
                end
                result["P(bar)"][floor(Int64, n)] = table["P(bar)"][n_idx][1]
            end
            result["T(K)"] = unique(table["T(K)"])
            result["node"] = nodes

        else # isobar or geotherm
            pt_steps = unique([table["P(bar)"] table["T(K)"]], dims=1)
            p_steps, t_steps = pt_steps[:,1], pt_steps[:,2]
            result = Dict{String, Vector{Float64}}(i => zeros(size(pt_steps,1)) for i in phase_names)

            id = 1
            # Loop through table
            for i in axes(pt_steps,1)
                # Index table 
                idx = (table["P(bar)"] .== pt_steps[i,1]) .& (table["T(K)"] .== pt_steps[i,2])
                # Index phase name and weight(kg) 
                name = table["Name"][idx]
                kg = table["phase,kg"][idx]
                # Calculate wt% and add to results dictionary
                for i in zip(name, kg)
                    result[i[1]][id] = (i[2]/nansum(kg)) * 100
                end
                id+=1
            end
            result["P(bar)"] = pt_steps[:,1]
            result["T(K)"] = pt_steps[:,2]
        end
    else 
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary.
        # Perplex sometimes returns duplicates of a single solution model, sum them.
        result = elementify(data, sumduplicates=true, verbose=false, importas=importas)
    end
    return result
end
"""
```julia
perplex_query_modes(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1, 
    \tnpoints::Integer=200, 
    \tinclude_fluid="y", 
    \timportas=:Dict,
)
```
```julia
perplex_query_modes(scratchdir::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, 
    \tinclude_fluid="y", 
    \timportas=:Dict,
)
```

Query modal mineralogy (mass proportions) along a specified P-T path using a
pre-computed pseudosection. Results are returned as a dictionary.
"""
function perplex_query_modes(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
        index::Integer=1, 
        npoints::Integer=200, 
        include_fluid="y", 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8/7.1.6 pseudosection
    write(fp,"$index\n3\nn\n$(first(T))\n$(first(P))\n$(last(T))\n$(last(P))\n$npoints\n25\nn\n$include_fluid\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary
        result = elementify(data, sumduplicates=true, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
function perplex_query_modes(scratchdir::String, P::AbstractArray, T::AbstractArray;
        index::Integer=1, 
        include_fluid="y", 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Write TP data to file
    fp = open(prefix*"TP.tsv", "w")
    for i in eachindex(T,P)
        write(fp,"$(T[i])\t$(P[i])\n")
    end
    close(fp)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    write(fp,"$index\n4\n2\nTP.tsv\n1\n25\nn\n$include_fluid\n0\n")  # v6.7.8/v7.1.6 pseudosection
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary
        result = elementify(data, sumduplicates=true, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
export perplex_query_modes


"""
```julia
perplex_query_system(scratchdir::String;
    \tindex::Integer=1, 
    \tdof::Integer=1, 
    \tinclude_fluid="y", 
    \tclean_units::Bool=true,
    \tnpoints::Integer=0,
    \tmanual_grid::Bool=npoints>0, 
    \timportas=:Dict,
)
```

Query all perplex-calculated properties for the system (with or without fluid)
along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d
grid / pseudosection (dof=2). Results are returned as a dictionary.
Set include_fluid="n" to return solid+melt only.
"""
function perplex_query_system(scratchdir::String;
        index::Integer=1, 
        dof::Integer=1, 
        include_fluid="y", 
        clean_units::Bool=true,
        npoints::Integer=0,
        manual_grid::Bool=npoints>0, 
        importas=:Dict,
    )
    # Query a pre-defined path (isobar or geotherm)
    werami = joinpath(Perple_X_jll.PATH[], "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")

    if manual_grid
        # Edit perplex_option.dat to specify a manual grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   F |/\" -i.backup $(prefix)perplex_option.dat")

        # v 7.1.8+ (same for both dof=1 and dof=2)
        write(fp,"$index\n3\n36\n1\n$include_fluid\nn\n$npoints\n0\n")
    else 
        # Edit perplex_option.dat to specify an automatic grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   T |/\" -i.backup $(prefix)perplex_option.dat")

        # v 7.1.6+ (same for both dof=1 and dof=2)
        write(fp,"$index\n3\n36\n1\n$include_fluid\nn\n1\n0\n") #uses coarsest grid option
    end
    
    close(fp)

    # Make sure there isn't already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        elements = data[1,:]

        # Renormalize weight percentages
        t = contains.(elements,"wt%")
        total_weight = nansum(Float64.(data[2:end,t]),dim=2)
        hasdata = count(x->x>0, total_weight) 
        if (hasdata > 0) && !(50 < nansum(total_weight)/hasdata < 150)
            @warn "Perple_X may be reporting incorrect or unnormalized system compositions"
        end
        data[2:end,t] .*= 100 ./ total_weight

        # Clean up element names
        if clean_units
            elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
            elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
        end

        # Convert to a dictionary
        result = elementify(data, elements, skipstart=1, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
"""
```julia
perplex_query_system(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1, 
    \tnpoints::Integer=200, 
    \tinclude_fluid="y", 
    \tclean_units::Bool=true, 
    \timportas=:Dict,
)
```
```julia
perplex_query_system(scratchdir::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, 
    \tinclude_fluid="y", 
    \tclean_units::Bool=true, 
    \timportas=:Dict,
)
```

Query all perplex-calculated properties for the system (with or without fluid)
along a specified P-T path using a pre-computed pseudosection. Results are
returned as a dictionary. Set include_fluid="n" to return solid+melt only.
"""
function perplex_query_system(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
        index::Integer=1, 
        npoints::Integer=200, 
        include_fluid="y", 
        clean_units::Bool=true, 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8+ pseudosection
    write(fp,"$index\n3\nn\n$(first(T))\n$(first(P))\n$(last(T))\n$(last(P))\n$npoints\n36\n1\n$include_fluid\n0\n")
    close(fp)

    # Make sure there isn't already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        elements = data[1,:]

        # Renormalize weight percentages
        t = contains.(elements,"wt%")
        total_weight = nansum(Float64.(data[2:end,t]),dim=2)
        hasdata = count(x->x>0, total_weight) 
        if (hasdata > 0) && !(50 < nansum(total_weight)/hasdata < 150)
            @warn "Perple_X may be reporting incorrect or unnormalized system compositions"
        end
        data[2:end,t] .*= 100 ./ total_weight

        # Clean up element names
        if clean_units
            elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
            elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
        end

        # Convert to a dictionary
        result = elementify(data, elements, skipstart=1, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
function perplex_query_system(scratchdir::String, P::AbstractArray, T::AbstractArray;
        index::Integer=1, 
        include_fluid="y", 
        clean_units::Bool=true, 
        importas=:Dict,
    )
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files
    include_fluid = parse_bool_to_yn(include_fluid)

    # Write TP data to file
    fp = open(prefix*"TP.tsv", "w")
    for i in eachindex(T,P)
        write(fp,"$(T[i])\t$(P[i])\n")
    end
    close(fp)

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8/7.1.6 pseudosection
    write(fp,"$index\n4\n2\nTP.tsv\n1\n36\n1\n$include_fluid\n0\n")
    close(fp)

    # Make sure there isn't already an output
    system("rm -f $(prefix)$(index)_1.tab*")

    # Add any needed library paths
    ldpathadjust = "DYLD_LIBRARY_PATH=$(first(Perple_X_jll.LIBPATH_list)):\$DYLD_LIBRARY_PATH"

    # Extract Perplex results with werami
    system("export $ldpathadjust; cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        elements = data[1,:]

        # Renormalize weight percentages
        t = contains.(elements,"wt%")
        total_weight = nansum(Float64.(data[2:end,t]),dim=2)
        hasdata = count(x->x>0, total_weight) 
        if (hasdata > 0) && !(50 < nansum(total_weight)/hasdata < 150)
            @warn "Perple_X may be reporting incorrect or unnormalized system compositions"
        end
        data[2:end,t] .*= 100 ./ total_weight

        # Clean up element names
        if clean_units
            elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
            elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
        end

        # Convert to a dictionary
        result = elementify(data, elements, skipstart=1, verbose=false, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
export perplex_query_system


## --- PerpleX interface: 3. PerplexTrace

    # Split binary and ternary feldspar solution models into An-Ab-Or endmembers
    function enumerate_feldspar_modes!(modes::Dict, fsp::Dict, fsp_model::AbstractString; mode_basis="wt")
        An_Ca = (238.12507+40.0784) / (15.999+40.0784)
        Ab_Na = (239.22853+22.98977*2) / (15.999+22.98977*2)
        Or_K  = (239.22853+39.09831*2) / (15.999+39.09831*2)
        CaO = haskey(fsp, "CaO") ? fsp["CaO"] : fsp["CAO"]
        Na2O = haskey(fsp, "Na2O") ? fsp["Na2O"] : fsp["NA2O"]
        K2O = fsp["K2O"]
        if mode_basis === "wt"
            AnAbOr = [An_Ca*CaO  Ab_Na*Na2O  Or_K*K2O] # Modes by weight
        elseif mode_basis === "vol"
            AnAbOr = [An_Ca*CaO/2.73  Ab_Na*Na2O/2.63  Or_K*K2O/2.59] # Modes by weight
        else
            @error "Invalid mode basis $mode_basis. Try \"wt\" or \"vol\" instead."
        end
        renormalize!(AnAbOr, dim=2)
        modes[fsp_model*"_anorthite"] = AnAbOr[:,1] .*  modes[fsp_model]
        modes[fsp_model*"_albite"] = AnAbOr[:,2] .*  modes[fsp_model]
        modes[fsp_model*"_orthoclase"] = AnAbOr[:,3] .*  modes[fsp_model]
        return modes
    end

    # Calculate GERM solid/melt Kds for the bulk solid assemblage
    function get_germ_bulk_kds(modes::Dict, trace_elements, si_index::Vector{<:Number})

        d = Dict{String,Vector{Float64}}()
        for e in trace_elements
            # Calculate bulk partition coeff.
            d[e] = zeros(size(modes["all_solids"]))
            for m in germ_kd["minerals"]
                for k in filter(x -> (germ_perplex_name_matches(m, x) | containsi(x, m)), keys(modes))
                    # Note that minerals that we don't have data for end up being
                    # treated like all elements are incompatible in them.
                    # Note, geometric mean = log average
                    nanadd!(d[e], modes[k]./modes["all_solids"] .* (10.0.^germ_kd[m][e][si_index]))
                end
            end
        end
        return d
    end

    # Ensure that more REEs are not partitioned into monazite than are structurally possible
    function germ_monazite_kd_corr(calculated::Dict, lree, si_index::Vector{<:Number})
        partitioned_REEt_in_monazite = zeros(length(si_index))
        for e in lree
            if haskey(calculated, e)
                nanadd!(partitioned_REEt_in_monazite, (calculated[e].* (10.0.^germ_kd["Monazite"][e][si_index]))./molarmass[e])
            end
        end
        return min.((596020.9/140.1161)./partitioned_REEt_in_monazite, 1)
    end

    # Update GERM kds in the presence of accessory phases
    function update_kds!(d::Dict, modes::Dict, trace_elements, si_index::Vector{<:Number}, monazite_kd_corr::Vector{<:Number})
        for e in trace_elements
            fill!(d[e], 0.)
            for m in germ_kd["minerals"]
                if containsi(m,"zircon")
                    zircon_kd_c = claiborne_zircon_kd.(e, modes["T(C)"])
                    zircon_kd_g = 10.0.^germ_kd["Zircon"][e][si_index]
                    zircon_kd = sqrt.(zircon_kd_c .* zircon_kd_g)
                    nanadd!(d[e], modes["zircon"]./modes["all_solids"] .* zircon_kd)
                elseif containsi(m,"monazite")
                    nanadd!(d[e], modes["monazite"]./modes["all_solids"] .* (10.0.^germ_kd[m][e][si_index])) .* monazite_kd_corr
                else
                    for k in filter(x -> (germ_perplex_name_matches(m, x) | containsi(x, m)), keys(modes))
                        # Note that minerals that we don't have data for end up being
                        # treated like all elements are incompatible in them.
                        # Note, geometric mean = log average
                        nanadd!(d[e], modes[k]./modes["all_solids"] .* (10.0.^germ_kd[m][e][si_index]))
                    end
                end
            end
        end
    end

    function perplextrace_query(scratchdir, composition::LinearTraceComposition, args...; 
            apatite = :Harrison,
            zircon = :Boehnke,
            sphene = :Ayers,
            monazite = :Montel,
            export_bulk_kds::Bool = false,
            export_mineral_kds::Bool = false,
            export_mineral_compositions::Bool = false,
            require_phase_for_export = "",
            kwargs...,
        )
        major_elements = collect(String.(filter(e->!isnan(composition[e]), majorelements(composition))))
        trace_elements = collect(String.(traceelements(composition)))

        bulk = perplex_query_system(scratchdir, args...; kwargs...)
        modes = perplex_query_modes(scratchdir, args...; kwargs...)

        @assert count(perplex_phase_is_melt, keys(modes)) > 0 "No melt model found"
        @assert count(perplex_phase_is_melt, keys(modes)) < 2 "Multiple melt models found"
        melt_model = only(filter(perplex_phase_is_melt, keys(modes)))
        melt = perplex_query_phase(scratchdir, melt_model, args...; kwargs...)
        dataset_uppercase = haskey(melt, "SIO2")

        # Add extra composite modes
        colshape = size(modes["T(K)"])
        modes["all_solids"] = fill(NaN, colshape)
        for m in filter(perplex_phase_is_solid,  keys(modes))
            nanadd!(modes["all_solids"], modes[m])
        end
        modes["all_fluids"] = fill(NaN, colshape)
        for m in filter(perplex_phase_is_fluid,  keys(modes))
            nanadd!(modes["all_fluids"], modes[m])
        end
        modes["all_melts"] = fill(NaN, colshape)
        for m in filter(perplex_phase_is_melt,  keys(modes))
            nanadd!(modes["all_melts"], modes[m])
        end

        for fsp_model in filter(x->contains(perplex_common_name(x), "feldspar"), keys(modes)) 
            fsp = perplex_query_phase(scratchdir, fsp_model, args...; kwargs...)
            enumerate_feldspar_modes!(modes, fsp, fsp_model)
        end

        # Add temperatures in C as well as K
        bulk["T(C)"] = bulk["T(K)"] .- 273.15
        modes["T(C)"] = modes["T(K)"] .- 273.15
        melt["T(C)"] = melt["T(K)"] .- 273.15
        meltrows = length(melt["T(K)"])

        # Create dictionary to hold solid composition and fill it using what we know from system and melt
        solid = Dict{String,Vector{Float64}}()
        solid["wt_pct"] = modes["all_solids"]
        solid_elements = setdiff(major_elements, ("CO2", "H2O", "O2",))
        for e in solid_elements
            solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
        end
        renormalize!(solid, solid_elements, total=100)

        ## --- First-pass calculation of trace elements

        # Since the mineral/melt partition coefficients we're using are averaged from GERM
        # as a function of silica (from 40-80), first figure out which coefficients we want
        si_index = round.(Int, melt[dataset_uppercase ? "SIO2" : "SiO2"] .|> x-> isnan(x) ? 0 : x) .- 39
        si_index[si_index.<1] .= 1
        si_index[si_index.>40] .= 40

        # Melt fraction
        F = modes[melt_model] ./ (modes[melt_model] .+ modes["all_solids"])

        # Use GERM partition coeffs to start
        d = get_germ_bulk_kds(modes, trace_elements, si_index)

        # Calculate trace elements as a function of melt fraction
        calculated = Dict{String, Union{Vector{Float64}, Vector{String}}}()
        calculated["elements"] = trace_elements
        for e in calculated["elements"]
            calculated[e] = composition[e] ./ (d[e].*(1.0.-F) + F)
        end

        # Prepare to add accessory minerals (zircon, apatite, sphene)
        # Note that all Fe is _already_ as FeO, since perplex melt has no Fe2O3
        # (redox balance is with O2 instead)

        if apatite === :Harrison
            # Add apatite using the apatite saturation temperature of Harrison and Watson
            e = "SiO2", "Al2O3", "CaO", "Na2O", "K2O", "T(C)"
            dataset_uppercase && (e = uppercase.(e))
            PSat = StatGeochem.Harrison_tapatiteP.((e .|> x -> melt[x])...)

            # Calculate mass of P in apatite. Note that you have to account for
            # not only the P in the melt, but also the P that would
            # reequilibrate into the melt from the solids at P=PSat. Note that mass
            # of melt actually cancels when you write out the equation for this.
            P_in_apatite = max.(calculated["P"] - PSat, 0) .* composition["P"]./calculated["P"]
            modes["apatite"] = P_in_apatite * 4.36008264/10_000
            nanadd!(modes["all_solids"], modes["apatite"])
            map!(x -> (x>0 ? x : NaN), modes["apatite"], modes["apatite"])  # NaN-out zero masses, for printing

            # Record ApSat in melt composition dictionary
            calculated["ApSat_P"] = PSat
            calculated["elements"] = ["ApSat_P"; calculated["elements"]]
        end

        if zircon === :Boehnke
            # Add zircon using the zircon saturation temperature of Boehnke et al. 2013
            # (Earth and Planetary Science Letters 453, 267–275. 10.1016/j.epsl.2016.07.014)
            e = "SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5"
            dataset_uppercase && (e = uppercase.(e))
            M = StatGeochem.tzircM((e .|> x -> haskey(melt,x) ? melt[x] : zeros(meltrows))...)
            ZrSat = StatGeochem.tzircZr((e .|> x -> haskey(melt,x) ? melt[x] : zeros(meltrows))..., melt["T(C)"])
            # Calculate mass of zirconium in zircon. Note that you have to account for
            # not only the zirconium in the melt, but also the zirconium that would
            # reequilibrate into the melt from the solids at Zr=ZrSat. Note that mass
            # of melt actually cancels when you write out the equation for this.
            # For now we'll also treat Hf as equivalent to Zr (b/c zircon-hafnon ssn)
            ZrHf = calculated["Zr"]+calculated["Hf"]
            Zr_in_zircon = max.(ZrHf - ZrSat, 0) .* (composition["Zr"]+composition["Hf"])./ZrHf
            # Convert from zirconium mass to zircon mass and from ppm to wt. %
            modes["zircon"] = Zr_in_zircon * 2.009/10_000
            nanadd!(modes["all_solids"], modes["zircon"])
            map!(x -> (x>0 ? x : NaN), modes["zircon"], modes["zircon"]) # NaN-out zero masses, for printing

            # Record M and ZrSat in melt composition dictionary
            calculated["M"] = M
            calculated["ZrnSat_Zr"] = ZrSat
            calculated["elements"] = ["M"; "ZrnSat_Zr"; calculated["elements"]]
        end

        if sphene === :Ayers
            # Add sphene using the saturation equation of Ayers et al. 2018.
            # (GSA abstract).
            e = "SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "T(C)"
            dataset_uppercase && (e = uppercase.(e))
            TiO2Sat = StatGeochem.Ayers_tspheneTiO2.((e .|> x -> haskey(melt,x) ? melt[x] : zeros(meltrows))...)
            modes["sphene"] = modes[melt_model]/100 .* max.(melt["TiO2"] - TiO2Sat, 0)*2.4545
            nanadd!(modes["all_solids"], modes["sphene"])
            nanadd!(modes[melt_model], -modes["sphene"])
            map!(x -> (x>0 ? x : NaN), modes["sphene"], modes["sphene"]) # NaN-out zero masses, for printing
            melt["TiO2"] = min.(melt["TiO2"], TiO2Sat)
        end

        if monazite === :Montel
            # Add monazite using the saturation equation of Montel et al. 1993.
            # (Chemical Geology 110, 127–146. 10.1016/0009-2541(93)90250-M)
            e = "SiO2", "TiO2", "Al2O3", "FeOT", "MgO", "CaO", "Na2O", "K2O", "Li2O", "H2O", "T(C)"
            dataset_uppercase && (e = uppercase.(e))
            REEtSat = StatGeochem.Montel_tmonaziteREE.((e .|> x -> haskey(melt,x) ? melt[x] : zeros(meltrows))...)
            lree = "La", "Ce", "Pr", "Nd", "Sm", "Gd"
            melt_lree = lree .|> x -> haskey(calculated, x) ? calculated[x] : zeros(meltrows)
            melt_REEt = StatGeochem.LREEt.(melt_lree...) # In PPM/mol.wt.
            bulk_REEt = StatGeochem.LREEt((lree .|> x -> haskey(composition, x) ? composition[x] : 0.)...)
            # REEt_in_monazite = modes[melt_model]/100 .* max.(melt_REEt - REEtSat, 0) # naive
            REEt_in_monazite = max.(melt_REEt - REEtSat, 0.) .* bulk_REEt./melt_REEt # accounting for back-equilibration
            modes["monazite"] = REEt_in_monazite .* (StatGeochem.LREEmolwt.(melt_lree...) .+ 94.969762)/10_000
            nanadd!(modes["all_solids"], modes["monazite"])
            map!(x -> (x>0 ? x : NaN), modes["monazite"], modes["monazite"]) # NaN-out zero masses, for printing
        end

        # In contrast to other phases, monazite partition coefficients may
        # actually be limited by mass balance!
        monazite_kd_corr = germ_monazite_kd_corr(calculated, lree, si_index)

        ## Recalculate bulk partition coeff. in presence of accessory phases
        update_kds!(d, modes, trace_elements, si_index, monazite_kd_corr)

        # Recalculate trace elements in melt as a function of melt fraction (equilibrium)
        for e in trace_elements
            calculated[e] .= composition[e] ./ (d[e].*(1.0.-F) + F)
        end

        #Re-Recalculate trace elements with monazite
        if haskey(modes, "monazite") && any(x->x>0, modes["monazite"])
            # Recalculate monazite mass-balance Kd adjustment
            all_ee3 = "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Y", "Th"
            monazite_kd_corr .= germ_monazite_kd_corr(calculated, all_ee3, si_index)

            # Re-recalculate bulk Kd
            update_kds!(d, modes, trace_elements, si_index, monazite_kd_corr)

            # Re-Recalculate trace elements in melt as a function of melt fraction (equilibrium)
            for e in trace_elements
                calculated[e] .= composition[e] ./ (d[e].*(1.0.-F) .+ F)
            end
        end

        # Adjust melt Zr and P contents to reflect zircon and apatite saturation,
        # respectively
        if haskey(modes, "zircon") && any(x->x>0, modes["zircon"])
            calculated["Zr"] .= min.(calculated["Zr"], ZrSat)
        end
        if haskey(modes, "apatite") && any(x->x>0, modes["apatite"])
            calculated["P"] .= min.(calculated["P"], PSat)
        end


        ## --- # # # # # # # # # # #  Export Results  # # # # # # # # # # # # #

        # Which rows to plot and export (only those with melt, for example?)
        exportrows = isempty(require_phase_for_export) ? trues(colshape) : modes[require_phase_for_export] .> 0

        ## Collect and print main output
        result = Dict{String, Union{Vector{String}, Vector{Float64}}}()
        result["elements"] = String[]

        # Modes
        for e in ("P(bar)", "T(C)", "all_fluids", "all_melts", "all_solids", melt_model)
            result[e] = modes[e]
            result["elements"] = [result["elements"]; e]
        end
        minerals = collect(setdiff(keys(modes), ("elements", "T(C)", "T(K)", "P(bar)", "all_fluids", "all_melts", "all_solids", melt_model)))
        t = sortperm(lowercase.(minerals))
        for e in minerals[t]
            # Don't print empty columns
            if any(x->x>0, modes[e][exportrows])
                result[e] = modes[e]
                result["elements"] = [result["elements"]; e]
            end
        end

        # Liquid composition
        for e in melt["elements"]
            # Don't print free energies, or empty columns
            if any(x->x>0, melt[e][exportrows])
                result["Melt_$(e)"] = melt[e]
                result["elements"] = [result["elements"]; "Melt_$(e)"]
            end
        end
        # Trace elements
        for e in calculated["elements"]
            result["Melt_$(e)"] = calculated[e]
            result["elements"] = [result["elements"]; "Melt_$(e)"]
        end

        # [optionally] Bulk germ_kds
        if export_bulk_kds
            for e in trace_elements
                result["D_$(e)"] = d[e]
                result["elements"] = [result["elements"]; "D_$(e)"]
            end
        end

        # [optionally] Individual mineral kds
        if export_mineral_kds
            for m in setdiff(germ_kd["minerals"], ("Zircon",)) # Everything but zircon
                for k in filter(x -> (germ_perplex_name_matches(m, x) | containsi(x, m)), keys(modes))
                    if any(x->x>0, modes[k][exportrows])
                        for e in trace_elements
                            if containsi(m, "monazite")
                                result["D_$(k)_$(e)"] = (10.0 .^ germ_kd[m][e][si_index]) .* monazite_kd_corr
                            else
                                result["D_$(k)_$(e)"] = 10.0 .^ germ_kd[m][e][si_index]
                            end
                            result["elements"] = [result["elements"]; "D_$(k)_$(e)"]
                        end
                    end
                end
            end
            # Zircon kds
            if haskey(modes, "zircon") && any(x->x>0, modes["zircon"][exportrows])
                for e in trace_elements
                    zircon_kd_c = claiborne_zircon_kd.(e, modes["T(C)"])
                    zircon_kd_g = 10.0.^germ_kd["Zircon"][e][si_index]
                    zircon_kd = sqrt.(zircon_kd_c .* zircon_kd_g)
                    result["D_Zircon_$(e)"] = zircon_kd
                    result["elements"] = [result["elements"]; "D_Zircon_$(e)"]
                end
            end
        end

        # [optionally] Individual mineral trace element compositions
        if export_mineral_compositions
            for m in setdiff(germ_kd["minerals"], ("Zircon",)) # Everything but zircon
                for k in filter(x -> (germ_perplex_name_matches(m, x) | containsi(x, m)), keys(modes))
                    if any(x->x>0, modes[k][exportrows])
                        for e in trace_elements
                            if containsi(m, "monazite")
                                result["$(k)_$(e)"] = (calculated[e] .* (10 .^ germ_kd[m][e][si_index]) .* monazite_kd_corr) .+ NaN .* .!(modes[k] .> 0)
                            else
                                result["$(k)_$(e)"] = (calculated[e] .* 10 .^ germ_kd[m][e][si_index]) .+ NaN .* .!(modes[k] .> 0)
                            end
                            result["elements"] = [result["elements"]; "$(k)_$(e)"]
                        end
                    end
                end
            end
            # Zircon composition
            if haskey(modes, "zircon") && any(x->x>0, modes["zircon"][exportrows])
                zircon_Zr = 496000.0*(modes["zircon"] .> 0)
                for e in trace_elements
                    zircon_kd_c = claiborne_zircon_kd.(e, modes["T(C)"])
                    zircon_kd_g = 10.0.^germ_kd["Zircon"][e][si_index]
                    zircon_kd = sqrt.(zircon_kd_c .* zircon_kd_g)
                    result["Zircon_$(e)"] = calculated[e] .* zircon_kd .+ (NaN .* .!(modes["zircon"] .> 0))
                    result["elements"] = [result["elements"]; "Zircon_$(e)"]
                    nanadd!(zircon_Zr, -result["Zircon_$(e)"])
                end
                # Add Ti in zircon
                aSiO2 = 1.0
                meltcomp = melt["SiO2"], melt["TiO2"], melt["Al2O3"], melt["FeO"], melt["MgO"], melt["CaO"], melt["Na2O"], melt["K2O"], calculated["P"].*70.9723/30.974/10000
                aTiO2 = min.(melt["TiO2"] ./ StatGeochem.Hayden_trutileTiO2.(meltcomp..., melt["T(C)"]), 1.0) # Titanium activity as a fraction of the TiO2 neeeded for rutile saturation
                result["Zircon_Ti"] = StatGeochem.Ferry_Ti_in_zircon.(melt["T(C)"], aSiO2, aTiO2)
                result["Zircon_Ti"] .+= (NaN .* .!(modes["zircon"] .> 0)) # NaN out if no zircon
                nanadd!(zircon_Zr, -result["Zircon_Ti"])

                # Zircon Zr concentration in ppm, calculated as 496000 less other trace elements
                result["Zircon_Zr"] = zircon_Zr .+ (NaN .* (zircon_Zr .<= 0))
            end
        end

        return result
    end
    export perplextrace_query

## --- Perplex name/string-related utilities

    # Parse falsey/truey intput into "y" or "n"
    function parse_bool_to_yn(x)
        if x=="y" || x =="yes" || (x isa Number && x > 0)
            return "y"
        else
            return "n"
        end
    end

    # Attempt to determine whether or not a given perplex dataset expects uppercase input
    function perplex_dataset_uppercase(dataset::AbstractString)
        if contains(dataset, "hp6") || contains(dataset, "hpha6")  || contains(dataset, "hpAQ") || contains(dataset, "DEW") || contains(dataset, "HKF")
            false
        else
            true
        end
    end

    # Translate between perplex names and germ names
    function germ_perplex_name_matches(germ_name, perplex_name)
        # Feldspar
        if germ_name == "Albite"
            any(perplex_name .== ["ab", "abh"]) ||
            any(contains.(perplex_name, ["albite",]))
        elseif germ_name == "Anorthite"
            perplex_name == "an" ||
            any(contains.(perplex_name, ["anorthite",]))
        elseif germ_name == "Orthoclase"
            any(perplex_name .== ["mic", "Kf", "San", "San(TH)"]) ||
            any(contains.(perplex_name, ["orthoclase",]))
        # Amphibole
        elseif germ_name == "Amphibole"
            any(lowercase(perplex_name) .== ["gl", "fgl", "rieb", "anth", "fanth", "cumm", "grun", "tr", "ftr", "ged", "parg", "ts"]) ||
            any(contains.(perplex_name, ["Amph", "GlTrTs", "Act(", "Anth"]))
        # Mica
        elseif germ_name == "Biotite"
            any(perplex_name .== ["ann"]) ||
            any(contains.(perplex_name, ["Bi(", "Bio("]))
        elseif germ_name == "Phlogopite"
            any(lowercase(perplex_name) .== ["naph", "phl", "fphl"])
        elseif germ_name == "Muscovite"
            any(lowercase(perplex_name) .== ["pheng(hp)", "mapa", "mu",]) ||
            any(contains.(perplex_name, ["Mica", "KN-Phen"]))
        # Pyroxene
        elseif germ_name == "Clinopyroxene"
            any(lowercase(perplex_name) .== ["di", "hed", "acm", "jd", "lcendi", "lcfshd",]) ||
            any(contains.(perplex_name, ["Augite", "Cpx", "Omph"])) ||
            any(lowercase(perplex_name) .== ["mont",]) # Assume monticellite has similar partition coefficients to cpx
        elseif germ_name == "Orthopyroxene"
            any(lowercase(perplex_name) .== ["en", "fs"]) ||
            contains(perplex_name, "Opx")
        # Cordierite
        elseif germ_name == "Cordierite"
            any(lowercase(perplex_name) .== ["crd", "fcrd", "hcrd", "mncrd"]) ||
            contains(perplex_name, "Crd")
        # Garnet
        elseif germ_name == "Garnet"
            any(lowercase(perplex_name) .== ["py", "spss", "alm", "andr", "gr"]) ||
            any(contains.(perplex_name, ["Grt", "Gt(", "Maj", "GrPyAlSp"]))
        elseif germ_name == "Zoisite"
            any(lowercase(perplex_name) .== ["zo", "cz", "ep", "fep"]) ||
            any(contains.(perplex_name, ["Ep(",]))
        # Oxides
        elseif germ_name == "Ilmenite"
            perplex_name == "ilm" || 
            any(contains.(perplex_name, ["Ilm", "IlHm", "IlGk"]))
        elseif germ_name == "Magnetite"
            perplex_name == "mt"
        elseif germ_name == "Rutile"
            perplex_name == "ru"
        # Feldspathoids
        elseif germ_name == "Leucite"
            perplex_name == "lc"
        elseif germ_name == "Nepheline"
            perplex_name == "ne" || contains(perplex_name, "Neph")
        # Olivine
        elseif germ_name == "Olivine"
            any(lowercase(perplex_name) .== ["fo", "fa"]) ||
            any(contains.(perplex_name, ["O(", "Ol("]))
        # Spinel
        elseif germ_name == "Spinel"
            any(lowercase(perplex_name) .== ["sp", "usp", "gahcsp"]) ||
            contains(perplex_name, "Sp(")
        # Accessories
        elseif germ_name == "Sphene"
            perplex_name == "sph"
        elseif germ_name == "Zircon"
            perplex_name == "zrc"
        elseif germ_name == "Baddeleyite"
            perplex_name == "bdy"
        else
            false
        end
    end
    export germ_perplex_name_matches

    function perplex_phase_is_fluid(phase_name)
        any(phase_name .== ["F", "WADDAH", "H2O"]) ||
        any(contains.(phase_name, ["Aq_", "F(", "Fluid"]))
    end
    export perplex_phase_is_fluid

    function perplex_phase_is_melt(phase_name)
        any(phase_name .== ["h2oL", "abL", "anL", "diL", "enL", "faL", "kspL", "qL", "silL"]) ||
        any(contains.(phase_name, ["liq", "melt", "LIQ", "MELTS"]))
    end
    export perplex_phase_is_melt

    function perplex_phase_is_solid(phase_name)
        !perplex_phase_is_fluid(phase_name) && !perplex_phase_is_melt(phase_name) &&
        !any(contains.(phase_name, ["P(", "T(", "Pressure", "Temperature", "elements", "minerals", "CO2", "Missing", "system", "O2"]))
    end
    export perplex_phase_is_solid

    function perplex_expand_name(name)
        abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm","oilm","lime", "mft", "mt", "mang", "bunsn", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi",)
        full_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote(ordered)", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite-bredigite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermaks pyroxene", "Diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "mg-tschermak", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite(Na-free)", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "fe-carpholite", "fe-sapphirine(793)", "mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite(14Ang)", "clinochlore(ordered)", "daphnite", "Fe-sudoite", "Mn-chlorite", "Sudoite", "antigorite", "chrysotile", "Fe-talc", "Kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)","lime", "magnesioferrite", "magnetite", "manganosite", "nickel oxide", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "Forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica",)
        t = name .== abbreviations
        if any(t)
            full_names[findfirst(t)]
        else
            name
        end
    end
    export perplex_expand_name

    function perplex_abbreviate_name(name)
        abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm", "oilm", "lime", "mft", "mt", "mang", "bunsn", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi",)
        full_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote(ordered)", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite-bredigite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermaks pyroxene", "Diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "mg-tschermak", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite(Na-free)", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "fe-carpholite", "fe-sapphirine(793)", "mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite(14Ang)", "clinochlore(ordered)", "daphnite", "Fe-sudoite", "Mn-chlorite", "Sudoite", "antigorite", "chrysotile", "Fe-talc", "Kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)", "lime", "magnesioferrite", "magnetite", "manganosite", "nickel oxide", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "Forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica",)
        t = name .== full_names
        if any(t)
            abbreviations[findfirst(t)]
        else
            name
        end
    end
    export perplex_abbreviate_name

    function perplex_common_name(name)
        abbreviations = ("cg", "cgh", "kls", "lrn", "esk", "syv", "hlt", "pyr", "tro", "mil", "cv", "bn", "ccp", "any", "mal", "azr", "frw", "mrw", "fwd", "mwd", "cpv", "apv", "npv", "mpv", "fpv", "ppv", "fper", "maj", "wu", "nagt", "ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm", "oilm", "lime", "mft", "mt", "mang", "bunsn", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi", "Aqfl(HGP)", "Cpx(HGP)", "Augite(G)", "Cpx(JH)", "Cpx(l)", "Cpx(h)", "Cpx(stx)", "Cpx(stx7)", "Omph(HP)", "Cpx(HP)", "Cpx(m)", "Cpx(stx8)", "Cps(HGP)", "Omph(GHP)", "cAmph(G)", "Cumm", "Gl", "Tr", "GlTrTsPg", "Amph(DHP)", "Amph(DPW)", "Ca-Amph(D)", "Na-Amph(D)", "Act(M)", "GlTrTsMr", "cAmph(DP)", "melt(HGPH)", "melt(G)", "melt(W)", "melt(HP)", "melt(HGP)", "pMELTS(G)", "mMELTS(G)", "LIQ(NK)", "LIQ(EF)", "Chl(W)", "Chl(HP)", "Chl(LWV)", "O(HGP)","O(JH)", "O(SG)", "O(HP)", "O(HPK)", "O(stx)", "O(stx7)", "Ol(m)", "O(stx8)", "Sp(HGP)", "Sp(JH)", "GaHcSp", "Sp(JR)", "Sp(GS)", "Sp(HP)", "Sp(stx)", "CrSp", "Sp(stx7)", "Sp(WPC)", "Sp(stx8)", "Pl(JH)", "Pl(h)", "Pl(stx8)", "Kf", "San", "San(TH)", "Gt(HGP)", "Grt(JH)", "Gt(W)", "CrGt", "Gt(MPF)", "Gt(B)", "Gt(GCT)", "Gt(HP)", "Gt(EWHP)", "Gt(WPH)", "Gt(stx)", "Gt(stx8)", "Gt(WPPH)", "ZrGt(KP)", "Maj", "Opx(HGP)", "Opx(JH)", "Opx(W)", "Opx(HP)", "CrOpx(HP)", "Opx(stx)", "Opx(stx8)", "Mica(W)", "Pheng(HP)", "MaPa", "Mica(CF)", "Mica(CHA1)", "Mica(CHA)", "Mica+(CHA)", "Mica(M)", "Mica(SGH)", "Ctd(W)", "Ctd(HP)", "Ctd(SGH)", "St(W)", "St(HP)", "Bi(HGP)", "Bi(W)", "Bio(TCC)", "Bio(WPH)", "Bio(HP)", "Crd(W)", "hCrd", "Sa(WP)", "Sapp(HP)", "Sapp(KWP)", "Sapp(TP)", "Osm(HP)", "F", "F(salt)", "COH-Fluid", "Aq_solven0", "WADDAH", "T", "Scap", "Carp", "Carp(M)", "Carp(SGH)", "Sud(Livi)", "Sud", "Sud(M)", "Anth", "o-Amph", "oAmph(DP)", "feldspar", "feldspar_B", "Pl(I1,HP)", "Fsp(C1)", "Do(HP)", "M(HP)", "Do(AE)", "Cc(AE)", "oCcM(HP)", "Carb(M)", "oCcM(EF)", "dis(EF)", "IlHm(A)", "IlGkPy", "Ilm(WPH)", "Ilm(WPH0)", "Neph(FB)", "Chum", "Atg(PN)", "B", "Pu(M)", "Stlp(M)", "Wus",)
        common_names = ("nepheline", "nepheline", "kalsilite", "larnite", "eskolaite", "sylvite", "halite", "pyrite", "troilite", "millerite", "covellite", "bornite", "chalcopyrite", "anhydrite", "malachite", "azurite", "Fe-ringwoodite", "Mg-ringwoodite", "Fe-wadsleyite", "Mg-wadsleyite", "perovskite", "Al-perovskite", "Na-perovskite", "Mg-bridgmanite", "Fe-bridgmanite", "post-perovskite", "ferropericlase", "majorite", "wustite", "Na-garnet", "akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermakite", "diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "Mg-tschermakite", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "Fe-carpholite", "Fe-sapphirine(793)", "Mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite", "clinochlore", "daphnite", "Fe-sudoite", "Mn-chlorite", "sudoite", "antigorite", "chrysotile", "Fe-talc", "kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)", "lime", "magnesioferrite", "magnetite", "manganosite", "nickel oxide", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica", "Aqueous fluid", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "melt", "melt", "melt", "melt", "melt", "melt", "melt", "melt", "melt", "chlorite", "chlorite", "chlorite", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "plagioclase", "plagioclase", "plagioclase", "k-feldspar", "k-feldspar", "k-feldspar", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "chloritoid", "chloritoid", "chloritoid", "staurolite", "staurolite", "biotite", "biotite", "biotite", "biotite", "biotite", "cordierite", "cordierite", "sapphirine", "sapphirine", "sapphirine", "sapphirine", "osumilite", "fluid", "fluid", "fluid", "fluid", "fluid", "talc", "scapolite", "carpholite", "carpholite", "carpholite", "sudoite", "sudoite", "sudoite", "orthoamphibole", "orthoamphibole", "orthoamphibole", "ternary feldspar", "ternary feldspar", "ternary feldspar", "ternary feldspar", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "ilmenite", "ilmenite", "ilmenite", "ilmenite", "nepheline", "clinohumite", "serpentine", "brucite", "pumpellyite", "stilpnomelane", "wüstite",)
        t = name .== abbreviations
        if any(t)
            common_names[findfirst(t)]
        elseif contains(name,"anorthite")
            "anorthite"
        elseif contains(name,"albite")
            "albite"
        elseif contains(name,"orthoclase")
            "orthoclase"
        else
            name
        end
    end
    export perplex_common_name

## --- End of File
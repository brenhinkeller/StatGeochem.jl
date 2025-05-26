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
    \tmode_basis::String="vol",  #["vol", "wt", "mol"]
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
        mode_basis::String="vol",
        composition_basis::String="wt",
        fluid_eos::Integer=5
    )

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
    for i in eachindex(composition, elements)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(first(P_range))\n$(last(P_range))\ny\n") # v7.1.6/7.1.8
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n5\n3\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(first(P_range))\n$(last(P_range))\ny\n") # v6.8.1

    # Whole-rock composition
    for i ∈ eachindex(composition, elements)
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

    println("Built problem definition")

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
    \tmode_basis::String="vol",  #["vol", "wt", "mol"]
    \tcomposition_basis::String="wt",  #["vol", "wt", "mol"]
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
    for i in eachindex(composition, elements)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\nn\nn\n2\n$(first(T_range))\n$(last(T_range))\n$P\ny\n") # v7.1.8
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n$fluid_eos\nn\nn\n2\n$(first(T_range))\n$(last(T_range))\n$P\ny\n") # v7.1.6
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n$fluid_eos\n3\nn\nn\n2\n$(first(T_range))\n$(last(T_range))\n$P\ny\n") # v6.8.1

    # Whole-rock composition
    for i ∈ eachindex(composition, elements)
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
    \tcomposition_basis::String="wt",  #["vol", "wt", "mol"]
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
    composition_basis::String="wt",  #["vol", "wt", "mol"]
    nonlinear_subdivision::Bool=false,
    fluid_eos::Integer=5,
    fractionate::Integer=0,
    )

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
    for i in eachindex(composition, elements)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\ny\n$PTfilename\ny\n") #7.1.8
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n5\ny\n$PTfilename\ny\n") #7.1.6
    # write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n$fluid_eos\ny\n$PTfilename\n2\ny\n") #6.8.7

    # Whole-rock composition
    for i ∈ eachindex(composition, elements)
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
    \tmode_basis::String="vol", #["vol", "wt", "mol"]
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
        mode_basis::String="vol",
        composition_basis::String="wt",
        fluid_eos::Number=5
    )        

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
    for i in eachindex(composition, elements)
        if !isnan(composition[i])
            elementstring *= elements[i] * "\n"
        end
    end
    write(fp,"$index\n$dataset\nperplex_option.dat\nn\n2\nn\nn\nn\n$elementstring\nn\n2\n$(first(T))\n$(last(T))\n$(first(P))\n$(last(P))\ny\n") # v6.8.7

    # Whole-rock composition
    for i ∈ eachindex(composition, elements)
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

## -- Perplex interface: 2. 0d queries

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

## --- Perplex interface: 3. 1d queries

# # We'll need this for when perplex messes up
# molarmass = Dict("SIO2"=>60.083, "TIO2"=>79.8651, "AL2O3"=>101.96007714, "FE2O3"=>159.6874, "FEO"=>71.8442, "MGO"=>40.304, "CAO"=>56.0774, "MNO"=>70.9370443, "NA2O"=>61.978538564, "K2O"=>94.19562, "H2O"=>18.015, "CO2"=>44.009, "P2O5"=>141.942523997)

"""
```julia
perplex_query_seismic(scratchdir::String;
    \tdof::Integer=1, index::Integer=1, include_fluid="n",
    \tmanual_grid::Bool=false, npoints::Integer=100)
```

Query perplex seismic results along a previously configured 1-d path (dof=1,
isobar or geotherm) or 2-d grid / pseudosection (dof=2).
Results are returned as a dictionary.
"""
function perplex_query_seismic(scratchdir::String;
    dof::Integer=1, index::Integer=1, include_fluid::String="n", importas=:Dict,
    manual_grid::Bool=false, npoints::Integer=100)
    # Query a pre-defined path (isobar or geotherm)

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
    \tindex::Integer=1, npoints::Integer=200, include_fluid="n")
```
```julia
perplex_query_seismic(scratchdir::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, npoints::Integer=200, include_fluid="n")
```
Query perplex seismic results along a specified P-T path using a pre-computed
pseudosection. Results are returned as a dictionary.
"""
function perplex_query_seismic(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    index::Integer=1, npoints::Integer=200, include_fluid="n", importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
    index::Integer=1, include_fluid="n", importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
    \tdof::Integer=1, index::Integer=1, include_fluid="y", clean_units::Bool=true
    \tmanual_grid::Bool=false, npoints::Integer=100)
```

Query all perplex-calculated properties for a specified phase (e.g. "Melt(G)")
along a previously configured 1-d path (dof=1, isobar, geotherm, or P–T path) or 2-d
grid / pseudosection (dof=2). Results are returned as a dictionary.
"""
function perplex_query_phase(scratchdir::String, phase::String;
    dof::Integer=1, index::Integer=1, include_fluid="y", clean_units::Bool=true, importas=:Dict,
    manual_grid::Bool=false, npoints::Integer=100)
    # Query a pre-defined path (isobar or geotherm)

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")

    if manual_grid
        # Edit perplex_option.dat to specify a manual grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   F |/\" -i.backup $(prefix)perplex_option.dat")
        # if dof == 1
            #v7.1.8+, 1d path
            write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\nn\n$npoints\n5\n0\n") #5 is for the case of immiscible phases
    
            # v7.1.6, 1d path
            # write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\nn\n1\n0\n")
    
            # v6.7.8, 1d path
            # write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\n5\n0\n")
            # If a named phase (e.g. feldspar) has multiple immiscible phases, average them (5)
        # elseif dof == 2
            # 2d grid
            # write(fp,"$index\n2\n36\n2\n$phase\n$include_fluid\nn\n1\n0\n") # v6.7.8/v7.1.6
        # else
            # error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        # end
    else 
        # Edit perplex_option.dat to specify an automatic grid size
        system("sed -e \"s/sample_on_grid .*|/sample_on_grid                   T |/\" -i.backup $(prefix)perplex_option.dat")

        #v7.1.8+, 1d path
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
        result = elementify(data,elements, skipstart=1,importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
"""
```julia
perplex_query_phase(scratchdir::String, phase::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1, npoints::Integer=200, include_fluid="y", clean_units::Bool=true, importas=:Dict)
```
```julia
perplex_query_phase(scratchdir::String, phase::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, npoints::Integer=200, include_fluid="y", clean_units::Bool=true, importas=:Dict)
```

Query all perplex-calculated properties for a specified phase (e.g. "Melt(G)")
along a specified P-T path using a pre-computed pseudosection. Results are
returned as a dictionary.
"""
function perplex_query_phase(scratchdir::String, phase::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    index::Integer=1, npoints::Integer=200, include_fluid="y", clean_units::Bool=true, importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        result = elementify(data,elements, skipstart=1,importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
function perplex_query_phase(scratchdir::String, phase::String, P::AbstractArray, T::AbstractArray;
    index::Integer=1, include_fluid="y", clean_units::Bool=true, importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        result = elementify(data,elements, skipstart=1,importas=importas)
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
    \tdof::Integer=1, index::Integer=1, include_fluid="y",
    \tmanualgrid::Bool=false, npoints::Integer=100)
```

Query modal mineralogy (mass proportions) along a previously configured 1-d
path (dof=1, isobar, geotherm, or P–T path) or 2-d grid / pseudosection (dof=2).
Results are returned as a dictionary.

Currently returns wt% 
"""
function perplex_query_modes(scratchdir::String;
    dof::Integer=1, index::Integer=1, include_fluid="y", importas=:Dict,
    manual_grid::Bool=false, npoints::Integer=100)
    # Query a pre-defined path (isobar or geotherm)

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
            t_steps = unique(table["T(K)"])
            result = Dict{String, Vector{Float64}}(i => zeros(length(t_steps)) for i in phase_names)
            id = 1
            # Loop through table
            for t in t_steps 
                # Index table 
                t_idx = table["T(K)"] .== t
                # Index phase name and weight(kg) 
                name = table["Name"][t_idx]
                kg = table["phase,kg"][t_idx]
                # Calculate wt% and add to results dictionary
                for i in zip(name, kg)
                    result[i[1]][id] = (i[2]/nansum(kg)) * 100
                end
                id+=1
            end
            result["T(K)"] = t_steps
        end
    else 
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary.
        # Perplex sometimes returns duplicates of a single solution model, sum them.
        result = elementify(data, sumduplicates=true, importas=importas)
    end
    return result
end
"""
```julia
perplex_query_modes(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1, npoints::Integer=200, include_fluid="y", manual_grid::Bool=false)
```
```julia
perplex_query_modes(scratchdir::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, npoints::Integer=200, include_fluid="y", manual_grid::Bool=false)
```

Query modal mineralogy (mass proportions) along a specified P-T path using a
pre-computed pseudosection. Results are returned as a dictionary.
"""
function perplex_query_modes(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    index::Integer=1, npoints::Integer=200, include_fluid="y", importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        result = elementify(data, sumduplicates=true, importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
function perplex_query_modes(scratchdir::String, P::AbstractArray, T::AbstractArray;
    index::Integer=1, include_fluid="y", importas=:Dict)
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        result = elementify(data, sumduplicates=true, importas=importas)
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
    \tindex::Integer=1, include_fluid="y", clean_units::Bool=true, dof::Integer=1, importas=:Dict,
    \tmanual_grid::Bool=false, npoints::Integer=100)
```?

Query all perplex-calculated properties for the system (with or without fluid)
along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d
grid / pseudosection (dof=2). Results are returned as a dictionary.
Set include_fluid="n" to return solid+melt only.
"""
function perplex_query_system(scratchdir::String;
    index::Integer=1, include_fluid="y", clean_units::Bool=true, dof::Integer=1, importas=:Dict,
    manual_grid::Bool=false, npoints::Integer=100)
    # Query a pre-defined path (isobar or geotherm)
    werami = joinpath(Perple_X_jll.PATH[], "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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

    # if dof == 1
    #     # v6.7.8, 1d path
    #     write(fp,"$index\n3\n36\n1\n$include_fluid\n0\n")
    # elseif dof == 2
    #     # v6.7.8, 2d grid
    #     write(fp,"$index\n2\n36\n1\n$include_fluid\nn\n1\n0\n")
    # else
    #     error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
    # end
    
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
        result = elementify(data,elements, skipstart=1,importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
"""
```julia
function perplex_query_system(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    \tindex::Integer=1, npoints::Integer=200, include_fluid="y",clean_units::Bool=true)
```
```julia
function perplex_query_system(scratchdir::String, P::AbstractArray, T::AbstractArray;
    \tindex::Integer=1, npoints::Integer=200, include_fluid="y",clean_units::Bool=true)
```

Query all perplex-calculated properties for the system (with or without fluid)
along a specified P-T path using a pre-computed pseudosection. Results are
returned as a dictionary. Set include_fluid="n" to return solid+melt only.
"""
function perplex_query_system(scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    index::Integer=1, npoints::Integer=200, include_fluid="y", clean_units::Bool=true, importas=:Dict)
    
    # Query a new path from a pseudosection
    werami = joinpath(Perple_X_jll.PATH[], "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        result = elementify(data,elements, skipstart=1,importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
function perplex_query_system(scratchdir::String, P::AbstractArray, T::AbstractArray;
    index::Integer=1, include_fluid="y", clean_units::Bool=true, importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(Perple_X_jll.PATH[], "werami") # path to PerpleX werami
    prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        result = elementify(data,elements, skipstart=1,importas=importas)
    catch
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end
export perplex_query_system
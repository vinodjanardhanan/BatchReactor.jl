module BatchReactor

using LightXML, Printf
using DifferentialEquations, Sundials
using IdealGas, GasphaseReactions, SurfaceReactions, ReactionCommons, RxnHelperUtils

#include("Constants.jl")


export batch_reactor

global o_streams

struct ConstantParams{T1}
    Asv::T1
    T::T1
end

struct UsrMech <: MechanismDefinition
end


"""
A common function is written for the reading input data
for all methods that utilizes batch reactor model. The input
    struct defines the common Parameters    
"""
struct InputData
    T::Float64
    p_initial::Float64
    Asv::Float64    
    tf::Float64
    gasphase::Array{String,1}
    mole_fracs::Array{Float64,1}
    thermo_obj::SurfaceReactions.IdealGas.SpeciesThermoObj
    gmd::Union{MechanismDefinition, Nothing}
    smd::Union{MechanismDefinition, Nothing} 
    umd::Union{MechanismDefinition, Nothing} 
end


"""
This is the calling function for executing the batch reactor with user defined rate calculation
#   Usage
    batch_reactor(input_file::AbstractString, lib_dir::AbstractString, user_defined::Function; sens= false)    
-   input_file: the xml input file for batch reactor
-   lib_dir: the direcrtory in which the data files are present. It must be the relative path
-   user_defined : a function that calculates the species source terms. Must be supplied by the user
-   sens : Boolean to be set as true when used along with the sensitivity code
"""
function batch_reactor(input_file::AbstractString, lib_dir::AbstractString, user_defined::Function; sens= false)    
    chem = Chemistry(false, false, true, user_defined)
    batch_reactor(input_file, lib_dir, sens, chem)
end


"""
This is the calling function for executing the batch reactor with chemistry input files 
#   Usage
batch_reactor(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
-   input_file: the xml input file for batch reactor
-   lib_dir: the direcrtory in which the data files are present. It must be the relative path
-   sens : Boolean to be set as true when used along with the sensitivity code
-   surfchem : Boolean to specify the calculation of surface reaction rates 
-   gaschem : Boolean to specify the calculation of gasphase reaction rates 
"""
function batch_reactor(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
    chem = Chemistry(surfchem, gaschem, false, f->())
    batch_reactor(input_file, lib_dir, sens, chem)
end



"""
A function for call from other packages, mainly intended for reactor network modeling
# Usage
batch_reactor(inlet_comp, T, p, time; Asv=1.0, chem, thermo_obj, md)
-   inlet_comp : A dictionary of species and its mole fractions at the reactor inlet 
-   T : operating temperature (K)
-   p : operating pressure (Pa) 
-   time : integration time (s)   
-   Asv : Surface area to volume ratio (important in the case of surface reactions. 1 in the case of gasphase chemistry)
-   thermo_obj : Species thermo object (Please refer IdealGas package documentation)
-   md : MechanismDefinition (Please refer to ReactionCommons documentation)
"""
function batch_reactor(inlet_comp, T, p, time; Asv=1.0, chem, thermo_obj, md)

    # constant parameters 
    cp = ConstantParams(Asv,T) 


    function  get_mole_fracs(species, inlet_comp)
        mole_fracs = zeros(length(species))
        for k in eachindex(species)
            if in(species[k], collect(keys(inlet_comp)))
                mole_fracs[k] = inlet_comp[species[k]]
            end
        end
        mole_fracs
    end

    if chem.surfchem
        species = collect(keys(inlet_comp))
        mole_fracs = get_mole_fracs(species, inlet_comp)
        covg = md.sm.si.ini_covg
        #Create the state object    
        surf_conc = similar(covg)
        rxn_rate = zeros(length(md.sm.reactions))
        n_species = length(mole_fracs) + length(covg)
        source = zeros(n_species)
        all_conc = zeros(n_species)
        sr_state = SurfaceRxnState(T, p, mole_fracs, covg, surf_conc, rxn_rate, source, all_conc)        
        params = (s_state = sr_state, thermo = thermo_obj, smd = md, cp = cp, chem = chem)        
        
    end

    if chem.gaschem
        species = md.gm.species
        mole_fracs = get_mole_fracs(md.gm.species, inlet_comp)
        conc = zeros(length(mole_fracs))
        source = zeros(length(conc))
        g_all = zeros(length(mole_fracs))
        Kp = zeros(length(md.gm.reactions))
        rxn_rate = zeros(length(Kp))
        gs_state = GasphaseState(T, p, mole_fracs, conc, rxn_rate, source, g_all, Kp)
        params = (g_state = gs_state, thermo = thermo_obj, gmd =  md, cp =  cp,  chem = chem)
    end

    #create the solution vector
    soln = zeros(length(mole_fracs))
    molefrac_to_massfrac!(soln, mole_fracs, thermo_obj.molwt)
    soln .*= density(mole_fracs, thermo_obj.molwt, T, p)
    if chem.surfchem
        append!(soln, covg)
    end

    
    alg = CVODE_BDF()
    t_span = (0,time)
    prob = ODEProblem(residual!,soln,t_span,params)
    sol = solve(prob, alg , reltol=1e-6, abstol=1e-10, save_everystep=false);        
    mass_fracs = sol.u[end][1:length(species)] ./ sum(sol.u[end][1:length(species)])
    mole_fracs = zeros(length(mass_fracs))
    massfrac_to_molefrac!(mole_fracs, mass_fracs, thermo_obj.molwt)
    return sol.t, Dict(species .=> mole_fracs)
    
end

#=
This is the general interface for staring the batch reactor code by reading the input files
=#
function batch_reactor(input_file::AbstractString, lib_dir::AbstractString, sens::Bool, chem::Chemistry)    
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)
    # read the input data 
    id = input_data(xmlroot, lib_dir, chem)    
    # create solution vector 
    soln = get_solution_vector(id, chem)
    # calculate the total number of species    
    n_species = length(id.gasphase) 
    if chem.surfchem
        n_species += length(id.smd.sm.species)
    end
    #storage for species source terms
    source = zeros(n_species)
    #space allocation for concentration 
    all_conc = zeros(n_species)
    #create output files for saving the data
    # g_stream = open("gas_profile.dat","w")
    g_stream = open(output_file(input_file, "gas_profile.dat"),"w")
    s_stream = open(output_file(input_file, "surface_covg.dat"),"w")  
    csv_g_stream = open(output_file(input_file, "gas_profile.csv"),"w")
    csv_s_stream = open(output_file(input_file, "surface_covg.csv"),"w")  
    global o_streams = (g_stream, s_stream, csv_g_stream, csv_s_stream)  
    create_header(g_stream,["t","T","p","rho"],id.gasphase)
    write_csv(csv_g_stream,["t","T","p","rho"],id.gasphase)
    if chem.surfchem
        create_header(s_stream,"t", "T" ,id.smd.sm.species)
        write_csv(csv_s_stream,"t", "T" ,id.smd.sm.species)
    end

    #define the Parameters
    cp = ConstantParams(id.Asv, id.T) 
    sr_state = gr_state = ur_state = nothing
    
    if chem.surfchem 
        surf_conc = similar(id.smd.sm.si.ini_covg)
        rxn_rate = zeros(length(id.smd.sm.reactions))        
        sr_state = SurfaceRxnState(id.T, id.p_initial, id.mole_fracs, id.smd.sm.si.ini_covg, surf_conc, rxn_rate, source, all_conc)
    end
    if chem.gaschem
        Kp = zeros(length(id.gmd.gm.reactions)) # Non allocating memory or the calculation of equilibrium constant 
        rxn_rate = zeros(length(Kp)) # Non allocating memory for the calculation of individual reactions 
        g_all = similar(id.mole_fracs) # Non allocating memory for the calculation of Gibb's free energy 
        gr_state = GasphaseState(id.T, id.p_initial, id.mole_fracs, all_conc, rxn_rate, source, g_all, Kp)
    end
    if chem.userchem && ! chem.surfchem && !chem.gaschem
        source = zeros(length(id.mole_fracs))         
        ur_state = UserDefinedState(id.T, id.p_initial, id.mole_fracs, id.thermo_obj.molwt, id.gasphase, source)
    end
    t_span = (0,id.tf)    
    
    params = (s_state = sr_state, g_state = gr_state, u_state = ur_state, thermo = id.thermo_obj, smd = id.smd, gmd = id.gmd, cp = cp, chem= chem)
    prob = ODEProblem(residual!,soln,t_span,params)
    if sens == true
        return (params, prob,t_span)
    end
    cb = FunctionCallingCallback(save_data)
    # sol = solve(prob, reltol=1e-6, abstol = 1e-8, save_everystep=false, callback=cb)
    sol = solve(prob, CVODE_BDF(), reltol=1e-6, abstol=1e-10, save_everystep=false,callback=cb);        
    
    close(g_stream)
    close(s_stream)    
    close(csv_g_stream)
    close(csv_s_stream)
    return Symbol(sol.retcode)
end


#=
Function for creating the solition vector. 
    The solution vector contains mass density of the species ie. ρ × Y_k    
=#
function get_solution_vector(id::InputData, chem)
    soln = zeros(length(id.gasphase))
    molefrac_to_massfrac!(soln,id.mole_fracs,id.thermo_obj.molwt)
    soln .*= density(id.mole_fracs,id.thermo_obj.molwt,id.T,id.p_initial)
    if chem.surfchem
        append!(soln,id.smd.sm.si.ini_covg)    
    end
    return soln
end


#=
Function for reading the common input parameters
=#
function input_data(xmlroot::XMLElement, lib_dir::AbstractString, chem)


    # locate the lib directory    
    thermo_file = "therm.dat"
    thermo_file = get_path(lib_dir, thermo_file)


    """ In case the gasphase mechanism is present then the gasphase species are read from the 
     gasphase mechanism file and not from the xml input. So this must be the first action before any
     other input parameters are read 
     """
    if chem.gaschem
        mech_file = get_text_from_xml(xmlroot,"gas_mech")
        # mech_file = lib_dir*"/"*mech_file      
        mech_file = get_path(lib_dir, mech_file)  
        gmd = compile_gaschemistry(mech_file)
        gasphase = gmd.gm.species
    else # If gasphase chemistry is not present, then get the gasphase species from xml
        #get the gasphase present
        gasphase = Array{String,1}
        gasphase = get_collection_from_xml(xmlroot,"gasphase")    
        gmd = nothing
    end 

    
    # create the thermo object 
    thermo_obj = IdealGas.create_thermo(gasphase, thermo_file)

    #Get the mole fractions from xml
    mole_fracs = get_molefraction_from_xml(xmlroot, thermo_obj.molwt, gasphase)
    
    #Get the temperature
    local T = get_value_from_xml(xmlroot,"T")

    #get the initial pressure 
    local p_initial = get_value_from_xml(xmlroot,"p")

    #get the surface area to volume ratio 
    local Asv = get_value_from_xml(xmlroot,"Asv")

    #get the integration time 
    local tf = get_value_from_xml(xmlroot,"time")

    #get the mechanism file if surface chemistry is involked 
    if chem.surfchem
        mech_file = get_text_from_xml(xmlroot,"surface_mech")
	    mech_file = get_path(lib_dir, mech_file) 
        #create the mechanism definition
        smd = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)
        # id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj,gmd,smd)
    elseif !chem.surfchem
        smd = nothing
    end
    

    # create the input data if Gasphase chemistry is invoked 
    if chem.gaschem || chem.surfchem
        id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj, gmd, smd, nothing)
    end

    if chem.userchem
        um = UsrMech()
        id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj,gmd, smd, um)
    end

    return id

end


#=
Residual function definition for the batch reactor
=#
function residual!(du,u,p,t)        

    if p[:chem].surfchem && !p[:chem].gaschem
        ng = length(p[:s_state].mole_frac)
    elseif p[:chem].gaschem && !p[:chem].surfchem
        ng = length(p[:g_state].mole_frac)
    elseif p[:chem].gaschem && p[:chem].surfchem
        ng = length(p[:g_state].mole_frac)
    elseif p[:chem].userchem
        ng = length(p[:u_state].mole_frac)
    end
    
    
    #density 
    ρ = sum(u[1:ng])
    #mass fractions 
    mass_fracs = u[1:ng]/ρ
        
    
    #calculate the molar production rates in case of surface chemistry 
    if p[:chem].surfchem 
        #mole fractions
        massfrac_to_molefrac!(p[:s_state].mole_frac,mass_fracs,p[:thermo].molwt)
        #average molecular weight
        mlwt_avg = average_molwt(p[:s_state].mole_frac,p[:thermo].molwt)
        #pressure update
        p[:s_state].p = ρ*RxnHelperUtils.R*p[:cp].T/mlwt_avg


        ns = length(p[:smd].sm.species)    
        #coverage update
        p[:s_state].covg = u[ng+1:ng+ns]
        SurfaceReactions.calculate_molar_production_rates!(p[:s_state],p[:thermo],p[:smd])
        p[:s_state].source *= p[:cp].Asv
    end    
    if p[:chem].gaschem
        #mole fractions
        massfrac_to_molefrac!(p[:g_state].mole_frac,mass_fracs,p[:thermo].molwt)
        #average molecular weight
        mlwt_avg = average_molwt(p[:g_state].mole_frac,p[:thermo].molwt)
        #pressure update
        p[:g_state].p = ρ*RxnHelperUtils.R*p[:cp].T/mlwt_avg

        GasphaseReactions.calculate_molar_production_rates!(p[:g_state], p[:gmd], p[:thermo])
    end
    
    if p[:chem].userchem
        p[:chem].udf(p[:u_state])
    end

    #gasphase residual
    if p[:chem].gaschem && !p[:chem].surfchem
        du[1:ng] = (p[:g_state].source[1:ng] .* p[:thermo].molwt)     
    elseif !p[:chem].gaschem && p[:chem].surfchem
        du[1:ng] = (p[:s_state].source[1:ng] .* p[:thermo].molwt)     
        du[ng+1:ng+ns] = (p[:s_state].source[ng+1:ng+ns] .* p[:smd].sm.si.site_coordination)/(p[:smd].sm.si.density*1e4)
    elseif p[:chem].gaschem && p[:chem].surfchem
        du[1:ng] = ( (p[:s_state].source[1:ng] .+ p[:g_state].source[1:ng]) .* p[:thermo].molwt)     
        du[ng+1:ng+ns] = (p[:s_state].source[ng+1:ng+ns] .* p[:smd].sm.si.site_coordination)/(p[:smd].sm.si.density*1e4)
    elseif p[:chem].userchem
        du[1:ng] = (p[:u_state].source[1:ng] .* p[:thermo].molwt)     
    end

    
end



#=
Function to save the variables into output file
=#
function save_data(u,t,integrator)
    
    g_stream, s_stream, csv_g_stream, csv_s_stream = o_streams
    #density 
    if integrator.p[:chem].surfchem 
        state = integrator.p[:s_state]        
    elseif integrator.p[:chem].gaschem
        state = integrator.p[:g_state]
    else
        state = integrator.p[:u_state]
    end
    ρ = sum(u[1:length(state.mole_frac)])
    write_to_file(g_stream,t,state.T,state.p,ρ,state.mole_frac)
    write_csv(csv_g_stream,t,state.T,state.p,ρ,state.mole_frac)
    if integrator.p[:chem].surfchem
        write_to_file(s_stream,t,integrator.p[:s_state].T,integrator.p[:s_state].covg)   
        write_csv(csv_s_stream,t,integrator.p[:s_state].T,integrator.p[:s_state].covg)   
    end
    @printf("%4e\n",t) 
end



end

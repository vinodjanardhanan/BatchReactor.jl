module BatchReactor

using LightXML, Printf
using Sundials, DifferentialEquations
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
    md::MechanismDefinition
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
        params = (sr_state, thermo_obj, md, cp, chem)        
        
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
        params = (gs_state, thermo_obj, md, cp,  chem)
    end

    #create the solution vector
    soln = zeros(length(mole_fracs))
    molefrac_to_massfrac!(soln, mole_fracs, thermo_obj.molwt)
    soln .*= density(mole_fracs, thermo_obj.molwt, T, p)
    if chem.surfchem
        append!(soln, covg)
    end

    t_span = (0,time)
    prob = ODEProblem(residual!,soln,t_span,params)
    sol = solve(prob, CVODE_BDF(), reltol=1e-6, abstol=1e-10, save_everystep=false);        
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
        n_species += length(id.md.sm.species)
    end
    #storage for species source terms
    source = zeros(n_species)
    #space allocation for concentration 
    all_conc = zeros(n_species)
    #create output files for saving the data
    # g_stream = open("gas_profile.dat","w")
    g_stream = open(output_file(input_file, "gas_profile.dat"),"w")
    s_stream = open(output_file(input_file, "surface_covg.dat"),"w")  
    global o_streams = (g_stream, s_stream)  
    create_header(g_stream,["t","T","p","rho"],id.gasphase)
    if chem.surfchem
        create_header(s_stream,"t", "T" ,id.md.sm.species)
    end

    #define the Parameters
    cp = ConstantParams(id.Asv, id.T) 
    
    if chem.surfchem
        surf_conc = similar(id.md.sm.si.ini_covg)
        rxn_rate = zeros(length(id.md.sm.reactions))        
        state = SurfaceRxnState(id.T, id.p_initial, id.mole_fracs, id.md.sm.si.ini_covg, surf_conc, rxn_rate, source, all_conc)
    end
    if chem.gaschem
        Kp = zeros(length(id.md.gm.reactions)) # Non allocating memory or the calculation of equilibrium constant 
        rxn_rate = zeros(length(Kp)) # Non allocating memory for the calculation of individual reactions 
        g_all = similar(id.mole_fracs) # Non allocating memory for the calculation of Gibb's free energy 
        state = GasphaseState(id.T, id.p_initial, id.mole_fracs, all_conc, rxn_rate, source, g_all, Kp)
    end

    if chem.userchem && ! chem.surfchem && !chem.gaschem
        source = zeros(length(id.mole_fracs))         
        state = UserDefinedState(id.T, id.p_initial, id.mole_fracs, id.thermo_obj.molwt, id.gasphase, source)
    end

    t_span = (0,id.tf)
    params = (state, id.thermo_obj, id.md, cp, chem)
    prob = ODEProblem(residual!,soln,t_span,params)
    if sens == true
        return (params, prob,t_span)
    end
    cb = FunctionCallingCallback(save_data)
    # sol = solve(prob,alg_hints=[:stiff] , reltol=1e-6, abstol = 1e-8, save_everystep=false, callback=cb)
    sol = solve(prob, CVODE_BDF(), reltol=1e-6, abstol=1e-10, save_everystep=false,callback=cb);        
    
    close(g_stream)
    close(s_stream)    

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
        append!(soln,id.md.sm.si.ini_covg)    
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
    end 

    # If gasphase chemistry is not present, then get the gasphase species from xml
    if !chem.gaschem
        #get the gasphase present
        gasphase = Array{String,1}
        gasphase = get_collection_from_xml(xmlroot,"gasphase")
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
        md = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)
        id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj,md)
    end

    # create the input data if Gasphase chemistry is invoked 
    if chem.gaschem
        id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj,gmd)
    end

    if chem.userchem
        um = UsrMech()
        id = InputData(T,p_initial,Asv,tf,gasphase,mole_fracs,thermo_obj,um)
    end

    return id

end


#=
Residual function definition for the batch reactor
=#
function residual!(du,u,p,t)        
    state = 1
    thermo_obj = 2
    md = 3
    cp = 4
    ng = length(p[state].mole_frac)
    

    #density 
    ρ = sum(u[1:ng])
    #mass fractions 
    mass_fracs = u[1:ng]/ρ
    #mole fractions
    massfrac_to_molefrac!(p[state].mole_frac,mass_fracs,p[thermo_obj].molwt)
    #average molecular weight
    mlwt_avg = average_molwt(p[state].mole_frac,p[thermo_obj].molwt)
    #pressure update
    p[state].p = ρ*RxnHelperUtils.R*p[cp].T/mlwt_avg
        
    
    #calculate the molar production rates in case of surface chemistry 
    if p[end].surfchem
        ns = length(p[md].sm.species)    
        #coverage update
        p[state].covg = u[ng+1:ng+ns]
        SurfaceReactions.calculate_molar_production_rates!(p[state],p[thermo_obj],p[md])
        p[state].source *= p[cp].Asv
    end
    if p[end].gaschem && !p[end].surfchem
        GasphaseReactions.calculate_molar_production_rates!(p[state], p[md], p[thermo_obj])
    end
    
    if p[end].userchem
        p[end].udf(p[state])
    end
    #gasphase residual
    du[1:ng] = (p[state].source[1:ng] .* p[thermo_obj].molwt)     
    #coverage residuals in case of surface chemistry 
    if p[end].surfchem
        du[ng+1:ng+ns] = (p[state].source[ng+1:ng+ns] .* p[md].sm.si.site_coordination)/(p[md].sm.si.density*1e4)
    end
    
end



#=
Function to save the variables into output file
=#
function save_data(u,t,integrator)
    state = 1
    g_stream, s_stream = o_streams
    #density 
    ρ = sum(u[1:length(integrator.p[state].mole_frac)])
    write_to_file(g_stream,t,integrator.p[state].T,integrator.p[state].p,ρ,integrator.p[state].mole_frac)
    if integrator.p[end].surfchem
        write_to_file(s_stream,t,integrator.p[state].T,integrator.p[state].covg)   
    end
    @printf("%4e\n",t) 
end



end

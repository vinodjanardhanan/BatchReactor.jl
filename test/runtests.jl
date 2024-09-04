using BatchReactor
using Test
using IdealGas, RxnHelperUtils, SurfaceReactions, GasphaseReactions, ReactionCommons

@testset "BatchReactor.jl" begin
        
    if Sys.isapple()  || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end

     @testset "Testing surface chemistry" begin        
        input_file = joinpath("batch_surf", "batch.xml")
        retcode = batch_reactor(input_file, lib_dir, surfchem=true)
        @test retcode == Symbol("Success")
    end

    @testset "Testing gas chemistry of h2 + o2" begin
        input_file = joinpath("batch_h2o2", "batch.xml")
        retcode = batch_reactor(input_file, lib_dir, gaschem=true)
        @test retcode == Symbol("Success")
    end

    @testset "Testing gas chemistry using grimech" begin
        input_file = joinpath("batch_ch4", "batch.xml")
        retcode = batch_reactor(input_file, lib_dir, gaschem=true)
        @test retcode == Symbol("Success")
    end
 
    @testset "Testing surface chemistry with interface call " begin
        inlet_comp = Dict("CH4"=>0.25,"H2O"=>0.0, "H2"=>0.0, "CO"=>0.0, "CO2"=>0.25, "O2"=>0.0, "N2"=>0.5)
        T = 1073.15
        p = 1e5
        t = 10
        gasphase = collect(keys(inlet_comp))
        thermo_obj = IdealGas.create_thermo(gasphase, get_path(lib_dir,"therm.dat"))
        md = SurfaceReactions.compile_mech(get_path(lib_dir,"ch4ni.xml"),thermo_obj,gasphase)
        chem = Chemistry(true, false, false, f->())
        geometry = (dia=0.005,cat_geom=1000.0)        
        retcodes = batch_reactor(inlet_comp, T, p, t; Asv = 10.0, chem=chem, thermo_obj=thermo_obj, md=md)                       
        @test retcodes[1][end] == t
    end

    @testset "Testing gasphase chemistry with interface call " begin
        
        mech_file = get_path(lib_dir, "h2o2.dat")
        gmd = compile_gaschemistry(mech_file)        
        gasphase = gmd.gm.species
        thermo_obj = IdealGas.create_thermo(gasphase, get_path(lib_dir,"therm.dat"))        
        inlet_comp = Dict("O2" => 0.25,"N2" => 0.5, "H2" => 0.25)
        
        T = 1073.15
        p = 1e5
        t = 10.0
                
        chem = Chemistry(false, true, false, f->())      
        geometry = (dia=0.005,cat_geom=1.0)               
        retcodes = batch_reactor(inlet_comp, T, p, t; chem=chem, thermo_obj=thermo_obj, md=gmd)                                
        @test retcodes[1][end] == t
    end


    @testset "Testing user defined chemistry " begin        
        function udf(state)
            state.source[1:end] .= 0.0            
        end
        input_file = joinpath("batch_udf", "batch.xml")
        retcode = batch_reactor(input_file, lib_dir, udf)
        @test retcode == Symbol("Success")        
    end 
end
using DrWatson, Test
@quickactivate :ImmuneBoostingODEs

@testset "ImmuneBoostingODEs tests" begin
    
@testset "Equilibrium calculations" begin 
    # parameters with no endemic equilibrium
    @testset let p1 = SirnsParameters(.01, 52., .01, .3, .8)
        @test equils(.01, 52, .01) == equils(p1)
        @test equils(p1) == 1.
        @test equili(p1) == .0
        u1 = equil(p1)
        @test u1 == [1., 0, 0, 0, 0]
    end
    # parameters with an endemic equilibrium
    @testset let p2 = SirnsParameters(88., 52., .01, .3, .8)
        p2b = @ntuple β0 = 88 γ = 52 μ = .01 ω = .8 ψ = .3
        @test equils(88, 52, .01) == equils(p2)
        S2 = equils(p2) 
        @test S2 < 1 
        @test S2 > 0
        @test equils(p2) == equils(p2b)
        I2 = equili(p2)
        @test I2 > 0
        @test I2 < 1
        @test equili(p2) == equili(p2b)
        u2 = equil(p2)
        u2b = equil(p2b)
        @test u2 != [1., 0, 0, 0, 0]
        @test sum(u2) ≈ 1
        @test minimum(u2) >= 0
        @test u2 == u2b
        I3 = equili(p2)
        @test I3 > 0
        @test I3 < 1
    end 
    # parameters with an extremely large birth and mortality rate
    @testset let p3 = SirnsParameters(880., 52., 66., .3, .8)
        S3 = equils(p3) 
        @test S3 < 1 
        @test S3 > 0
        u3 = equil(p3)
        @test u3 != [1., 0, 0, 0, 0]
        @test minimum(u3) >= 0
        @test sum(u3) ≈ 1
    end 
end 
@testset "Equilibrium eigenvalues" begin 
    import ImmuneBoostingODEs: maxequileigen
    # test equileigen for the disease-free equilibrium, where we know the solution 
    R0 = .2 
    γ = 50.
    μ = .001
    ω = .5
    β = R0 * (γ + μ)
    p1 = SirnsParameters(β, γ, μ, .0, ω)
    Λ1 = equileigen(p1)
    Λ1calc = [ [ β - γ - μ, -μ ]; ones(3) .* (-3 * ω - μ) ]
    sort!(Λ1calc) # to put calculated values in same order as equileigen results
    @testset "Compare function-calculated eigenvalues to known results" begin
        @testset for i ∈ eachindex(Λ1calc) @test Λ1[i] == Λ1calc[i] end   
    end
    # test that disease-free equileigen is not influenced by ψ
    p2 = SirnsParameters(β, γ, μ, .0, ω)
    Λ2 = equileigen(p2)
    @testset "Compare disease=free eigenvalues with different ψ" begin
        @testset for i ∈ eachindex(Λ2) @test Λ1[i] == Λ2[i] end  
    end  
    # check that other functions do what they are meant to 
    @testset "Additional eigenvalue functions" begin
        @test maxequileigen(p1) == maximum(Λ1calc)
        @test realmaxequileigen(p1) == maximum(Λ1calc)
        @test realmaxequileigen(β, γ, μ, 100., ω) == maximum(Λ1calc)
    end
end
@testset "Model functions" begin
    import Base: oneto
    import ImmuneBoostingODEs: compartmentinds
    p1 = SirnsParameters(.01, 52., .01, .3, .8)
    @testset "No endemic equilibrium, times starts at 0" begin
        tspan = ( 0., 10. )
        saveat = [ 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p1, equalrs = false, t0 = 0)
        sol = run_sirns(u0, p1, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t)
        end
        inds = compartmentinds(sol)
        @test inds == oneto(6)
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    @testset "No endemic equilibrium, times starts at -1.3" begin
        tspan = ( -1.3, 10.1 )
        saveat = [ -1.3, -1, -.27, 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p1, equalrs = false, t0 = -1.3)
        sol = run_sirns(u0, p1, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t)
        end
        inds = compartmentinds(sol)
        @test inds == 4:9
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    p2 = SirnsParameters(88., 52., .01, .3, .8)
        @testset "Endemic equilibrium, times starts at 0" begin
        tspan = ( 0., 10. )
        saveat = [ 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p2, equalrs = false, t0 = 0)
        sol = run_sirns(u0, p2, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t)
        end
        inds = compartmentinds(sol)
        @test inds == oneto(6)
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    @testset "Endemic equilibrium, times starts at -1.3" begin
        tspan = ( -1.3, 10.1 )
        saveat = [ -1.3, -1, -.27, 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p2, equalrs = false, t0 = -1.3)
        sol = run_sirns(u0, p2, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t)
        end
        inds = compartmentinds(sol)
        @test inds == 4:9
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    p3 = SirnsParameters(88., 52., .01, .3, .8)
    @testset "Large mortality rate, times starts at 0" begin
        tspan = ( 0., 10. )
        saveat = [ 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p3, equalrs = false, t0 = 0)
        sol = run_sirns(u0, p3, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t)
        end
        inds = compartmentinds(sol)
        @test inds == oneto(6)
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    @testset "Large mortality rate, times starts at -1.3" begin
        tspan = ( -1.3, 10.1 )
        saveat = [ -1.3, -1, -.27, 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p3, equalrs = false, t0 = -1.3)
        sol = run_sirns(u0, p3, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t)
        end
        inds = compartmentinds(sol)
        @test inds == 4:9
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    ϕ = .3
    p4 = SirnsParameters(.01, .1, ϕ, 52., .01, .3, .8)
    @testset "Seasonal forcing (1), times starts at 0" begin
        tspan = ( 0., 10. )
        saveat = [ 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p4, equalrs = false, t0 = 0)
        sol = run_sirns(u0, p4, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t - ϕ)
        end
        inds = compartmentinds(sol)
        @test inds == oneto(6)
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    @testset "Seasonal forcing (1), times starts at -1.3" begin
        tspan = ( -1.3, 10.1 )
        saveat = [ -1.3, -1, -.27, 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p4, equalrs = false, t0 = -1.3)
        sol = run_sirns(u0, p4, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t - ϕ)
        end
        inds = compartmentinds(sol)
        @test inds == 4:9
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    ϕ = .4
    p5 = SirnsParameters(88., .15, ϕ, 52., 7.2, .3, .8)
        @testset "Seasonal forcing (2), times starts at 0" begin
        tspan = ( 0., 10. )
        saveat = [ 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p5, equalrs = false, t0 = 0)
        sol = run_sirns(u0, p5, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t - ϕ)
        end
        inds = compartmentinds(sol)
        @test inds == oneto(6)
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
    @testset "Seasonal forcing (2), times starts at -1.3" begin
        tspan = ( -1.3, 10.1 )
        saveat = [ -1.3, -1, -.27, 0, .2, .5, 3.3, 8.65, 10 ]
        S0 = .5 
        I0 = .05
        u0 = sirns_u0(S0, I0; p = p5, equalrs = false, t0 = -1.3)
        sol = run_sirns(u0, p5, tspan; saveat)
        @testset for (i, t) ∈ enumerate(saveat)
            @test sum(sol[i][1:5]) ≈ 1 
            @test minimum(sol[i][1:5]) >= 0
            @test maximum(sol[i][1:5]) <= 1
            # check that "x1" is doing what it is meant to
            @test sol[i][6] ≈ cos(2π * t - ϕ)
        end
        inds = compartmentinds(sol)
        @test inds == 4:9
        cc = modelcompartments(sol, :cc)
        @test length(cc) == length(inds)
        cases = casespertimeblock(cc)
        @test length(cases) == length(inds) - 1
    end
end

end # @testset "ImmuneBoostingODEs tests" 

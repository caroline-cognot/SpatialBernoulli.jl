using SpatialBernoulli
using Test
using Random
rng = MersenneTwister(1234)

#################### D = 36 ######################################
	# define locations in the unit square
	my_locations = vcat(([x y] for x in 0:0.2:1 for y in 0:0.2:1)...)
	nlocs = length(my_locations[:, 1])


	my_distance = [sqrt(sum(abs2, my_locations[i, :] - my_locations[j, :])) for i in axes(my_locations, 1), j in axes(my_locations, 1)]

	# randomly generate a SpatialBernoulli
	my_λ = rand(nlocs)# [0.5+ 0.1*i for i in 1:nlocs] 
	my_range = 0.3
	my_sill = 1.0
	my_order = 1 / 2
	d = SB(my_range, my_sill, my_order, my_λ, my_distance)
	# generate data
	n = 2000
	y=Bool.(rand(rng,d, n))

	# fit a model using initial values
	init_range = 0.5
	init_order = 1.5
	init_lambda = fill(0.4, nlocs)
	init_d = SB(init_range, 1.0, init_order, init_lambda, my_distance)

	# using full likelihood is a bad idea.
	# @timed sol = fit_mle(init_d,y; m = 100*length(init_d), return_sol = false)
	# #maximise without the order !
	#  @timed sol = fit_mle(init_d,y; m = 10*length(init_d), return_sol = false, order = my_order)

	# pairwise maximisation is a lot better
	tdist = maximum(my_distance) / 1
	wp = 1.0 .* (my_distance .< tdist)
	 
    @testset "n= $n, D=$(length(d)), fit_mle, high m, check Bernoulli proba , then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol1 = fit_mle(init_d, y, wp; order = my_order, m = 100 * 2, return_sol = true, maxiters = 2000)
        @test isapprox(my_λ, sol1[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol1[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol1[1].λ; rtol=0.02)
        @test isapprox(my_range, sol1[1].range; rtol=0.1)
        @test isapprox(my_range, sol1[1].range; rtol=0.05)
        @test isapprox(my_range, sol1[1].range; rtol=0.02)
    end 
   


	
    @testset "n= $n,  D=$(length(d)),fit_mle, low m, check Bernoulli proba, then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol2 = fit_mle(init_d, y, wp; order = my_order, m = 30 * 2, return_sol = true, maxiters = 2000)
        @test isapprox(my_λ, sol2[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol2[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol2[1].λ; rtol=0.02)
        @test isapprox(my_range, sol2[1].range; rtol=0.1)
        @test isapprox(my_range, sol2[1].range; rtol=0.05)
        @test isapprox(my_range, sol2[1].range; rtol=0.02)
    end 
   


    @testset "n= $n,  D=$(length(d)),fit_mle_vfast, check Bernoulli proba, then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol3 = fit_mle_vfast(init_d, y, wp; order = my_order, return_sol = true, maxiters = 20000)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.02)
        @test isapprox(my_range, sol3[1].range; rtol=0.1)
        @test isapprox(my_range, sol3[1].range; rtol=0.05)
        @test isapprox(my_range, sol3[1].range; rtol=0.02)
    end 
   


#################### D = 121 ##################################################"""
	# define locations in the unit square
	my_locations = vcat(([x y] for x in 0:0.1:1 for y in 0:0.1:1)...)
	nlocs = length(my_locations[:, 1])


	my_distance = [sqrt(sum(abs2, my_locations[i, :] - my_locations[j, :])) for i in axes(my_locations, 1), j in axes(my_locations, 1)]

	# randomly generate a SpatialBernoulli
	my_λ = rand(nlocs)# [0.5+ 0.1*i for i in 1:nlocs] 
	my_range = 0.3
	my_sill = 1.0
	my_order = 1 / 2
	d = SB(my_range, my_sill, my_order, my_λ, my_distance)
	# generate data
	n = 5000
	y=Bool.(rand(rng,d, n))

	# fit a model using initial values
	init_range = 0.2
	init_order = 1.5
	init_lambda = fill(0.4, nlocs)
	init_d = SB(init_range, 1.0, init_order, init_lambda, my_distance)

	# pairwise maximisation is a lot better
	tdist = maximum(my_distance) / 1
	wp = 1.0 .* (my_distance .< tdist)

    @testset "n= $n,  D=$(length(d)),fit_mle_vfast, check Bernoulli proba, then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol3 = fit_mle_vfast(init_d, y, wp; order = my_order, return_sol = true, maxiters = 20000)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.02)
        @test isapprox(my_range, sol3[1].range; rtol=0.1)
        @test isapprox(my_range, sol3[1].range; rtol=0.05)
        @test isapprox(my_range, sol3[1].range; rtol=0.02)
    end 
   
####################### D = 16 #####################################



	# define locations in the unit square
	my_locations = vcat(([x y] for x in 0:0.3:1 for y in 0:0.3:1)...)
	nlocs = length(my_locations[:, 1])


	my_distance = [sqrt(sum(abs2, my_locations[i, :] - my_locations[j, :])) for i in axes(my_locations, 1), j in axes(my_locations, 1)]

	# randomly generate a SpatialBernoulli
	my_λ = rand(nlocs)# [0.5+ 0.1*i for i in 1:nlocs] 
	my_range = 0.3
	my_sill = 1.0
	my_order = 1 / 2
	d = SB(my_range, my_sill, my_order, my_λ, my_distance)
	# generate data
	n = 2000
	y=Bool.(rand(rng,d, n))

	# fit a model using initial values
	init_range = 0.2
	init_order = 1.5
	init_lambda = fill(0.4, nlocs)
	init_d = SB(init_range, 1.0, init_order, init_lambda, my_distance)

	# using full likelihood is a bad idea.
	# @timed sol = fit_mle(init_d,y; m = 100*length(init_d), return_sol = false)
	# #maximise without the order !
	#  @timed sol = fit_mle(init_d,y; m = 10*length(init_d), return_sol = false, order = my_order)

	# pairwise maximisation is a lot better
	tdist = maximum(my_distance) / 1
	wp = 1.0 .* (my_distance .< tdist)
	 
    @testset "n= $n, D=$(length(d)), fit_mle, high m, check Bernoulli proba , then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol1 = fit_mle(init_d, y, wp; order = my_order, m = 100 * 2, return_sol = true, maxiters = 20000)
        @test isapprox(my_λ, sol1[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol1[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol1[1].λ; rtol=0.02)
        @test isapprox(my_range, sol1[1].range; rtol=0.1)
        @test isapprox(my_range, sol1[1].range; rtol=0.05)
        @test_broken isapprox(my_range, sol1[1].range; rtol=0.02)
    end 
   


	
    @testset "n= $n, D=$(length(d)),fit_mle, low m, check Bernoulli proba, then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol2 = fit_mle(init_d, y, wp; order = my_order, m = 30 * 2, return_sol = true, maxiters = 20000)
        @test isapprox(my_λ, sol2[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol2[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol2[1].λ; rtol=0.02)
        @test isapprox(my_range, sol2[1].range; rtol=0.1)
        @test isapprox(my_range, sol2[1].range; rtol=0.05)
        @test_broken isapprox(my_range, sol2[1].range; rtol=0.02)
    end 
   


    @testset "n= $n,  D=$(length(d)),fit_mle_vfast, check Bernoulli proba, then range (rtol = 0.1, then 0.05, then 0.02)" begin
        sol3 = fit_mle_vfast(init_d, y, wp; order = my_order, return_sol = true, maxiters = 20000)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.1)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.05)
        @test isapprox(my_λ, sol3[1].λ; rtol=0.02)
        @test isapprox(my_range, sol3[1].range; rtol=0.1)
        @test isapprox(my_range, sol3[1].range; rtol=0.05)
        @test_broken isapprox(my_range, sol3[1].range; rtol=0.02)
    end 
   
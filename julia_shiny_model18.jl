#making a model in julia 1.6.1
#
try
	#using these packages
	using Pkg, StatsBase, Compat, Random, Distributions, Plots, LinearAlgebra
	using Combinatorics, Random, LinearAlgebra, Base.Iterators, SHA
	using JLD2, CSV, RData, DataFrames, StatsPlots, Distributed, FastaIO
	using Measurements, LightXML, Statistics, PyPlot, ArgParse
	using Pkg
catch
	#add packages that are needed
	Pkg.add.(["RCall", "RData", "StatsBase", "Compat", "Random", "Distributions", "Plots", "LinearAlgebra"])
	Pkg.add.(["Random", "DataFrames", "CSV", "Combinatorics", "SHA", "JLD2", "StatsPlots", "Distributed", "FastaIO"])
	Pkg.add.(["Measurements", "LightXML", "Statistics", "PyPlot", "ArgParse"])
	using Pkg, StatsBase, Compat, Random, Distributions, Plots, LinearAlgebra
	using Combinatorics, Random, LinearAlgebra, Base.Iterators, SHA
	using JLD2, CSV, RData, DataFrames, StatsPlots, Distributed, FastaIO
	using Measurements, LightXML, Statistics, PyPlot, ArgParse

end

#########################################################################################
#Function list
#########################################################################################

#############################
#initialise function
#############################
popi = 100
popsi = 4
migi = 0.0
conv = 0.8
indiv = 5
evni = 0
length_i = 10
refresh = function(;popi = 100, popsi = 4, migi = 0.0, conv = 0.8, indiv = 5, evni = 0, length_i = 10, max_popi = 10000)

	#number of populations considered
	global pops
	pops = popsi

	#allele length
	global allele
	allele = length_i

	#starting alleles
	global init_allele
	init_allele = indiv

	#sequence of random haplotypes
	global pop1
	pop1 = convert(Array{Int8, 2}, rand(0:3, (allele, init_allele)));

	#frequency of haplotypes
	global freq1
	freq1 = round.(rand(init_allele, 1, pops), digits = 6);
		#normalise over 1
		freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

		#have to define this outside the function for some reason
	# #global frequency vector
	# global all1
	# all1 = deepcopy(freq1);

	#define genomic type
	global gene_type
	gene_type = repeat(["genomic"], size(pop1, 1));

	#define chromosome
	global chrm
	chrm = ones(size(pop1, 1));

	#define position on chromosome
	global pos
	pos = collect(1:size(pop1,2));

	#define population
	global pop
	pop = ones(size(pop1,2));

	#population size
	global popsize1
	popsize1 = reshape(repeat(popi:popi, pops), (1,1,pops));

	#maximum popsize
	global max_pop
	max_pop = max_popi;

	#clutch size, affects maximum reproductive capacity of individual
	global clutch
	clutch = 100;

	#conversion effeciency of drive
	global conversion
	conversion = conv

	#define substitution matrix
	global sub_matrix
	#A = 0, T = 1, G = 2, C = 3, gap = 4
	sub_matrix = Dict(0 => [0.25, 0.25, 0.25, 0.25, 0], 1 => [0.25, 0.25, 0.25, 0.25, 0],
	 2 => [0.25, 0.25, 0.25, 0.25, 0], 3 => [0.25, 0.25, 0.25, 0.25, 0], 4 => [0, 0, 0, 0, 1]);

	#sex bias dictionary
	global sex_dict
	sex_dict = Dict(0 => 0.5, 1 => conversion, 2 => 1)


	#the mutation rate per generation
	global mu_rate
	mu_rate = 1e-5;

	#storage for all freqs in simulation
	global freq_all
	freq_all = freq1;

	#selection pressure 0 <= pressure <= 1
	global pressure
	pressure = zeros(1, 1, pops)

	#resistant haplotype that all others are compared to
	global res_haplo
	res_haplo = convert(Array{Int8, 2}, rand(0:3, (size(pop1, 1),1)))

	# #defines fitness of each haplotype against the reference
	# global haplo_fitness
	# haplo_fitness = mapslices(haplo_fit, pop1, dims = 1)

	#track total popsize
	global total_pop1
	total_pop1 = deepcopy(popsize1)

	#which haplotypes are drives
	global drivers1
	drivers1 = zeros(1, size(pop1, 2))
	drivers1[1] = 1

	#sex_bias
	global sexes
	sexes = fill(0.5, 1, 1, pops);

	#sex bias
	global sex_bias
	sex_bias = fill(0.5, 1, 1, pops);

	global total_sex
	total_sex = fill(0.5, 1, 1, pops);

	global trait_distro
	trait_distro = rand(Beta(2,100), size(pop1, 1))

	#multi fitness traits
	# global traits1
	# # traits1 = [rand(Normal(0.5, 0.1), 1, init_allele); zeros(2, init_allele)] ;
	# traits1 = [seq_traits(pop1, trait_distro); zeros(2, size(pop1, 2))]

	#migration matrix
	global mig_rates
	mig_rates = fill(migi, pops, pops)
	mig_rates[diagind(mig_rates)] .= 0

	#migration input
	global mig_input
	mig_input = migi

	#environmental factor
	global env
	env = zeros(1, 1, pops)

	# global mig_ref
	# mig_ref = (mig_a = mig_a, mig_b = mig_b, mig_c = mig_c, mig_d = mig_d, mig_e = mig_e)

	global exposure
	exposure = fill(1.0, 1, 1, pops)

	# global m_conver
	# m_conver = m_conver_make(pop1, drivers1)

	global ud_tox
	ud_tox = zeros(indiv, indiv)

	global Fis
	Fis = zeros(1,1,pops)

	global gamma
	gamma = ones(size(pop1, 1))

	global crash
	crash = false

	global dominance
	dominance = fill(1, 3)

end

##########################
#function end
##########################


#############################
#drift function - new frequecies with binomial
#############################

#define binomial distribution for change in frequency
#more complexity will be a unique distro for each haplotype give:
#popsize, clutch size

new_freq = function(x, popu)
	#need to for loop over popu
	#success of binomial should be eqaul to frequency
	#need to dot this function
	success = x

	distro = Binomial(2 * popsize1[:, :, popu][1], success)

	freq = rand(distro, 1)[1] / (2 * popsize1[:, :, popu][1])

end


##########################
#function end
##########################



#############################
#mutation function
#############################
mutate = function(haplo, position, sub_matrix, popu)

#haplo = the haplotypes that will mutate
#position = the position within the haplotype that will mutate
#sub_matrix = dictionary with mutation rates per base in order ATGC

#reset mu_haplo
	mu_haplo = Array{Int8}(undef, size(pop1, 1), 0)

	#input a haplotype array and a position array to mutate a base
	for i in 1:size(haplo, 2)

		#mutate single base in one haplo of input
		new_base = sample(0:3,
		Weights(sub_matrix[haplo[:,i][position[i]]]), 1)

		#check that the mutation is new, if so, add it
		if new_base != haplo[:,i][position[i]]
			temp_haplo = copy(haplo[:,i])
			temp_haplo[position[i]] = new_base[1]
			mu_haplo = hcat(mu_haplo, temp_haplo)
		end
	end

	return (mu_haplo)
end

##########################
#function end
##########################


##########################
#use mutate function
##########################

do_mutate = function(population, popsize, mu_rate)

	#number of mutations
	mu_number = trunc(Int, ceil(size(population, 1) * popsize[:, :, popu] * mu_rate))

	#column and row to be mutated in pop1
	mu_col = sample(1:size(population, 2), mu_number)
	mu_row = sample(1:size(population, 1), mu_number)
	#find haplotypes that mutated, copy them with mutated base
	#and add a new frequency of 1/2N

	#find the column where mutation occured, div divides without remainder
	#add copy of mutated haplotype to population array
	#get columns to mutate
	new_haplo = population[:, mu_col, popu]
	sites = mu_row

	#mutate the desired bases and return an array of new haplos
	mu_haplo = mutate(new_haplo, sites, sub_matrix)

end


#improved mutation module

new_mute = function(pop1, popsize1, mu_rate, sub_matrix, drivers1, gamma)

	#new alleles
	# global new_haplos
	new_haplos = Array{Int8}(undef, size(pop1,1), 0)
	new_freqs = Array{Float64}(undef, size(popsize1, 3), 0)
	new_drivers = Array{Float64}(undef, 1, 0)
	new_traits = Array{Float64}(undef, 3, 0)
	#mutations per population
	for popu in 1:pops

		#number of mutations in population
		mu_number = rand(Binomial(popsize1[popu], mu_rate), 1)[1]

		# #location of mutations - doesnt account for allele freqs
		# mu_locations = sample(1:prod(size(pop1)), mu_number)

		#alleles where mutations occur
		mu_alleles = sample(1:size(pop1, 2), Weights(vec(freq1[:,:,popu])), mu_number)

		#position of mutations along allele
		mu_position = sample(1:size(pop1, 1), Weights(gamma), mu_number)

		#combine position and allele
		location = hcat.(mu_position, mu_alleles)

		for i in location

			#copy the column allele and mutate
			new_allele = copy(pop1[:,i[2]])
			old_base = new_allele[i[1]]
			new_base = sample(0:4, Weights(sub_matrix[old_base]), 1)[1]

			if old_base != new_base
				new_allele[i[1]] = new_base

				#global new_haplos, new_freqs, new_drivers #needed if working inside function
				new_haplos = hcat(new_haplos, new_allele)

				#add new frequency but only to the right population
				temp_freq = zeros(size(popsize1,3))
				temp_freq[popu] = 1 / (2 * popsize1[popu])
				new_freqs = hcat(new_freqs, temp_freq)

				#add new haplos to drivers1 array
				new_drivers = hcat(new_drivers, drivers1[i[2]])

				# #add new trait values deviating from original
				# temp_traits = traits1[:, i[2]] .+ [rand(Normal(0, 0.05), 1); 0; 0]
				# new_traits = hcat(new_traits, temp_traits)

				#add new trait values based on sequence
				temp_traits = seq_traits(new_allele, trait_distro)[1]
				former_traits = copy(traits1[:, i[2]])
				former_traits[1] = temp_traits
				new_traits = hcat(new_traits, former_traits)

			end

		end

	end
	return(new_haplos, reshape(transpose(new_freqs), :, 1, size(popsize1,3)), new_drivers, new_traits)
end


#generate parameters for new sequences

##########################
#generate new trait values - same for all pops
##########################

new_trait = function(new_haplos, popsize1)

	#generate random trait values for new haplos
	new_traits = rand(Normal(.035, .01), 3, size(new_haplos, 2))
	#new_traits = fill(0.05, 3, size(new_haplos, 2), size(popsize1, 3))
end

seq_traits = function(pop1, trait_distro)

traits_1 = []

	for i in 1:size(pop1, 2)
		len = size(pop1, 1)
		push!(traits_1, sum((pop1[1:end,i] .== zeros(len, 1)) .* trait_distro))
	end
	return(transpose(traits_1))
end






##########################
#function end
##########################


##########################
#concatenate frequencies with preceding zeros
##########################

#old function
function cat_freq(freq_new, freq_old)

	#find number of new haplos compared to previous
	no_haplos = length(freq_new) - size(freq_old, 1)

	#fill an array of zeros for the number of generations and number of new haplos
	gap_fill = fill(0.0, no_haplos, size(freq_old, 2))

	#add to to the bottom of old frequency
	freq_old = vcat(freq_old, gap_fill)

	#cat the two freqs together
	freq_old = hcat(freq_old, freq_new)

end

#new function combine all1 and freq1

combine_all1 = function(all1, freq1, pops)

	#find how many new alleles there are
	newbies = size(freq1, 1) - size(all1, 1)

	#add zeros to all1 to accomodate new alleles
	all1_new = cat(dims = 1, all1, repeat(zeros(newbies, size(all1, 2)), outer = [1,1,pops]))

	all1_new = cat(dims = 2, all1_new, freq1)

	all1_new
end


##########################
#function end
##########################


##########################
#fitness module
##########################

haplo_fit = function(population)

	#input is population haplotypes, output is fitness of haplotypes
	#define fitness for each haplotype = likihood to survive
	#currently ratio of similarity to res_haplo is fitness
	#need to mapslices function

	fitness = sum(res_haplo .== population) / length(res_haplo)

end

##########################
#end
##########################

##########################
#selection diploids
##########################

#fitness of each pair of haploids
#selection is additive of each haplotype
selection_dips_add = function(fitness)

	#length of fitness
	len = length(fitness)

	#repeated fitness
	fit = repeat(fitness, len)

	#sum haplo fits for each combination of haplos
	dip_fit = fit + rotr90(fit)

	#fitness > 1 should be returned to 1
	dip_fit[dip_fit .> 1] .= 1

	return(LowerTriangular(dip_fit))

end


#fitness of each pair of haploids
#selection is dominant of each haplotype
selection_dips_dom = function(fitness)

	#length of fitness
	len = length(fitness)

	#repeated fitness
	fit = repeat(fitness, len)

	#cat the two arrays then take the max of each position
	temp = cat(fit, rotr90(fit), dims = 3)
	dip_fit = reshape(maximum(temp, dims = 3), (len, len))

	#fitness > 1 should be returned to 1
	dip_fit[dip_fit .> 1] .= 1

	return(LowerTriangular(dip_fit))

end

##########################
#multi-trait selection
##########################

# #how many alleles
# len = length(freq1)
#
# #index position of each diploid
# dips = collect.(collect(product(1:len, 1:len)))

#select each diploid and calculate fitness for many traits
#index the traits by the dips index and sum the columns for now
	trait_calc_add = function(dips, traits1, popu, dominance)
		#need to dot this function
		trait_score = Array{Float64}(undef, size(dips))
		for i in vec(dips)
			#base fitness
			base = mean(traits1[1, i])

			#fitness cost homos trait2
			if traits1[2, i][1] == traits1[2, i][2]
				cost1 = traits1[2, i][1]
			end

			#fitness cost heteros trait2
			if traits1[2, i][1] != traits1[2, i][2]
				cost1 = minimum(traits1[2, i]) * dominance[2]
			end

			#fitness cost homos trait3
			if traits1[3, i][1] == traits1[2, i][2]
				cost2 = traits1[3, i][1]
			end

			#fitness cost heteros trait3
			if traits1[3, i][1] != traits1[2, i][2]
				cost2 = minimum(traits1[3, i]) * dominance[3]
			end


			trait_score[i[1], i[2]] = sum([base, cost1, cost2])

			# trait_score[i[1], i[2]] = sum(traits1[:, i]) + (randn(1) * env[:,:,popu][1])[1]
		end
	trait_score
	end

	trait_calc_dom = function(dips, traits1, popu)
		#need to dot this function
		trait_score = Array{Float64}(undef, size(dips))
		for i in vec(dips)
		trait_score[i[1], i[2]] = maximum(sum(traits1[:, i], dims = 1)) + (randn(1) * env[:,:,popu][1])[1]
		end
	trait_score
	end



##########################
#end
##########################


#sample how many survive after selection
survival = function(dip_freq, dip_fit, popu, exposure)

	#survival(frequency of diploids, fitness of diploids)
	#sample how many expect to survive due to selection
	#binomial with trials = number of diploid, success = fitness
	#need to dot this function
	#survival is now how many die
	trials = trunc(Int, floor(exposure * dip_freq * popsize1[:,:,popu][1]))
	success = (1 - (dip_fit)) * pressure[:,:,popu][1]
	distro = Binomial(trials, success)
	surviving = rand(distro, 1)[1] / (popsize1[:,:,popu][1])

end


##########################
#make haploids from diploids
##########################

dip2hap = function(diploids)
	#how to return from diploid freqs to haploid freqs
	hap_freqs = []
	for i in 1:size(diploids, 1)

		#sum along the rows and columns and subtract the one counted twice
		new_freq = (sum(diploids[i,:]) / 2) + (sum(diploids[:,i]) / 2)
		 			+ diploids[i,i]

		if new_freq < 1e-20
			new_freq = 0
		end

		push!(hap_freqs, new_freq)

	end

	return(hap_freqs)

end


##########################
#function end
##########################


##########################
#find drivers
##########################

	driver_split = function(drivers)

		#separate each 1 in drivers into a unique array / new row
		len = length(drivers)
		base = zeros(Int8, len, len)

		#add the one for each driver in a new row
		split_drive = base .+ Diagonal(vec(drivers))

		#only keep the rows where drivers is 1
		split_drive = split_drive[Bool.(vec(drivers)), :]

	end



	driver_inter = function(drivers)
	#number of haplotypes

	drivers .+ transpose(drivers)

	end

##########################
#function end
##########################


##########################
#gene drive conversion - change frequencies
##########################

	converter = function(diploids, het_drivers)

		#frequency of hets should decrease and freq of
		#homo drive should increase

		heteros = map(x -> 1 - (x - 1)^2, het_drivers)
		homos = map(x -> x * (x - 1) * (x - 1.5), het_drivers)

		#find heterozygotes and record the change in their freqs
		change = diploids .* heteros .* conversion

		#add change to freq of homo drive
		diploids = diploids + (homos .* sum(change)) - change


		return(diploids)

	end



	converter1 = function(diploids, split_drive, m_conver, popu, drivers1)

		#get hetero/homo for gene drive
		het_drivers = LowerTriangular(split_drive .+ transpose(split_drive))
		heteros = map(x -> 1 - (x - 1)^2, het_drivers)

		position = findall(split_drive .==1)[1]
		drive_number = sum(drivers1[:, 1:position])
		conversions = transpose(m_conver[trunc(Int64, drive_number), :])

		#find opposite of drivers1 so can make
		#it that drivers cant convert drives
		not_drivers = (map(x -> -x + 1, drivers1))

		#make it relative conversion????
		conversions = conversions ./ maximum(conversions)
		conversions = conversions .* not_drivers #.* 2.5

		#need to remove the position of the homozygote gene drive
		homo_drive = findall(split_drive .== 0)
		conversions = conversions[:, homo_drive]

		targets = findall(heteros .== 1)
		convertees = diploids[targets]
		convert_rates = conversions
		conv_change = zeros(length(convertees))

		for i in 1:length(convertees)

			trials = trunc(Int64,  popsize1[popu])
			success = convertees[i] * convert_rates[i]
			distro = Binomial(trials, success)
			conv_change[i] = rand(distro, 1)[1] / popsize1[popu]
		end

		heteros[findall(heteros .== 1)] = heteros[findall(heteros .== 1)] .* conv_change
		return(heteros)
	end


converter2 = function(diploids, split_drive, m_conver, popu, drivers1)

	#get hetero/homo for gene drive
	het_drivers = LowerTriangular(split_drive .+ transpose(split_drive))
	heteros = map(x -> 1 - (x - 1)^2, het_drivers)
	targets = het_drivers .> 0

	#position of gene drive allele
	homo_pos = findall(split_drive .== 1)[1]
	drive_num = trunc(Int64, sum(drivers1[1:homo_pos]))

	#opposite of drivers1
	wilds = (-1 .* drivers1) .+ 1

	#conversion rate per allele
	#removing chance to convert other drives
	conversion = transpose(m_conver[drive_num, :]) .* wilds

	convertees = diploids .* targets
	convertees1 = convertees[targets]

	#conversion process
	conv_change = zeros(length(convertees1))
	for i in 1:length(convertees1)

		trials = trunc(Int64,  popsize1[popu] )
		success = conversion[i] * convertees1[i]
		distro = Binomial(trials, success)
		conv_change[i] = rand(distro, 1)[1] / popsize1[popu]
	end

	convertees[targets] .= conv_change

	return(convertees)

end






m_conver_make = function(pop1, drivers1, priming, g_start, g_end, t_start, t_end)

	#make empty array
	m_conver = Array{Float64}(undef, sum(drivers1 .== 1), size(pop1, 2))
	len = length(pop1[g_start:g_end,1])

	drivers = findall(drivers1 .== 1)
	for i in 1:sum(drivers1 .== 1)

		#only pick the drive alleles
		this_one = drivers[i][2]

		m_conver[i, :] = sum(pop1[g_start:g_end, this_one] .== pop1[t_start:t_end,:], dims = 1) ./ len .* priming

	end
	return(m_conver)
end


m_conver_scan = function(pop1, drivers1, priming, grna_seq)

	#make empty array
	m_conver = Array{Float64}(undef, sum(drivers1 .== 1), size(pop1, 2))
	len = length(grna_seq)

	drivers = findall(drivers1 .== 1)


	for i in 1:sum(drivers1 .== 1)

		#only pick the drive alleles
		this_one = drivers[i][2]

		grna = pop1[1:len, this_one]

		#iterate over sequence
		for k in 1:size(pop1, 2)

			mean_conv = []

				for j in 1:(size(pop1, 1) - len)

					push!(mean_conv, sum(grna .== pop1[j:(j + len - 1), k]) / len * priming)

				end

			m_conver[i, k] = maximum(mean_conv)

		end



	end
	return(m_conver)
end




	under_d = function(diploids, ud_tox)

		toxins = diploids .* ud_tox

		toxins = mapslices(x -> x ./ sum(x), toxins, dims = [1,2])

		toxins
	end






##########################
#function end
##########################


##########################
#sex bias conversion - male bias
#needs work!
##########################

#expected ratio of M/F offspring

	sex_conv = function(diploids, het_drivers)

		#diploids is the frequency of diploids
		#drivers is the homo/hetero state of drives

		#number of trials in binomial
		trials = trunc(Int, diploids * popsize)

		#success rate of binomial = ratio of M/F
		success = sex_dict[het_drivers]

		#Binomial distribution to see how many males born
		distro = Binomial(trials, success)

		sex = rand(distro, 1)[1] / trials



	end

##########################
#function end
##########################


##########################
#make diploids
##########################
make_dips = function(freq1, Fis)
	#need to for loop over pops
	#make a square matrix of freqs
	#make diploids with frequency of each diploid
	diploids = freq1 * transpose(freq1)

	#multiple everything except the diagonal by 2 to compensate for lower
	heteros = map(x -> -x + 1, Diagonal(ones(length(freq1))))
	homos = abs.(heteros .-1)

	#adjust for inbreeding
	diploids_Fis = Diagonal(vec(sum(diploids .* heteros .* (Fis), dims = 2)))
	diploids_homo = diploids .* homos .+ diploids_Fis
	diploids_hets = LowerTriangular(2 .* diploids .* heteros .* (1 - Fis))

	diploids = diploids_homo .+ diploids_hets
	diploids = diploids ./ sum(diploids)

end



##########################
#function end
##########################




##########################
#migration sampling
##########################

mig_sample = function(freq1, rate)


	#migration(frequency of haploids, popsize)
	#sample how many expect to survive due to selection
	#binomial with trials = number of diploid, success = fitness
	#need to dot this function
	trials = trunc(Int, round(from * fpop))
	success = rate
	distro = Binomial(trials, success)
	migrating = rand(distro, 1)[1] / fpop

end


##########################
#function end
##########################



##########################
#death rate per group
##########################
#returns the number of survivors

death = function(d_rate, pops)

    #function to see how many of each group die based on the
    #death rate ~ popsize ~ carrying capacity
    #need to dot this function

    #distribution for chance of survival (1 - death)
    ded = Binomial(pops, (1 - d_rate))

    survivors = rand(ded, 1)[1]

end


##########################
#gene drive total function
##########################

gene_drive_post = function(freq1, pop1, popsize1, drivers1, traits1, m_conver, Fis, sexes)

    # global freq1
    # global all
    # global pop1
    # global haplo_fitness
    # global total_pop
    # global popsize
    # global sexes
    # global drivers1
    # global traits
	# global m_conver

###################################################################
#mutation
###################################################################

if mutation == true

	#mutate sequences to make new alleles
	new_haplos, new_freqs, new_drivers, new_traits = new_mute(pop1, popsize1, mu_rate, sub_matrix, drivers1, gamma)

	#add sequences to pop1
	pop1 = hcat(pop1, new_haplos)

	#add new frequencies to freq1
	freq1 = cat(dims = 1, freq1, new_freqs)
	freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

	#adjust the drivers1 array
	drivers1 = hcat(drivers1, new_drivers)

    #add new traits to array - multi-trait fitness
    #traits1 = cat(dims = 2, traits1, new_trait(new_haplos, popsize1))
	# traits1 = hcat(traits1, new_traits)
	temp = seq_traits(pop1, trait_distro)
	traits1 = hcat(traits1, new_traits)
	traits1[1,:] = temp ./ maximum(temp) ./ 2


	# # #add new conversion efficiency (by similarity)
	# m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)

	#scanning m_conver
	m_conver = m_conver_scan(pop1, drivers1, priming, grna)
end



###################################################################
#diploids/selection
###################################################################

    #make diploids with frequency of each diploid, sums to 1
    diploids = zeros(size(freq1,1), size(freq1,1), pops)
    for popu in 1:pops
        diploids[:,:,popu] = make_dips(freq1[:,:,popu], Fis[:,:,popu][1])
    end
#######################################################################
#######################################################################

###################################################################
#gene drive action - conversion
###################################################################
# drivers1 = zeros(size(freq1,1))
# drivers1[1] = 1
    #split the drivers array so each gene drive is solo
    #then iterate over each driver array
    #add the one for each driver in a new row
    split_drive = driver_split(drivers1)


#repeat for each row of split_drive the conversion process
for popu in 1:pops
	for i in randperm(size(split_drive,1))

        #convert het drives to homo drives
        diplo_change = LowerTriangular(converter2(diploids[:,:,popu], split_drive[i,:], m_conver, popu, drivers1))

		diploids[:,:,popu] = diploids[:,:,popu] .- diplo_change
        diploids[:,:,popu] = diploids[:,:,popu] .+ Diagonal(vec(sum(diplo_change) .* split_drive[i,:]))

		#make sure no values below 0 or above 1
		diploids[:,:,popu] = abs.((diploids[:,:,popu] .> 0) .* diploids[:,:,popu])
		diploids[:,:,popu] = diploids[:,:,popu] ./ sum(diploids[:,:,popu])

    end
end

###################################################################
#selection
###################################################################

    # #matrix of fitness for each diploid
    # fitness = selection_dips_add(haplo_fitness)

    #how many alleles
    len = length(freq1[:,:,1])

    #index position of each diploid
    dips = collect.(collect(product(1:len, 1:len)))

    #matrix of diploid fitness for multi-traits
	fitness = zeros(len, len, pops)
	for popu in 1:pops
		fitness[:,:,popu] = LowerTriangular(trait_calc_add(dips, traits1, popu, dominance))
	end

	fitness = fitness ./ maximum(fitness)


    # #set max fitness to 1 and min to 0
    # fitness[fitness .> 1] .= 1
	# fitness[fitness .< 0] .= 0

    #survival - see how many survive by fitness
	#survival is now how many die
    survivors = zeros(size(freq1,1), size(freq1,1), pops)
    for popu in 1:pops
        survivors[:,:,popu] = survival.(diploids[:,:,popu], fitness[:,:,popu], popu, exposure[:,:,popu])
    end

	survivors = diploids .- survivors
	survivors[survivors .< 0] .= 0


    #modify population size by how many died
	for popu in 1:pops
		popsize1[popu] = trunc.(Int64, popsize1[:,:,popu][1] * sum(survivors[:,:,popu]))
	end

	#issues arise if popsize = 0, make it equal to 1?
	popsize1[popsize1 .<= 0] .= 1

    #diploids frequency is survivors
    diploids = survivors

    #return diploids to sum = 1?
    diploids = mapslices(x -> x ./ sum(x), diploids, dims = [1,2])


###################################################################
#gene drive action - conversion sex bias
###################################################################

if crash == true
	for popu in 1:pops
	    #ratio of sexes
	    ratios = LowerTriangular(fill(0.5, size(diploids[:,:,popu])))

	    #change driver to be 1 (all males)
		het_drivers = drivers1 .+ transpose(drivers1)
		ratios[Diagonal(het_drivers .== 2)] .= 1

	    #convert sex ratio to average sex bias
	    sex_bias[:,:,popu] .= sum(diploids[:,:,popu] .* ratios)


	end

	#add sex bias to plotting vector
	sexes = vcat(sexes, sex_bias)
end


###################################################################
#Return to haploids for next generation
###################################################################

    #return the survivor freqs to haploid freqs
    for popu in 1:pops
        freq1[:,:,popu] = dip2hap(diploids[:,:,popu])
    end


    # #return new_freq to sum of 1
    # freq1 = freq1 ./ sum(freq1)

    #increase in population size by clutch and ratio of females
    for popu in 1:pops
        popsize1[:,:,popu] = [trunc(Int, popsize1[:,:,popu][1] * 2 *
                                        (1 - sex_bias[:,:,popu][1]) * clutch)]
    end


    #limit popsize to some maximum
    popsize1[popsize1 .> max_pop] .= max_pop

	#script crashes if pop = 0
	popsize1[popsize1 .<= 0] .= 1
    # #record the new popsize for plotting
    # total_pop = vcat(total_pop, popsize)


###################################################################
#Drift
###################################################################
    #drift
    for popu in 1:pops
        freq1[:,:,popu] = new_freq.(freq1[:,:,popu], popu)
		if sum(freq1[:,:,popu]) == 0
			freq1[1,:,popu] .= 1
		end
    end


    #sum freq to 1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])



###################################################################
#adjust arrays so that they dont include freq == 0
###################################################################

    # #find all haplotypes with freq != 0 which are retained
    # keep = freq1 .!= 0
    #
    # #only keep freqs with keep
    # freq1 = freq1[keep]
    #
    # #only keep the haplotypes with keep
    # pop1 = pop1[:,keep]
    #
    # # #only keep haplo fitness with keep
    # # haplo_fitness = haplo_fitness[:,keep]
    #
    # #only keep traits with keep
    # traits = traits[:,keep]
    #
    # #only keep drivers with keep
    # drivers = drivers[:,keep]
    #
    # #remove haplotypes that died in all
    # all = all[keep,:]


return(freq1, pop1, popsize1, drivers1, traits1, m_conver, sexes)



end


###################################################################
###################################################################
###################################################################
#running model as a single function
###################################################################
###################################################################
###################################################################
















#testing page #works on model 18
#################################################################
#set-up parameters
#################################################################

run_model_post = function(; pop_input = 100000, pops_input = 1, indiv_input = 3, gens_input = 100, fis_input = 0.0, migration_input = false, exposure1 = 1.0, exposure2 = 1.0,
                       dominance_input = 0.0, fitness_input = 0.0, conversion_input = 1.0, resistance_input = 0.0, drive_freq = 0.001, resfreq_input = 0.1,
                       migrates_input = 1e-3, iterations_input = 1)

pyplot()
ps = repeat([Plots.plot(1, label = "")], pops_input)

global h1, ps
global freq1, pop1, all1, total_pop1, total_traits, sexes, sex_bias, pressure, dominance, mutation
global popsize1, drivers1, traits1, m_conver, total_conv, total_sex, mig_rates, exposure, pops


for iter in 1:iterations_input

    #initialise some key variables
    #popsi is number of populations
    refresh(popi = pop_input, popsi = pops_input, indiv = indiv_input, length_i = 50, max_popi = pop_input)

    #################################################################
    #population parameters
    #################################################################

    #number of generations to cycle over
    gens = gens_input

    # #carrying capacity of area
    # max_pop = 100

    #inbreeding coefficient
    Fis = fill(fis_input, 1, 1, pops)

    #number of offspring per mating pair
    clutch = 10

    #################################################################
    #variables that change code base
    #################################################################
    migration = migration_input
    dynamic = false
    mutation = false
    cleaning = false
    crash = false
    plotting = false
    histograms = false #plot conversion efficiency distributions
    sequence = "random" #"fasta" or "random"


    #################################################################
    #Selection variables
    #################################################################
    #whether or not selection pressure is present
    #in each population
    pressure = fill(1, 1, 1, pops) #1 is on 0 is off
    # pressure = reshape([1, 1], 1, 1, pops) #1 is on 0 is off

    #The percentage of a population exposed to selection pressure
    exposure = fill(1.0, 1, 1, pops)
    exposure[1] = exposure1

    if pops > 1
        exposure[2] = exposure2
    end


    #dominance of fitness cost
    #second value is dominance of fitness cost
    dominance = fill(1.0, 3) #1 = dom, .5 = add, 0 = recess #1 is general fitness
    dominance[2] = dominance_input #fitness cost of drive
    dominance[3] = 1 #fitness cost of resistance

    #random environmental affect on fitness
    env .= 0.0

    #fitness values based of sequence, arbitary full fitness sequence
    #relative fitness
    fitness = "set" #"seq" "relative"
    traits1 = [seq_traits(pop1, trait_distro); zeros(2, size(pop1, 2))]
    traits1[1,:] .= (traits1[1,:] ./ maximum(traits1[1,:]))

    #fitness cost of gene drive
    traits1[2,1] = -(fitness_input)
    traits1[2,2] = -0.0 #fitness cost of resistance to drive

    if fitness == "set"
        traits1[1,:] .= 1
    end

    if fitness == "relative" #untested#
        #trade of between resistance and fitness
        #as resistance increases, fitness decreases
        traits1[3,:] = (m_conver .- 1) ./ 5
        traits1[3,1] = 0
    end



    #################################################################
    #Sequence variables
    #################################################################
    #mutation rate of sequences
    mu_rate = 0

    #distribution of mutations
    gamma .= 1

    #################################################################
    #Conversion variables
    #################################################################

    conversion_distribution = "set" #"set" "bipolar" "beta" "normal" "sequence"

    if conversion_distribution == "bipolar"
        percent_peak1 = 0.5
        peak1 = 0.7
        se1 = 0.2
        peak2 = 0.2
        se2 = 0.2
        maxi = 1
        mini = 0.05
    end

    if conversion_distribution == "set"
        conv_set = conversion_input
        res = resistance_input
    end

    if conversion_distribution == "beta"
        shape1 = 3
        shape2 = 0.6
    end

    if conversion_distribution == "normal"
        nmean = 0.8
        nstd = 0.1
        maxi = 1
        mini = 0.01
    end

    if conversion_distribution == "sequence"
        grna = "ATGC"
        grna = split(grna, "")
        #affect of gRNA interaction
        priming = 1

    end


    if conversion_distribution == "bipolar"
        #bimodal conversion distribution
        #binomial to say where each is sampled from
        binom1 = rand(Binomial(1, percent_peak1), size(freq1,1))'
        binom2 = abs.(binom1 .- 1)

        norma = Normal(peak1, se1)
        m_conver1 = rand(norma, 1, size(pop1, 2)) .* binom1

        normb = Normal(peak2, se2)
        m_conver2 = rand(normb, 1, size(pop1, 2)) .* binom2

        m_conver = m_conver1 .+ m_conver2

        m_conver[m_conver .< 0] .= mini
        m_conver[m_conver .> 1] .= maxi
    end

    if conversion_distribution == "set"
        m_conver = fill(conv_set, 1, size(pop1, 2))
        m_conver[2] = res
    end

    if conversion_distribution == "beta"
        # artificial m_conver values
        conv_beta = Beta(shape1, shape2)
        m_conver = rand(conv_beta, 1, size(pop1, 2))
    end

    if conversion_distribution == "normal"
        # artificial m_conver values
        conv_norm = Normal(nmean, nstd)
        m_conver = rand(conv_norm, 1, size(pop1, 2))
        m_conver[m_conver .< 0] .= mini
        m_conver[m_conver .> 1] .= maxi
    end

    if conversion_distribution == "sequence"
        m_conver = m_conver_scan(pop1, drivers1, priming, grna)
    end

    if histograms == true
        histogram(m_conver', bins = 50, label = "", title = "initial distribution", xaxis = "conversion efficiency", yaxis = "frequency")
    end


    #################################################################
    #Population variation
    #################################################################

    #only works for 2 populations
    pop_variation = "none" #"subtle" "none"

    #the conversion efficiency at which the populations are separated
    cutoff = 0.5

    #the severity of the difference in populations
    magnitude = 10

    if pops == 2

        if pop_variation == "strict"
            freq1[:,:,1] .= freq1[:,:,1] .* transpose(m_conver .>  cutoff)
            freq1[:,:,2] .= freq1[:,:,2] .* transpose(m_conver .<= cutoff)
        end

        if pop_variation == "subtle"
            adjust1 = ones(size(freq1,1), 1)
            adjust1[m_conver' .>= cutoff] .*= magnitude
            adjust2 = ones(size(freq1,1), 1)
            adjust2[m_conver' .< cutoff] .*= magnitude

            freq1[:,:,1] .*= adjust1
            freq1[:,:,2] .*= adjust2
        end

        if pop_variation == "none"
            #do nothing
        end

    end

    if histograms == true
        ps = repeat([Plots.plot(1)], pops)
        ps[1] = scatter(m_conver',freq1[:,:,1], label = "", title = "Population 1", xaxis = "conversion efficiency", yaxis = "frequency")
        ps[2] = scatter(m_conver',freq1[:,:,2], label = "", title = "Population 2", xaxis = "conversion efficiency", yaxis = "frequency")
        plot(ps...)
    end


    #################################################################
    #Gene Drive starting frequency
    #################################################################

    #starting frequency of gene drive
    starting_freq = drive_freq
    res_freq = resfreq_input
    # val = sum(freq1[2:end,1,1]) + starting_freq
    # freq1[1,1,1] = val * starting_freq
    # freq1 ./= val

    # if pops > 1
    #     freq1[1,:,:] = [0.001, 0.0] #need array to be same size as pops
    # end

    freq1[1] = starting_freq
    freq1[2] = res_freq
    freq1[3] = 1 - freq1[1] - freq1[2]

    if pops > 1
        freq1[1,1,2] = 0
        freq1[2,1,2] = .4
        freq1[3,1,2] = 1 - freq1[1,1,2] - freq1[2,1,2]
    end



    #make all pop allele freqs sum to 1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
    all1 = deepcopy(freq1)

    #################################################################
    #Migration rates
    #################################################################

    mig_style = "set" #"variable"

    if mig_style == "set"
        mig_rates .= migrates_input
    end

    if mig_style == "variable"
        mig_rates = reshape([0, .1, .1, 0], pops, pops) #vector needs to be pops squared length
    end


    #average trait value
    total_traits = Array{Float64}(undef, 0, pops)
    total_conv = Array{Float64}(undef, 0, pops)


    all1 = deepcopy(freq1)

# #change seed for simulation but keep starting scenario
# Random.seed!(50)

    for q in 1:gens
        # global freq1, pop1, all1, total_pop1, total_traits, sexes, sex_bias, pressure, dominance
        # global popsize1, drivers1, traits1, m_conver, total_conv, total_sex, mig_rates, exposure


        #1 generation of gene drive
        cycle1 = gene_drive_post(freq1, pop1, popsize1, drivers1, traits1, m_conver, Fis, sexes)
        #store new values for next cyle
        freq1 = cycle1[1]
        pop1 = cycle1[2]
        popsize1 = cycle1[3]
        drivers1 = cycle1[4]
        traits1 = cycle1[5]
        m_conver = cycle1[6]


		if migration == true
		    #migration
		    #simulates how much of each allele travels to each pop
		    migrants = repeat([freq1[:,:,1]],pops, pops)
		    for k in 1:pops
		        for j in 1:pops

		            migrants[j,k] = freq1[:,:,j] .* mig_rates[j,k] .* popsize1[j]

		        end
		    end

		    # #adjust each pops frequencies based on migration
		    # #each row of migrants represents all the leavers
		    # #each column represents all the arrivers
		    for popu in 1:pops

		        freq1[:,:,popu] = freq1[:,:,popu] .+ sum(hcat(migrants[:,popu]...) ./ repeat(fill(popsize1[popu][1], 1, size(popsize1, pops)), outer = size(pop1, 2)), dims = 2)

		    end

		    #adjust popsizes
		    popsize1 = round.(popsize1 .+ reshape(sum(sum.(migrants), dims = 1), 1,1,pops) .- reshape(sum(sum.(migrants), dims = 2), 1,1,pops))

		    #normalise freq1
		    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
		end



	        #store new freq in all for plotting purposes
	        all1 = combine_all1(all1, freq1, pops)


	        #store total_pop
	        total_pop1 = vcat(total_pop1, popsize1)
	        total_sex = vcat(total_sex, sexes)


	        #population cleaning
	        #if allele is present for less than 10 generations,
	        #remove it from:
	        #pop1, freq1, drivers1, traits1, m_conver

	            temp_conv = []
	        for j in 1:pops
	            push!(temp_conv, sum(m_conver .* repeat(transpose(freq1[:,:,pops]), size(m_conver, 1)) / size(m_conver, 1)))
	        end
	        total_conv = vcat(total_conv, transpose(temp_conv))


	    #cleaning dead alleles out
	    if cleaning == true
	             #to remove
	             to_keep = collect(1:size(all1,1))
	             for k in 1:size(all1,1)

	                if sum(all1[k, end, :] .< 1e-6) == pops &&
	                    sum(maximum(all1[k, :, :], dims = 1) .< 0.001) == pops

	                    deleteat!(to_keep, findall(to_keep .== k))

	                end
	            end



	            #only keep the alleles specified in to_keep
	            pop1 = pop1[:, to_keep]
	            freq1 = freq1[to_keep, :, :]
	            drivers1 = drivers1[: ,to_keep]
	            traits1 = traits1[:, to_keep]
	            m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)
	            all1 = all1[to_keep, :, :]
	    end



	            temp_traits = []
	        for j in 1:pops
	            push!(temp_traits, 2 * sum(transpose(freq1[:,:,j]) .* sum(traits1, dims = 1)))
	        end
	        temp_traits[temp_traits .> 1] .= 1
	        total_traits = vcat(total_traits, transpose(temp_traits))



	        #dynamic bar plotting
	        #plot data - make drives negative values
	        #pdrives = transpose(map(x -> (-2 * x) + 1, drivers1))
	        #pfreq = sort(freq1 .* pdrives, dims = 2)
	    if dynamic == true
	        dat = transpose(freq1[:,1,:]) .* repeat(vec(popsize1), 1, size(freq1,1)) ./ max_pop
	        nam = repeat(["Bana", "Pala"], outer = size(freq1,1))
	        key = sort(unique(nam))
	        val = collect(1:length(key)) .- 0.5
	        dict = Dict(zip(key, val))
	        ctg = repeat(1:size(freq1,1), inner = size(freq1, 3))
	        p = groupedbar(nam, dat, bar_position = :stack,
	            group = ctg,
	            label = "", title = "Gene Drive Simulation \n Generation $(q)", linecolour = :match,
	            ylim = (0,1.3), xlab = "Populations", ylab = "Frequency",
	            ann = [ (dict["Bana"], 1.05, popsize1[1]),
	                    (dict["Pala"], 1.05, popsize1[2]),
	                    ])

	        display(p)
	        # savefig("frame$(q).png")
	    end


	        # #dynamic pie charts
	        # ps = repeat([Plots.plot(1)], pops)
	        #
	        # for popu in 1:pops
	        #     ps[popu] = pie(freq1[:,:,popu])
	        # end
	        #
	        # display(Plots.plot(ps...))


    # println(q)
    end






    # elapsed = round(time() - start, digits = 3)
    # println(elapsed, " seconds")
    ###################################################################
    #Plotting
    ###################################################################
    #make empty array to store plots
    # pyplot()
	# ps = repeat([Plots.plot(1)], pops)

    #make gene drives solid and others dashed
    styles = fill(:dash, length(drivers1))
    styles[vec(drivers1 .== 1)] .= :solid
    styles = reshape(styles, 1, length(styles))

    pop_names = ["Target", "Neighbour", "Neighbour", "Ghana"]

    for popu in 1:pops

        ps[popu] = Plots.plot(ps[popu], transpose(all1[:,:,popu]),
                #annotate = [(30, 0.3, text("Pressure = $(pressure[popu])", 8)),
                #            (30, 0.1, text("Env = $(env[popu])", 8))],
                legend = :topleft,
                xlab = "Generation",
                ylab = "Frequency",
                guidefontsize= 9,
                titlefontsize = 12,
                title = pop_names[popu],
                # left_margin = [-5mm 0mm],
                layout = 1,
                ylim = (0,1),
                linestyle = styles,
                width = 1.5,
				color = [1,2,3]',
                label = "")

            #     #plotting population size
            # Plots.plot!(total_pop1[:,:,popu] ./ max_pop,
            #     label = "",
            #     colour = "black",
            #     ls = :dash,
            #     width = 2)

            #     #plotting average fitness
            # Plots.plot!(total_traits[:, popu],
            # label = "",
            # colour = "black",
            # ls = :dot,
            # width = 2)

            #     #plotting average conversion
            # Plots.plot!(total_conv[:, popu],
            # label = "",
            # colour = "black",
            # ls = :dot,
            # width = 1.5)

            # Plots.vline!([200], color = "black", linewidth = 5, linestyle = :dot, label = "")
    end


	#plot all plots
	h1 = Plots.plot(ps..., xlim = (0,gens), ylim = (0,1.05), reuse = true, label = "", size = ((600 + ((pops > 1) * 600)),(400 + ((pops > 2) * 400))))
	# Plots.annotate!(50, 0.5, popsize1[1])




end



return(h1)

end




###################################################################
#Save file
###################################################################

# #save to file
# all1 = reshape(all1[1,:,:], 1, :, pops)
#
# @save "/home/student.unimelb.edu.au/bcamm/Documents/Uni/PhD/Modelling/Output/Test3/Test_run_c$(conversion)_m$(n)_t$(t)_rep$(m).jld2" all1 total_pop1
#
# end
# end
# end
# end
#
#



# #how to make unique identifiers for each sequence
# #works for very large bit strings
# mapslices(join, pop1, dims = 1)
# bytes2hex.(sha256.(mapslices(join, pop1, dims = 1)))
#







####




####
##########################
#gene drive total function
##########################

gene_drive_pre = function(freq1, pop1, popsize1, drivers1, traits1, m_conver, Fis, sexes)

    # global freq1
    # global all
    # global pop1
    # global haplo_fitness
    # global total_pop
    # global popsize
    # global sexes
    # global drivers1
    # global traits
	# global m_conver

###################################################################
#mutation
###################################################################

if mutation == true

	#mutate sequences to make new alleles
	new_haplos, new_freqs, new_drivers, new_traits = new_mute(pop1, popsize1, mu_rate, sub_matrix, drivers1, gamma)

	#add sequences to pop1
	pop1 = hcat(pop1, new_haplos)

	#add new frequencies to freq1
	freq1 = cat(dims = 1, freq1, new_freqs)
	freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

	#adjust the drivers1 array
	drivers1 = hcat(drivers1, new_drivers)

    #add new traits to array - multi-trait fitness
    #traits1 = cat(dims = 2, traits1, new_trait(new_haplos, popsize1))
	# traits1 = hcat(traits1, new_traits)
	temp = seq_traits(pop1, trait_distro)
	traits1 = hcat(traits1, new_traits)
	traits1[1,:] = temp ./ maximum(temp) ./ 2


	# # #add new conversion efficiency (by similarity)
	# m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)

	#scanning m_conver
	m_conver = m_conver_scan(pop1, drivers1, priming, grna)
end



###################################################################
#diploids/selection
###################################################################

    #make diploids with frequency of each diploid, sums to 1
    diploids = zeros(size(freq1,1), size(freq1,1), pops)
    for popu in 1:pops
        diploids[:,:,popu] = make_dips(freq1[:,:,popu], Fis[:,:,popu][1])
    end
#######################################################################
#######################################################################

###################################################################
#selection
###################################################################

    # #matrix of fitness for each diploid
    # fitness = selection_dips_add(haplo_fitness)

    #how many alleles
    len = length(freq1[:,:,1])

    #index position of each diploid
    dips = collect.(collect(product(1:len, 1:len)))

    #matrix of diploid fitness for multi-traits
	fitness = zeros(len, len, pops)
	for popu in 1:pops
		fitness[:,:,popu] = LowerTriangular(trait_calc_add(dips, traits1, popu, dominance))
	end

	fitness = fitness ./ maximum(fitness)


    # #set max fitness to 1 and min to 0
    # fitness[fitness .> 1] .= 1
	# fitness[fitness .< 0] .= 0

    #survival - see how many survive by fitness
	#survival is now how many die
    survivors = zeros(size(freq1,1), size(freq1,1), pops)
    for popu in 1:pops
        survivors[:,:,popu] = survival.(diploids[:,:,popu], fitness[:,:,popu], popu, exposure[:,:,popu])
    end

	survivors = diploids .- survivors
	survivors[survivors .< 0] .= 0


    #modify population size by how many died
	for popu in 1:pops
		popsize1[popu] = trunc.(Int64, popsize1[:,:,popu][1] * sum(survivors[:,:,popu]))
	end

	#issues arise if popsize = 0, make it equal to 1?
	popsize1[popsize1 .<= 0] .= 1

    #diploids frequency is survivors
    diploids = survivors

    #return diploids to sum = 1?
    diploids = mapslices(x -> x ./ sum(x), diploids, dims = [1,2])


###################################################################
#gene drive action - conversion
###################################################################
# drivers1 = zeros(size(freq1,1))
# drivers1[1] = 1
    #split the drivers array so each gene drive is solo
    #then iterate over each driver array
    #add the one for each driver in a new row
    split_drive = driver_split(drivers1)


#repeat for each row of split_drive the conversion process
for popu in 1:pops
	for i in randperm(size(split_drive,1))

        #convert het drives to homo drives
        diplo_change = LowerTriangular(converter2(diploids[:,:,popu], split_drive[i,:], m_conver, popu, drivers1))

		diploids[:,:,popu] = diploids[:,:,popu] .- diplo_change
        diploids[:,:,popu] = diploids[:,:,popu] .+ Diagonal(vec(sum(diplo_change) .* split_drive[i,:]))

		#make sure no values below 0 or above 1
		diploids[:,:,popu] = abs.((diploids[:,:,popu] .> 0) .* diploids[:,:,popu])
		diploids[:,:,popu] = diploids[:,:,popu] ./ sum(diploids[:,:,popu])

    end
end


###################################################################
#gene drive action - conversion sex bias
###################################################################

if crash == true
	for popu in 1:pops
	    #ratio of sexes
	    ratios = LowerTriangular(fill(0.5, size(diploids[:,:,popu])))

	    #change driver to be 1 (all males)
		het_drivers = drivers1 .+ transpose(drivers1)
		ratios[Diagonal(het_drivers .== 2)] .= 1

	    #convert sex ratio to average sex bias
	    sex_bias[:,:,popu] .= sum(diploids[:,:,popu] .* ratios)


	end

	#add sex bias to plotting vector
	sexes = vcat(sexes, sex_bias)
end


###################################################################
#Return to haploids for next generation
###################################################################

    #return the survivor freqs to haploid freqs
    for popu in 1:pops
        freq1[:,:,popu] = dip2hap(diploids[:,:,popu])
    end


    # #return new_freq to sum of 1
    # freq1 = freq1 ./ sum(freq1)

    #increase in population size by clutch and ratio of females
    for popu in 1:pops
        popsize1[:,:,popu] = [trunc(Int, popsize1[:,:,popu][1] * 2 *
                                        (1 - sex_bias[:,:,popu][1]) * clutch)]
    end


    #limit popsize to some maximum
    popsize1[popsize1 .> max_pop] .= max_pop

	#script crashes if pop = 0
	popsize1[popsize1 .<= 0] .= 1
    # #record the new popsize for plotting
    # total_pop = vcat(total_pop, popsize)


###################################################################
#Drift
###################################################################
    #drift
    for popu in 1:pops
        freq1[:,:,popu] = new_freq.(freq1[:,:,popu], popu)
		if sum(freq1[:,:,popu]) == 0
			freq1[1,:,popu] .= 1
		end
    end


    #sum freq to 1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])



###################################################################
#adjust arrays so that they dont include freq == 0
###################################################################

    # #find all haplotypes with freq != 0 which are retained
    # keep = freq1 .!= 0
    #
    # #only keep freqs with keep
    # freq1 = freq1[keep]
    #
    # #only keep the haplotypes with keep
    # pop1 = pop1[:,keep]
    #
    # # #only keep haplo fitness with keep
    # # haplo_fitness = haplo_fitness[:,keep]
    #
    # #only keep traits with keep
    # traits = traits[:,keep]
    #
    # #only keep drivers with keep
    # drivers = drivers[:,keep]
    #
    # #remove haplotypes that died in all
    # all = all[keep,:]


return(freq1, pop1, popsize1, drivers1, traits1, m_conver, sexes)



end


###################################################################
###################################################################
###################################################################
#running model as a single function
###################################################################
###################################################################
###################################################################
















#testing page #works on model 18
#################################################################
#set-up parameters
#################################################################

run_model_pre = function(; pop_input = 100000, pops_input = 1, indiv_input = 3, gens_input = 100, fis_input = 0.0, migration_input = false, exposure1 = 1.0, exposure2 = 1.0,
                       dominance_input = 0.0, fitness_input = 0.0, conversion_input = 1.0, resistance_input = 0.0, drive_freq = 0.001, resfreq_input = 0.1,
                       migrates_input = 1e-3, iterations_input = 1)

global h1
global freq1, pop1, all1, total_pop1, total_traits, sexes, sex_bias, pressure, dominance, mutation
global popsize1, drivers1, traits1, m_conver, total_conv, total_sex, mig_rates, exposure, pops

gr()
ps = repeat([Plots.plot(1, label = "")], pops_input)

for iter in 1:iterations_input

    #initialise some key variables
    #popsi is number of populations
    refresh(popi = pop_input, popsi = pops_input, indiv = indiv_input, length_i = 50, max_popi = pop_input)

    #################################################################
    #population parameters
    #################################################################

    #number of generations to cycle over
    gens = gens_input

    # #carrying capacity of area
    # max_pop = 100

    #inbreeding coefficient
    Fis = fill(fis_input, 1, 1, pops)

    #number of offspring per mating pair
    clutch = 10

    #################################################################
    #variables that change code base
    #################################################################
    migration = migration_input
    dynamic = false
    mutation = false
    cleaning = false
    crash = false
    plotting = false
    histograms = false #plot conversion efficiency distributions
    sequence = "random" #"fasta" or "random"


    #################################################################
    #Selection variables
    #################################################################
    #whether or not selection pressure is present
    #in each population
    pressure = fill(1, 1, 1, pops) #1 is on 0 is off
    # pressure = reshape([1, 1], 1, 1, pops) #1 is on 0 is off

    #The percentage of a population exposed to selection pressure
    exposure = fill(1.0, 1, 1, pops)
    exposure[1] = exposure1

    if pops > 1
        exposure[2] = exposure2
    end


    #dominance of fitness cost
    #second value is dominance of fitness cost
    dominance = fill(1.0, 3) #1 = dom, .5 = add, 0 = recess #1 is general fitness
    dominance[2] = dominance_input #fitness cost of drive
    dominance[3] = 1 #fitness cost of resistance

    #random environmental affect on fitness
    env .= 0.0

    #fitness values based of sequence, arbitary full fitness sequence
    #relative fitness
    fitness = "set" #"seq" "relative"
    traits1 = [seq_traits(pop1, trait_distro); zeros(2, size(pop1, 2))]
    traits1[1,:] .= (traits1[1,:] ./ maximum(traits1[1,:]))

    #fitness cost of gene drive
    traits1[2,1] = -(fitness_input)
    traits1[2,2] = -0.0 #fitness cost of resistance to drive

    if fitness == "set"
        traits1[1,:] .= 1
    end

    if fitness == "relative" #untested#
        #trade of between resistance and fitness
        #as resistance increases, fitness decreases
        traits1[3,:] = (m_conver .- 1) ./ 5
        traits1[3,1] = 0
    end



    #################################################################
    #Sequence variables
    #################################################################
    #mutation rate of sequences
    mu_rate = 0

    #distribution of mutations
    gamma .= 1

    #################################################################
    #Conversion variables
    #################################################################

    conversion_distribution = "set" #"set" "bipolar" "beta" "normal" "sequence"

    if conversion_distribution == "bipolar"
        percent_peak1 = 0.5
        peak1 = 0.7
        se1 = 0.2
        peak2 = 0.2
        se2 = 0.2
        maxi = 1
        mini = 0.05
    end

    if conversion_distribution == "set"
        conv_set = conversion_input
        res = resistance_input
    end

    if conversion_distribution == "beta"
        shape1 = 3
        shape2 = 0.6
    end

    if conversion_distribution == "normal"
        nmean = 0.8
        nstd = 0.1
        maxi = 1
        mini = 0.01
    end

    if conversion_distribution == "sequence"
        grna = "ATGC"
        grna = split(grna, "")
        #affect of gRNA interaction
        priming = 1

    end


    if conversion_distribution == "bipolar"
        #bimodal conversion distribution
        #binomial to say where each is sampled from
        binom1 = rand(Binomial(1, percent_peak1), size(freq1,1))'
        binom2 = abs.(binom1 .- 1)

        norma = Normal(peak1, se1)
        m_conver1 = rand(norma, 1, size(pop1, 2)) .* binom1

        normb = Normal(peak2, se2)
        m_conver2 = rand(normb, 1, size(pop1, 2)) .* binom2

        m_conver = m_conver1 .+ m_conver2

        m_conver[m_conver .< 0] .= mini
        m_conver[m_conver .> 1] .= maxi
    end

    if conversion_distribution == "set"
        m_conver = fill(conv_set, 1, size(pop1, 2))
        m_conver[2] = res
    end

    if conversion_distribution == "beta"
        # artificial m_conver values
        conv_beta = Beta(shape1, shape2)
        m_conver = rand(conv_beta, 1, size(pop1, 2))
    end

    if conversion_distribution == "normal"
        # artificial m_conver values
        conv_norm = Normal(nmean, nstd)
        m_conver = rand(conv_norm, 1, size(pop1, 2))
        m_conver[m_conver .< 0] .= mini
        m_conver[m_conver .> 1] .= maxi
    end

    if conversion_distribution == "sequence"
        m_conver = m_conver_scan(pop1, drivers1, priming, grna)
    end

    if histograms == true
        histogram(m_conver', bins = 50, label = "", title = "initial distribution", xaxis = "conversion efficiency", yaxis = "frequency")
    end


    #################################################################
    #Population variation
    #################################################################

    #only works for 2 populations
    pop_variation = "none" #"subtle" "none"

    #the conversion efficiency at which the populations are separated
    cutoff = 0.5

    #the severity of the difference in populations
    magnitude = 10

    if pops == 2

        if pop_variation == "strict"
            freq1[:,:,1] .= freq1[:,:,1] .* transpose(m_conver .>  cutoff)
            freq1[:,:,2] .= freq1[:,:,2] .* transpose(m_conver .<= cutoff)
        end

        if pop_variation == "subtle"
            adjust1 = ones(size(freq1,1), 1)
            adjust1[m_conver' .>= cutoff] .*= magnitude
            adjust2 = ones(size(freq1,1), 1)
            adjust2[m_conver' .< cutoff] .*= magnitude

            freq1[:,:,1] .*= adjust1
            freq1[:,:,2] .*= adjust2
        end

        if pop_variation == "none"
            #do nothing
        end

    end

    if histograms == true
        ps = repeat([Plots.plot(1)], pops)
        ps[1] = scatter(m_conver',freq1[:,:,1], label = "", title = "Population 1", xaxis = "conversion efficiency", yaxis = "frequency")
        ps[2] = scatter(m_conver',freq1[:,:,2], label = "", title = "Population 2", xaxis = "conversion efficiency", yaxis = "frequency")
        plot(ps...)
    end


    #################################################################
    #Gene Drive starting frequency
    #################################################################

    #starting frequency of gene drive
    starting_freq = drive_freq
    res_freq = resfreq_input
    # val = sum(freq1[2:end,1,1]) + starting_freq
    # freq1[1,1,1] = val * starting_freq
    # freq1 ./= val

    # if pops > 1
    #     freq1[1,:,:] = [0.001, 0.0] #need array to be same size as pops
    # end

    freq1[1] = starting_freq
    freq1[2] = res_freq
    freq1[3] = 1 - freq1[1] - freq1[2]

    if pops > 1
        freq1[1,1,2] = 0
        freq1[2,1,2] = .4
        freq1[3,1,2] = 1 - freq1[1,1,2] - freq1[2,1,2]
    end



    #make all pop allele freqs sum to 1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
    all1 = deepcopy(freq1)

    #################################################################
    #Migration rates
    #################################################################

    mig_style = "set" #"variable"

    if mig_style == "set"
        mig_rates .= migrates_input
    end

    if mig_style == "variable"
        mig_rates = reshape([0, .1, .1, 0], pops, pops) #vector needs to be pops squared length
    end


    #average trait value
    total_traits = Array{Float64}(undef, 0, pops)
    total_conv = Array{Float64}(undef, 0, pops)


    all1 = deepcopy(freq1)

# #change seed for simulation but keep starting scenario
# Random.seed!(50)

    for q in 1:gens
        # global freq1, pop1, all1, total_pop1, total_traits, sexes, sex_bias, pressure, dominance
        # global popsize1, drivers1, traits1, m_conver, total_conv, total_sex, mig_rates, exposure


        #1 generation of gene drive
        cycle1 = gene_drive_pre(freq1, pop1, popsize1, drivers1, traits1, m_conver, Fis, sexes)
        #store new values for next cyle
        freq1 = cycle1[1]
        pop1 = cycle1[2]
        popsize1 = cycle1[3]
        drivers1 = cycle1[4]
        traits1 = cycle1[5]
        m_conver = cycle1[6]


		if migration == true
		    #migration
		    #simulates how much of each allele travels to each pop
		    migrants = repeat([freq1[:,:,1]],pops, pops)
		    for k in 1:pops
		        for j in 1:pops

		            migrants[j,k] = freq1[:,:,j] .* mig_rates[j,k] .* popsize1[j]

		        end
		    end

		    # #adjust each pops frequencies based on migration
		    # #each row of migrants represents all the leavers
		    # #each column represents all the arrivers
		    for popu in 1:pops

		        freq1[:,:,popu] = freq1[:,:,popu] .+ sum(hcat(migrants[:,popu]...) ./ repeat(fill(popsize1[popu][1], 1, size(popsize1, pops)), outer = size(pop1, 2)), dims = 2)

		    end

		    #adjust popsizes
		    popsize1 = round.(popsize1 .+ reshape(sum(sum.(migrants), dims = 1), 1,1,pops) .- reshape(sum(sum.(migrants), dims = 2), 1,1,pops))

		    #normalise freq1
		    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
		end



	        #store new freq in all for plotting purposes
	        all1 = combine_all1(all1, freq1, pops)


	        #store total_pop
	        total_pop1 = vcat(total_pop1, popsize1)
	        total_sex = vcat(total_sex, sexes)


	        #population cleaning
	        #if allele is present for less than 10 generations,
	        #remove it from:
	        #pop1, freq1, drivers1, traits1, m_conver

	            temp_conv = []
	        for j in 1:pops
	            push!(temp_conv, sum(m_conver .* repeat(transpose(freq1[:,:,pops]), size(m_conver, 1)) / size(m_conver, 1)))
	        end
	        total_conv = vcat(total_conv, transpose(temp_conv))


	    #cleaning dead alleles out
	    if cleaning == true
	             #to remove
	             to_keep = collect(1:size(all1,1))
	             for k in 1:size(all1,1)

	                if sum(all1[k, end, :] .< 1e-6) == pops &&
	                    sum(maximum(all1[k, :, :], dims = 1) .< 0.001) == pops

	                    deleteat!(to_keep, findall(to_keep .== k))

	                end
	            end



	            #only keep the alleles specified in to_keep
	            pop1 = pop1[:, to_keep]
	            freq1 = freq1[to_keep, :, :]
	            drivers1 = drivers1[: ,to_keep]
	            traits1 = traits1[:, to_keep]
	            m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)
	            all1 = all1[to_keep, :, :]
	    end



	            temp_traits = []
	        for j in 1:pops
	            push!(temp_traits, 2 * sum(transpose(freq1[:,:,j]) .* sum(traits1, dims = 1)))
	        end
	        temp_traits[temp_traits .> 1] .= 1
	        total_traits = vcat(total_traits, transpose(temp_traits))



	        #dynamic bar plotting
	        #plot data - make drives negative values
	        #pdrives = transpose(map(x -> (-2 * x) + 1, drivers1))
	        #pfreq = sort(freq1 .* pdrives, dims = 2)
	    if dynamic == true
	        dat = transpose(freq1[:,1,:]) .* repeat(vec(popsize1), 1, size(freq1,1)) ./ max_pop
	        nam = repeat(["Bana", "Pala"], outer = size(freq1,1))
	        key = sort(unique(nam))
	        val = collect(1:length(key)) .- 0.5
	        dict = Dict(zip(key, val))
	        ctg = repeat(1:size(freq1,1), inner = size(freq1, 3))
	        p = groupedbar(nam, dat, bar_position = :stack,
	            group = ctg,
	            label = "", title = "Gene Drive Simulation \n Generation $(q)", linecolour = :match,
	            ylim = (0,1.3), xlab = "Populations", ylab = "Frequency",
	            ann = [ (dict["Bana"], 1.05, popsize1[1]),
	                    (dict["Pala"], 1.05, popsize1[2]),
	                    ])

	        display(p)
	        # savefig("frame$(q).png")
	    end


	        # #dynamic pie charts
	        # ps = repeat([Plots.plot(1)], pops)
	        #
	        # for popu in 1:pops
	        #     ps[popu] = pie(freq1[:,:,popu])
	        # end
	        #
	        # display(Plots.plot(ps...))


    # println(q)
    end






    # elapsed = round(time() - start, digits = 3)
    # println(elapsed, " seconds")
    ###################################################################
    #Plotting
    ###################################################################
    #make empty array to store plots

    #make gene drives solid and others dashed
	styles = fill(:dash, length(drivers1))
    styles[vec(drivers1 .== 1)] .= :solid
    styles = reshape(styles, 1, length(styles))

    pop_names = ["Target", "Neighbour", "Neighbour", "Ghana"]

    for popu in 1:pops

        ps[popu] = Plots.plot(ps[popu], transpose(all1[:,:,popu]),
                #annotate = [(30, 0.3, text("Pressure = $(pressure[popu])", 8)),
                #            (30, 0.1, text("Env = $(env[popu])", 8))],
                legend = :topleft,
                xlab = "Generation",
                ylab = "Frequency",
                guidefontsize= 9,
                titlefontsize = 12,
                title = pop_names[popu],
                # left_margin = [-5mm 0mm],
                layout = 1,
                ylim = (0,1),
                linestyle = styles,
                width = 1.5,
				color = [1,2,3]',
                label = "")

            #     #plotting population size
            # Plots.plot!(total_pop1[:,:,popu] ./ max_pop,
            #     label = "",
            #     colour = "black",
            #     ls = :dash,
            #     width = 2)

            #     #plotting average fitness
            # Plots.plot!(total_traits[:, popu],
            # label = "",
            # colour = "black",
            # ls = :dot,
            # width = 2)

            #     #plotting average conversion
            # Plots.plot!(total_conv[:, popu],
            # label = "",
            # colour = "black",
            # ls = :dot,
            # width = 1.5)

            # Plots.vline!([200], color = "black", linewidth = 5, linestyle = :dot, label = "")
    end




    #plot all plots
	h1 = Plots.plot(ps..., xlim = (0,gens), ylim = (0,1.05), reuse = true, label = "", size = ((600 + ((pops > 1) * 600)),(400 + ((pops > 2) * 400))))
    # Plots.annotate!(50, 0.5, popsize1[1])


end

return(h1)

end




using PyPlot
using JLD2
using DataFrames
using Plots
using Interact, Mux
# using KernelDensity, StatsBase

#need to run model in the repl

gr(size = (500,400)) #local
# gr(size = (200,100)) #VM

mp = @manipulate throttle= 0.1 for  zygotic = widget(Dict("Post-Zygotic" => 1, "Pre-Zygotic" => 2), label="Dataset"),

                                    popsi = Widgets.slider((1:1:3), value = 1, label = "Number of Populations"),
                                    popi = Widgets.slider((1:1:7), value = 0.0, label = "Population size (1e-value-)"),
                                    gensi = Widgets.slider((10:10:500), value = 100, label = "Number of Generations"),
                                    convi = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Conversion Efficiency"),
                                    fitnesscost = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Fitness Cost"),
                                    exposurerate = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Exposure Rate"),
                                    dominance = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Degree of Dominance"),
                                    resistance = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Resistance Level"),
                                    inbreeding = Widgets.slider((0:0.01:1.0), value = 0.0, label = "Inbreeding"),
                                    drivefreq = Widgets.slider((0:0.01:0.3), value = 0.0, label = "Drive Frequency"),
                                    resfreq = Widgets.slider((0:0.01:0.7), value = 0.0, label = "Resistance Frequency"),
                                    migrants = Widgets.slider((-1.0:-1.0:-5.0), value = 0.0, label = "Migration Rate"),
                                    migon = Widgets.dropdown(([true, false]), value = false, label = "Migration?"),
                                    iterations = Widgets.slider((0:1:1), value = 0.0, label = "10 ^ Iterations")



        # plt.close("all")

        if zygotic == 1

            run_model_post(pop_input = 10 ^ (popi),
                      pops_input = popsi,
                      indiv_input = 3,
                      gens_input = gensi,
                      fis_input = inbreeding,
                      migration_input = migon,
                      exposure1 = exposurerate,
                      exposure2 = 1.0,
                      dominance_input = dominance,
                      fitness_input = fitnesscost,
                      conversion_input = convi,
                      resistance_input = resistance,
                      drive_freq = drivefreq,
                      resfreq_input = resfreq,
                      migrates_input = 10 ^ migrants,
                      iterations_input = 10 ^ iterations)

      elseif zygotic == 2

          run_model_pre(pop_input = 10 ^ (popi),
                    pops_input = popsi,
                    indiv_input = 3,
                    gens_input = gensi,
                    fis_input = inbreeding,
                    migration_input = migon,
                    exposure1 = exposurerate,
                    exposure2 = 1.0,
                    dominance_input = dominance,
                    fitness_input = fitnesscost,
                    conversion_input = convi,
                    resistance_input = resistance,
                    drive_freq = drivefreq,
                    resfreq_input = resfreq,
                    migrates_input = 10 ^ migrants,
                    iterations_input = 10 ^ iterations)
    end





end

@layout! mp vbox(hbox(pad(1em, :zygotic)),
                 hbox(pad(1em, :popsi), pad(1em, :popi), pad(1em, :gensi)),
                 hbox(pad(1em, :convi), pad(1em, :inbreeding), pad(1em, :migon)),
                 hbox(pad(1em, :fitnesscost), pad(1em, :exposurerate), pad(1em, :migrants)),
                 hbox(pad(1em, :resistance), pad(1em, :resfreq)),
                 hbox(pad(1em, :dominance), pad(1em, :drivefreq), pad(1em, :iterations)),

                 observe(_))

ui = dom"div"(mp);
WebIO.webio_serve(page("/", req -> ui), 8000);

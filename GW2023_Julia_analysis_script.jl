# GW2023_Julia_analysis_script.jl

# Started by Darren Irwin on 14June2023.
# Based on R code written for Greenish Warbler analysis (Irwin et al. 2016), 
# and then the North American warbler analyses (Irwin et al. 2019),
# and then the 2019 Greenish Warbler analysis.
# This one is with the new 2022 Biozeron genome assembly.
# This adapted from the R script called GW2022_R_analysis_script.R
# as an experiment to see how well I can do switch over the analysis to Julia.
# Also basing much of this on IrwinLabGenomicsAnalysisScript.jl 
# which I developed earlier to help Libby and then Rashika 
# process scripts faster.


# If needing to load packages, run what is in this comment section. Will take some time to install the packages:
#=
import Pkg; Pkg.add("CSV") # took less than a minute
Pkg.add("DataFrames") # took about a minute
Pkg.add("Plots") # seems to install and working more simply than Makie (but less powerful)
Pkg.add("Haversine") # for great circle (Haversine) distances
Pkg.add("Distributions") # this seemed to fix a problem installing GLMakie
Pkg.add("MultivariateStats")
Pkg.add("StatsBase")
Pkg.add("Impute")
Pkg.add("JLD2")
Pkg.add("CairoMakie")
Pkg.add("PrettyTables") # for printing nice tables to REPL
 =#

using CSV # for reading in delimited files
using DataFrames # for storing data as type DataFrame
# using Plots # for plotting (uses GR backend by default)
using Haversine # for calculating Great Circle (haversine) distances between sites
using MultivariateStats # for Principal Coordinates Analysis (multidimensional scaling)
using DelimitedFiles # for reading delimited files (the genotypic data)
using Impute # for imputing missing genotypes
using JLD2 # for saving data
using CairoMakie # for plots

CairoMakie.activate!()  # this makes CairoMakie the main package for figures (in case another loaded)

include(".SNPlots.jl") # load file containing custom-built functions

using .SNPlots # actually make SNPlots module available with SNPlots.functionName(),
# or if functions are exported from SNPlots then they are available.

# choose working directory
cd("/Users/darrenirwin_office_iMac/Dropbox/Darren's current work/")

# Load LatLongs of locations ----

latlong_filepath = "GW manuscript 2022_2023/GW_locations_LatLong_2023_forR.txt"

latlongs = DataFrame(CSV.File(latlong_filepath))

print(latlongs)

# make a quick plot to inspect latlong data:
scatter(latlongs.long_E, latlongs.lat_N)

# remove nitidus
latlongs2 = latlongs[Not(latlongs.subspecies .== "nitidus"), :]

# make a matrix of great circle distances (Haversine, assuming spherical Earth which is really close)
geoPoints = GeoLocation.(latlongs2.long_E, latlongs2.lat_N)
# this next line is so neat--uses list comprehension to make a matrix of pairwise calculations
distances = [(HaversineDistance(geoPoints[i], geoPoints[j])/1000) for i in eachindex(geoPoints), j in eachindex(geoPoints)]

# Construct a matrix of distances around the ring ----

# Now construct the matrix for around the ring:
# to check names of sites: latlongs2.Location_name
# get some key distances
function getIndex(name, nameVector = latlongs2.Location_name)
    findfirst(isequal(name), nameVector)
end

index_AA = getIndex("Ala_Archa")
index_PK = getIndex("Naran_Pakistan")
index_LN = getIndex("Langtang")
index_EM = getIndex("Emeishan")
index_XN = getIndex("Xining")
index_BJ = getIndex("Beijing")
index_last = nrow(latlongs2)

dist_PK_to_LN = distances[index_PK, index_LN]
dist_LN_to_EM = distances[index_LN, index_EM]
dist_EM_to_BJ = distances[index_EM, index_BJ]

# This next part will assume locations in the input file are arranged in order around ring:
distsAroundRing = Matrix{Float32}(undef, size(distances)[1], size(distances)[2])
# accept all distances within viridanus:

# function for accepting straight-line great circle dists as distances between sets of sites
acceptDists = function(straightGreatCircleDists, start, finish, distsAroundRing)
    distsAroundRing[start:finish, start:finish] = straightGreatCircleDists[start:finish, start:finish]
    return(distsAroundRing)
end

# accept all distances within viridanus:
distsAroundRing = acceptDists(distances, 1, index_AA, distsAroundRing)

# accept dist from AA to PK:
distsAroundRing = acceptDists(distances, index_AA, index_PK, distsAroundRing)

# accept all distances from PK to LN:
distsAroundRing = acceptDists(distances, index_PK, index_LN, distsAroundRing)

# accept dist from LN to EM:
distsAroundRing = acceptDists(distances, index_LN, index_EM, distsAroundRing)

# accept dists between EM, XN, BJ:
distsAroundRing = acceptDists(distances, index_EM, index_BJ, distsAroundRing)

# accept all distances within plumbeitarsus:
distsAroundRing = acceptDists(distances, index_BJ, index_last, distsAroundRing)

display(distsAroundRing)

# function for adding up distances measured through certain sites:
addDists = function(set1start, set1end, set2start, set2end, distsAroundRing)
    firstDists = repeat(distsAroundRing[set1start:(set1end-1), set1end], 1, set2end-set2start+1)
    secondDists = repeat(transpose(distsAroundRing[set1end, set2start:set2end]), set1end-set1start, 1)
    totalDists = firstDists + secondDists
    distsAroundRing[set1start:(set1end-1), set2start:set2end] = totalDists
    distsAroundRing[set2start:set2end, set1start:(set1end-1)] = transpose(totalDists)
    return(distsAroundRing)
end

# dists from viridanus to PK are sum of dists to AA plus AA to PK:
distsAroundRing = addDists(1, index_AA, index_PK, index_PK, distsAroundRing)

# dists from "northwest of PK" to Himalayas are sum of ringdists to PK plus PK to locations up to LN:
distsAroundRing = addDists(1, index_PK, index_PK+1, index_LN, distsAroundRing)

# dists from "west / northwest of LN" to EM are sum of dists to LN plus LN to EM:
distsAroundRing = addDists(1, index_LN, index_EM, index_EM, distsAroundRing)

# dists from "west / northwest of EM" to China are sum of dists to EM plus EM to (XN, BJ):
distsAroundRing = addDists(1, index_EM, index_XN, index_BJ, distsAroundRing)

# dists from "west of BJ" to east Siberia are sum of dists to BJ plus BJ to other plumbeitarsus:
distsAroundRing = addDists(1, index_BJ, index_BJ+1, index_last, distsAroundRing)

# conduct Principal Coordinates Analysis on the distances around the ring,
# to make a single location around ring axis:
PCO_around_ring = fit(MDS, distsAroundRing; distances=true, maxoutdim=1)
# add this as a column to the data frame:
latlongs2[:, :LocationAroundRing] = vec(-predict(PCO_around_ring))
# another way: 
#latlongs2.LocationAroundRing = vec(-predict(PCO_around_ring))


# Before running below, changed 012NA file back into 012minus1 file, 
# using commands like below in Terminal, so can be read as integer:
# cd cd /Users/darrenirwin/Dropbox/Darren\'s\ current\ work/GW_genomics_2022_with_new_genome/GW2022_GBS_012NA_files
# cat GW2022_all4plates.genotypes.SNPs_only.whole_genome.max2allele_noindel.vcf.maxmiss60.MQ20.lowHet.tab.012NA | sed 's/NA/-1/g' > GW2022_all4plates.genotypes.SNPs_only.whole_genome.max2allele_noindel.vcf.maxmiss60.MQ20.lowHet.tab.012minus1

# PCA whole-genome ----
# Load 012NA file containing only variable sites throughout genome;
# construct a PCA based on all sites passing an Fst threshold between the "groups" below;
# and all individuals in "groups.to.plot.PCA" according to colors in "group.colors.PCA"
groups = ["vir","troch_LN","plumb"]  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
groups_to_plot_PCA = ["vir","vir_misID","vir_S","nit", "lud_PK", "lud_KS", "lud_central", "lud_Sath", "lud_ML","troch_west","troch_LN","troch_EM","obs","plumb_BJ","plumb","plumb_vir"]
group_colors_PCA = ["blue","blue","turquoise1","grey","seagreen4","seagreen3","seagreen2","olivedrab3","olivedrab2","olivedrab1","yellow","gold","orange","pink","red","purple"]

# choose path and filename for the 012NA files
baseName = "GW_genomics_2022_with_new_genome/GW2022_GBS_012NA_files/GW2022_all4plates.genotypes.SNPs_only.whole_genome"
filenameTextMiddle = ".max2allele_noindel.vcf.maxmiss"
# indicate percent threshold for missing genotypes for each SNP--
# this was set by earlier filtering, and is just a record-keeper for the filenames:
missingGenotypeThreshold = 60 
filenameTextEnd = ".MQ20.lowHet.tab"

tagName = ".Aug2023."   # choose a tag name for this analysis
# indicate name of metadata file, a text file with these column headings:
# ID	location	group	Fst_group	plot_order
metadataFile = "GW_genomics_2022_with_new_genome/GW_all4plates.Fst_groups.txt"
# load metadata
metadata = DataFrame(CSV.File(metadataFile)) # the CSV.File function interprets the correct delimiter
num_metadata_cols = ncol(metadata)
num_individuals = nrow(metadata) 
# specify window size (number of bp with info) and step size
#windowSize = 10000   
#stepSize = windowSize  # could change if wanting overlapping windows (not built in yet to Julia version)

# read in individual names for this dataset
individuals_file_name = string(baseName, filenameTextMiddle, missingGenotypeThreshold, filenameTextEnd, ".012.indv")
ind = DataFrame(CSV.File(individuals_file_name; header=["ind"], types=[String])) 
indNum = size(ind, 1) # number of individuals
if num_individuals != indNum
    println("WARNING: number of rows in metadata file different than number of individuals in .indv file")
end
# read in position data for this dataset
position_file_name = string(baseName, filenameTextMiddle, missingGenotypeThreshold, filenameTextEnd, ".012.pos")
pos_whole_genome = DataFrame(CSV.File(position_file_name; header=["chrom", "position"], types=[String, Int]))
# read in genotype data
column_names = ["null"; string.("c.", pos_whole_genome.chrom, ".", pos_whole_genome.position)]    
genotype_file_name = string(baseName, filenameTextMiddle, missingGenotypeThreshold, filenameTextEnd, ".012minus1") 
@time if 1 <= indNum <= 127   
    geno = readdlm(genotype_file_name, '\t', Int8, '\n'); # this has been sped up dramatically, by first coverting "NA" to -1
elseif 128 <= indNum <= 32767
    geno = readdlm(genotype_file_name, '\t', Int16, '\n'); # this needed for first column, which is number of individual; Int16 not much slower on import than Int8
else
    print("Error: Number of individuals in .indv appears outside of range from 1 to 32767")
end
loci_count = size(geno, 2) - 1   # because the first column is not a SNP (just a count from zero)
print(string("Read in genotypic data at ", loci_count," loci for ", indNum, " individuals. \n"))
# took 324 seconds to: Read in genotypic data at 2431709 loci for 310 individuals.
# In Julia 1.9, now took only 52 seconds!

# Check that individuals are same in genotype data and metadata ----
ind_with_metadata = hcat(ind, metadata)
print(ind_with_metadata)
print("\n")  # prints a line break 
if isequal(ind_with_metadata.ind, ind_with_metadata.ID)
    println("Good news: names of individuals in metadata file and genotype ind file match perfectly.")
else
    println("WARNING: names of individuals in metadata file and genotype ind file do not completely match.")
end

# Filter individuals ----
# If need to filter out particular individuals:
filter = false
filter_out_inds = [] # if filtering out individuals, specify their row number here, e.g. "[20, 103]"
if filter
    # Specify individuals to filter out:
    ind_with_metadata_indFiltered = ind_with_metadata[Not(filter_out_inds), :]
    geno_indFiltered = geno[Not(filter_out_inds), :]
    println("Specific individuals filtered out as requested")
else
    ind_with_metadata_indFiltered = ind_with_metadata
    geno_indFiltered = geno
    println("No specific individuals requested to be filtered out")
end

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
SNPmissing_percent_allowed_per_ind = 40   # this is the percentage threshold
threshold_missing = loci_count * SNPmissing_percent_allowed_per_ind/100
numMissings = sum(geno .== -1, dims=2)
selection = vec(numMissings .<= threshold_missing) # the vec command converts to BitVector rather than BitMatrix--important below
geno_indFiltered = geno_indFiltered[selection, :]
# print filtered out individuals:
ind_with_metadata_indFiltered.ind[selection.==false]
println(ind_with_metadata_indFiltered.ind[selection.==false])
# ["GW_Armando_plate1_JG08G02", "GW_Armando_plate1_JG10G01", "GW_Armando_plate1_NO_BC_TTGW05", 
# "GW_Armando_plate1_NO_DNA", "GW_Armando_plate1_TTGW21", "GW_Armando_plate1_TTGW71", 
# "GW_Armando_plate2_NO_BC_TTGW05", "GW_Armando_plate2_NO_DNA", "GW_Armando_plate2_TTGW15", 
# "GW_Lane5_AA10", "GW_Lane5_DA6", "GW_Lane5_LN11", 
# "GW_Liz_GBS_Liz5101", "GW_Liz_GBS_Liz5101_R", "GW_Liz_GBS_Liz5118", 
# "GW_Liz_GBS_Liz5139", "GW_Liz_GBS_Liz5142", "GW_Liz_GBS_Liz5150", 
# "GW_Liz_GBS_Liz5159", "GW_Liz_GBS_Liz5162", "GW_Liz_GBS_Liz5169", 
# "GW_Liz_GBS_Liz5171", "GW_Liz_GBS_Liz5172", "GW_Liz_GBS_Liz5174", 
# "GW_Liz_GBS_Liz5176", "GW_Liz_GBS_Liz5177", "GW_Liz_GBS_Liz5180", 
# "GW_Liz_GBS_Liz5186", "GW_Liz_GBS_Liz5187", "GW_Liz_GBS_Liz5192", 
# "GW_Liz_GBS_Liz5195", "GW_Liz_GBS_Liz6012", "GW_Liz_GBS_Liz6203", 
# "GW_Liz_GBS_Liz6766", "GW_Liz_GBS_P_fusc", "GW_Liz_GBS_P_h_man", 
# "GW_Liz_GBS_P_humei", "GW_Liz_GBS_P_inor", "GW_Liz_GBS_S_burk"]
# at a glance looks identical to results R code filter (also 39 removed, leaving 271)
ind_with_metadata_indFiltered = ind_with_metadata_indFiltered[selection, :]


# Filter SNPs ----

# filter out SNPs with too many missing genotypes:
# (remember that first column is arbitrary row number in input file)
missing_genotypes_per_SNP = sum(geno_indFiltered .== -1, dims=1)
missing_genotypes_percent_allowed_per_site = 5   # this is the percentage threshold
threshold_genotypes_missing = size(geno_indFiltered)[1] * missing_genotypes_percent_allowed_per_site/100
selection = vec(missing_genotypes_per_SNP .<= threshold_genotypes_missing)
geno_ind_SNP_filtered = geno_indFiltered[:, selection] 
pos_SNP_filtered = pos_whole_genome[selection[Not(1)],:]  # the Not(1) is needed because first column in geno is arbitrary row number

# the above leaves 1003400, the same number as the equivalent R code for filtering at 5% maximum missing genotypes

#### INSERTING 22 AUGUST 2023: ADD SECOND ROUND OF FILTERING OF INDIVIDUALS BASED 
# ON MISSING GENOTYPES:

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
SNPmissing_percent_allowed_per_ind_round2 = 10   # this is the percentage threshold
threshold_missing = (size(geno_ind_SNP_filtered, 2) - 1) * SNPmissing_percent_allowed_per_ind_round2/100
numMissings = sum(geno_ind_SNP_filtered .== -1, dims=2)
selection = vec(numMissings .<= threshold_missing) # the vec command converts to BitVector rather than BitMatrix--important below
geno_ind_SNP_ind_filtered = geno_ind_SNP_filtered[selection, :]
# print filtered out individuals:
ind_with_metadata_indFiltered.ind[selection.==false]
println(ind_with_metadata_indFiltered.ind[selection.==false])
# filtered out: ["GW_Armando_plate1_TTGW74", "GW_Armando_plate2_TTGW54", "GW_Lane5_AA8", "GW_Lane5_YK1"]
ind_with_metadata_indFiltered = ind_with_metadata_indFiltered[selection, :]
# this leaves 267 individuals and 1003400 loci, with no individuals missing more than 10% of genotypes
# and no loci missing in more than 5% of individuals.

# CONVERT -1 to MISSING
genos_with_missing = Matrix{Union{Missing, Int32}}(geno_ind_SNP_ind_filtered)
genos_with_missing[genos_with_missing .== -1] .= missing


# remove first column of genotype data (which is arbitrary)
genosOnly = genos_with_missing[:,Not(1)]

# conduct a simple Principal Components Analysis on the individuals around the ring
# (can expand this later to all the bells and whistles of the R version)

# I haven't found a great way to impute missing genotypes in Julia,
# but this has caused me to rethink the PCA methodology. Why impute genotypes at all?
# Why not just generate a distance matrix, and do PCO on that?
# Make a pairwise genetic distance matrix using the shared genotypes,
# with distances based on 2 changes for both alleles different (different homozygotes)
# and 1 change for heterozygote vs. homozygote.
# THIS IS NOW DONE FAR BELOW--BUT IN THE END I DON'T THINK IT IS VERY USEFUL
# BECAUSE IT LOSES TOO MUCH INFO. WITHIN-CLUSTER DISTANCES ARE EXPANDED A BIT,
# COMPARED TO THE PCA BASED ON GENOTYPES DIRECTLY. SO STICKING WITH THE IMPUTATION
# AND PCA AS DONE HERE:


#######
# cycle through chromosomes, imputing and saving genotypes in each: 
chromosomes_to_process = ["gw2",
                            "gw1",
                            "gw3",
                            "gwZ",
                            "gw1A",
                            "gw4",
                            "gw5",
                            "gw7",
                            "gw6",
                            "gw8",
                            "gw9",
                            "gw11",
                            "gw12",
                            "gw10",
                            "gw13",
                            "gw14",
                            "gw18",
                            "gw20",
                            "gw15",
                            "gw1B",
                            "gws100",
                            "gw17",
                            "gw19",
                            "gws101",
                            "gw4A",
                            "gw21",
                            "gw26",
                            "gws102",
                            "gw23",
                            "gw25",
                            "gws103",
                            "gw22",
                            "gws104",
                            "gw28",
                            "gw27",
                            "gw24",
                            "gws105",
                            "gws106",
                            "gws107",
                            "gws108",
                            "gws109",
                            "gws110",
                            "gws112"]

for i in eachindex(chromosomes_to_process)
    chrom = chromosomes_to_process[i]
    regionText = string("chr", chrom)
    loci_selection = (pos_SNP_filtered.chrom .== chrom)
    pos_SNP_filtered_region = pos_SNP_filtered[loci_selection,:]
    genosOnly_region_for_imputing = Matrix{Union{Missing, Float32}}(genosOnly[:,loci_selection])
    @time imputed_genos = Impute.svd(genosOnly_region_for_imputing)
    filename = string(baseName, tagName, regionText, ".imputedMissing.jld2")
    jldsave(filename; imputed_genos, ind_with_metadata_indFiltered, pos_SNP_filtered_region)
    println(string("Chromosome ", chrom, ": Saved real and imputed genotypes for ", size(pos_SNP_filtered_region, 1)," SNPs and ", size(genosOnly_region_for_imputing, 1)," filtered individuals."))
end
# most chromosomes took less than 200 seconds, some far less, but chr 17 and 27 took more than 1000
# length of time only loosely based on size.

for i in eachindex(chromosomes_to_process)
    chrom = chromosomes_to_process[i]
    regionText = string("chr", chrom)
    filename = string(baseName, tagName, regionText, ".imputedMissing.jld2")
    imputed_genos = load(filename, "imputed_genos")
    ind_with_metadata_indFiltered = load(filename, "ind_with_metadata_indFiltered")
    pos_SNP_filtered_region = load(filename, "pos_SNP_filtered_region")
    println(string("Loaded ",filename))
    println(string(regionText, ": ", size(imputed_genos,2), " SNPs from ", size(imputed_genos,1), " individuals"))
    plotPCA(imputed_genos, ind_with_metadata_indFiltered, groups_to_plot_PCA; regionText)
end


# this now works well--appears to NOT be scaling the variables to variance 1, which is good
@time PCA_indGenos = fit(PCA, matrixForPCA; method = :svd, maxoutdim=3); # good to suppress output of this--otherwise many lines
PCA_values = predict(PCA_indGenos, matrixForPCA)


f = CairoMakie.Figure()
ax = Axis(f[1, 1], 
    autolimitaspect = 1,  # this sets the axis scaling to be equal in terms of the units
    title = string("PCA of greenish warblers, ", regionText),
    xlabel = "PC1",
    ylabel = "PC2")
for i in eachindex(groups_to_plot_PCA) 
    selection = ind_with_metadata_indFiltered.Fst_group .== groups_to_plot_PCA[i]
    CairoMakie.scatter!(ax, PC1[selection], PC2[selection], marker = :diamond, color=group_colors_PCA[i], markersize=10, strokewidth=0.5)
end
f


#######

genosOnly_for_imputing = Matrix{Union{Missing, Float32}}(genosOnly)
#genosOnly_for_imputing = genosOnly_for_imputing[:,1:100000]

regionText = "wholeGenome"
filename = string(baseName, tagName, regionText, ".imputedMissing.jld2")
# to do the imputing, do this by setting to true, but TAKES A LONG TIME:
do_imputing = false
if do_imputing
    @time imputed_genosOnly = Impute.svd(genosOnly_for_imputing)
    # took almost 2 hours!
    # to speed up, maybe could set convergence tolerance to higher value, currently: 
    # tol::Float64: convergence tolerance (default: 1e-10) (But tried this on smaller chr, and did not help)
    jldsave(filename; imputed_genosOnly, ind_with_metadata_indFiltered, pos_SNP_filtered)
    print("Saved matrix of real and imputed genotypes for filtered individuals. \n")
else # load the already saved imputing
    imputed_genosOnly = load(filename, "imputed_genosOnly")
end

# including the :rows vs. :cols doesn't seem to matter
#@time imputed_genosOnly = Impute.svd(genosOnly_for_imputing; dims=:rows)

# sum(imputed_genosOnly .!= imputed_genosOnlyB)

#findfirst(x -> 0.1 < x < 0.9, imputed_genosOnly  )

#imputed_genosOnly_noMissing = Matrix{Float32}(imputed_genosOnly)



plotPCA()
# looks virtually identical to the R plot used in my powerpoints. Great!


# option to filter out all but selected chromosome (or set of them):
chooseChrom = "one_chr"      # "one_chr" or "whole_genome" or "all_autosomes"
if chooseChrom == "one_chr"
    chrom = "gw26"
    regionText = string("chr", chrom)
    loci_selection = (pos_SNP_filtered.chrom .== chrom)
    pos_SNP_filtered_region = pos_SNP_filtered[loci_selection,:]
    genosOnly_region = genosOnly[:,loci_selection]
elseif chooseChrom == "all_autosomes"
    regionText = "all_autosomes"	
    loci_selection = (pos_SNP_filtered.chrom .!= "gwZ")
    pos_SNP_filtered_region = pos_SNP_filtered[loci_selection,:]
    genosOnly_region = genosOnly[:,loci_selection]
elseif chooseChrom == "whole_genome"
    regionText = "whole_genome"
    pos_SNP_filtered_region = pos_SNP_filtered
    genosOnly_region = genosOnly
else 
    print("Warning: Missing a correct chooseChrom setting")
end
# tested the above options and they work

genosOnly_region_for_imputing = Matrix{Union{Missing, Float32}}(genosOnly_region)


@time imputed_genos = Impute.svd(genosOnly_region_for_imputing)
# On Mac Studio, takes only 27 sec for chr gw26 14101 SNPs on chr gw26
# Below times are on my laptop:
# 47 sec for default tol = 1e-10 and 14101 SNPs on chr gw26 
# 44 sec for tol = 1e-8 and 14101 SNPs on chr gw26 (looks same as above)
# 45 sec for tol = 1e-6 and 14101 SNPs on chr gw26 (looks same as above)
# 46 sec for tol = 1e-4 and 14101 SNPs on chr gw26 (looks same as above)
# 44 sec for tol = 1e-3 and 14101 SNPs on chr gw26 (looks same as above)
# 2.5 sec for tol = 1e-2 and 14101 SNPs on chr gw26 (slightly different from above)
# I guess should just leave as default

# 37 seconds for tol = 1e-4 and 11000 SNPs on chr gw28
# 35 seconds for tol = 1e-5
# 36 s for tol = 1e-6
# when I cut maxiter down to 50, time was halved (compared to 100), but repeatability went down.
# 1.5 s for tol = 1e-2
# 32 s for tol = 1e-3
# I guess should just stick with default--it has best repeatability and saving time seems to require large loss in repeatability.



matrixForPCA = Matrix{Float32}(transpose(imputed_genos))

# this now works well--appears to NOT be scaling the variables to variance 1, which is good
@time PCA_indGenos = fit(PCA, matrixForPCA; method = :svd, maxoutdim=3); # good to suppress output of this--otherwise many lines
PCA_values = predict(PCA_indGenos, matrixForPCA)

# option to flip axes by putting minus sign in these lines:
PC1 = -PCA_values[1,:]
PC2 = -PCA_values[2,:]

using PrettyTables
PCtable = DataFrame(ind = ind_with_metadata_indFiltered.ind,
                    Fst_group = ind_with_metadata_indFiltered.Fst_group,
                    PC1 = PC1, 
                    PC2 = PC2)

pretty_table(PCtable, crop = :none)
PCtable.PC1[PCtable.PC1 .> 5]
# find individuals according to PC values:
PCtable[findall(20 .<  PCtable.PC1 .< 30 .&& -10 .< PCtable.PC2 .< -5), :]


f = CairoMakie.Figure(dpi=1000)
ax = Axis(f[1, 1], 
    autolimitaspect = 1,  # this sets the axis scaling to be equal in terms of the units
    title = string("PCA of greenish warblers, ", regionText),
    xlabel = "PC1",
    ylabel = "PC2")
for i in eachindex(groups_to_plot_PCA) 
    selection = ind_with_metadata_indFiltered.Fst_group .== groups_to_plot_PCA[i]
    CairoMakie.scatter!(ax, PC1[selection], PC2[selection], marker = :diamond, color=group_colors_PCA[i], markersize=10, strokewidth=0.5)
end
f
# save as high-res bitmap image:
save("test_figure.png", f, px_per_unit=2)
# save as vector graphics image:
save("test_figure.pdf", f)


# Genotype-by-individual plot ----
# Graph genotypes of high-Fst loci along a chromosome
# choose only loci that are variable in the dataset (SNPs), and above an Fst threshhold
# groups.to.compare =  "Fst_among"   #"Fst_among"  #"vir_troch"   # can choose a pair or Fst_among for multi-group Fst_among
set = "west_side_of_ring"  #"east_side_of_ring"    #"67_inds_around_ring"  # "west_side_of_ring"

if set == "67_inds_around_ring"
    groups = ["vir","troch_LN","plumb"] # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)   
    plotGroups = ["vir","vir_S","lud_PK","lud_KS","lud_central","troch_LN","troch_EM","obs", "plumb_BJ","plumb"]
    plotGroupColors = ["blue","turquoise1","seagreen4","seagreen3","seagreen2","yellow","gold","orange", "pink","red"]
    numIndsToPlot = [10, 5, 6, 2, 7, 15, 15, 15, 15, 15] # maximum number of individuals to plot from each group
    group1 = "vir"   # these groups will determine the color used in the graph
    group2 = "plumb"
    groupsToCompare = "vir_plumb"   #"Fst_among"  #"vir_troch_LN"       #"vir_plumb"      #"troch_LN_plumb"      #"vir_troch_LN"
    Fst_cutoff = 0.7
    missingFractionAllowed = 0.2  # only show SNPs with less than this fraction of missing data among individuals
elseif set == "37_inds_around_ring_plusAllVirPlumb"
    groups = ["vir","troch_LN","plumb"] # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
    plotGroups = ["vir","lud","troch_LN","troch_EM","obs", "obs_plumb","plumb"]
    plotGroupColors = ["blue","seagreen4","yellow","gold","orange", "pink","red"]
    numIndsToPlot = [100, 15, 15, 15, 15, 15, 100] # maximum number of individuals to plot from each group
    group1 = "vir"   # these groups will determine the color used in the graph
    group2 = "plumb"
    groupsToCompare = "Fst_among"
    Fst_cutoff = 0.7
    missingFractionAllowed = 0.2  # only show SNPs with less than this fraction of missing data among individuals
elseif set == "west_side_of_ring"
    groups = ["vir","troch_LN"] # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
    plotGroups = ["vir","vir_misID","vir_S","nit", "lud_PK", "lud_KS", "lud_central", "lud_Sath", "lud_ML","troch_west","troch_LN"]
    plotGroupColors = ["blue","blue","turquoise1","grey","seagreen4","seagreen3","seagreen2","olivedrab3","olivedrab2","olivedrab1","yellow"]
    numIndsToPlot = [15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15] # maximum number of individuals to plot from each group
    group1 = "vir"   # these groups will determine the color used in the graph
    group2 = "troch_LN"
    groupsToCompare = "vir_troch_LN" # "Fst_among"
    Fst_cutoff = 0.6
    missingFractionAllowed = 0.2  # only show SNPs with less than this fraction of missing data among individuals
elseif set == "all_ludlowi_plus_a_few_other"
    groups = ["vir","troch_LN","plumb"] # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
    plotGroups = ["vir","vir_S","nit", "lud_PK", "lud_KS", "lud_central", "lud_Sath", "lud_ML","troch_west","troch_LN","plumb"]
    plotGroupColors = ["blue","turquoise1","grey","seagreen4","seagreen3","seagreen2","olivedrab3","olivedrab2","olivedrab1","yellow","red"]
    numIndsToPlot = [4, 4, 4, 1000, 1000, 1000, 1000, 1000, 1000, 4, 4] # maximum number of individuals to plot from each group
    group1 = "vir"   # these groups will determine the color used in the graph
    group2 = "troch_LN"
    groupsToCompare = "vir_troch_LN" # "Fst_among"
    Fst_cutoff = 0.6
    missingFractionAllowed = 0.2  # only show SNPs with less than this fraction of missing data among individuals
elseif set == "east_side_of_ring"
    groups = ["troch_LN","obs","plumb"] # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
    plotGroups = ["troch_LN","troch_EM","obs","obs_plumb","plumb"]
    plotGroupColors = ["yellow","gold","orange","pink","red"]
    numIndsToPlot = [15, 15, 15, 15, 15] # maximum number of individuals to plot from each group
    group1 = "troch_LN"   # these groups will determine the color used in the graph
    group2 = "plumb"
    groupsToCompare = "troch_LN_plumb"
    Fst_cutoff = 0.7
    missingFractionAllowed = 0.2  # only show SNPs with less than this fraction of missing data among individuals
elseif set == "vir_plumb"
    groups = ["vir","plumb"]
    plotGroups = ["vir","plumb_vir","plumb"]
    plotGroupColors = ["blue","purple","red"]
    numIndsToPlot = [100,100,100] # maximum number of individuals to plot from each group
    group1 = "vir"   # these groups will determine the color used in the graph
    group2 = "plumb"
    groupsToCompare = "vir_plumb"
    Fst_cutoff =  0.7
    missingFractionAllowed = 0.2  # only show SNPs with less than this fraction of missing data among individuals
end

# Calculate allele freqs and sample sizes (use column Fst_group)
freqs, sampleSizes, freqsRowNames = getFreqsAndSampleSizes(genosOnly, ind_with_metadata_indFiltered.Fst_group, groups)
println("Calculated population allele frequencies and sample sizes")

# calculate WC84_Fst 
Fst, FstNumerator, FstDenominator, pairwiseNamesFst = getFst(freqs, sampleSizes, groups; among=true)  # set among to FALSE if no among Fst wanted (some things won't work without it) 
println("Calculated Fst values")

genosOnly_included, ind_with_metadata_included = limitIndsToPlot(plotGroups, numIndsToPlot, genosOnly, ind_with_metadata_indFiltered)

chr = "gw26"

#regionInfo = chooseChrRegion(pos_SNP_filtered, chr; positionMin=1, positionMax=1000000)

regionInfo = chooseChrRegion(pos_SNP_filtered, chr; positionMin=1, positionMax=NaN) # this gets the maximum position for the chromosome

# NOTE FOR LATER: SHOULD REALLY GET CHROMOSOME LENGTH FOR positionMax

#Fst_cutoff = 0.6
#missingFractionAllowed = 0.2

plotInfo = plotGenotypeByIndividual(groupsToCompare, Fst_cutoff, missingFractionAllowed,
    regionInfo, pos_SNP_filtered, Fst, pairwiseNamesFst, 
    genosOnly_included, ind_with_metadata_included, freqs, plotGroups, plotGroupColors)
# plotInfo contains a tuple with: (f, plottedGenotype, locations, plottedMetadata)
  

# for debugging the above function, set these first:
#= 
pos = pos_SNP_filtered
genoData = genosOnly_included
indMetadata = ind_with_metadata_included
 =#





########################











############

# function to calculate pairwise genetic distance between two individuals,
# as a proportion of maximum possible given all non-missing genotype pairs:
getGenoDistLinear = function(genotypes1, genotypes2)
    diffs = abs.(genotypes1 .- genotypes2)
    sum_diffs = sum(skipmissing(diffs))
    count_missing = sum(ismissing.(diffs))
    return sum_diffs / 2(length(diffs) - count_missing)
end

# will try genotype dists in terms of different allele pairings. So:
# 2 identical homozgyotes = 0 dist
# homozygote vs. heterozygote = 2 dist
# heterozygote vs. heterozygote = 2 dist
# homozygote vs. different homozygote = 4 dist
getGenoDist = function(genotypes1, genotypes2)
    # below, note use of isequal--it essentially ignores the missings, returning false for them
    # homozygote vs. different homozygote = 4 dist
    homozygotes_diff = sum(isequal.(genotypes1, 0) .&& isequal.(genotypes2, 2) .|| isequal.(genotypes2, 0) .&& isequal.(genotypes1, 2))
    homozygote_heterozygote_diff = sum(isequal.(genotypes1, 1) .&& (isequal.(genotypes2, 0) .|| isequal.(genotypes2, 2)) .|| 
                                    isequal.(genotypes2, 1) .&& (isequal.(genotypes1, 0) .|| isequal.(genotypes1, 2)))
    heterozygotes_diff = sum(isequal.(genotypes1, 1) .&& isequal.(genotypes2, 1))
    sum_diffs = 4*homozygotes_diff + 2*homozygote_heterozygote_diff + 2*heterozygotes_diff
    not_missing = sum(.!ismissing.(genotypes1) .|| .!ismissing.(genotypes2))
    return sum_diffs / 4not_missing
end

#test
genotypes1 = [1, 0, 2, missing, 0, 2]
genotypes2 = [1, 2, missing, missing, 0, 1]

@time getGenoDist(genotypes1, genotypes2)

# getGenoDist = function(genotypes1, genotypes2)
#     @time diffs = abs.(genotypes1 .- genotypes2)
#     @time sum_diffs = sum(skipmissing(diffs))
#     @time count_missing = sum(ismissing.(diffs))
#     return @time sum_diffs / 2(length(diffs) - count_missing)
# end

@time getGenoDist(genosOnly[1,:], genosOnly[2,:])

@time getGenoDist(view(genosOnly, 1, :), view(genosOnly, 2, :))
# this speeds things up: 0.019 sec compared to 0.04 sec without using view

#= getGenoDistMatrix = function(genotypesAll, numInds) # sped this up dramatically by including numInds to streamline precompiling
    return [getGenoDist(genosOnly[i,:], genosOnly[j,:]) for i in 1:numInds, j in 1:numInds]
end =#

getGenoDistMatrix = function(genotypes, numInds) # sped this up dramatically by including numInds to streamline precompiling
    return [getGenoDist(view(genotypes, i, :), view(genotypes, j, :)) for i in 1:numInds, j in 1:numInds]
end

@time genoDistMatrix = getGenoDistMatrix(genosOnly_region, 271)
# took 48 seconds with 40 individuals (for whole genome of 1003400 loci)
# took 0.2 sec with 40 inds for just chr 26 (14101 loci)
# took 1.2 sec with 100 inds for just chr 26 (14101 loci)
# took 10.6 sec with 271 inds for just chr 26 (14101 loci)
# took 106 sec with 271 inds for just chr Z (52592 loci)

# new allele pairings approach:
# took 51 sec with 100 inds for just chr 26 (14101 loci)
# took 116 sec with 150 inds for just chr 26 (14101 loci)
# took 376 sec with 271 inds for just chr 26 (14101 loci)

@time PCO_indGenos = fit(MDS, genoDistMatrix; distances=true, maxoutdim=3)
PC_values = predict(PCO_indGenos)

scatter(PC_values[1,:], PC_values[2,:])



# function for doing PCO using Euclidean distance 
# (but this doesn't work with missing data):

getEuclidDist = function(genotypes1, genotypes2)
    squares = (genotypes1 .- genotypes2).^2
    sumOfSquares = sum(squares)
    return sqrt(sumOfSquares)
end


getEuclidDistMatrix = function(genotypesAll, numInds) # sped this up dramatically by including numInds to streamline precompiling
    return [getEuclidDist(genosOnly[i,:], genosOnly[j,:]) for i in 1:numInds, j in 1:numInds]
end

genotypes1 = [0,1,2,2]
genotypes2 = [0,0,0,2]

@time EuclidDistMatrix = getEuclidDistMatrix(genosOnly, 271)
#took 255 seconds

@time PCO_Euclid = fit(MDS, EuclidDistMatrix; distances=true, maxoutdim=3)
PC_values = predict(PCO_Euclid)
# super fast

scatter(PC_values[1,:], PC_values[2,:])

# The Z chromosome PCA looks no good using this approach--I think it is separating out males and females.

################

@time getGenoDistMatrix(genosOnly)

@time geno_distances = [getGenoDist(genosOnly[i,:], genosOnly[j,:]) for i in 1:2, j in 1:2]


@time genotypes1 = geno_ind_SNP_filtered[1,Not(1)]

@time genotypes1 = geno_ind_SNP_filtered[1,:]

genotypes2 = geno_ind_SNP_filtered[2,Not(1)]

@time diffs = abs.(genotypes1 .- genotypes2)
    @time sum_diffs = sum(skipmissing(diffs))
    @time count_missing = sum(ismissing.(diffs))




geno_distances = [getGenoDist(geno_ind_SNP_filtered[i,Not(1)], geno_ind_SNP_filtered[j,Not(1)]) for i in eachindex(geno_ind_SNP_filtered[:,1]), j in eachindex(geno_ind_SNP_filtered[:,1])]








PCA_object = fit(PCA, transpose(temp2[:,Not(1)]); method = :cov, maxoutdim=3)

predict(PCA_object)

## NEED TO CONVERT -1 TO MISSING!!
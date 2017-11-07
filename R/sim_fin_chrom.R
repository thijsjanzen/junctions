sim_fin_chrom <-
function (pop_size, initial_heterozygosity, total_runtime, morgan, 
    seed, R) 
{
    .Call("_junctions_sim_fin_chrom", PACKAGE = "junctions", pop_size, 
        initial_heterozygosity, total_runtime, morgan, 
        seed, R)
}


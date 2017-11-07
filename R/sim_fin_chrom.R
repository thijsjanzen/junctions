sim_fin_chrom <-
function (pop_size, initial_heterozygosity, total_runtime, morgan, 
    seed, R) 
{
    .Call("_junctions_sim_inf_chrom", PACKAGE = "junctions", pop_size, 
        R, initial_heterozygosity, total_runtime, morgan, 
        seed)
}


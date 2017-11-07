sim_inf_chrom <-
function (pop_size, initial_heterozygosity, total_runtime, morgan, 
    markers, seed) 
{
    .Call("_junctions_sim_inf_chrom", PACKAGE = "junctions", pop_size, 
        initial_heterozygosity, total_runtime, morgan, markers, 
        seed)
}

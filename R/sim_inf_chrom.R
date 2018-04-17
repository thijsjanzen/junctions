sim_inf_chrom <- function(pop_size = 100,
                          initial_heterozygosity = 0.5,
                          total_runtime = 100,
                          morgan = 1,
                          markers = -1,
                          seed = 42) {
    set.seed(seed)
    .Call("_junctions_sim_inf_chrom", PACKAGE = "junctions",
        pop_size, initial_heterozygosity, total_runtime,
        morgan, markers, seed)
}
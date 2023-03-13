

vx <- junctions::sim_multiple_chrom(pop_size = 1000,
                                    freq_ancestor_1 = 0.5,
                                    total_runtime = 100,
                                    size_in_morgan = c(0.1),
                                    markers = 1000,
                                    num_indiv_sampled = 3,
                                    time_points = c(1, 10, 50, 100, 500, 999),
                                    num_threads = 6,
                                    verbose = TRUE,
                                    record_true_junctions = TRUE)
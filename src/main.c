#include "mc_def.h"
#include "memory.h"

/*
 * Overlap between current and initial (t=0) spin config:
 *   overlap = (1/N) sum_i S_i(t) · S_i(0)
 * Lower overlap indicates config has evolved away from initial state.
 */
static double compute_overlap(struct BindStruct *X, int int_T) {
    int all_i, All_N;
    double sum;

    All_N = X->Def.All_N;
    sum = 0.0;
    for (all_i = 0; all_i < All_N; all_i++) {
        sum += X->Def.sx[int_T][all_i] * X->Def.prev_sx[int_T][all_i] +
               X->Def.sy[int_T][all_i] * X->Def.prev_sy[int_T][all_i] +
               X->Def.sz[int_T][all_i] * X->Def.prev_sz[int_T][all_i];
    }
    return sum / (double)All_N;
}

/* Copy current config to init (t=0) for overlap. Called once at measurement start. */
static void copy_spins_to_init(struct BindStruct *X, int int_T) {
    int All_N = X->Def.All_N;
    memcpy(X->Def.prev_sx[int_T], X->Def.sx[int_T], (size_t)All_N * sizeof(double));
    memcpy(X->Def.prev_sy[int_T], X->Def.sy[int_T], (size_t)All_N * sizeof(double));
    memcpy(X->Def.prev_sz[int_T], X->Def.sz[int_T], (size_t)All_N * sizeof(double));
}

static double magnetization_sq(struct BindStruct *X, int int_T) {
    int all_i, All_N;
    double mx, my, mz;

    mx = 0.0;
    my = 0.0;
    mz = 0.0;
    All_N = X->Def.All_N;

    for (all_i = 0; all_i < All_N; all_i++) {
        mx += X->Def.sx[int_T][all_i];
        my += X->Def.sy[int_T][all_i];
        mz += X->Def.sz[int_T][all_i];
    }

    return (mx * mx + my * my + mz * mz) / ((double)All_N * (double)All_N);
}

/* Order parameter magnitude: |M|/N. For Ising: |sum sz|/N. */
static double magnetization_mag(struct BindStruct *X, int int_T) {
    int all_i, All_N;
    double mx, my, mz;

    mx = 0.0;
    my = 0.0;
    mz = 0.0;
    All_N = X->Def.All_N;

    for (all_i = 0; all_i < All_N; all_i++) {
        mx += X->Def.sx[int_T][all_i];
        my += X->Def.sy[int_T][all_i];
        mz += X->Def.sz[int_T][all_i];
    }

    return sqrt(mx * mx + my * my + mz * mz) / (double)All_N;
}

/* NER overlap: q(t) = (1/N) sum_i S_i(t) · S_i(0) using Ini_sx/Ini_sy/Ini_sz. */
static double compute_ner_overlap(struct BindStruct *X, int int_T) {
    int all_i, All_N;
    double sum;

    All_N = X->Def.All_N;
    sum = 0.0;
    for (all_i = 0; all_i < All_N; all_i++) {
        sum += X->Def.sx[int_T][all_i] * X->Def.Ini_sx[all_i] +
               X->Def.sy[int_T][all_i] * X->Def.Ini_sy[all_i] +
               X->Def.sz[int_T][all_i] * X->Def.Ini_sz[all_i];
    }
    return sum / (double)All_N;
}

/*
 * Accepted CLI forms:
 *   1) ./MC_simple
 *   2) ./MC_simple param.def
 *   3) ./MC_simple param.def lattice.def interaction.def
 */
static int parse_input_files(int argc, char **argv, const char **param_file,
                             const char **lattice_file,
                             const char **interaction_file) {
    *param_file = "param.def";
    *lattice_file = "lattice.def";
    *interaction_file = "interaction.def";

    if (argc == 2) {
        *param_file = argv[1];
        return 0;
    }

    if (argc == 4) {
        *param_file = argv[1];
        *lattice_file = argv[2];
        *interaction_file = argv[3];
        return 0;
    }

    if (argc == 1)
        return 0;

    return -1;
}

#define MPI_SEED_STRIDE 1000003UL

int main(int argc, char **argv) {
    struct MCMainCalStruct X;
    struct SimpleWorkArrays W;
    dsfmt_t dsfmt;
    const char *param_file;
    const char *lattice_file;
    const char *interaction_file;
    const char *seed_env;
    char *endptr;
    FILE *fp;

    int num_temp, All_N, ni_max;
    int Burn_in, Total_Step, Sample;
    int int_samp, int_T, step;
    int n_exchange_accept, n_exchange_try;
    long total_exchange_accept, total_exchange_try;
    unsigned long base_seed;
    unsigned long seed;
    double Ini_T, Delta_T;
    int rc;

#ifndef NO_MPI
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

    memset(&X, 0, sizeof(X));
    memset(&W, 0, sizeof(W));
    rc = 1;

    if (parse_input_files(argc, argv, &param_file, &lattice_file,
                          &interaction_file) != 0) {
        fprintf(stderr,
                "Usage: %s [param.def] or %s [param.def lattice.def "
                "interaction.def]\n",
                argv[0], argv[0]);
#ifndef NO_MPI
        MPI_Finalize();
#endif
        return 1;
    }

#ifndef NO_MPI
    X.Bind.Def.myrank = myrank;
    X.Bind.Def.total_proc = nprocs;
#else
    X.Bind.Def.myrank = MASTER;
    X.Bind.Def.total_proc = 1;
#endif

    /* Read model/simulation settings from input files. */
    if (read_param(param_file, &X.Bind.Def) != 0)
        return 1;
    if (read_lattice(lattice_file, &X.Bind.Def) != 0)
        return 1;
    if (read_interaction(interaction_file, &X.Bind.Def) < 0)
        return 1;

    num_temp = X.Bind.Def.num_temp;
    All_N = X.Bind.Def.All_N;
    ni_max = X.Bind.Def.ni_max;
    Burn_in = X.Bind.Def.Burn_in;
    Total_Step = X.Bind.Def.Total_Step;
    Sample = X.Bind.Def.Sample;
    Ini_T = X.Bind.Def.Ini_T;
    Delta_T = X.Bind.Def.Delta_T;

    if (num_temp <= 0 || All_N <= 0 || ni_max <= 0) {
        fprintf(stderr, "Error: invalid lattice/temperature setting\n");
        goto cleanup;
    }
    if (Burn_in < 0 || Total_Step <= 0 || Sample <= 0) {
        fprintf(
            stderr,
            "Error: Burn_in >= 0, Total_Step > 0, Sample > 0 are required\n");
        goto cleanup;
    }

    if (allocate_minimal_mc_arrays(&X) != 0) {
        fprintf(stderr, "Error: memory allocation failed\n");
        goto cleanup;
    }

    if (allocate_work_arrays(num_temp, &W) != 0) {
        fprintf(stderr, "Error: memory allocation failed\n");
        goto cleanup;
    }

    /* Build neighbor table (Nr, J) from interaction file. Required before
     * initial() or initial_ner(). */
    if (build_lattice(interaction_file, &X.Bind) != 0) {
        fprintf(stderr, "Error: build_lattice failed\n");
        goto cleanup;
    }

    /* NER mode: use single temperature for validation. */
    if (X.Bind.Def.enable_ner) {
        num_temp = 1;
        X.Bind.Def.num_temp = 1;
    }

    base_seed = 11456UL;
    seed_env = getenv("MC_SIMPLE_SEED");
    if (seed_env != NULL && seed_env[0] != '\0') {
        base_seed = strtoul(seed_env, &endptr, 10);
        if (endptr == seed_env || *endptr != '\0') {
            if (X.Bind.Def.myrank == MASTER) {
                fprintf(stderr,
                        "Warning: invalid MC_SIMPLE_SEED='%s'; using default %lu\n",
                        seed_env, base_seed);
            }
            base_seed = 11456UL;
        }
    }

    total_exchange_accept = 0;
    total_exchange_try = 0;
    if (X.Bind.Def.myrank == MASTER) {
        printf("MC_simple stage1 (%s, %s)\n",
#ifndef NO_MPI
               "MPI",
#else
               "no MPI",
#endif
               X.Bind.Def.enable_ner ? "NER mode"
               : (X.Bind.Def.enable_exchange ? "with exchange MC" : "no exchange MC"));
        printf("L=(%d,%d,%d) orb=%d spin_dim=%d All_N=%d ni_max=%d\n",
               X.Bind.Def.L_x, X.Bind.Def.L_y, X.Bind.Def.L_z, X.Bind.Def.orb_num,
               X.Bind.Def.spin_dim, All_N, ni_max);
        printf("Burn_in=%d Total_Step=%d Sample=%d num_temp=%d\n", Burn_in,
               Total_Step, Sample, num_temp);
    }

    if (X.Bind.Def.enable_ner) {
        /*
         * NER (Non-Equilibrium Relaxation) mode:
         * No burn-in. Quench from controlled initial state, track m(t), m2(t),
         * e(t), q(t) as time series. Average over samples, output to NER_result.dat.
         */
        double *ner_m, *ner_m2, *ner_e, *ner_q;

        ner_m = (double *)calloc((size_t)Total_Step, sizeof(double));
        ner_m2 = (double *)calloc((size_t)Total_Step, sizeof(double));
        ner_e = (double *)calloc((size_t)Total_Step, sizeof(double));
        ner_q = (double *)calloc((size_t)Total_Step, sizeof(double));
        if (ner_m == NULL || ner_m2 == NULL || ner_e == NULL || ner_q == NULL) {
            fprintf(stderr, "Error: NER accumulator allocation failed\n");
            free(ner_m);
            free(ner_m2);
            free(ner_e);
            free(ner_q);
            goto cleanup;
        }

        for (int_samp = 0; int_samp < Sample; int_samp++) {
#ifndef NO_MPI
            seed = base_seed + MPI_SEED_STRIDE * (unsigned long)X.Bind.Def.myrank
                   + 8945UL * (unsigned long)int_samp;
            if (int_samp == 0 && X.Bind.Def.myrank == MASTER) {
                int r;
                for (r = 0; r < X.Bind.Def.total_proc; r++) {
                    unsigned long s = base_seed + MPI_SEED_STRIDE * (unsigned long)r;
                    printf("(rank=%d, seed=%lu)\n", r, s);
                }
            }
#else
            seed = base_seed + 8945UL * (unsigned long)int_samp;
#endif
            dsfmt_init_gen_rand(&dsfmt, seed);

            initial_ner(&dsfmt, &(X.Bind));

            memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));

            /* Time-series loop: MC update then measure at each step. */
            for (step = 0; step < Total_Step; step++) {
                X.Bind.Def.int_T = 0;
                X.Bind.Phys.T = Ini_T;
                MC(&dsfmt, &(X.Bind));

                ner_m[step] += magnetization_mag(&(X.Bind), 0);
                ner_m2[step] += magnetization_sq(&(X.Bind), 0);
                ner_e[step] += X.Bind.Phys.Energy[0] / (double)All_N;
                ner_q[step] += compute_ner_overlap(&(X.Bind), 0);
            }
            if (X.Bind.Def.myrank == MASTER) {
                printf("NER sample %d/%d done (seed=%lu)\n", int_samp + 1, Sample, seed);
            }
        }

#ifndef NO_MPI
        {
            int n_total = Sample * X.Bind.Def.total_proc;
            double inv_total = 1.0 / (double)n_total;
            MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : ner_m,
                       ner_m, Total_Step, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : ner_m2,
                       ner_m2, Total_Step, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : ner_e,
                       ner_e, Total_Step, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : ner_q,
                       ner_q, Total_Step, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            if (X.Bind.Def.myrank == MASTER) {
                for (step = 0; step < Total_Step; step++) {
                    ner_m[step] *= inv_total;
                    ner_m2[step] *= inv_total;
                    ner_e[step] *= inv_total;
                    ner_q[step] *= inv_total;
                }
            }
        }
#endif

        fp = NULL;
        if (X.Bind.Def.myrank == MASTER) {
            fp = fopen("NER_result.dat", "w");
        }
        if (fp != NULL) {
#ifdef NO_MPI
            {
                double inv_sample = 1.0 / (double)Sample;
                fprintf(fp, "# t  T  m  m2  e_per_site  q\n");
                for (step = 0; step < Total_Step; step++) {
                    fprintf(fp, "%d  %.8f  %.10f  %.10f  %.10f  %.10f\n",
                            step, Ini_T,
                            ner_m[step] * inv_sample,
                            ner_m2[step] * inv_sample,
                            ner_e[step] * inv_sample,
                            ner_q[step] * inv_sample);
                }
            }
#else
            /* MPI: already divided by Sample*total_proc in reduce block */
            fprintf(fp, "# t  T  m  m2  e_per_site  q\n");
            for (step = 0; step < Total_Step; step++) {
                fprintf(fp, "%d  %.8f  %.10f  %.10f  %.10f  %.10f\n",
                        step, Ini_T,
                        ner_m[step], ner_m2[step], ner_e[step], ner_q[step]);
            }
#endif
            fclose(fp);
        } else if (X.Bind.Def.myrank == MASTER) {
            fprintf(stderr, "Warning: cannot open NER_result.dat for write\n");
        }

        free(ner_m);
        free(ner_m2);
        free(ner_e);
        free(ner_q);

        rc = 0;
        goto cleanup;
    }

    /*
     * Equilibrium mode: outer loop over independent samples (different RNG seeds).
     * We average observables over these samples at the end.
     */
    for (int_samp = 0; int_samp < Sample; int_samp++) {
#ifndef NO_MPI
        seed = base_seed + MPI_SEED_STRIDE * (unsigned long)X.Bind.Def.myrank
               + 8945UL * (unsigned long)int_samp;
        if (int_samp == 0 && X.Bind.Def.myrank == MASTER) {
            int r;
            for (r = 0; r < X.Bind.Def.total_proc; r++) {
                unsigned long s = base_seed + MPI_SEED_STRIDE * (unsigned long)r;
                printf("(rank=%d, seed=%lu)\n", r, s);
            }
        }
#else
        seed = base_seed + 8945UL * (unsigned long)int_samp;
#endif
        dsfmt_init_gen_rand(&dsfmt, seed);

        /* Initialize spins and neighbor-based effective fields. */
        initial(&dsfmt, &(X.Bind));

        memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));

        /* Burn-in: evolve the system without collecting measurements. */
        for (step = 0; step < Burn_in; step++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                X.Bind.Def.int_T = int_T;
                X.Bind.Phys.T = Ini_T + Delta_T * (double)int_T;
                MC(&dsfmt, &(X.Bind));
            }
            /* Exchange MC: attempt swap between adjacent temperatures */
            if (X.Bind.Def.enable_exchange) {
                for (int_T = 0; int_T < num_temp - 1; int_T++) {
                    (void)attempt_exchange(&dsfmt, &(X.Bind), int_T);
                }
            }
        }

        memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));
        memset(W.sample_E, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_E2, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_M2, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_overlap, 0, (size_t)num_temp * sizeof(double));
        n_exchange_accept = 0;
        n_exchange_try = 0;

        /* Store t=0 config for overlap = S(t)·S(0) (after burn-in, before measurement) */
        for (int_T = 0; int_T < num_temp; int_T++) {
            copy_spins_to_init(&(X.Bind), int_T);
        }

        /* Measurement phase: update then record E/N, E^2, M^2, overlap each step. */
        for (step = 0; step < Total_Step; step++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                X.Bind.Def.int_T = int_T;
                X.Bind.Phys.T = Ini_T + Delta_T * (double)int_T;
                MC(&dsfmt, &(X.Bind));
            }
            /* Exchange MC: attempt swap between adjacent temperatures */
            if (X.Bind.Def.enable_exchange) {
                for (int_T = 0; int_T < num_temp - 1; int_T++) {
                    n_exchange_try++;
                    n_exchange_accept +=
                        attempt_exchange(&dsfmt, &(X.Bind), int_T);
                }
            }
            for (int_T = 0; int_T < num_temp; int_T++) {
                double e_total = X.Bind.Phys.Energy[int_T];
                W.sample_E[int_T] += e_total / (double)All_N;
                W.sample_E2[int_T] += e_total * e_total;
                W.sample_M2[int_T] += magnetization_sq(&(X.Bind), int_T);
                W.sample_overlap[int_T] += compute_overlap(&(X.Bind), int_T);
            }
        }

        /*
         * Convert step sums to per-sample averages.
         * Specific heat per site:
         *   C = ( <E^2> - <E>^2 ) / (N * T^2)
         */
        for (int_T = 0; int_T < num_temp; int_T++) {
            double T_cur = Ini_T + Delta_T * (double)int_T;
            double e_avg = W.sample_E[int_T] / (double)Total_Step;
            double e2_avg = W.sample_E2[int_T] / (double)Total_Step;
            double m2_avg = W.sample_M2[int_T] / (double)Total_Step;
            double a_avg = (double)X.Bind.Phys.ratio_1[int_T] /
                           ((double)All_N * (double)Total_Step);
            double e_total_avg = e_avg * (double)All_N;
            double c_avg = 0.0;
            if (T_cur > 0.0) {
                c_avg = (e2_avg - e_total_avg * e_total_avg) /
                        ((double)All_N * T_cur * T_cur);
            }
            W.accum_E[int_T] += e_avg;
            W.accum_C[int_T] += c_avg;
            W.accum_M2[int_T] += m2_avg;
            W.accum_A[int_T] += a_avg;
            W.accum_overlap[int_T] += W.sample_overlap[int_T] / (double)Total_Step;
        }
        total_exchange_accept += n_exchange_accept;
        total_exchange_try += n_exchange_try;
        if (X.Bind.Def.myrank == MASTER) {
            printf("sample %d/%d done (seed=%lu)\n", int_samp + 1, Sample, seed);
        }
    }

#ifndef NO_MPI
    {
        int n_total = Sample * X.Bind.Def.total_proc;
        double inv_total = 1.0 / (double)n_total;
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : W.accum_E,
                   W.accum_E, num_temp, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : W.accum_C,
                   W.accum_C, num_temp, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : W.accum_M2,
                   W.accum_M2, num_temp, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : W.accum_A,
                   W.accum_A, num_temp, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : W.accum_overlap,
                   W.accum_overlap, num_temp, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : &total_exchange_accept,
                   &total_exchange_accept, 1, MPI_LONG, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(X.Bind.Def.myrank == MASTER ? MPI_IN_PLACE : &total_exchange_try,
                   &total_exchange_try, 1, MPI_LONG, MPI_SUM, MASTER, MPI_COMM_WORLD);
        if (X.Bind.Def.myrank == MASTER) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                W.accum_E[int_T] *= inv_total;
                W.accum_C[int_T] *= inv_total;
                W.accum_M2[int_T] *= inv_total;
                W.accum_A[int_T] *= inv_total;
                W.accum_overlap[int_T] *= inv_total;
            }
        }
    }
#endif

    fp = NULL;
    if (X.Bind.Def.myrank == MASTER) {
        fp = fopen("MC_simple_result.dat", "w");
    }
    if (fp != NULL) {
        fprintf(fp, "# T  E_per_site  C_per_site  M2  acceptance  overlap\n");
        if (total_exchange_try > 0) {
            fprintf(fp, "# exchange_acceptance = %.6f\n",
                    (double)total_exchange_accept / (double)total_exchange_try);
        }
    } else if (X.Bind.Def.myrank == MASTER) {
        fprintf(stderr,
                "Warning: cannot open MC_simple_result.dat for write\n");
    }

    if (X.Bind.Def.myrank == MASTER) {
        printf("\n# T  E_per_site  C_per_site  M2  acceptance  overlap\n");
        if (total_exchange_try > 0) {
            printf("# exchange_acceptance = %.6f\n",
                   (double)total_exchange_accept / (double)total_exchange_try);
        }
        /* Final average: for MPI, already divided by Sample*total_proc in reduce */
        for (int_T = 0; int_T < num_temp; int_T++) {
            double T = Ini_T + Delta_T * (double)int_T;
#ifdef NO_MPI
            double E = W.accum_E[int_T] / (double)Sample;
            double C = W.accum_C[int_T] / (double)Sample;
            double M2 = W.accum_M2[int_T] / (double)Sample;
            double A = W.accum_A[int_T] / (double)Sample;
            double Ov = W.accum_overlap[int_T] / (double)Sample;
#else
            double E = W.accum_E[int_T];
            double C = W.accum_C[int_T];
            double M2 = W.accum_M2[int_T];
            double A = W.accum_A[int_T];
            double Ov = W.accum_overlap[int_T];
#endif
            printf("%.8f  %.10f  %.10f  %.10f  %.10f  %.10f\n", T, E, C, M2, A, Ov);
            if (fp != NULL) {
                fprintf(fp, "%.8f  %.10f  %.10f  %.10f  %.10f  %.10f\n", T, E, C, M2, A, Ov);
            }
        }
    }

    if (fp != NULL)
        fclose(fp);

    rc = 0;

cleanup:
    free_work_arrays(&W);
    free_minimal_mc_arrays(&X);
#ifndef NO_MPI
    MPI_Finalize();
#endif
    return rc;
}

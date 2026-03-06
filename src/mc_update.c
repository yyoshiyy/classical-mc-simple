#include "mc_def.h"

/*
 * Attempt Exchange MC swap between adjacent replicas int_T and int_T+1.
 * Acceptance: P = min(1, exp((1/T_i - 1/T_j) * (E_i - E_j)))
 * On accept: swap sx/sy/sz, env_sx/env_sy/env_sz, Energy, ratio_1.
 * Returns 1 if accepted, 0 otherwise.
 */
int attempt_exchange(dsfmt_t *dsfmt, struct BindStruct *X, int int_T) {
    int j = int_T + 1;
    int All_N = X->Def.All_N;
    double T_i, T_j, E_i, E_j;
    double prob, r;
    int all_i;
    double tmp;

    if (j >= X->Def.num_temp)
        return 0;

    T_i = X->Def.Ini_T + X->Def.Delta_T * (double)int_T;
    T_j = X->Def.Ini_T + X->Def.Delta_T * (double)j;
    E_i = X->Phys.Energy[int_T];
    E_j = X->Phys.Energy[j];

    prob = exp((1.0 / T_i - 1.0 / T_j) * (E_i - E_j));
    if (prob > 1.0)
        prob = 1.0;

    r = dsfmt_genrand_close_open(dsfmt);
    if (r >= prob)
        return 0;

    /* Swap spins */
    for (all_i = 0; all_i < All_N; all_i++) {
        tmp = X->Def.sx[int_T][all_i];
        X->Def.sx[int_T][all_i] = X->Def.sx[j][all_i];
        X->Def.sx[j][all_i] = tmp;

        tmp = X->Def.sy[int_T][all_i];
        X->Def.sy[int_T][all_i] = X->Def.sy[j][all_i];
        X->Def.sy[j][all_i] = tmp;

        tmp = X->Def.sz[int_T][all_i];
        X->Def.sz[int_T][all_i] = X->Def.sz[j][all_i];
        X->Def.sz[j][all_i] = tmp;
    }

    /* Swap effective fields */
    for (all_i = 0; all_i < All_N; all_i++) {
        tmp = X->Def.env_sx[int_T][all_i];
        X->Def.env_sx[int_T][all_i] = X->Def.env_sx[j][all_i];
        X->Def.env_sx[j][all_i] = tmp;

        tmp = X->Def.env_sy[int_T][all_i];
        X->Def.env_sy[int_T][all_i] = X->Def.env_sy[j][all_i];
        X->Def.env_sy[j][all_i] = tmp;

        tmp = X->Def.env_sz[int_T][all_i];
        X->Def.env_sz[int_T][all_i] = X->Def.env_sz[j][all_i];
        X->Def.env_sz[j][all_i] = tmp;
    }

    /* Swap Energy and ratio_1 */
    tmp = X->Phys.Energy[int_T];
    X->Phys.Energy[int_T] = X->Phys.Energy[j];
    X->Phys.Energy[j] = tmp;

    {
        int rtmp = X->Phys.ratio_1[int_T];
        X->Phys.ratio_1[int_T] = X->Phys.ratio_1[j];
        X->Phys.ratio_1[j] = rtmp;
    }

    return 1;
}

/*
 * One full Metropolis sweep for a single temperature slot (int_T):
 *   - Visit all spins once.
 *   - Propose a new spin orientation (distribution depends on spin_dim).
 *   - Accept/reject with Metropolis rule.
 *
 * Performance note:
 *   The code stores local effective fields env_* and updates them
 *   incrementally after each accepted move. This avoids recalculating
 *   neighbor sums from scratch for every proposal.
 */
void MC(dsfmt_t *dsfmt, struct BindStruct *X) {
    int all_i, all_j;
    int n_i;
    int int_T, ni_max, spin_dim;
    double delta_E, T, lambda;
    double tmp, tmp_2, theta;
    double x, y, z;
    double tmp_prob, rand_prob;

    int_T = X->Def.int_T;
    T = X->Phys.T;
    lambda = X->Def.lambda;
    ni_max = X->Def.ni_max;
    spin_dim = X->Def.spin_dim;

    for (all_i = 0; all_i < X->Def.All_N; all_i++) {
        /*
         * Proposal step:
         *   Heisenberg -> random unit vector on S^2,
         *   XY         -> random unit vector on S^1,
         *   Ising      -> deterministic flip sz -> -sz.
         */
        if (spin_dim == 3) {
            /* Heisenberg: Marsaglia method for uniform S^2 sampling */
            x = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
            y = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
            tmp = x * x + y * y;
            while (tmp > 1.0) {
                x = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                y = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                tmp = x * x + y * y;
            }
            tmp_2 = sqrt(1 - tmp);
            x = 2 * x * tmp_2;
            y = 2 * y * tmp_2;
            z = 1 - 2 * tmp;
        } else if (spin_dim == 2) {
            /* XY: uniform angle on S^1 */
            theta = 2.0 * PI * dsfmt_genrand_close_open(dsfmt);
            x = cos(theta);
            y = sin(theta);
            z = 0.0;
        } else {
            /* Ising: flip sz -> -sz */
            x = 0.0;
            y = 0.0;
            z = -X->Def.sz[int_T][all_i];
        }

        /*
         * Energy difference from old spin S_old to proposed S_new:
         *   dE = (S_new - S_old) · env_i
         * where env_i already includes J-neighbor sum (+field terms).
         */
        delta_E = (x - X->Def.sx[int_T][all_i]) * X->Def.env_sx[int_T][all_i] +
                  (y - X->Def.sy[int_T][all_i]) * X->Def.env_sy[int_T][all_i] +
                  (z - X->Def.sz[int_T][all_i]) * X->Def.env_sz[int_T][all_i];

        if (delta_E < 0.0) {
            /* Always accept downhill moves. */
            for (n_i = 0; n_i < ni_max; n_i++) {
                all_j = X->Def.Nr[n_i][all_i];
                /*
                 * Incremental env update on neighbors:
                 * env_j += J_ji * (S_new - S_old)
                 * (with anisotropy factor lambda on z component).
                 */
                X->Def.env_sx[int_T][all_j] +=
                    X->Def.J[n_i][all_i] * (x - X->Def.sx[int_T][all_i]);
                X->Def.env_sy[int_T][all_j] +=
                    X->Def.J[n_i][all_i] * (y - X->Def.sy[int_T][all_i]);
                X->Def.env_sz[int_T][all_j] += lambda * X->Def.J[n_i][all_i] *
                                               (z - X->Def.sz[int_T][all_i]);
            }
            X->Def.sx[int_T][all_i] = x;
            X->Def.sy[int_T][all_i] = y;
            X->Def.sz[int_T][all_i] = z;
            X->Phys.Energy[int_T] += delta_E;
            X->Phys.ratio_1[int_T] += 1;
        } else {
            /* Metropolis acceptance for uphill moves. */
            tmp_prob = exp(-delta_E / T);
            rand_prob = dsfmt_genrand_close_open(dsfmt);

            if (tmp_prob > rand_prob) {
                /* Accepted uphill move: apply the same incremental updates. */
                for (n_i = 0; n_i < ni_max; n_i++) {
                    all_j = X->Def.Nr[n_i][all_i];
                    X->Def.env_sx[int_T][all_j] +=
                        X->Def.J[n_i][all_i] * (x - X->Def.sx[int_T][all_i]);
                    X->Def.env_sy[int_T][all_j] +=
                        X->Def.J[n_i][all_i] * (y - X->Def.sy[int_T][all_i]);
                    X->Def.env_sz[int_T][all_j] +=
                        lambda * X->Def.J[n_i][all_i] *
                        (z - X->Def.sz[int_T][all_i]);
                }
                X->Def.sx[int_T][all_i] = x;
                X->Def.sy[int_T][all_i] = y;
                X->Def.sz[int_T][all_i] = z;
                X->Phys.Energy[int_T] += delta_E;
                X->Phys.ratio_1[int_T] += 1;
            }
        }
    }
}

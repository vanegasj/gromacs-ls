/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/* When calculating RF or Ewald interactions we calculate the electrostatic
 * forces and energies on excluded atom pairs here in the non-bonded loops.
 */
#if defined CHECK_EXCLS && (defined CALC_COULOMB || defined LJ_EWALD)
#define EXCL_FORCES
#endif

{
    int cj;
#ifdef ENERGY_GROUPS
    int egp_cj;
#endif
    int i;

    cj = l_cj[cjind].cj;

#ifdef ENERGY_GROUPS
    egp_cj = nbat->energrp[cj];
#endif
    // locals needs the values of rvdw and rcut for the impulse corrections
    real rvdw = std::sqrt(rvdw2);
    real rcut = std::sqrt(rcut2);
    // end locals
    for (i = 0; i < UNROLLI; i++)
    {
        int ai;
        int type_i_off;
        int j;

        ai = ci*UNROLLI + i;

        type_i_off = type[ai]*ntype2;

        for (j = 0; j < UNROLLJ; j++)
        {
            int             aj;
            real            dx, dy, dz;
            real            rsq, rinv;
            real            rinvsq, rinvsix;
            real            c6, c12;
            real            FrLJ6 = 0, FrLJ12 = 0, frLJ = 0;
            real            VLJ gmx_unused;
            // Locals constants needed to compute the elasticity tensor from vdw and coul interactions
            real            phi_coul = 0, kappa_coul = 0, phi_lj = 0, kappa_lj = 0, phi_coul_ic = 0, kappa_coul_ic = 0, phi_lj_ic = 0, kappa_lj_ic = 0, skipmask_rvdw;
            real            dfw = locals_grid->GetSpacing(), dfwsq = dfw*dfw, deltavdw = 0, deltavdwsq = 0, deltacoul = 0, deltacoulsq = 0, rinvl = 0, rinvsql = 0, rinvsixl = 0;
            bool            bCoulEwald = EEL_FULL(ic->eeltype); // Locals needs to know whether we are using plain coul or Ewald for elasticity calculations
            bool            bCoulCut = (ic->eeltype == eelCUT);
#if defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
            real            r, rsw;
#endif

#ifdef CALC_COULOMB
            real qq;
            real fcoul;
#ifdef CALC_COUL_TAB
            real rs, frac;
            int  ri;
            real fexcl;
#endif
#ifdef CALC_ENERGIES
            real vcoul;
#endif
#endif
            real fscal;
            real fx, fy, fz;

            /* A multiply mask used to zero an interaction
             * when either the distance cutoff is exceeded, or
             * (if appropriate) the i and j indices are
             * unsuitable for this kind of inner loop. */
            real skipmask;

#ifdef CHECK_EXCLS
            /* A multiply mask used to zero an interaction
             * when that interaction should be excluded
             * (e.g. because of bonding). */
            int interact;

            interact = ((l_cj[cjind].excl>>(i*UNROLLI + j)) & 1);
#ifndef EXCL_FORCES
            skipmask = interact;
#else
            skipmask = (cj == ci_sh && j <= i) ? 0.0 : 1.0;
#endif
#else
#define interact 1.0
            skipmask = 1.0;
#endif

            // cppcheck-suppress unreadVariable
            VLJ = 0;

            aj = cj*UNROLLJ + j;

            dx  = xi[i*XI_STRIDE+XX] - x[aj*X_STRIDE+XX];
            dy  = xi[i*XI_STRIDE+YY] - x[aj*X_STRIDE+YY];
            dz  = xi[i*XI_STRIDE+ZZ] - x[aj*X_STRIDE+ZZ];

            rsq = dx*dx + dy*dy + dz*dz;

            /* Prepare to enforce the cut-off. */
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            /* 9 flops for r^2 + cut-off check */

            // Ensure the distances do not fall below the limit where r^-12 overflows.
            // This should never happen for normal interactions.
            rsq = std::max(rsq, NBNXN_MIN_RSQ);

#ifdef COUNT_PAIRS
            npair++;
#endif

            rinv = gmx::invsqrt(rsq);
            rinvl = rinv; // locals rinv
            /* 5 flops for invsqrt */

            /* Partially enforce the cut-off (and perhaps
             * exclusions) to avoid possible overflow of
             * rinvsix when computing LJ, and/or overflowing
             * the Coulomb table during lookup. */
            rinv = rinv * skipmask;

            rinvsq  = rinv*rinv;

#ifdef HALF_LJ
            if (i < UNROLLI/2)
#endif
            {
                c6      = nbfp[type_i_off+type[aj]*2  ];
                c12     = nbfp[type_i_off+type[aj]*2+1];
                if (locals_grid != NULL)
                {
                    if (locals_grid->GetContribType() == mds_cou)
                    {
                        c6 = 0.0;
                        c12 = 0.0;
                    }
                }

#if defined LJ_CUT || defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6   = c6*rinvsix;
                FrLJ12  = c12*rinvsix*rinvsix;
                frLJ    = FrLJ12 - FrLJ6;
                /* 7 flops for r^-2 + LJ force */

                // begin locals compute the kappa and phi values needed to compute elasticity for plain LJ interactions
                phi_lj = -c12*rinvsix*rinvsix*rinv + c6*rinvsix*rinv;
                kappa_lj = 13*c12*rinvsix*rinvsix*rinvsq - 7*c6*rinvsix*rinvsq;

                // locals impulsive correction for particles near the cutoff
                deltavdw = (rvdw-1.0/rinvl);
                deltavdwsq = deltavdw*deltavdw;
                if (deltavdwsq < dfwsq)
                {
                    rinvsql = rinvl*rinvl;
                    rinvsixl = rinvsql*rinvsql*rinvsql;
                    phi_lj_ic = (c12*rinvsixl*rinvsixl/12.0 - c6*rinvsixl/6.0)/dfw;
                    kappa_lj_ic = (-c12*rinvsixl*rinvsixl*rinvl + c6*rinvsixl*rinvl)/dfw;
                }
                // end locals

#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                VLJ     = (FrLJ12 + c12*ic->repulsion_shift.cpot)/12 -
                    (FrLJ6 + c6*ic->dispersion_shift.cpot)/6;
                /* 7 flops for LJ energy */
#endif
#endif

#if defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
                /* Force or potential switching from ic->rvdw_switch */
                r       = rsq*rinv;
                rsw     = r - ic->rvdw_switch;
                rsw     = (rsw >= 0.0 ? rsw : 0.0);
#endif
#ifdef LJ_FORCE_SWITCH
                frLJ   +=
                    -c6*(ic->dispersion_shift.c2 + ic->dispersion_shift.c3*rsw)*rsw*rsw*r
                    + c12*(ic->repulsion_shift.c2 + ic->repulsion_shift.c3*rsw)*rsw*rsw*r;
                // begin locals adjustments to phi_lj and kappa_lj
                phi_lj += c6*(ic->dispersion_shift.c2 + ic->dispersion_shift.c3*rsw)*rsw*rsw
                 - c12*(ic->repulsion_shift.c2 + ic->repulsion_shift.c3*rsw)*rsw*rsw;
                kappa_lj += c6*(2.0*ic->dispersion_shift.c2 + 3.0*ic->dispersion_shift.c3*rsw)*rsw*rsw
                            - c12*(2.0*ic->repulsion_shift.c2 + 3.0*ic->repulsion_shift.c3*rsw)*rsw*rsw;
                // end locals

#if defined CALC_ENERGIES
                VLJ    +=
                    -c6*(-ic->dispersion_shift.c2/3 - ic->dispersion_shift.c3/4*rsw)*rsw*rsw*rsw
                    + c12*(-ic->repulsion_shift.c2/3 - ic->repulsion_shift.c3/4*rsw)*rsw*rsw*rsw;
#endif
#endif

#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                /* Masking should be done after force switching,
                 * but before potential switching.
                 */
                /* Need to zero the interaction if there should be exclusion. */
                VLJ     = VLJ * interact;
#endif

#ifdef LJ_POT_SWITCH
                {
                    real sw, dsw;

                    sw    = 1.0 + (swV3 + (swV4+ swV5*rsw)*rsw)*rsw*rsw*rsw;
                    dsw   = (swF2 + (swF3 + swF4*rsw)*rsw)*rsw*rsw;

                    // locals adjustments to phi_lj and kappa_lj
                    //Code for phi and kappa for switching
                    real ddsw =  (6*swV3 + (12*swV4 + 20*swV5*rsw)*rsw)*rsw;
                    real phi_lj0 = phi_lj;
                    real kappa_lj0 = kappa_lj;

                    phi_lj = phi_lj0*sw + dsw*VLJ;
                    kappa_lj = kappa_lj0*sw + 2.0*phi_lj0*dsw + ddsw*VLJ;
                    //end phi and kappa section

                    frLJ  = frLJ*sw - r*VLJ*dsw;
                    VLJ  *= sw;

                }
#endif

#ifdef LJ_EWALD
                {
                    real            c6grid, rinvsix_nm, cr2, expmcr2, poly;
#ifdef CALC_ENERGIES
                    real            sh_mask;
#endif

#ifdef LJ_EWALD_COMB_GEOM
                    c6grid       = ljc[type[ai]*2]*ljc[type[aj]*2];
                    if (locals_grid != NULL)
                    {
                        if (locals_grid->GetContribType() == mds_cou)
                        {
                            c6grid = 0.0;
                        }
                    }
#elif defined LJ_EWALD_COMB_LB
                    {
                        real sigma, sigma2, epsilon;

                        /* These sigma and epsilon are scaled to give 6*C6 */
                        sigma   = ljc[type[ai]*2] + ljc[type[aj]*2];
                        epsilon = ljc[type[ai]*2+1]*ljc[type[aj]*2+1];

                        sigma2  = sigma*sigma;
                        c6grid  = epsilon*sigma2*sigma2*sigma2;
                        if (locals_grid != NULL)
                        {
                            if (locals_grid->GetContribType() == mds_cou)
                            {
                                c6grid = 0.0;
                            }
                        }
                    }
#else
#error "No LJ Ewald combination rule defined"
#endif

#ifdef CHECK_EXCLS
                    /* Recalculate rinvsix without exclusion mask */
                    rinvsix_nm   = rinvsq*rinvsq*rinvsq;
#else
                    rinvsix_nm   = rinvsix;
#endif
                    cr2          = lje_coeff2*rsq;
#if GMX_DOUBLE
                    expmcr2      = exp(-cr2);
#else
                    expmcr2      = expf(-cr2);
#endif
                    poly         = 1 + cr2 + 0.5*cr2*cr2;

                    /* Subtract the grid force from the total LJ force */
                    frLJ        += c6grid*(rinvsix_nm - expmcr2*(rinvsix_nm*poly + lje_coeff6_6));
#ifdef CALC_ENERGIES
                    /* Shift should only be applied to real LJ pairs */
                    sh_mask      = lje_vc*interact;

                    VLJ         += c6grid/6*(rinvsix_nm*(1 - expmcr2*poly) + sh_mask);
#endif
                }
#endif          /* LJ_EWALD */

// Always do a rvdw != rcoul check for local stress calculations
//#ifdef VDW_CUTOFF_CHECK
                /* Mask for VdW cut-off shorter than Coulomb cut-off */
                {
                    skipmask_rvdw = (rsq < rvdw2);
                    frLJ         *= skipmask_rvdw;
                    phi_lj       *= skipmask_rvdw;
                    kappa_lj     *= skipmask_rvdw;
#ifdef CALC_ENERGIES
                    VLJ *= skipmask_rvdw;
#endif
                }
//#else
//#if defined CALC_ENERGIES
                /* Need to zero the interaction if r >= rcut */
//                VLJ     = VLJ * skipmask;
                /* 1 more flop for LJ energy */
//#endif
//#endif          /* VDW_CUTOFF_CHECK */


#ifdef CALC_ENERGIES
#ifdef ENERGY_GROUPS
                Vvdw[egp_sh_i[i]+((egp_cj>>(nbat->neg_2log*j)) & egp_mask)] += VLJ;
#else
                Vvdw_ci += VLJ;
                /* 1 flop for LJ energy addition */
#endif
#endif
            }

#ifdef CALC_COULOMB
            /* Enforce the cut-off and perhaps exclusions. In
             * those cases, rinv is zero because of skipmask,
             * but fcoul and vcoul will later be non-zero (in
             * both RF and table cases) because of the
             * contributions that do not depend on rinv. These
             * contributions cannot be allowed to accumulate
             * to the force and potential, and the easiest way
             * to do this is to zero the charges in
             * advance. */
            qq = skipmask * qi[i] * q[aj];
            // begin locals
            if (locals_grid != NULL)
            {
                if (locals_grid->GetContribType() == mds_vdw)
                {
                    qq = 0.0;
                }
            }
            // Locals Calculate Elasticity Constants using a plain cutoff when using PME
            if (bCoulEwald)
            {
                phi_coul = -qq*rinvsq*interact;
                kappa_coul = 2*qq*rinvsq*rinv*interact;
                // locals impulsive correction for particles near the cutoff
                deltacoul = (rcut-1.0/rinvl);
                deltacoulsq = deltacoul*deltacoul;
                if (deltacoulsq < dfwsq)
                {
                    rinvsql = rinvl*rinvl;
                    phi_coul_ic = (qq*rinvl)/dfw;
                    kappa_coul_ic = (-qq*rinvsql)/dfw;
                }
                // end locals
            }
#ifdef CALC_COUL_RF
            fcoul  = qq*(interact*rinv*rinvsq - k_rf2);
            /* 4 flops for RF force */
            // Locals Calculate Elasticity Constants for plain and Ewald coulomb potential
            if ((ic->eeltype == eelRF) || (ic->eeltype == eelRF_ZERO))
            {
                // Locals Calculate Elasticity Constants for reaction-field coulomb potential
                phi_coul = qq*(-rinvsq*interact + k_rf2/rinv);
                kappa_coul = qq*(2*rinvsq*rinv*interact + k_rf2);
            }
            else //use plain cutoff electrostatitcs for everything else ....
            {
                phi_coul = -qq*rinvsq*interact;
                kappa_coul = 2*qq*rinvsq*rinv*interact;
                // locals impulsive correction for particles near the cutoff
                deltacoul = (rcut-1.0/rinvl);
                deltacoulsq = deltacoul*deltacoul;
                if ((deltacoulsq < dfwsq) && bCoulCut)
                {
                    rinvsql = rinvl*rinvl;
                    phi_coul_ic = (qq*rinvl)/dfw;
                    kappa_coul_ic = (-qq*rinvsql)/dfw;
                }
                // end locals
            }
#ifdef CALC_ENERGIES
            vcoul  = qq*(interact*rinv + k_rf*rsq - c_rf);
            /* 4 flops for RF energy */
#endif
#endif

#ifdef CALC_COUL_TAB
            rs     = rsq*rinv*ic->tabq_scale;
            ri     = (int)rs;
            frac   = rs - ri;
#if !GMX_DOUBLE
            /* fexcl = F_i + frac * (F_(i+1)-F_i) */
            fexcl  = tab_coul_FDV0[ri*4] + frac*tab_coul_FDV0[ri*4+1];
#else
            /* fexcl = (1-frac) * F_i + frac * F_(i+1) */
            fexcl  = (1 - frac)*tab_coul_F[ri] + frac*tab_coul_F[ri+1];
#endif
            fcoul  = interact*rinvsq - fexcl;
            /* 7 flops for float 1/r-table force */
#ifdef CALC_ENERGIES
#if !GMX_DOUBLE
            vcoul  = qq*(interact*(rinv - ic->sh_ewald)
                         -(tab_coul_FDV0[ri*4+2]
                           -halfsp*frac*(tab_coul_FDV0[ri*4] + fexcl)));
            /* 7 flops for float 1/r-table energy (8 with excls) */
#else
            vcoul  = qq*(interact*(rinv - ic->sh_ewald)
                         -(tab_coul_V[ri]
                           -halfsp*frac*(tab_coul_F[ri] + fexcl)));
#endif
#endif
            fcoul *= qq*rinv;
#endif

#ifdef CALC_ENERGIES
#ifdef ENERGY_GROUPS
            Vc[egp_sh_i[i]+((egp_cj>>(nbat->neg_2log*j)) & egp_mask)] += vcoul;
#else
            Vc_ci += vcoul;
            /* 1 flop for Coulomb energy addition */
#endif
#endif
#endif

#ifdef CALC_COULOMB
#ifdef HALF_LJ
            if (i < UNROLLI/2)
#endif
            {
                fscal = frLJ*rinvsq + fcoul;
                /* 2 flops for scalar LJ+Coulomb force */
            }
#ifdef HALF_LJ
            else
            {
                fscal = fcoul;
            }
#endif
#else
            fscal = frLJ*rinvsq;
#endif
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;

            /* begin stress tensor */
            if (locals_grid != NULL)
            {
                int  lpatIDs[2];
                lpatIDs[0] = xi_id[i]; lpatIDs[1] = x_id[aj];

                // remove the 'far away' particles
                if (lpatIDs[0] != -1 && lpatIDs[1] != -1)
                {
                    int cont_type = locals_grid->GetContribType();
                    if (cont_type == mds_all || cont_type == mds_vdw || cont_type == mds_cou)
                    {
                        real ix = xi[i*XI_STRIDE+XX]; real jx = x[aj*X_STRIDE+XX];
                        real iy = xi[i*XI_STRIDE+YY]; real jy = x[aj*X_STRIDE+YY];
                        real iz = xi[i*XI_STRIDE+ZZ]; real jz = x[aj*X_STRIDE+ZZ];

                        rvec lpR[2], lpF[2];
                        lpR[0][0] = ix; lpR[0][1] = iy; lpR[0][2] = iz;
                        lpR[1][0] = jx; lpR[1][1] = jy; lpR[1][2] = jz;
                        lpF[0][0] = fx;  lpF[0][1] = fy;  lpF[0][2] = fz;
                        lpF[1][0] = -fx; lpF[1][1] = -fy; lpF[1][2] = -fz;
                        if (skipmask > 0)
                        {
                            locals_grid->DistributeInteraction(2, lpR, lpF, lpatIDs);
                            locals_grid->DistributeElasticity(lpR[0], lpR[1], lpR[0], lpR[1], phi_lj, kappa_lj);
#ifdef CALC_COULOMB
                            locals_grid->DistributeElasticity(lpR[0], lpR[1], lpR[0], lpR[1], phi_coul,kappa_coul);
#endif
                        }
#ifdef LJ_CUT
                        //if (deltavdwsq < dfwsq) // uncomment this line to include impulse correction from particles below and above the cutoff
                        if (deltavdwsq < dfwsq && skipmask_rvdw > 0) // uncomment this line to include impulse correction only from particles below the cutoff
                        {
                            locals_grid->DistributeElasticity(lpR[0], lpR[1], lpR[0], lpR[1], -phi_lj_ic, -kappa_lj_ic);
                            if (ic->vdwtype == evdwCUT && ic->vdw_modifier == eintmodNONE)
                            {
                                real lj_ic = phi_lj_ic*rinvl;
                                lpF[0][0] = lj_ic*dx;  lpF[0][1] = lj_ic*dy;  lpF[0][2] = lj_ic*dz;
                                lpF[1][0] = -lpF[0][0]; lpF[1][1] = -lpF[0][1]; lpF[1][2] = -lpF[0][2];
                                locals_grid->DistributeInteraction(2, lpR, lpF, lpatIDs);
                            }
                        }
#endif
#ifdef CALC_COULOMB
                        //if (deltacoulsq < dfwsq) // uncomment this line to include impulse correction from particles below and above the cutoff
                        if (deltacoulsq < dfwsq && skipmask > 0) // uncomment this line to include impulse correction only from particles below the cutoff
                        {
                            locals_grid->DistributeElasticity(lpR[0], lpR[1], lpR[0], lpR[1], -phi_coul_ic, -kappa_coul_ic);
                            if ((bCoulCut && ic->coulomb_modifier == eintmodNONE) || bCoulEwald)
                            {
                                real coul_ic = phi_coul_ic*rinvl;
                                lpF[0][0] = coul_ic*dx;  lpF[0][1] = coul_ic*dy;  lpF[0][2] = coul_ic*dz;
                                printf("icx = %e, icy = %e, icz = %e\n", lpF[0][0],  lpF[0][1],  lpF[0][2]);
                                lpF[1][0] = -lpF[0][0]; lpF[1][1] = -lpF[0][1]; lpF[1][2] = -lpF[0][2];
                                locals_grid->DistributeInteraction(2, lpR, lpF, lpatIDs);
                            }
                        }
#endif

                        //if (ic->eeltype == eelCUT)
                        //printf("ai = %d, aj = %d, fx = %e, fy = %e, fz = %e, px = %e, py = %e, pz = %e\n", xi_id[i], x_id[aj], fx, fy, fz, phi_lj*dx*rinv, phi_lj*dy*rinv, phi_lj*dz*rinv);
                        //printf("r = %e, phi_ic = %e, kappa_ic = %e\n", 1/rinv, phi_lj_ic, kappa_lj_ic);
                    }
                }
            }
            /* end stress tensor */

            /* Increment i-atom force */
            fi[i*FI_STRIDE+XX] += fx;
            fi[i*FI_STRIDE+YY] += fy;
            fi[i*FI_STRIDE+ZZ] += fz;
            /* Decrement j-atom force */
            f[aj*F_STRIDE+XX]  -= fx;
            f[aj*F_STRIDE+YY]  -= fy;
            f[aj*F_STRIDE+ZZ]  -= fz;
            /* 9 flops for force addition */
        }
    }
}

#undef interact
#undef EXCL_FORCES

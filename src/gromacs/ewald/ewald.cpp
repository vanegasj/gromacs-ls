/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file contains function definitions necessary for
 * computing energies and forces for the plain-Ewald long-ranged part,
 * and the correction for overall system charge for all Ewald-family
 * methods.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "ewald.h"

#include <math.h>
#include <stdio.h>

#include <cstdlib>

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

struct gmx_ewald_tab_t
{
    int        nx, ny, nz, kmax;
    real       fs; //ir->fourier_spacing
    cvec     **eir;
    t_complex *tab_xy, *tab_qxyz;
};

void init_ewald_tab(struct gmx_ewald_tab_t **et, const t_inputrec *ir, FILE *fp)
{
    snew(*et, 1);
    if (fp)
    {
        fprintf(fp, "Will do ordinary reciprocal space Ewald sum.\n");
    }

    (*et)->nx       = ir->nkx+1;
    (*et)->ny       = ir->nky+1;
    (*et)->nz       = ir->nkz+1;
    (*et)->kmax     = std::max((*et)->nx, std::max((*et)->ny, (*et)->nz));
    (*et)->fs       = ir->fourier_spacing;
    (*et)->eir      = NULL;
    (*et)->tab_xy   = NULL;
    (*et)->tab_qxyz = NULL;
}

//! Calculates wave vectors.
static void calc_lll(const rvec box, rvec lll)
{
    lll[XX] = 2.0*M_PI/box[XX];
    lll[YY] = 2.0*M_PI/box[YY];
    lll[ZZ] = 2.0*M_PI/box[ZZ];
}

//! Make tables for the structure factor parts
static void tabulateStructureFactors(int natom, rvec x[], int kmax, cvec **eir, rvec lll)
{
    int  i, j, m;

    if (kmax < 1)
    {
        printf("Go away! kmax = %d\n", kmax);
        exit(1);
    }

    for (i = 0; (i < natom); i++)
    {
        for (m = 0; (m < 3); m++)
        {
            eir[0][i][m].re = 1;
            eir[0][i][m].im = 0;
        }

        for (m = 0; (m < 3); m++)
        {
            eir[1][i][m].re = cos(x[i][m]*lll[m]);
            eir[1][i][m].im = sin(x[i][m]*lll[m]);
        }
        for (j = 2; (j < kmax); j++)
        {
            for (m = 0; (m < 3); m++)
            {
                eir[j][i][m] = cmul(eir[j-1][i][m], eir[1][i][m]);
            }
        }
    }
}

real do_ewald(t_inputrec *ir, t_blocka *excl,
              rvec x[],        rvec f[],
              real chargeA[],  real chargeB[],
              rvec box, matrix full_box,
              t_commrec *cr,   int natoms,
              matrix lrvir,    real ewaldcoeff,
              real lambda,     real *dvdlambda,
              struct gmx_ewald_tab_t *et,
              mds::StressGrid *locals_grid)
{
    real     factor     = -1.0/(4*ewaldcoeff*ewaldcoeff);
    real     scaleRecip = 4.0*M_PI/(box[XX]*box[YY]*box[ZZ])*ONE_4PI_EPS0/ir->epsilon_r; /* 1/(Vol*e0) */
    real    *charge, energy_AB[2], energy;
    rvec     lll;
    int      lowiy, lowiz, ix, iy, iz, n, q;
    real     tmp, cs, ss, ak, akv, mx, my, mz, m2, scale;
    gmx_bool bFreeEnergy;
    // local stress
    rvec fij, mvec, rij, fk;
    int ai, aj, ai1, ai2, ai3, d;
    real fscal, qq;
    // local stress

    if (cr != NULL)
    {
        if (PAR(cr))
        {
            gmx_fatal(FARGS, "No parallel Ewald. Use PME instead.\n");
        }
    }


    if (!et->eir) /* allocate if we need to */
    {
        snew(et->eir, et->kmax);
        for (n = 0; n < et->kmax; n++)
        {
            snew(et->eir[n], natoms);
        }
        snew(et->tab_xy, natoms);
        snew(et->tab_qxyz, natoms);
    }

    bFreeEnergy = (ir->efep != efepNO);

    clear_mat(lrvir);

    calc_lll(box, lll);
    tabulateStructureFactors(natoms, x, et->kmax, et->eir, lll);

    /* tracking forces on each particle as a test, need arrays for that */
    //rvec fsum_pairs_ls[natoms] = {{0.0,0.0,0.0,},};
    //rvec fsum_pairs_ew[natoms] = {{0.0,0.0,0.0,},};

    /* tracking pressure, need a 3x3 matrix */
    rvec Pxxz = {0.0,0.0,0.0,};
    rvec Pyxz = {0.0,0.0,0.0,};
    rvec Pzxz = {0.0,0.0,0.0,};

    for (q = 0; q < (bFreeEnergy ? 2 : 1); q++)
    {
        if (!bFreeEnergy)
        {
            charge = chargeA;
            scale  = 1.0;
        }
        else if (q == 0)
        {
            charge = chargeA;
            scale  = 1.0 - lambda;
        }
        else
        {
            charge = chargeB;
            scale  = lambda;
        }
        lowiy        = 0;
        lowiz        = 1;
        energy_AB[q] = 0;

        /* begin stress tensor */
        /*if (locals_grid != NULL &&
                (locals_grid->GetContribType() == mds_all ||
                 locals_grid->GetContribType() == mds_ewal))
        {
            printf("/ncalled the ewald function and am now summing\n");
            //printf("\nLx: %6.2f, Ly: %6.2f, Lz: %6.2f\n",box[0],box[1],box[2]);
            real ang_av = 0.0, ang;
            real mx[2*et->nx - 1];
            real my[2*et->ny - 1];
            real mz[2*et->nz - 1];
            real ak[2*et->nx - 1][2*et->ny - 1][2*et->nz - 1];
            int ixp, iyp, izp, counter=0;

            for (ix = -et->nx + 1; ix < et->nx; ix++)
            {
                ixp = ix + et->nx - 1;
                mx[ixp] = ix*lll[XX];
                for (iy = -et->ny + 1; iy < et->ny; iy++)
                {
                    iyp = iy + et->ny - 1;
                    my[iyp] = iy*lll[YY];
                    for (iz = -et->nz + 1; iz < et->nz; iz++)
                    {
                        izp = iz + et->nz - 1;
                        mz[izp]  = iz*lll[ZZ];
                        m2 = mx[ixp]*mx[ixp] + my[iyp]*my[iyp] + mz[izp]*mz[izp];
                        ak[ixp][iyp][izp] = scale*exp(m2*factor)/m2;
                        //printf("ak = %6.4e\n",ak[ix][iy][iz]);
                    }
                }
            }
            //natoms = 500;
            for (ai = 0; ai < natoms; ai++)
            {
                ai1 = excl->index[ai];
                ai2 = excl->index[ai+1];

                //printf("adding interaction of particle %i and", ai);
                for (aj = ai+1; aj < natoms; aj++)
                {
                    bool exclude = false;
                    for (ai3 = ai1; ai3 < ai2; ++ai3)
                    {
                        if (excl->a[ai3] == aj)
                        {
                            exclude = true;
                            break;
                        }
                    }

                    // need to exclude bonded pairs here
                    if (!exclude)
                    {
                        //printf(" %i:\n", aj);
                        qq = charge[ai]*charge[aj]*scaleRecip;

                        // place them in the box, position wise
                        rvec xx; rvec xy;
                        for (d = 0; d < DIM; d++)
                        {
                            // position xx
                            xx[d] = x[ai][d];
                            if (xx[d] < -box[d])
                                xx[d] += box[d];
                            if (xx[d] >= box[d])
                                xx[d] -= box[d];

                            // position xy
                            xy[d] = x[aj][d];
                            if (xy[d] < -box[d])
                                xy[d] += box[d];
                            if (xy[d] >= box[d])
                                xy[d] -= box[d];
                        }

                        // now take the minimal distance
                        clear_rvec(rij);
                        rvec_sub(xx, xy, rij);
                        for (d = 0; d < DIM; d++)
                        {
                            if (rij[d] > 0.5*box[d])
                            {
                                rij[d] -= box[d];
                            }
                            else if (rij[d] <= -0.5*box[d])
                            {
                                rij[d] += box[d];
                            }
                        }

                        // clear the old force calculate the new
                        clear_rvec(fij);
                        for (ix = -et->nx + 1; ix < et->nx; ix++)
                        {
                            ixp = ix + et->nx - 1;
                            //mx = 0.5*ix*lll[XX];
                            for (iy = -et->ny + 1; iy < et->ny; iy++)
                            {
                                iyp = iy + et->ny - 1;
                                //my = 0.5*iy*lll[YY];
                                for (iz = -et->nz + 1; iz < et->nz; iz++)
                                {
                                    izp = iz + et->nz - 1;
                                    //mz  = 0.5*iz*lll[ZZ];
                                    //m2  = mx*mx + my*my + mz*mz;
                                    if (abs(ix) + abs(iy) + abs(iz) > 0)
                                    {
                                        //ak = 2*exp(m2*factor)/m2;
                                        mvec[0] = mx[ixp]; mvec[1] = my[iyp]; mvec[2] = mz[izp];
                                        fscal = ak[ixp][iyp][izp]*sin(iprod(mvec, rij));
                                        svmul(fscal, mvec, fk);
                                        rvec_inc(fij, fk);
                                        //printf(" ix: %d, iy: %d, iz: %d, m2: %5.2e, ak: %5.2e, fscal: %5.2e\n",ix,iy,iz,m2,ak,fscal);
                                    }
                                }
                            }
                        }
                        svmul(qq, fij, fij);
                        ang = 360.0*acos(abs(iprod(fij,rij))/(sqrt(iprod(fij,fij))*sqrt(iprod(rij,rij))))/(2.0*M_PI);
                        ang_av += ang;
                        counter += 1;
                        printf(" %i on %i: %5.4f degrees; ",ai, aj, ang);
                        printf(" ri: %6.4f, %6.4f, %6.4f; rj: %6.4f, %6.4f, %6.4f; rij: %6.4f, %6.4f, %6.4f, |rij| = %6.4f\n",
                                xx[0],xx[1],xx[2],xy[0],xy[1],xy[2],rij[0]/box[XX],rij[1]/box[YY],rij[2]/box[ZZ], norm(rij));
                        rvec fij_unit, rij_unit;
                        unitv(fij, fij_unit);
                        unitv(rij, rij_unit);
                        //printf("\nai = %d, aj = %d", ai, aj);
                        //printf("rij = %6.4f %6.4f %6.4f, fij = %6.4f %6.4f %6.4f; |rij| = %6.4f; |k| = %6.4f; |rij|/|k| = %6.4f\n", rij_unit[0], rij_unit[1], rij_unit[2], fij_unit[0], fij_unit[1], fij_unit[2], norm(rij), norm(mvec), norm(rij)/norm(mvec));
                        // call mdstress library here
                        int lpatIDs[2];
                        lpatIDs[0] = ai; lpatIDs[1] = aj;

                        rvec lpR[2], lpF[2];
                        lpR[0][0] = x[ai][0]; lpR[0][1] = x[ai][1]; lpR[0][2] = x[ai][2];
                        lpR[1][0] = x[ai][0] - rij[0]; lpR[1][1] = x[ai][1] - rij[1]; lpR[1][2] = x[ai][2] - rij[2];
                        lpF[0][0] = fij[0];  lpF[0][1] = fij[1];  lpF[0][2] = fij[2];
                        lpF[1][0] = -fij[0]; lpF[1][1] = -fij[1]; lpF[1][2] = -fij[2];
                        locals_grid->DistributeInteraction(2, lpR, lpF, lpatIDs);

                        rvec_inc(fsum_pairs_ls[ai], fij);
                        rvec_dec(fsum_pairs_ls[aj], fij);

                        Pxxz[0] += 0.5*rij[0]*fij[0]; Pxxz[1] += 0.5*rij[0]*fij[1]; Pxxz[2] += 0.5*rij[0]*fij[2];
                        Pyxz[0] += 0.5*rij[1]*fij[0]; Pyxz[1] += 0.5*rij[1]*fij[1]; Pyxz[2] += 0.5*rij[1]*fij[2];
                        Pzxz[0] += 0.5*rij[2]*fij[0]; Pzxz[1] += 0.5*rij[2]*fij[1]; Pzxz[2] += 0.5*rij[2]*fij[2];
                    }
                }
                printf("\n");
                //printf("\n ai = %d ; Fls = %6.4f %6.4f %6.4f", ai, fcumul[XX], fcumul[YY], fcumul[ZZ]);
                //printf("ai = %d\n", ai);
            }
            printf("Ang_av = %6.4f\n\n",ang_av/counter);
        }*/

        /* end stress tensor */

        for (ix = 0; ix < et->nx; ix++)
        {
            mx = ix*lll[XX];
            for (iy = lowiy; iy < et->ny; iy++)
            {
                my = iy*lll[YY];
                if (iy >= 0)
                {
                    for (n = 0; n < natoms; n++)
                    {
                        et->tab_xy[n] = cmul(et->eir[ix][n][XX], et->eir[iy][n][YY]);
                    }
                }
                else
                {
                    for (n = 0; n < natoms; n++)
                    {
                        et->tab_xy[n] = cmul(et->eir[ix][n][XX], conjugate(et->eir[-iy][n][YY]));
                    }
                }
                for (iz = lowiz; iz < et->nz; iz++)
                {
                    mz  = iz*lll[ZZ];
                    m2  = mx*mx+my*my+mz*mz;
                    ak  = exp(m2*factor)/m2;
                    akv = 2.0*ak*(1.0/m2-factor);
                    if (iz >= 0)
                    {
                        for (n = 0; n < natoms; n++)
                        {
                            et->tab_qxyz[n] = rcmul(charge[n], cmul(et->tab_xy[n],
                                                                    et->eir[iz][n][ZZ]));
                        }
                    }
                    else
                    {
                        for (n = 0; n < natoms; n++)
                        {
                            et->tab_qxyz[n] = rcmul(charge[n], cmul(et->tab_xy[n],
                                                                    conjugate(et->eir[-iz][n][ZZ])));
                        }
                    }

                    cs = ss = 0;
                    for (n = 0; n < natoms; n++)
                    {
                        cs += et->tab_qxyz[n].re;
                        ss += et->tab_qxyz[n].im;
                    }
                    energy_AB[q]  += ak*(cs*cs+ss*ss);
                    tmp            = scale*akv*(cs*cs+ss*ss);
                    lrvir[XX][XX] -= tmp*mx*mx;
                    lrvir[XX][YY] -= tmp*mx*my;
                    lrvir[XX][ZZ] -= tmp*mx*mz;
                    lrvir[YY][YY] -= tmp*my*my;
                    lrvir[YY][ZZ] -= tmp*my*mz;
                    lrvir[ZZ][ZZ] -= tmp*mz*mz;
                    for (n = 0; n < natoms; n++)
                    {
                        /*tmp=scale*ak*(cs*tab_qxyz[n].im-ss*tab_qxyz[n].re);*/
                        tmp       = scale*ak*(cs*et->tab_qxyz[n].im-ss*et->tab_qxyz[n].re);
                        f[n][XX] += tmp*mx*2*scaleRecip;
                        f[n][YY] += tmp*my*2*scaleRecip;
                        f[n][ZZ] += tmp*mz*2*scaleRecip;
                        
                        //fsum_pairs_ew[n][XX] += tmp*mx*2*scaleRecip;
                        //fsum_pairs_ew[n][YY] += tmp*my*2*scaleRecip;
                        //fsum_pairs_ew[n][ZZ] += tmp*mz*2*scaleRecip;
#if 0
                        f[n][XX] += tmp*mx;
                        f[n][YY] += tmp*my;
                        f[n][ZZ] += tmp*mz;
#endif
                    }
                    lowiz = 1-et->nz;
                }
                lowiy = 1-et->ny;
            }
        }
    }

    /* print pressure tensor */
    /*printf("ls_pressure:\n");
    printf("Pxx: %18.12e, Pyx: %18.12e, Pzx: %18.12e\n", Pxxz[0], Pxxz[1], Pxxz[2]);
    printf("Pxy: %18.12e, Pyy: %18.12e, Pzy: %18.12e\n", Pyxz[0], Pyxz[1], Pyxz[2]);
    printf("Pxz: %18.12e, Pyz: %18.12e, Pzz: %18.12e\n\n", Pzxz[0], Pzxz[1], Pzxz[2]);
    for (ai = 0; ai < natoms; ai++)
    {
        printf("ls forces on %03i: %18.12e, %18.12e, %18.12e\n", ai, fsum_pairs_ls[ai][XX],
                fsum_pairs_ls[ai][YY], fsum_pairs_ls[ai][ZZ]);
        printf("ew forces on %03i: %18.12e, %18.12e, %18.12e\n\n", ai, fsum_pairs_ew[ai][XX],
                fsum_pairs_ew[ai][YY], fsum_pairs_ew[ai][ZZ]);
    }*/

    if (!bFreeEnergy)
    {
        energy = energy_AB[0];
    }
    else
    {
        energy      = (1.0 - lambda)*energy_AB[0] + lambda*energy_AB[1];
        *dvdlambda += scaleRecip*(energy_AB[1] - energy_AB[0]);
    }

    lrvir[XX][XX] = -0.5*scaleRecip*(lrvir[XX][XX]+energy);
    lrvir[XX][YY] = -0.5*scaleRecip*(lrvir[XX][YY]);
    lrvir[XX][ZZ] = -0.5*scaleRecip*(lrvir[XX][ZZ]);
    lrvir[YY][YY] = -0.5*scaleRecip*(lrvir[YY][YY]+energy);
    lrvir[YY][ZZ] = -0.5*scaleRecip*(lrvir[YY][ZZ]);
    lrvir[ZZ][ZZ] = -0.5*scaleRecip*(lrvir[ZZ][ZZ]+energy);

    lrvir[YY][XX] = lrvir[XX][YY];
    lrvir[ZZ][XX] = lrvir[XX][ZZ];
    lrvir[ZZ][YY] = lrvir[YY][ZZ];
    
    /*printf("ew_pressure:\n");
    printf("Pxx: %18.12e, Pyx: %18.12e, Pzx: %18.12e\n", lrvir[XX][XX], lrvir[YY][XX], lrvir[ZZ][XX]);
    printf("Pxy: %18.12e, Pyy: %18.12e, Pzy: %18.12e\n", lrvir[XX][YY], lrvir[YY][YY], lrvir[ZZ][YY]);
    printf("Pxz: %18.12e, Pyz: %18.12e, Pzz: %18.12e\n\n", lrvir[XX][ZZ], lrvir[YY][ZZ], lrvir[ZZ][ZZ]);*/

    energy *= scaleRecip;

    return energy;
}

real ewald_charge_correction(t_commrec *cr, t_forcerec *fr, real lambda,
                             matrix box,
                             real *dvdlambda, tensor vir)

{
    real vol, fac, qs2A, qs2B, vc, enercorr;
    int  d;

    if (MASTER(cr))
    {
        /* Apply charge correction */
        vol = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];

        fac = M_PI*ONE_4PI_EPS0/(fr->epsilon_r*2.0*vol*vol*gmx::square(fr->ewaldcoeff_q));

        qs2A = fr->qsum[0]*fr->qsum[0];
        qs2B = fr->qsum[1]*fr->qsum[1];

        vc = (qs2A*(1 - lambda) + qs2B*lambda)*fac;

        enercorr = -vol*vc;

        *dvdlambda += -vol*(qs2B - qs2A)*fac;

        for (d = 0; d < DIM; d++)
        {
            vir[d][d] += vc;
        }

        if (debug)
        {
            fprintf(debug, "Total charge correction: Vcharge=%g\n", enercorr);
        }
    }
    else
    {
        enercorr = 0;
    }

    return enercorr;
}

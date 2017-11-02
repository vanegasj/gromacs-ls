/*
 * gmx_density3D - Custom tool to calculate the 3D density from a
 * trajectory. Created by Alejandro Torres-Sanchez and
 * Juan M. Vanegas. Send any comments/questions/bugs to
 * juan.m.vanegas@gmail.com
 * Nov. 2017
 *
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

#include "gmxpre.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

// Modulo operation
int modulo(int a, int b)
{
    int ret;
    ret = a % b;
    if(ret < 0)
    ret+=b;
    return ret;
}

// Finds the indices on the grid for a given set of coordinates
void grid_givecoord(matrix box, int nx, int ny, int nz,rvec pt, int *i, int *j, int *k, int pbc)
{

    real rxx,ryx,ryy,rzx,rzy,rzz;

    /*Looks up indices for grid:
    fractional indices are invbox * coordinate;
    grid indices are then nx*f_ind[XX], etc.
    */
    real tmp=1.0/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
    rxx = box[YY][YY]*box[ZZ][ZZ]*tmp;
    ryy = box[XX][XX]*box[ZZ][ZZ]*tmp;
    rzz = box[XX][XX]*box[YY][YY]*tmp;

    *i = nx * pt[XX] * rxx;
    *j = ny * pt[YY] * ryy;
    *k = nz * pt[ZZ] * rzz;

    if(pbc) {
    *i = modulo(*i,nx);
    *j = modulo(*j,ny);
    *k = modulo(*k,nz);
    }

    if(pt[0] < 0) *i =*i-1;
    if(pt[1] < 0) *j =*j-1;
    if(pt[2] < 0) *k =*k-1;
}

void calc_density3D(const char *fn, int **index, int gnx[], double **slDensity,
                  t_topology *top, int ePBC, int nr_grps, int *gridx, int *gridy, int *gridz,
                  real *gridsp, matrix avgbox, const gmx_output_env_t *oenv)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  double invvol,mass,invgridsp,dummy1,dummy2,factor;
  int natoms;            /* nr. atoms in trj */
  int tx,ty,tz,sx,sy,sz,nx,ny,nz,d,dd;
  t_trxstatus *status;  
  int  **slCount,         /* nr. of atoms in one slice for a group */
      i,ii,iii,j,jj,jjj,k,kk,kkk,l,n,               /* loop indices */
      teller = 0,
      nr_frames = 0;     /* number of frames */
  real t,x,y,z;
  char *buf;             /* for tmp. keeping atomname */
  gmx_rmpbc_t  gpbc=NULL;
  rvec gridspv,pt;

  if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");

  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr);
  for (d = 0; d < 3; d++)
      for (dd = 0; dd < 3; dd++)
        avgbox[d][dd] = 0.0;
  if (*gridx == 0)
    *gridx = box[XX][XX]/(*gridsp);
  if (*gridy == 0)
    *gridy = box[YY][YY]/(*gridsp);
  if (*gridz == 0)
    *gridz = box[ZZ][ZZ]/(*gridsp);

  printf("gridx = %d; gridy = %d, gridz = %d \n", *gridx, *gridy, *gridz);
  nx = *gridx;
  ny = *gridy;
  nz = *gridz;

  snew(*slDensity, nx*ny*nz);
  /*********** Start processing trajectory ***********/
  do {
    gmx_rmpbc(gpbc,natoms,box,x0);

    for (d = 0; d < 3; d++)
      for (dd = 0; dd < 3; dd++)
        avgbox[d][dd] += box[d][dd];

    invvol = nx*ny*nz/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
    teller++;

    // Define the grid spacing in each direction
    gridspv[XX] = box[XX][XX]/nx;
    gridspv[YY] = box[YY][YY]/ny;
    gridspv[ZZ] = box[ZZ][ZZ]/nz;
    invgridsp = 1/(gridspv[XX]*gridspv[YY]*gridspv[ZZ]);

    for (n = 0; n < nr_grps; n++) {
      for (l = 0; l < gnx[n]; l++) {   /* loop over all atoms in index file */
        pt[0] = x0[index[n][l]][0];
        pt[1] = x0[index[n][l]][1];
        pt[2] = x0[index[n][l]][2];

        // Get the coordinates of the point in the grid
        grid_givecoord(box,nx,ny,nz,pt,&ii,&jj,&kk,0);
        iii=ii;
        jjj=jj;
        kkk=kk;

        // Spread it
        for(i=1;i>=-1;i-=2){
          iii+=i;
          dummy1 = i * invgridsp * (pt[0]-(ii+0.5*(1-i))*gridspv[0]);
          for(j=1;j>=-1;j-=2){
            jjj+=j;
            dummy2 = dummy1 * j * (pt[1]-(jj+0.5*(1-j))*gridspv[1]);
            for(k=1;k>=-1;k-=2){
              kkk+=k;
              factor = dummy2 * k * (pt[2]-(kk+0.5*(1-k))*gridspv[2]);
              mass = top->atoms.atom[index[n][l]].m;
              (*slDensity)[modulo(iii,nx)*nz*ny+modulo(jjj,ny)*nz+modulo(kkk,nz)] += factor*mass*invvol;
            }
          }
        }
      }
    }
    nr_frames++;
  } while (read_next_x(oenv,status,&t,x0,box));
  gmx_rmpbc_done(gpbc);

  for (d = 0; d < 3; d++)
    for (dd = 0; dd < 3; dd++)
      avgbox[d][dd] /= nr_frames;

  /*********** done with status file **********/
  close_trj(status);

  fprintf(stderr,"\nRead %d frames from trajectory. Calculating density\n", nr_frames);
  for (i = 0; i < nx*ny*nz; i++)
    (*slDensity)[i] /= nr_frames;

  sfree(x0);  /* free memory used by coordinate array */
}

int gmx_density3D(int argc,char *argv[])
{
  const char *desc[] = {
    "Computes the density in 3D over the simulation volume and stores it "
    "in a grid with a given spacing. Use the tensortools program (part "
    "of mdstresslib) to convert the output into a netcdf file "
    "that can be read by a number of programs such as paraview. "
    "This is a custom utility part of GROMACS-LS. Please send any "
    "questions/comments/bugs to Juan M. Vanegas at juan.m.vanegas@gmail.com "
  };

  gmx_output_env_t *oenv;
  static const char *dens_opt[] = 
    { NULL, "mass", "number", "charge", NULL };
  static int  ngrps   = 1;       /* nr. of groups              */
  real gridsp;           /* grid spacing in nm         */
  int gridx,gridy,gridz;  /* grid size in x,y, and z    */
  gridx = 0;
  gridy = 0;
  gridz = 0;
  gridsp = 0.05;
  t_pargs pa[] = {
    { "-gridsp",  FALSE, etREAL, {&gridsp},
      "Divide the box in 3D using a grid with a spacing given by gridsp in nm" },
    { "-gridx",  FALSE, etINT, {&gridx},
      "Override the size of the grid in the x direction" },
    { "-gridy",  FALSE, etINT, {&gridy},
      "Override the size of the grid in the y direction" },
    { "-gridz",  FALSE, etINT, {&gridz},
      "Override the size of the grid in the z direction" },
    { "-dens",    FALSE, etENUM, {dens_opt},
      "Density"}
  };

  const char        *bugs[] = {"",};

  double *density;      /* density array         */
  
  char **grpname;        /* groupnames                 */
  int  *ngx;             /* sizes of groups            */
  t_topology *top;       /* topology                   */
  int  ePBC, d, dd;
  int   **index;     /* indices for all groups     */
  int  i,bdouble;
  bdouble = 3;
  double dbox[3][3];
  matrix avgbox;
  FILE *outfp;

  t_filenm  fnm[] = {    /* files for g_density        */
    { efTRX, "-f", NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffOPTRD },
    { efTPR, NULL, NULL,  ffREAD },
    { efDAT,"-o","density",ffWRITE }
  };

#define NFILE asize(fnm)

  if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
  {
    return 0;
  }

  top = read_top(ftp2fn(efTPR,NFILE,fnm),&ePBC);     /* read topology file */

  if (dens_opt[0][0] == 'm') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = top->atoms.atom[i].m*AMU/(NANO*NANO*NANO);
  } else if (dens_opt[0][0] == 'n') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = 1;
  } else if (dens_opt[0][0] == 'c') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = top->atoms.atom[i].q;
  }

  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);

  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,ngx,index,grpname);
  calc_density3D(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, top, ePBC, ngrps, &gridx, &gridy, &gridz, &gridsp, avgbox, oenv);
  
  for (d = 0; d < 3; d++)
      for (dd = 0; dd < 3; dd++)
        dbox[d][dd] = avgbox[d][dd];
  outfp = fopen(ftp2fn(efDAT,NFILE,fnm),"w");
  fwrite(&bdouble,sizeof(int),1,outfp);
  fwrite(dbox,sizeof(double),9,outfp);
  fwrite(&gridx,sizeof(int),1,outfp);
  fwrite(&gridy,sizeof(int),1,outfp);
  fwrite(&gridz,sizeof(int),1,outfp);
  fwrite(density,sizeof(double),gridx*gridy*gridz,outfp);
  fclose(outfp);
  do_view(oenv,opt2fn("-o",NFILE,fnm), "-nxy");       /* view xvgr file */
  return 0;
}

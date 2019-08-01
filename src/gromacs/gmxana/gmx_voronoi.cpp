/* 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

/* Tool to compute voronoi tesselations in 2D and 3D using 
   the voro++ library http://math.lbl.gov/voro++/
   If you have any comments/questions/bugs please contact
   Juan M. Vanegas at juan.m.vanegas@gmail.com
   Jan. 2014 - Version 8 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gmxpre.h"

#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/tpxio.h"
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

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/forcerec.h"

#include "voro++/voro++.hh"

#define xNone 0
#define xTs 1
#define xAll 2
#define xAll_Ts 3

typedef int atom_id;

void put_atom_in_box(matrix box,rvec x)
{
  int i,m,d;
  for (d=0; d<DIM; d++){
    if (x[d] < 0.0)
      x[d] += box[d][d];
    else if (x[d] >= box[d][d])
      x[d] -= box[d][d];
  }
}

using namespace voro;

void avg_areas(int nr_grps, int grpsize[], double areas[], double grp_areas[], double grp_areas_var[])
{
  int n,t,i;
  double sum, diff;

  t = 0;
  for (n = 0; n < nr_grps; n++){
    sum = 0.0;
    for (i = 0; i < grpsize[n]; i++){
      sum += areas[t];
      t += 1;
    }
    grp_areas[n] = sum/grpsize[n];
  }
  t = 0;
  for (n = 0; n < nr_grps; n++){
    sum = 0.0;
    for (i = 0; i < grpsize[n]; i++){
      diff = areas[t] - grp_areas[n];
      sum += diff*diff;
      t += 1;
    }
    grp_areas_var[n] = sum/grpsize[n];
  }
}

/* Compute the number of atoms in each molecule and
   the number of molecules when using the -mol option
   and given and index file with all the atoms in the group */
void get_mol_size(atom_id **index, int nr_grps, int grpsize[], char *grpname[], t_block *mols, int nr_atoms_mol[], int nr_mols[])
{
  int i,n,mol;
  i = 0;
  for (n=0; n<nr_grps; n++){
    while (mols->index[i] < index[n][0]){
      i += 1;
    }
    if (i == 0 && n > 0)
      printf("\nSomething is wrong with your index file or the molecules are not whole!\n");
    mol = 0;
    while (mols->index[mol+i] < index[n][grpsize[n]-1])
      mol += 1;
    nr_atoms_mol[n] = grpsize[n]/mol;
    if (nr_atoms_mol[n] == 1)
      nr_mols[n] = grpsize[n];
    else
      nr_mols[n] = mol;
    //printf("\nGroup Size = %f \n",nr_atoms_mol[n]);
    printf("\nIdentified %d molecules in group %s with %d atoms per molecule\n\n",nr_mols[n],grpname[n],nr_atoms_mol[n]);
  }
}

/* Since the voro++ library can only compute 3D tesselations,
   for 2d tesselations we set the coordinate of the normal axis for all points to
   the same value = 0.5, since the voronoi box in the axis direction
   extends from 0 to 1.0 */
void copy_rvec2d(int axis, matrix gmx_box, rvec a, rvec b)
{
  copy_rvec(a,b);
  //put_atom_in_box(gmx_box, b);
  if (axis == XX)
    b[XX] = 0.5;
  else if (axis == YY)
    b[YY] = 0.5;
  else
    b[ZZ] = 0.5;
}

void copy_rvec3d(matrix gmx_box, rvec a, rvec b)
{
  copy_rvec(a,b);
  //put_atom_in_box(gmx_box,b);
}

/* Initialize the parameters needed to start the voronoi tesselation */
void init_voro_par2d(int nr_ndx, int axis, matrix gmx_box, double vor_box[], int gridn[])
{
  /* gridn is an int vector that tells the voro++ library how many boxes
     to divide the system for the tesselation. Optimally you want 3 to 8 particles
     per box  */
  double gridf; /* gridf = factor to adjust the number of boxes in each dim accordingly to the system size */
  int i;
  for(i=0;i<DIM;i++){
    vor_box[i] = 1.0;
    gridn[i] = 1;
  }
  if (axis == XX){
    vor_box[YY] = gmx_box[YY][YY];
    vor_box[ZZ] = gmx_box[ZZ][ZZ];
    gridf = vor_box[ZZ]/vor_box[YY];
    gridn[1] = sqrt(nr_ndx/(3*gridf));
    gridn[2] = gridn[1]*gridf;
  } else if (axis == YY){
    vor_box[XX] = gmx_box[XX][XX];
    vor_box[ZZ] = gmx_box[ZZ][ZZ];
    gridf = vor_box[ZZ]/vor_box[XX];
    gridn[0] = sqrt(nr_ndx/(3*gridf));
    gridn[2] = gridn[0]*gridf;
  } else{
    vor_box[XX] = gmx_box[XX][XX];
    vor_box[YY] = gmx_box[YY][YY];
    gridf = vor_box[YY]/vor_box[XX];
    gridn[0] = sqrt(nr_ndx/(3*gridf));
    gridn[1] = gridn[0]*gridf;
  }
}

void init_voro_par3d(int nr_ndx, matrix gmx_box, double vor_box[], int gridn[])
{
  /* gridn is an int vector that tells the voro++ library how many boxes
     to divide the system for the tesselation. Optimally you want 3 to 8 particles
     per box  */
  double gfxy,gfxz; /* gfx* = factor to adjust the number of boxes in each dim accordingly to the system size */
  vor_box[XX] = gmx_box[XX][XX];
  vor_box[YY] = gmx_box[YY][YY];
  vor_box[ZZ] = gmx_box[ZZ][ZZ];
  gfxy = vor_box[YY]/vor_box[XX];
  gfxz = vor_box[ZZ]/vor_box[XX];
  gridn[0] = pow(nr_ndx/(3*gfxy*gfxz), 1/3.0);
  gridn[1] = pow(nr_ndx/(3*gfxy*gfxz), 1/3.0)*gfxy;
  gridn[2] = pow(nr_ndx/(3*gfxy*gfxz), 1/3.0)*gfxz;
  //printf("\nGridx %d, gridy %d, gridz %d \n",gridn[0],gridn[1],gridn[2]);
  //printf("\nBox x %f, Box y %f, Box z %f \n",gmx_box[XX][XX],gmx_box[YY][YY],gmx_box[ZZ][ZZ]);
}

/* calculate the COM of each molecule in x0 and save the COM coord into xmol */
void calc_mol_com(int nr_grps, int nr_mols[], int nr_atoms_mol[], atom_id **index,t_atoms *atoms, rvec *x0,rvec *xmol)
{
  int  ai,a0,mol,i,d,n;
  dvec xm;
  double m,mtot;

  mol = 0;
  for (n=0; n<nr_grps; n++){
    a0 = index[n][0];
    for (i=0; i<nr_mols[n]; i++){
      clear_dvec(xm);
      mtot = 0.0;
      for (ai=a0; ai<a0+nr_atoms_mol[n]; ai++){
        m = atoms->atom[ai].m;
        for (d=0; d<DIM;d++)
          xm[d] += m*x0[ai][d];
        mtot += m;
      }
      a0 += nr_atoms_mol[n];
      for (d=0; d<DIM; d++)
        xmol[mol][d] = xm[d]/mtot;
      mol += 1;
    }
  }
}

/* compute the voronoi tesselation and output selected quantities */
void compute_voronoi(const char *fn, atom_id **index, int grpsize[], t_topology *top, int ePBC, int axis, int nr_grps, char *grpname[],
                     gmx_bool bMol, gmx_bool vorPBC, gmx_bool vert, gmx_bool vorAA, gmx_bool vor3d, gmx_bool bRad, real tessRadius, int fullout,
                     t_forcerec *fr,
                     const char *avgfile, const char *tsfile, const char *allfile, const char *vfile, const char *pfile,
                     const gmx_output_env_t *oenv)
{
  rvec *x0;              /* coords without pbc */
  rvec *xmol, px;        /* coords of COM for each molecule; temp position vector */
  matrix gmx_box;        /* system box (3x3) */
  double vor_box[3];     /* box used to initialize voronoi tesselation */
  int natoms;            /* total nr. of atoms in trj */
  int nr_atoms_mol[nr_grps];
  int nr_mols[nr_grps];  /* nr. of mols in each selected group */
  int nr_ndx=0;          /* total number of particles in the tesselation */
  int nr_frames=0;       /* nr. atoms in trj */
  int gridn[3];          /* gridi size for voronoi tesselation (see init_voro_par) */
  double *areas, *frame_areas, *frame_areas_full, *areas_full;
  double grp_areas[nr_grps], grp_areas_var[nr_grps], avg_grp_areas_var[nr_grps], sum;
  double gmx_frame_area, gmx_tot_area = 0.0, vor_tot_area = 0.0;
  double cell_area, vor_frame_area;
  t_trxstatus *status;
  int nr_full, *fgrpsize;
  bool xper=FALSE, yper=FALSE, zper=FALSE;
  int i,j,k,n,t,pid,s;     /* number of frames */
  real tt;
  gmx_rmpbc_t  gpbc=NULL;
  const char *label[2] = {"volume","area"};
  const char *units[2] = {"nm^3","nm^2"};

  if (vor3d)
    s = 0;
  else
    s = 1;

  if (axis < 0 || axis >= DIM)
    gmx_fatal(FARGS,"Invalid axes. Terminating\n");

  if (bMol){
    get_mol_size(index, nr_grps, grpsize, grpname, &top->mols, nr_atoms_mol, nr_mols);
    for (i = 0; i < nr_grps; i++)
      grpsize[i] = nr_mols[i];
  }

  if (vorAA){
    get_mol_size(index, nr_grps, grpsize, grpname, &top->mols, nr_atoms_mol, nr_mols);
    nr_full = 0;
    snew(fgrpsize, nr_grps);
    for (i = 0; i < nr_grps; i++){
      nr_full += nr_mols[i];
      fgrpsize[i] = nr_mols[i];
      //printf("\n\n fgrpsize = %d\n\n",fgrpsize[i]);
    }
    snew(areas_full, nr_full);
    snew(frame_areas_full, nr_full);
    for (i = 0; i < nr_full; i++){
      areas_full[i] = 0.0;
      frame_areas_full[i] = 0.0;
   }
  }

  if ((natoms = read_first_x(oenv,&status,fn,&tt,&x0,gmx_box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");

  /* Compute the VDW radii for all the atoms from the C6 and C12 LJ parameters */
  int AA_index[natoms];
  double radii[natoms];
  if (bRad){
    double c6,c12;
    int ii;
    for (i = 0; i < natoms; i++){
      ii = top->atoms.atom[i].type;
      c6 = C6(fr->nbfp,fr->ntype,ii,ii)/6.0;
      c12 = C12(fr->nbfp,fr->ntype,ii,ii)/12.0;
      if (c6 > 0.0)
        radii[i] = (int)(1000000*pow(c12/c6,1/6.0)/2.0)/1000000.0; // keeping only 6 sig digits for radius to avoid problems with the tesselation
      else
        radii[i] = 0.0;
      //printf("Atom %i radius = %f nm\n",i,radii[i]);
      AA_index[i] = 0;
    }
  }else{
    for (i = 0; i < natoms; i++){
      radii[i] = tessRadius;
      AA_index[i] = 0;
      //printf("\nAtom %i radius = %f nm\n",i,radii[i]);
    }
  }

  for (i = 0; i < nr_grps; i++){
    nr_ndx += grpsize[i];
    avg_grp_areas_var[i] = 0.0;
  }

  if (bMol)
    snew(xmol, nr_ndx);

  snew(areas,nr_ndx);
  snew(frame_areas,nr_ndx);

  for (i = 0; i < nr_ndx; i++){
    areas[i] = 0.0;
    frame_areas[i] = 0.0;
  }

  //gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr);

  if (vorPBC){
    if (axis == XX){
      yper = TRUE;
      zper = TRUE;
    }
    if (axis == YY){
      xper = TRUE;
      zper = TRUE;
    }
    else{
      xper = TRUE;
      yper = TRUE;
    }
  }

  if (vorPBC && vor3d){
    xper = TRUE;
    yper = TRUE;
    zper = TRUE;
  }

  FILE * tsfp;
  FILE * allfp;

  if (fullout == xTs){
    tsfp = fopen(tsfile,"w");
    fprintf(tsfp,"# Average %s of molecules in the selected index groups\n",label[s]);
    fprintf(tsfp,"# All %s units in %s\n",label[s], units[s]);
    fprintf(tsfp,"# Time (ps)\t");
    for (n = 0; n < nr_grps; n++)
      fprintf(tsfp,"Per mol. Avg %s\t+/- S.D.\tTotal %s\t+/- S.D.\t", grpname[n],grpname[n]);
    fprintf(tsfp,"Total voronoi\n");
  } else if (fullout == xAll){
    allfp = fopen(allfile,"w");
  } else if (fullout == xAll_Ts){
    allfp = fopen(allfile,"w");
    tsfp = fopen(tsfile,"w");
    fprintf(tsfp,"# Average %s of molecules in the selected index groups\n",label[s]);
    fprintf(tsfp,"# All %s units in %s\n",label[s], units[s]);
    fprintf(tsfp,"# Time (ps)\t\t");
    for (n = 0; n < nr_grps; n++)
      fprintf(tsfp,"Per mol. Avg %s\t+/- S.D.\tTotal %s\t+/- S.D.\t", grpname[n],grpname[n]);
    fprintf(tsfp,"Total voronoi\n");
  }
  
  /*********** Start processing trajectory ***********/
  int compute_cells = 0;
  int maxcells, nparticles, id;
  double x,y,z,r;
  do {
    if (vorAA)
      nparticles = natoms;
    else
      nparticles = nr_ndx;

    for (i = 0; i < nr_ndx; i++){
      frame_areas[i] = 0.0;
    }
    particle_order vorpo(nparticles);

    /* Compute parameters and initialize voronoi tesselation */
    if (vor3d)
      init_voro_par3d(natoms, gmx_box, vor_box, gridn);
    else
      init_voro_par2d(nr_ndx, axis, gmx_box, vor_box, gridn);

    gmx_frame_area = vor_box[XX]*vor_box[YY]*vor_box[ZZ];
    gmx_tot_area += gmx_frame_area;
    container_poly vorcon(0.0,vor_box[XX],0.0,vor_box[YY],0.0,vor_box[ZZ],gridn[0],gridn[1],gridn[2],xper,yper,zper,8);

    //gmx_rmpbc(gpbc,natoms,gmx_box,x0);
    put_atoms_in_box(ePBC,gmx_box,natoms,x0);
    pid = 0;
    maxcells = 0;

    if (bMol){
      calc_mol_com(nr_grps, nr_mols, nr_atoms_mol, index, &top->atoms, x0, xmol);
      for (n = 0; n < nr_grps; n++) {
        for (i = 0; i < nr_mols[n]; i++) {   /* loop over all molecules in each group and add them to voronoi container*/
          copy_rvec2d(axis, gmx_box, xmol[pid], px);
          vorcon.put(vorpo,pid,px[XX],px[YY],px[ZZ],1.0);
          pid += 1;
          maxcells += 1;
        }
      }
    }
    else if (vor3d){
      for (n = 0; n < nr_grps; n++) {
        for (i = 0; i < grpsize[n]; i++) {   /* loop over all atoms in each group and add them to voronoi container*/
          copy_rvec3d(gmx_box, x0[index[n][i]], px);
          vorcon.put(vorpo,pid,px[XX],px[YY],px[ZZ],radii[index[n][i]]);
          pid += 1;
          maxcells += 1;
          AA_index[index[n][i]] = 1;
        }
      }
      for (i = 0; i < natoms; i++) {   // loop over the remaining atoms in the trajectory and add them to voronoi container although their volumes will not be computed
          if (AA_index[i] == 0) {
            copy_rvec3d(gmx_box, x0[i], px);
            vorcon.put(vorpo,pid,px[XX],px[YY],px[ZZ],radii[i]);
            pid += 1;
          }
      }
    }
    else {
      for (n = 0; n < nr_grps; n++) {
        for (i = 0; i < grpsize[n]; i++) {   /* loop over all atoms in each group and add them to voronoi container*/
          copy_rvec2d(axis, gmx_box, x0[index[n][i]], px);
          vorcon.put(vorpo,pid,px[XX],px[YY],px[ZZ],radii[index[n][i]]);
          pid += 1;
          maxcells += 1;
        }
      }
    }
    voronoicell c;
    c_loop_order vl(vorcon, vorpo);
    vor_frame_area = 0.0;
    i=0;
    if (vl.start()){   /* loop over all particles in voronoi container and get areas */
      do{
        if (vorcon.compute_cell(c,vl)){
          cell_area = c.volume();
          vl.pos(id,x,y,z,r);
          //printf("\ni= %d, vol=%6.12f, radius = %6.12f, xyz = %f %f %f",id, cell_area,r,x,y,z);
          frame_areas[id] = cell_area;
          areas[id] += cell_area;
          vor_frame_area += cell_area;
          compute_cells += 1;
        }
        i+=1;
      }while((vl.inc()) && (i < maxcells));
    }
    vor_tot_area += vor_frame_area;

    /* If computing a global tesselation over all atoms, sum the atoms corresponding to each molecule */
    if (vorAA){
      t = 0;
      k = 0;
      for (n = 0; n < nr_grps; n++){
        for (i = 0; i < nr_mols[n]; i++){
          sum = 0.0;
          for (j = 0; j < nr_atoms_mol[n]; j++){
            sum += frame_areas[t];
            t += 1;
          }
          frame_areas_full[k] = sum;
          areas_full[k] += sum;
          k += 1;
        }
      }
    }
    if (vorAA){
        avg_areas(nr_grps, fgrpsize, frame_areas_full, grp_areas, grp_areas_var);
    } else {
        avg_areas(nr_grps, grpsize, frame_areas, grp_areas, grp_areas_var);
    }

    for (n = 0; n < nr_grps; n++)
      avg_grp_areas_var[n] += grp_areas_var[n];

    if (fullout == xTs){
      fprintf(tsfp,"%6.3f\t\t",tt);
      for (n = 0; n < nr_grps; n++)
        if (vorAA)
          fprintf(tsfp,"%6.6f\t%6.6f\t%6.6f\t%6.6f\t", grp_areas[n], sqrt(grp_areas_var[n]), grp_areas[n]*fgrpsize[n], sqrt(grp_areas_var[n]*fgrpsize[n]));
        else
          fprintf(tsfp,"%6.6f\t%6.6f\t%6.6f\t%6.6f\t", grp_areas[n], sqrt(grp_areas_var[n]), grp_areas[n]*grpsize[n], sqrt(grp_areas_var[n]*grpsize[n]));
      fprintf(tsfp,"%6.6f\n", vor_frame_area);
    } else if (fullout == xAll){
      fprintf(allfp,"###############################\n");
      fprintf(allfp,"# Time frame = %6.4f ps\n",tt);
      fprintf(allfp,"# Box %s = %6.6f %s\n",label[s], gmx_frame_area, units[s]);
      fprintf(allfp,"# Total voronoi %s = %6.6f %s\n", label[s], vor_frame_area, units[s]);
      t = 0;
      for (n = 0; n < nr_grps; n++){
        fprintf(allfp,"# %s of each molecule in group %s:\n",label[s], grpname[n]);
        if (vorAA){
          for (i = 0; i < fgrpsize[n]; i++){
            fprintf(allfp,"%6.6f\t", frame_areas_full[t]);
            t += 1;
          }
        } else {
          for (i = 0; i < grpsize[n]; i++){
            fprintf(allfp,"%6.6f\t", frame_areas[t]);
            t += 1;
          }
        }
        fprintf(allfp,"\n");
      } 
    } else if (fullout == xAll_Ts){
      fprintf(allfp,"###############################\n");
      fprintf(allfp,"# Time frame = %6.4f ps\n",tt);
      fprintf(allfp,"# Box %s = %6.6f %s\n",label[s], gmx_frame_area, units[s]);
      fprintf(allfp,"# Total voronoi %s = %6.6f %s\n", label[s], vor_frame_area, units[s]);
      t = 0;
      for (n = 0; n < nr_grps; n++){
        fprintf(allfp,"# %s of each molecule in group %s:\n",label[s], grpname[n]);
        if (vorAA){
          for (i = 0; i < fgrpsize[n]; i++){
            fprintf(allfp,"%6.6f\t", frame_areas_full[t]);
            t += 1;
          }
        } else {
          for (i = 0; i < grpsize[n]; i++){
            fprintf(allfp,"%6.6f\t", frame_areas[t]);
            t += 1;
          }
        }
        fprintf(allfp,"\n");
      } 
      fprintf(tsfp,"%6.3f\t\t",tt);
      for (n = 0; n < nr_grps; n++)
        if (vorAA)
          fprintf(tsfp,"%6.6f\t%6.6f\t%6.6f\t%6.6f\t", grp_areas[n], sqrt(grp_areas_var[n]), grp_areas[n]*fgrpsize[n], sqrt(grp_areas_var[n]*fgrpsize[n]));
        else
          fprintf(tsfp,"%6.6f\t%6.6f\t%6.6f\t%6.6f\t", grp_areas[n], sqrt(grp_areas_var[n]), grp_areas[n]*grpsize[n], sqrt(grp_areas_var[n]*grpsize[n]));
      fprintf(tsfp,"%6.6f\n", vor_frame_area);
    }

    if (vert && nr_frames == 0){ /* output optional voronoi information */
      vorcon.draw_cells_gnuplot(vfile);
      vorcon.draw_particles(pfile);
    }
    nr_frames++;
  } while (read_next_x(oenv,status,&tt,x0,gmx_box));

  printf("Cells computed this frame: %i\n", compute_cells);
  fprintf(stderr,"\nRead %d frames from trajectory.\n", nr_frames);
  gmx_rmpbc_done(gpbc);

  /*********** done with status file **********/
  //close_trj(status);


  if (fullout == xTs){
    fclose(tsfp);
  } else if (fullout == xAll){
    fclose(allfp);
  } else if (fullout == xAll_Ts){
    fclose(allfp);
    fclose(tsfp);
  }

  if (vorAA){
    for (i = 0; i < nr_full; i++)
      areas_full[i] /= nr_frames;
    avg_areas(nr_grps, fgrpsize, areas_full, grp_areas, grp_areas_var);
  } else {
    for (n = 0; n < nr_ndx; n++)
      areas[n] /= nr_frames;
    avg_areas(nr_grps, grpsize, areas, grp_areas, grp_areas_var);
  }

  vor_tot_area /= nr_frames;
  gmx_tot_area /= nr_frames;
  for (n = 0; n < nr_grps; n++){
    avg_grp_areas_var[n] /= nr_frames;
  }
  /* write areas to file */

  FILE * avgfp;
  avgfp = fopen(avgfile,"w");
  fprintf(avgfp,"# Total number of frames analyzed = %i\n", nr_frames);
  fprintf(avgfp,"# Average box %s = %6.6f %s\n",label[s], gmx_tot_area, units[s]);
  fprintf(avgfp,"# Average total voronoi %s = %6.6f %s\n",label[s], vor_tot_area, units[s]);
  t = 0;
  for (n = 0; n < nr_grps; n++){
    fprintf(avgfp,"# Average %s per molecule in group %s = %6.6f +/- %6.6f %s\n", label[s], grpname[n], grp_areas[n], sqrt(avg_grp_areas_var[n]), units[s]);
    if (vorAA)
        fprintf(avgfp,"# Average %s for all molecules in group %s = %6.6f +/- %6.6f %s\n", label[s], grpname[n], grp_areas[n]*fgrpsize[n], sqrt(avg_grp_areas_var[n]*fgrpsize[n]), units[s]);
    else
        fprintf(avgfp,"# Average %s for all molecules in group %s = %6.6f +/- %6.6f %s\n", label[s], grpname[n], grp_areas[n]*grpsize[n], sqrt(avg_grp_areas_var[n]*grpsize[n]), units[s]);
    fprintf(avgfp,"# Average %s of each molecule in group %s:\n", label[s], grpname[n]);
    if (vorAA){
      for (i = 0; i < fgrpsize[n]; i++){
        fprintf(avgfp,"%6.6f\t", areas_full[t]);
        t += 1;
      }
    } else {
      for (i = 0; i < grpsize[n]; i++){
        fprintf(avgfp,"%6.6f\t", areas[t]);
        t += 1;
      }
    }
    fprintf(avgfp,"\n");
  }
  fclose(avgfp);
  //if (bMol)
    //sfree(xmol);
  //if (vorAA){
    //sfree(areas_full);
    //sfree(frame_areas_full);
    //sfree(fgrpsize);
  //}
  //sfree(areas);
  //sfree(frame_areas);
  //sfree(x0);  /* free memory used by coordinate array */
}

int gmx_voronoi(int argc,char *argv[])
{
  const char *desc[] = {
    "[THISMODULE] is a tool to compute voronoi tesselations and calculate volumes/areas from a structure (pdb, gro, etc.) ",
    "or a simulation trajectory. In addition to the positions of the particles, you will need a tpr and an index file.[PAR]",
    "The default behavior of [THISMODULE] is to compute a radical voronoi tesselation (option [TT]-radtess[tt]) in 3D (option [TT]-3d[tt]) using ",
    "all atoms (option [TT]-aa[tt]) in the system. In a radical tesselation, the particle radii (obtained from the VdW parameters) ",
    "are used to weight the location of the boundaries between them. Tesselations can be computed assumming periodic boundary conditions ",
    "(option [TT]-pbc[tt]) or hard walls ([TT]-nopbc[tt]). Conventional voronoi tesselation (non-radical) can be performed with the option ",
    "[TT]-noradtess[tt], and the default particle radii can be set with option [TT]-radius[tt]. Please note that if the default particle ",
    "radius is too large, the tesselation may not provide accurate results.[PAR]",
    "[THISMODULE] was created with molecular quantities in mind, therefore it assumes that index groups correspond to groups of whole molecules",
    "In addition to 3D volume calculations, [THISMODULE] can compute a 2D tesselation (option [TT]-no3d[tt]) of the particle's projection onto a plane defined by a ",
    "normal axis (option [TT]-normal[tt]). In this case, the output will be the areas "
    
    ""
  };
  t_tpxheader header;
  t_inputrec  *ir;
  t_state     *state;
  gmx_mtop_t  *mtop;
  t_forcerec  *fr;
  t_commrec   *cr;
  gmx_hw_info_t *hwinfo;
  matrix      box;
  gmx_output_env_t * oenv;
  static int fullout = xNone;
  static int axis = 2;          /* normal to memb. default z  */
  static const char *axtitle="Z";
  static const char *fullsel="ts";
  static int  ngrps   = 1;       /* nr. of groups              */
  gmx_bool vorPBC=TRUE;
  gmx_bool bMol=FALSE;
  gmx_bool vert=FALSE;
  gmx_bool vorAA=TRUE;
  gmx_bool vor3d=TRUE;
  gmx_bool bRad=TRUE;
  real tessRadius = 1.0;
  t_pargs pa[] = {
    { "-normal",    FALSE, etSTR, {&axtitle},
      "Take the normal on the membrane in direction X, Y or Z." },
    { "-ng",   FALSE, etINT, {&ngrps},
      "Number of groups to compute voronoi cells on" },
    { "-pbc",  FALSE, etBOOL, {&vorPBC},
      "Compute Voronoi cells using periodic boundary conditions" },
    { "-mol",  FALSE, etBOOL, {&bMol},
      "Assume that groups from index file contain whole molecules instead of individual atoms, and compute Voronoi cells based on the center of mass of each molecule" },
    { "-vert",  FALSE, etBOOL, {&vert},
      "Output the information needed to plot the voronoi tesselation. Use -ov to specify the file containing the vertices, and -op to specify the file containing the tesselation points" },
    { "-aa", FALSE, etBOOL, {&vorAA},
      "Compute a complete tesselation over all atoms and sum for each molecule" },
    { "-radtess", FALSE, etBOOL, {&bRad},
      "Compute Radical tesselation where the vdw radii are used to weight the voronoi cells" },
    { "-radius", FALSE, etREAL, {&tessRadius},
      "Default particle radius for non-radical tesselation"},
    { "-3d", FALSE, etBOOL, {&vor3d},
      "Compute tesselation in 3D, turns on -aa"},
    { "-full", FALSE, etSTR, {&fullsel},
      "Choose to output additional volume information: none, all (all volumes at each time frame), ts (avg volume at each time frame), all_ts (both all and ts)"}
  };

  const char *bugs[] = {    " ",  };
  char **grpname;        /* groupnames                 */
  int  *grpsize;         /* sizes of groups            */
  t_topology *top;       /* topology                   */
  int  ePBC;
  atom_id   **index;     /* indices for all groups     */
  int  i,ntopatoms;

  t_filenm  fnm[] = {    /* files for g_voronoi        */
    { efTRX, "-f", NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffREAD },
    { efTPR, NULL, NULL,  ffREAD },
    { efDAT,"-avg","avgvolumes", ffWRITE },
    { efDAT,"-ts","timeseries", ffOPTWR },
    { efDAT,"-all","allvolumes", ffOPTWR },
    { efDAT,"-ov","vertices", ffOPTWR },
    { efDAT,"-op","points", ffOPTWR }
  };

#define NFILE asize(fnm)
  cr = init_commrec();
  if (!parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs,
                    &oenv))
      return 0;

  /* Calculate axis */
  axis = toupper(axtitle[0]) - 'X';
  top = read_top(ftp2fn(efTPR,NFILE,fnm),&ePBC);
  snew(ir, 1);
  snew(state, 1);
  snew(mtop, 1);
  read_tpx_state(ftp2fn(efTPR,NFILE,fnm),ir,state,mtop);
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(grpsize,ngrps);
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,grpsize,index,grpname);

  if (vor3d)
    vorAA = TRUE;

  if (bRad && bMol)
    gmx_fatal(FARGS,"-rad is incompatible with -mol\n");

  if (vor3d && bMol)
    gmx_fatal(FARGS,"Cannot use -3d and -mol together.\n");

  // only need the C6 and C12 parameters, no need for a full forcerec
  // initialization
  fr = mk_forcerec();
  fr->bBHAM = (mtop->ffparams.functype[0] == F_BHAM);
  fr->ntype = mtop->ffparams.atnr;
  fr->nbfp  = mk_nbfp(&mtop->ffparams, fr->bBHAM);

  if (strncmp(fullsel,"none",6) == 0) {
    fullout = xNone;
  } else if (strncmp(fullsel,"all",6) == 0) {
    fullout = xAll;
  } else if(strncmp(fullsel,"ts",6) == 0){
    fullout = xTs;
  } else if(strncmp(fullsel,"all_ts",6) == 0){
    fullout = xAll_Ts;
  } else {
    gmx_fatal(FARGS,"Wrong setting of -full flag.\n");
  }
  
  compute_voronoi(ftp2fn(efTRX,NFILE,fnm), index, grpsize, top, ePBC, axis, ngrps, grpname, bMol, vorPBC, vert, vorAA, vor3d,
                  bRad, tessRadius, fullout, fr, opt2fn("-avg",NFILE,fnm), opt2fn("-ts",NFILE,fnm), opt2fn("-all",NFILE,fnm), opt2fn("-ov",NFILE,fnm), opt2fn("-op",NFILE,fnm), oenv);
  return 0;
}

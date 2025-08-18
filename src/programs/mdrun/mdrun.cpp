/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016, by the GROMACS development team, led by
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
/*! \defgroup module_mdrun Implementation of mdrun
 * \ingroup group_mdrun
 *
 * \brief This module contains code that implements mdrun.
 */
/*! \internal \file
 *
 * \brief This file implements mdrun
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "config.h"

#include <stdio.h>
#include <string.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"

#include "mdrun_main.h"
#include "runner.h"

extern mds::StressGrid locals_grid;

/*! \brief Return whether either of the command-line parameters that
 *  will trigger a multi-simulation is set */
static bool is_multisim_option_set(int argc, const char *const argv[])
{
    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], "-multi") == 0 || strcmp(argv[i], "-multidir") == 0)
        {
            return true;
        }
    }
    return false;
}

//! Implements C-style main function for mdrun
int gmx_mdrun(int argc, char *argv[])
{
    const char   *desc[] = {
        "[THISMODULE] is the main computational chemistry engine",
        "within GROMACS. Obviously, it performs Molecular Dynamics simulations,",
        "but it can also perform Stochastic Dynamics, Energy Minimization,",
        "test particle insertion or (re)calculation of energies.",
        "Normal mode analysis is another option. In this case [TT]mdrun[tt]",
        "builds a Hessian matrix from single conformation.",
        "For usual Normal Modes-like calculations, make sure that",
        "the structure provided is properly energy-minimized.",
        "The generated matrix can be diagonalized by [gmx-nmeig].[PAR]",
        "The [TT]mdrun[tt] program reads the run input file ([TT]-s[tt])",
        "and distributes the topology over ranks if needed.",
        "[TT]mdrun[tt] produces at least four output files.",
        "A single log file ([TT]-g[tt]) is written.",
        "The trajectory file ([TT]-o[tt]), contains coordinates, velocities and",
        "optionally forces.",
        "The structure file ([TT]-c[tt]) contains the coordinates and",
        "velocities of the last step.",
        "The energy file ([TT]-e[tt]) contains energies, the temperature,",
        "pressure, etc, a lot of these things are also printed in the log file.",
        "Optionally coordinates can be written to a compressed trajectory file",
        "([TT]-x[tt]).[PAR]",
        "The option [TT]-dhdl[tt] is only used when free energy calculation is",
        "turned on.[PAR]",
        "Running mdrun efficiently in parallel is a complex topic topic,",
        "many aspects of which are covered in the online User Guide. You",
        "should look there for practical advice on using many of the options",
        "available in mdrun.[PAR]",
        "ED (essential dynamics) sampling and/or additional flooding potentials",
        "are switched on by using the [TT]-ei[tt] flag followed by an [REF].edi[ref]",
        "file. The [REF].edi[ref] file can be produced with the [TT]make_edi[tt] tool",
        "or by using options in the essdyn menu of the WHAT IF program.",
        "[TT]mdrun[tt] produces a [REF].xvg[ref] output file that",
        "contains projections of positions, velocities and forces onto selected",
        "eigenvectors.[PAR]",
        "When user-defined potential functions have been selected in the",
        "[REF].mdp[ref] file the [TT]-table[tt] option is used to pass [TT]mdrun[tt]",
        "a formatted table with potential functions. The file is read from",
        "either the current directory or from the [TT]GMXLIB[tt] directory.",
        "A number of pre-formatted tables are presented in the [TT]GMXLIB[tt] dir,",
        "for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard-Jones potentials with",
        "normal Coulomb.",
        "When pair interactions are present, a separate table for pair interaction",
        "functions is read using the [TT]-tablep[tt] option.[PAR]",
        "When tabulated bonded functions are present in the topology,",
        "interaction functions are read using the [TT]-tableb[tt] option.",
        "For each different tabulated interaction type used, a table file name must",
        "be given. For the topology to work, a file name given here must match a",
        "character sequence before the file extension. That sequence is: an underscore,",
        "then a 'b' for bonds, an 'a' for angles or a 'd' for dihedrals,",
        "and finally the matching table number index used in the topology.[PAR]",
        "The options [TT]-px[tt] and [TT]-pf[tt] are used for writing pull COM",
        "coordinates and forces when pulling is selected",
        "in the [REF].mdp[ref] file.[PAR]",
        "Finally some experimental algorithms can be tested when the",
        "appropriate options have been given. Currently under",
        "investigation are: polarizability.",
        "[PAR]",
        "The option [TT]-membed[tt] does what used to be g_membed, i.e. embed",
        "a protein into a membrane. This module requires a number of settings",
        "that are provided in a data file that is the argument of this option.",
        "For more details in membrane embedding, see the documentation in the",
        "user guide. The options [TT]-mn[tt] and [TT]-mp[tt] are used to provide",
        "the index and topology files used for the embedding.",
        "[PAR]",
        "The option [TT]-pforce[tt] is useful when you suspect a simulation",
        "crashes due to too large forces. With this option coordinates and",
        "forces of atoms with a force larger than a certain value will",
        "be printed to stderr.",
        "[PAR]",
        "Checkpoints containing the complete state of the system are written",
        "at regular intervals (option [TT]-cpt[tt]) to the file [TT]-cpo[tt],",
        "unless option [TT]-cpt[tt] is set to -1.",
        "The previous checkpoint is backed up to [TT]state_prev.cpt[tt] to",
        "make sure that a recent state of the system is always available,",
        "even when the simulation is terminated while writing a checkpoint.",
        "With [TT]-cpnum[tt] all checkpoint files are kept and appended",
        "with the step number.",
        "A simulation can be continued by reading the full state from file",
        "with option [TT]-cpi[tt]. This option is intelligent in the way that",
        "if no checkpoint file is found, GROMACS just assumes a normal run and",
        "starts from the first step of the [REF].tpr[ref] file. By default the output",
        "will be appending to the existing output files. The checkpoint file",
        "contains checksums of all output files, such that you will never",
        "loose data when some output files are modified, corrupt or removed.",
        "There are three scenarios with [TT]-cpi[tt]:[PAR]",
        "[TT]*[tt] no files with matching names are present: new output files are written[PAR]",
        "[TT]*[tt] all files are present with names and checksums matching those stored",
        "in the checkpoint file: files are appended[PAR]",
        "[TT]*[tt] otherwise no files are modified and a fatal error is generated[PAR]",
        "With [TT]-noappend[tt] new output files are opened and the simulation",
        "part number is added to all output file names.",
        "Note that in all cases the checkpoint file itself is not renamed",
        "and will be overwritten, unless its name does not match",
        "the [TT]-cpo[tt] option.",
        "[PAR]",
        "With checkpointing the output is appended to previously written",
        "output files, unless [TT]-noappend[tt] is used or none of the previous",
        "output files are present (except for the checkpoint file).",
        "The integrity of the files to be appended is verified using checksums",
        "which are stored in the checkpoint file. This ensures that output can",
        "not be mixed up or corrupted due to file appending. When only some",
        "of the previous output files are present, a fatal error is generated",
        "and no old output files are modified and no new output files are opened.",
        "The result with appending will be the same as from a single run.",
        "The contents will be binary identical, unless you use a different number",
        "of ranks or dynamic load balancing or the FFT library uses optimizations",
        "through timing.",
        "[PAR]",
        "With option [TT]-maxh[tt] a simulation is terminated and a checkpoint",
        "file is written at the first neighbor search step where the run time",
        "exceeds [TT]-maxh[tt]\\*0.99 hours. This option is particularly useful in",
        "combination with setting [TT]nsteps[tt] to -1 either in the mdp or using the",
        "similarly named command line option. This results in an infinite run,",
        "terminated only when the time limit set by [TT]-maxh[tt] is reached (if any)"
        "or upon receiving a signal."
        "[PAR]",
        "When [TT]mdrun[tt] receives a TERM signal, it will stop as soon as",
        "checkpoint file can be written, i.e. after the next global communication step.",
        "When [TT]mdrun[tt] receives an INT signal (e.g. when ctrl+C is",
        "pressed), it will stop at the next neighbor search step or at the",
        "second global communication step, whichever happens later.",
        "In both cases all the usual output will be written to file.",
        "When running with MPI, a signal to one of the [TT]mdrun[tt] ranks",
        "is sufficient, this signal should not be sent to mpirun or",
        "the [TT]mdrun[tt] process that is the parent of the others.",
        "[PAR]",
        "Interactive molecular dynamics (IMD) can be activated by using at least one",
        "of the three IMD switches: The [TT]-imdterm[tt] switch allows one to terminate",
        "the simulation from the molecular viewer (e.g. VMD). With [TT]-imdwait[tt],",
        "[TT]mdrun[tt] pauses whenever no IMD client is connected. Pulling from the",
        "IMD remote can be turned on by [TT]-imdpull[tt].",
        "The port [TT]mdrun[tt] listens to can be altered by [TT]-imdport[tt].The",
        "file pointed to by [TT]-if[tt] contains atom indices and forces if IMD",
        "pulling is used."
        "[PAR]",
        "When [TT]mdrun[tt] is started with MPI, it does not run niced by default."
    };
    t_commrec    *cr;
    t_filenm      fnm[] = {
        { efTPR, NULL,      NULL,       ffREAD },
        { efTRN, "-o",      NULL,       ffWRITE },
        { efCOMPRESSED, "-x", NULL,     ffOPTWR },
        { efCPT, "-cpi",    NULL,       ffOPTRD | ffALLOW_MISSING },
        { efCPT, "-cpo",    NULL,       ffOPTWR },
        { efSTO, "-c",      "confout",  ffWRITE },
        { efEDR, "-e",      "ener",     ffWRITE },
        { efLOG, "-g",      "md",       ffWRITE },
        { efXVG, "-dhdl",   "dhdl",     ffOPTWR },
        { efXVG, "-field",  "field",    ffOPTWR },
        { efXVG, "-table",  "table",    ffOPTRD },
        { efXVG, "-tablep", "tablep",   ffOPTRD },
        { efXVG, "-tableb", "table",    ffOPTRDMULT },
        { efTRX, "-rerun",  "rerun",    ffOPTRD },
        { efXVG, "-tpi",    "tpi",      ffOPTWR },
        { efXVG, "-tpid",   "tpidist",  ffOPTWR },
        { efEDI, "-ei",     "sam",      ffOPTRD },
        { efXVG, "-eo",     "edsam",    ffOPTWR },
        { efXVG, "-devout", "deviatie", ffOPTWR },
        { efXVG, "-runav",  "runaver",  ffOPTWR },
        { efXVG, "-px",     "pullx",    ffOPTWR },
        { efXVG, "-pf",     "pullf",    ffOPTWR },
        { efXVG, "-ro",     "rotation", ffOPTWR },
        { efLOG, "-ra",     "rotangles", ffOPTWR },
        { efLOG, "-rs",     "rotslabs", ffOPTWR },
        { efLOG, "-rt",     "rottorque", ffOPTWR },
        { efMTX, "-mtx",    "nm",       ffOPTWR },
        { efRND, "-multidir", NULL,      ffOPTRDMULT},
        { efDAT, "-membed", "membed",   ffOPTRD },
        { efTOP, "-mp",     "membed",   ffOPTRD },
        { efNDX, "-mn",     "membed",   ffOPTRD },
        { efDAT, "-ols",    "localstress", ffWRITE },
        { efXVG, "-if",     "imdforces", ffOPTWR },
        { efXVG, "-swap",   "swapions", ffOPTWR }
    };
    const int     NFILE = asize(fnm);

    /* Command line options ! */
    gmx_bool          bDDBondCheck  = TRUE;
    gmx_bool          bDDBondComm   = TRUE;
    gmx_bool          bTunePME      = TRUE;
    gmx_bool          bVerbose      = FALSE;
    gmx_bool          bRerunVSite   = FALSE;
    gmx_bool          bConfout      = TRUE;
    gmx_bool          bReproducible = FALSE;
    gmx_bool          bIMDwait      = FALSE;
    gmx_bool          bIMDterm      = FALSE;
    gmx_bool          bIMDpull      = FALSE;

    int               npme          = 0;
    int               nstlist       = 0;
    int               nmultisim     = 0;
    int               nstglobalcomm = -1;
    int               repl_ex_nst   = 0;
    int               repl_ex_seed  = -1;
    int               repl_ex_nex   = 0;
    int               nstepout      = 100;
    int               resetstep     = -1;
    gmx_int64_t       nsteps        = -2;   /* the value -2 means that the mdp option will be used */
    int               imdport       = 8888; /* can be almost anything, 8888 is easy to remember */

    rvec              realddxyz                   = {0, 0, 0};
    const char       *ddrank_opt[ddrankorderNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };
    const char       *dddlb_opt[] =
    { NULL, "auto", "no", "yes", NULL };
    const char       *thread_aff_opt[threadaffNR+1] =
    { NULL, "auto", "on", "off", NULL };
    const char       *nbpu_opt[] =
    { NULL, "auto", "cpu", "gpu", "gpu_cpu", NULL };
    real              rdd                   = 0.0, rconstr = 0.0, dlb_scale = 0.8, pforce = -1;
    char             *ddcsx                 = NULL, *ddcsy = NULL, *ddcsz = NULL;
    real              cpt_period            = 15.0, max_hours = -1;
    gmx_bool          bTryToAppendFiles     = TRUE;
    gmx_bool          bKeepAndNumCPT        = FALSE;
    gmx_bool          bResetCountersHalfWay = FALSE;
    gmx_output_env_t *oenv                  = NULL;
    
    /* Local Stress parameters
     */
    real localsmindihangle=0.0;
    int nstlocals=0;
    real localsgridspacing=0.1;
    real localsimpulsewidth=0.0001;
    int localsgridx=0;
    int localsgridy=0;
    int localsgridz=0;
    int localsskip=1;
    int localsdebugprint=0;
    const char * localsenum = "all";
    const char * localsfdenum = "ccfd";
    gmx_bool localsdispcor = TRUE;
    gmx_bool localscuda = FALSE;
    gmx_bool localspbc = FALSE;

    /* Non transparent initialization of a complex gmx_hw_opt_t struct.
     * But unfortunately we are not allowed to call a function here,
     * since declarations follow below.
     */
    gmx_hw_opt_t    hw_opt = {
        0, 0, 0, 0, threadaffSEL, 0, 0,
        { NULL, FALSE, 0, NULL }
    };

    t_pargs         pa[] = {

        { "-dd",      FALSE, etRVEC, {&realddxyz},
          "Domain decomposition grid, 0 is optimize" },
        { "-ddorder", FALSE, etENUM, {ddrank_opt},
          "DD rank order" },
        { "-npme",    FALSE, etINT, {&npme},
          "Number of separate ranks to be used for PME, -1 is guess" },
        { "-nt",      FALSE, etINT, {&hw_opt.nthreads_tot},
          "Total number of threads to start (0 is guess)" },
        { "-ntmpi",   FALSE, etINT, {&hw_opt.nthreads_tmpi},
          "Number of thread-MPI threads to start (0 is guess)" },
        { "-ntomp",   FALSE, etINT, {&hw_opt.nthreads_omp},
          "Number of OpenMP threads per MPI rank to start (0 is guess)" },
        { "-ntomp_pme", FALSE, etINT, {&hw_opt.nthreads_omp_pme},
          "Number of OpenMP threads per MPI rank to start (0 is -ntomp)" },
        { "-pin",     FALSE, etENUM, {thread_aff_opt},
          "Whether mdrun should try to set thread affinities" },
        { "-pinoffset", FALSE, etINT, {&hw_opt.core_pinning_offset},
          "The lowest logical core number to which mdrun should pin the first thread" },
        { "-pinstride", FALSE, etINT, {&hw_opt.core_pinning_stride},
          "Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core" },
        { "-gpu_id",  FALSE, etSTR, {&hw_opt.gpu_opt.gpu_id},
          "List of GPU device id-s to use, specifies the per-node PP rank to GPU mapping" },
        { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
          "Check for all bonded interactions with DD" },
        { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
          "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
        { "-rdd",     FALSE, etREAL, {&rdd},
          "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
        { "-rcon",    FALSE, etREAL, {&rconstr},
          "Maximum distance for P-LINCS (nm), 0 is estimate" },
        { "-dlb",     FALSE, etENUM, {dddlb_opt},
          "Dynamic load balancing (with DD)" },
        { "-dds",     FALSE, etREAL, {&dlb_scale},
          "Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in order to "
          "provide a margin in which dynamic load balancing can act while preserving the minimum cell size." },
        { "-ddcsx",   FALSE, etSTR, {&ddcsx},
          "HIDDENA string containing a vector of the relative sizes in the x "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsy",   FALSE, etSTR, {&ddcsy},
          "HIDDENA string containing a vector of the relative sizes in the y "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-ddcsz",   FALSE, etSTR, {&ddcsz},
          "HIDDENA string containing a vector of the relative sizes in the z "
          "direction of the corresponding DD cells. Only effective with static "
          "load balancing." },
        { "-gcom",    FALSE, etINT, {&nstglobalcomm},
          "Global communication frequency" },
        { "-nb",      FALSE, etENUM, {&nbpu_opt},
          "Calculate non-bonded interactions on" },
        { "-nstlist", FALSE, etINT, {&nstlist},
          "Set nstlist when using a Verlet buffer tolerance (0 is guess)" },
        { "-tunepme", FALSE, etBOOL, {&bTunePME},
          "Optimize PME load between PP/PME ranks or GPU/CPU" },
        { "-v",       FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy" },
        { "-pforce",  FALSE, etREAL, {&pforce},
          "Print all forces larger than this (kJ/mol nm)" },
        { "-reprod",  FALSE, etBOOL, {&bReproducible},
          "Try to avoid optimizations that affect binary reproducibility" },
        { "-cpt",     FALSE, etREAL, {&cpt_period},
          "Checkpoint interval (minutes)" },
        { "-cpnum",   FALSE, etBOOL, {&bKeepAndNumCPT},
          "Keep and number checkpoint files" },
        { "-append",  FALSE, etBOOL, {&bTryToAppendFiles},
          "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names" },
        { "-nsteps",  FALSE, etINT64, {&nsteps},
          "Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)" },
        { "-maxh",   FALSE, etREAL, {&max_hours},
          "Terminate after 0.99 times this time (hours)" },
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do multiple simulations in parallel" },
        { "-replex",  FALSE, etINT, {&repl_ex_nst},
          "Attempt replica exchange periodically with this period (steps)" },
        { "-nex",  FALSE, etINT, {&repl_ex_nex},
          "Number of random exchanges to carry out each exchange interval (N^3 is one suggestion).  -nex zero or not specified gives neighbor replica exchange." },
        { "-reseed",  FALSE, etINT, {&repl_ex_seed},
          "Seed for replica exchange, -1 is generate a seed" },
        { "-localsgrid",  FALSE, etREAL, {&localsgridspacing},
          "Spacing for local stress grid (default = 0.1 nm)" },
        { "-localsiw",  FALSE, etREAL, {&localsimpulsewidth},
          "Apply an impulsive correction for plain cutoff potentials (VDW or elec) using a delta function with a finite width (default = 0.0001 nm)" },
        { "-nstlp",  FALSE, etINT, {&nstlocals},
          "HIDDENFrequency of writing local stress grid to file (default = 0)" },
        { "-lsgridx", FALSE, etINT, {&localsgridx},
          "Set the local stress grid size in the x direction (default use box[XX][XX]/localsgrid)"},
        { "-lsgridy", FALSE, etINT, {&localsgridy},
          "Set the local stress grid size in the y direction (default use box[YY][YY]/localsgrid)"},
        { "-lsgridz", FALSE, etINT, {&localsgridz},
          "Set the local stress grid size in the z direction (default use box[ZZ][ZZ]/localsgrid)"},
        { "-lscont", FALSE, etSTR, {&localsenum},
          "Select which contribution to write to output (default = all): all, vdw, coul, angles, bonds, dihp, dihi, dihrb, lincs, settle, shake, cmap, vel, none"},
        { "-lsfd", FALSE, etSTR, {&localsfdenum},
          "Select the type of force decomposition to be used: ccfd (covariant central force decomposition, default), ncfd (non-covariant central force decomposition), or gld (Goetz-Lipowsky decomposition)"},
        { "-lsdispcor",  FALSE, etBOOL, {&localsdispcor},
          "Include contribution from dispersion correction." },
        { "-lspbc",  FALSE, etBOOL, {&localspbc},
          "Correct periodic boundary conditions in mdstress library." },
        { "-lsmindihang",  FALSE, etREAL, {&localsmindihangle},
          "Don't include dihedral local stress contributions if the sin(|phi|) is less than this factor. Use this flag if there is a dihedral potential (e.g. CHARMM36 lipid FF) that has been parametrized with a min/max that is not 0/Pi and the stress profiles show large noise that does not converge with additional frames. A -lsmindihang value of 0.0005 is typically sufficient to fix this problem." },
        { "-lsskip",  FALSE, etINT, {&localsskip},
          "Only compute the local stress every nth frame" },
        { "-lscuda",  FALSE, etBOOL, {&localscuda},
          "HIDDENEnable CUDA operations when calculating local stress contributions" },
        { "-lsdebugprint",  FALSE, etINT, {&localsdebugprint},
          "HIDDENPrint various elements of the elasticity tensor to stdout at a given frame interval for debugging purposes" },
        { "-imdport",    FALSE, etINT, {&imdport},
          "HIDDENIMD listening port" },
        { "-imdwait",  FALSE, etBOOL, {&bIMDwait},
          "HIDDENPause the simulation while no IMD client is connected" },
        { "-imdterm",  FALSE, etBOOL, {&bIMDterm},
          "HIDDENAllow termination of the simulation from IMD client" },
        { "-imdpull",  FALSE, etBOOL, {&bIMDpull},
          "HIDDENAllow pulling in the simulation from IMD client" },
        { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
          "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
        { "-confout", FALSE, etBOOL, {&bConfout},
          "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
        { "-stepout", FALSE, etINT, {&nstepout},
          "HIDDENFrequency of writing the remaining wall clock time for the run" },
        { "-resetstep", FALSE, etINT, {&resetstep},
          "HIDDENReset cycle counters after these many time steps" },
        { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
          "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" }
    };
    unsigned long   Flags;
    ivec            ddxyz;
    int             dd_rank_order;
    int             localscontrib;
    int             localscontribc;
    int             localsfdecomp;
    gmx_bool        bDoAppendFiles, bStartFromCpt;
    FILE           *fplog;
    int             rc;
    char          **multidir = NULL;

    cr = init_commrec();

    unsigned long PCA_Flags = PCA_CAN_SET_DEFFNM;
    // With -multi or -multidir, the file names are going to get processed
    // further (or the working directory changed), so we can't check for their
    // existence during parsing.  It isn't useful to do any completion based on
    // file system contents, either.
    if (is_multisim_option_set(argc, argv))
    {
        PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
    }

    /* Comment this in to do fexist calls only on master
     * works not with rerun or tables at the moment
     * also comment out the version of init_forcerec in md.c
     * with NULL instead of opt2fn
     */
    /*
       if (!MASTER(cr))
       {
       PCA_Flags |= PCA_NOT_READ_NODE;
       }
     */

    if (!parse_common_args(&argc, argv, PCA_Flags, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    dd_rank_order = nenum(ddrank_opt);
    hw_opt.thread_affinity = nenum(thread_aff_opt);

    if (strcmp(localsfdenum,"ccfd") == 0) {
      localsfdecomp = mds_ccfd;
      printf("\nSelected force decomposition: %s\n", localsfdenum);
    }else if (strcmp(localsfdenum,"ncfd") == 0) {
      localsfdecomp = mds_ncfd;
      printf("\nSelected force decomposition: %s\n", localsfdenum);
    }else{
      printf("\nOption not recognized, will use covariant central force decomposition\n");
      localsfdecomp = mds_ccfd;
    }

    printf("\nSelected contribution: %s\n",localsenum);
    if (strcmp(localsenum,"all") == 0) {
      printf("\nWill write all contributions to the local stress\n");
      localscontrib = mds_all;
    }else if(strcmp(localsenum,"vdw") == 0){
      printf("\nWill only write vdw contributions to the local stress\n");
      localscontrib = mds_vdw;
    }else if(strcmp(localsenum,"coul") == 0){
      printf("\nWill only write coulomb contributions to the local stress\n");
      localscontrib = mds_cou;
    }else if(strcmp(localsenum,"angles") == 0){
      printf("\nWill only write angle contributions to the local stress\n");
      localscontrib = mds_ang;
    }else if(strcmp(localsenum,"bonds") == 0){
      printf("\nWill only write bonding contributions to the local stress\n");
      localscontrib = mds_bnd;
    }else if(strcmp(localsenum,"dihp") == 0){
      printf("\nWill only write proper dihedral contributions to the local stress\n");
      localscontrib = mds_dip;
    }else if(strcmp(localsenum,"dihi") == 0){
      printf("\nWill only write inproper dihedral contributions to the local stress\n");
      localscontrib = mds_dii;
    }else if(strcmp(localsenum,"diho") == 0){
      printf("\nWill only write other dihedral contributions to the local stress\n");
      localscontrib = mds_dio;
    }else if(strcmp(localsenum,"dihrb") == 0){
      printf("\nWill only write RB dihedral contributions to the local stress\n");
      localscontrib = mds_drb;
    }else if(strcmp(localsenum,"lincs") == 0){
      printf("\nWill only write LINCS constraints contributions to the local stress\n");
      localscontrib = mds_lin;
    }else if(strcmp(localsenum,"settle") == 0){
      printf("\nWill only write SETTLE water constraints contributions to the local stress\n");
      localscontrib = mds_set;
    }else if(strcmp(localsenum,"shake") == 0){
      printf("\nWill only write SHAKE constraints contributions to the local stress\n");
      localscontrib = mds_sha;
    }else if(strcmp(localsenum,"vel") == 0){
      printf("\nWill only write velocity contributions to the local stress\n");
      localscontrib = mds_kin;
    }else if(strcmp(localsenum,"cmap") == 0){
      printf("\nWill only write CMAP contributions to the local stress\n");
      localscontrib = mds_cmp;
    }else if(strcmp(localsenum,"none") == 0){
      printf("\nWill not write any contributions to the local stress\n");
      localscontrib = mds_none;
    }else{
      printf("\nOption not recognized, will write all contributions to the local stress\n");
      localscontrib = mds_all;
    }

    /* always set the maximum number of threads */
    locals_grid.SetMaxThreads(hw_opt.nthreads_tmpi);
    
    /* initialize what we can of locals_grid */
    if (false == locals_grid.settings.initialized) {
        if (localsdispcor == FALSE)
            locals_grid.DisableDispersionCorrection();
        if (localscuda == TRUE)
            locals_grid.EnableCuda();
        locals_grid.SetContribType(localscontrib);
        locals_grid.SetForceDecomposition(localsfdecomp);
        locals_grid.SetMinDihAngle(localsmindihangle);

        if(localsgridspacing<=0)
        {
            gmx_fatal(FARGS,"Cannot do local stress with spacing (-localsgrid) <= 0.0\n");
        }

        locals_grid.SetGridSpacing(localsgridspacing);
        locals_grid.SetImpulseWidth(localsimpulsewidth);
        locals_grid.SetGridNCells(
                localsgridx,
                localsgridy,
                localsgridz);
        locals_grid.SetDebugPrint(localsdebugprint);
    }

    /* now check the -multi and -multidir option */
    if (opt2bSet("-multidir", NFILE, fnm))
    {
        if (nmultisim > 0)
        {
            gmx_fatal(FARGS, "mdrun -multi and -multidir options are mutually exclusive.");
        }
        nmultisim = opt2fns(&multidir, "-multidir", NFILE, fnm);
    }

    if (repl_ex_nst != 0 && nmultisim < 2)
    {
        gmx_fatal(FARGS, "Need at least two replicas for replica exchange (option -multi)");
    }

    if (repl_ex_nex < 0)
    {
        gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
    }

    if (nmultisim >= 1)
    {
#if !GMX_THREAD_MPI
        init_multisystem(cr, nmultisim, multidir, NFILE, fnm);
#else
        gmx_fatal(FARGS, "mdrun -multi or -multidir are not supported with the thread-MPI library. "
                  "Please compile GROMACS with a proper external MPI library.");
#endif
    }

    if (!opt2bSet("-cpi", NFILE, fnm))
    {
        // If we are not starting from a checkpoint we never allow files to be appended
        // to, since that has caused a ton of strange behaviour and bugs in the past.
        if (opt2parg_bSet("-append", asize(pa), pa))
        {
            // If the user explicitly used the -append option, explain that it is not possible.
            gmx_fatal(FARGS, "GROMACS can only append to files when restarting from a checkpoint.");
        }
        else
        {
            // If the user did not say anything explicit, just disable appending.
            bTryToAppendFiles = FALSE;
        }
    }

    handleRestart(cr, bTryToAppendFiles, NFILE, fnm, &bDoAppendFiles, &bStartFromCpt);

    Flags = opt2bSet("-rerun", NFILE, fnm) ? MD_RERUN : 0;
    Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
    Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
    Flags = Flags | (bTunePME      ? MD_TUNEPME      : 0);
    Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
    Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
    Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
    Flags = Flags | (bDoAppendFiles  ? MD_APPENDFILES  : 0);
    Flags = Flags | (opt2parg_bSet("-append", asize(pa), pa) ? MD_APPENDFILESSET : 0);
    Flags = Flags | (bKeepAndNumCPT ? MD_KEEPANDNUMCPT : 0);
    Flags = Flags | (bStartFromCpt ? MD_STARTFROMCPT : 0);
    Flags = Flags | (bResetCountersHalfWay ? MD_RESETCOUNTERSHALFWAY : 0);
    Flags = Flags | (opt2parg_bSet("-ntomp", asize(pa), pa) ? MD_NTOMPSET : 0);
    Flags = Flags | (bIMDwait      ? MD_IMDWAIT      : 0);
    Flags = Flags | (bIMDterm      ? MD_IMDTERM      : 0);
    Flags = Flags | (bIMDpull      ? MD_IMDPULL      : 0);

    /* We postpone opening the log file if we are appending, so we can
       first truncate the old log file and append to the correct position
       there instead.  */
    if (MASTER(cr) && !bDoAppendFiles)
    {
        gmx_log_open(ftp2fn(efLOG, NFILE, fnm), cr,
                     Flags & MD_APPENDFILES, &fplog);
    }
    else
    {
        fplog = NULL;
    }

    ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
    ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
    ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);
    /* Disable MDStress here (there is no enable)
     */

    rc = gmx::mdrunner(&hw_opt, fplog, cr, NFILE, fnm, oenv, bVerbose,
                       nstglobalcomm, ddxyz, dd_rank_order, npme, rdd, rconstr,
                       dddlb_opt[0], dlb_scale, ddcsx, ddcsy, ddcsz,
                       nbpu_opt[0], nstlist,
                       nsteps, nstepout, resetstep,
                       nmultisim, repl_ex_nst, repl_ex_nex, repl_ex_seed,
                       pforce, cpt_period, max_hours, imdport, nstlocals,
                       localscontrib, localspbc, localsskip, Flags);

    /* Log file has to be closed in mdrunner if we are appending to it
       (fplog not set here) */
    if (MASTER(cr) && !bDoAppendFiles)
    {
        gmx_log_close(fplog);
    }

    return rc;
}

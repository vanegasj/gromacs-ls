#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# FindMDStress
#
# - Locate headers and libraries for the Molecular Dynamics Stress
#   (mdstress) library.
#
# Usage: find_package(MDStress)
#
# FindMDStress defines the following variables:
#
# MDSTRESS_FOUND          - True if both headers and libraries were located
# MDSTRESS_INCLUDE_DIRS   - Where to find mdstress headers
# MDSTRESS_LIBRARIES      - Libraries to link with to enable mdstress usage
#

message(STATUS "Looking for mds_stressgrid.h")
find_path(MDSTRESS_INCLUDE_DIRS "mdstress/mds_stressgrid.h")
if(MDSTRESS_INCLUDE_DIRS)
    message(STATUS "mds_stressgrid.h - found in ${MDSTRESS_INCLUDE_DIRS}")
else()
    message(STATUS "mds_stressgrid.h - not found")
endif()

message(STATUS "Looking for mdstress library")
find_library(MDSTRESS_LIBRARIES "mdstress")

find_package_handle_standard_args(MDStress
    REQUIRED_VARS MDSTRESS_LIBRARIES MDSTRESS_INCLUDE_DIRS)

mark_as_advanced(MDSTRESS_INCLUDE_DIRS MDSTRESS_LIBRARIES)


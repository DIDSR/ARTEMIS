ARTEMIS
=======

ARTEMIS (pArticle transport, Recombination, and Trapping in sEMiconductor Imaging Simulations) is a Monte Carlo package for X-ray, electron and electron-hole pair (EHP) transport. X-ray and secondary electron interactions in the presence of an electric field are modeled by PENELOPE 2006, and the locations of inelastic electron interactions are coupled in space and time to the transport routine for EHP simulation. ARTEMIS was developed at the U. S. Food and Drug Administration (FDA), Center for Devices and Radiological Health, Office of Science and Engineering Laboratories, Division of Imaging and Applied Mathematics. Please report to the authors any issue/bug that you may encounter.

The source code and documentation of ARTEMIS are openly distributed at the website: http://code.google.com/p/artemis/ . The following disclaimer notice applies to the code and documentation developed exclusively at the FDA (this disclaimer is provided at the beginning of each file developed at the FDA):

Code disclaimer
---------------

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 

Documentation and reference
---------------------------

The ARTEMIS code was first introduced in the paper listed below, which should be referenced by researchers using this code.

-Yuan Fang, Andreu Badal, Nicholas Allec, Karim S. Karim, and Aldo Badano, "Spatiotemporal Monte Carlo transport methods in x-ray semiconductor detectors: Application to pulse-height spectroscopy in a-Se", Medical Physics 39, pp. 308–319 (2012). 

ABSTRACT:

Purpose: The authors describe a detailed Monte Carlo (MC) method for the coupled transport of ionizing particles and charge carriers in amorphous selenium (a-Se) semiconductor x-ray detectors, and model the effect of statistical variations on the detected signal.

Methods: A detailed transport code was developed for modeling the signal formation process in semiconductor x-ray detectors. The charge transport routines include three dimensional spatial and temporal models of electron-hole pair transport taking into account recombination and trapping. Many electron-hole pairs are created simultaneously in bursts from energy deposition events. Carrier transport processes include drift due to external field and Coulombic interactions, and diffusion due to Brownian motion.

Results: Pulse-height spectra (PHS) have been simulated with different transport conditions for a range of monoenergetic incident x-ray energies and mammography radiation beam qualities. Two methods for calculating Swank factors from simulated PHS are shown, one using the entire PHS distribution, and the other using the photopeak. The latter ignores contributions from Compton scattering and K-fluorescence. Comparisons differ by approximately 2% between experimental measurements and simulations.

Conclusions: The a-Se x-ray detector PHS responses simulated in this work include three dimensional spatial and temporal transport of electron-hole pairs. These PHS were used to calculate the Swank factor and compare it with experimental measurements. The Swank factor was shown to be a function of x-ray energy and applied electric field. Trapping and recombination models are all shown to affect the Swank factor.
External software used by ARTEMIS: PENELOPE and penEasy

PENELOPE (version 2006) is a general purpose code that performs Monte Carlo simulation of coupled electron-photon transport in arbitrary materials and in the energy range from 50 eV to 1 GeV. The standard geometry model used by PENELOPE (PENGEOM) is based on defining objects as the volume limited by a set of quadric surfaces. Despite this model can be used to describe complex geometries, it is not adequate to represent biological structures with arbitrary shapes. The PENELOPE subroutines are copyrighted by the Universitat de Barcelona and can be obtained for free at http://www.nea.fr/abs/html/nea-1525.html or at http://www-rsicc.ornl.gov/ . A comprehensive description of the physical and geometrical models implemented in PENELOPE 2006 can be found at:

-F. Salvat, J. M. Fernandez-Varea and J. Sempau, "PENELOPE, A Code System for Monte Carlo Simulation of Electron and Photon Transport", Workshop Proceeding Issy-les-Moulineaux OECD/NEA 2006. ISBN: 92-64-02145-0. Document available at http://www.oecd-nea.org/dbprog/penelope.pdf . 

    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    C  PENELOPE/PENGEOM (version 2006)                                     C
    C  Copyright (c) 2001-2006                                             C
    C  Universitat de Barcelona                                            C
    C                                                                      C
    C  Permission to use, copy, modify, distribute and sell this software  C
    C  and its documentation for any purpose is hereby granted without     C
    C  fee, provided that the above copyright notice appears in all        C
    C  copies and that both that copyright notice and this permission      C
    C  notice appear in all supporting documentation. The Universitat de   C
    C  Barcelona makes no representations about the suitability of this    C
    C  software for any purpose. It is provided "as is" without express    C
    C  or implied warranty.                                                C
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

PenEasy is a general-purpose simulation package for PENELOPE developed by Josep Sempau. The package contains a modular main program and several tally options and source models that facilitate the simulation of medical physics applications. The modular structure of the code makes it easy to develop additional tools that extend the applicability of PENELOPE to new fields of study. The main program of ARTEMIS is essentially a custom version of the penEasy's main program. The penEasy package is copyrighted by the Universitat Politecnica de Catalunya and distributed at http://www.upc.es/inte/downloads/penEasy.htm .

-J. Sempau, A. Badal and L. Brualla, "A PENELOPE-based system for the automated Monte Carlo simulation of clinacs and voxelized geometries—application to far-from-axis fields," Medical Physics 38, pp. 5887 (2011). 

    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    c  penEasy                                                                     c
    c  Copyright (c) 2004-2008                                                     c
    c  Universitat Politecnica de Catalunya                                        c
    c                                                                              c
    c  Permission to use, copy, modify and re-distribute copies of this software   c
    c  or parts of it and its documentation for any purpose is hereby granted      c
    c  without fee, provided that this copyright notice appears in all copies.     c
    c  The Universitat Politecnica de Catalunya makes no representations about     c
    c  the suitability of this software for any purpose. It is provided "as is"    c
    c  without express or implied warranty.                                        c
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Code compilation and usage
--------------------------

To compile ARTEMIS, the user needs a current FORTRAN and C compiler. For the executable in this package, Intel FORTRAN and C compilers version 11.1 are used, but the program can be compiled also with GCC (g++ and gfortran). For plotting ARTEMIS results, the distribution includes example scripts for GNUPLOT, a command-driven plotting program.

If you are not using the executable provided with the distribution, ARTEMIS can be compiled using the provide compile.sh script file in the root directory of the distribution.

Note: PENELOPE 2006 source code files penelope.f and pengeom.f are not distributed with the ARTEMIS package, but are needed for compiling. If these files are needed, the contact information for the code developer is provided at the bottom of the page.

To begin an ARTEMIS simulation, the input file specifying the geometry, material data and code options has to be re-directed to the executable's standard input at execution time.

   $./artemis.x  < input.in > output.out &

Questions, bug reports, feature suggestions, etc, can be posted in the Issues section in this website. For further information, the users can also directly contact the code developer at the address: Yuan.Fang@fda.hhs.gov. 

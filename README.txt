##########################################################################################################################################################
#
# ****Disclaimer****
#  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
#  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
#  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
#  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
#  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
#  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
#  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
#  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are
#  derived from it, and any modified versions bear some notice that they have been modified.
#
#	@file    README.txt
#       @author  Yuan Fang (Yuan.Fang@fda.hhs.gov)
#       @date    Jan 02, 2014
#
##########################################################################################################################################################


*************
INTRODUCTION
*************

ARTEMIS (pArticle transport, Recombination, and Trapping in sEMiconductor Imaging Simulations) is a Monte Carlo package for X-ray, electron and 
electron-hole pair (EHP) transport. X-ray and secondary electron interactions in the presence of an electric field are modeled by PENELOPE 2006, 
and the locations of inelastic electron interactions are coupled in space and time to the transport routine for EHP simulation. ARTEMIS was developed 
at the U. S. Food and Drug Administration (FDA), Center for Devices and Radiological Health, Office of Science and Engineering Laboratories, Division 
of Imaging and Applied Mathematics. Please report to the authors any issue/bug that you may encounter. 

The source code and documentation of ARTEMIS are openly distributed at the website: http://code.google.com/p/artemis/ . The above disclaimer notice 
applies to the code and documentation developed exclusively at the FDA. 


*******************************
CODE COMPILATION AND EXECUTION
*******************************
To compile ARTEMIS the user needs a current FORTRAN and C compiler. 

To run ARTEMIS
	./ARTEMIS_v1.0.x < input.in 
	
Sample simulation inputs and outputs are included under the 'Demo' folder. For more details, read 'MANUAL_ARTEMIS.pdf'.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																//
//    NOTE: PENELOPE 2006 SOURCE CODE FILES PENELOPE.F AND PENGEOM.F ARE NOT DISTRIBUTED WITH THE ARTEMIS PACKAGE,	 	//
//    BUT ARE NEEDED FOR COMPILING. IF THESE FILES ARE NEEDED, PLEASE CONTACT YUAN FANG AT Yuan.Fang(at)fda.hhs.gov.	        //
//																//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


**********************************
artemis v1.0 PACKAGE CONTENTS
**********************************

	* EHP Transport SOURCE CODES: 	'transEHP.c' and 'transEHP.h' files. 
	
	* penEasy SOURCE CODES: 	All fortran files. Files 'penEasy_EMfield_EDEtally.f' and 'tallyEnergyDepositionEvents.f' were derived from original penEasy files. These files have been modified to output energy deposition events information.
	
	* Demo:				This folder contains sample input and output files for running simulations using artemis, for Mo and W input spectra and a range of histories.
	
	* GNUPLOT_script:		artemis gnuplot scripts to plot pulse-height spectrum.
	
	* input.in:			ARTEMIS input files.
	
	* MANUAL_ARTEMIS.pdf: 		Reference manual for artemis v1.0.
	
	* ARTEMIS_v1.0.x:		ARTEMIS v1.0 executable.
	
	* compile_v1.0.sh:		Script for compiling artemis_v1.0. The user will require two more files PENELOPE.F, PENGEOM.F to compile this package successfully. These files are not distributed with this package. If these files are needed, please contact Yuan Fang at Yuan.Fang(at) fda.hhs.gov.
	
	* README.txt:			this file.
	

For more details, read 'MANUAL_ARTEMIS.pdf'.
	

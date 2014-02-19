//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ****Disclaimer****
//  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
//  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
//  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
//  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
//  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
//  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
//  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
//  decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are
//  derived from it, and any modified versions bear some notice that they have been modified.
//
//	     @file    		transEHP.c
//       @Description	MC EHP transport main program
//       @author  		Yuan Fang (Yuan.Fang@fda.hhs.gov)
//       @date    		Jan 02, 2014
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "transEHP.h"

void init(double time){
        delT = time;
		Ldiff_e = sqrt(6*u_e*kb*T/q*delT); //maximum diffusion distance
        Ldiff_h = sqrt(6*u_h*kb*T/q*delT);
        du_e = delT*u_e; //diffusion parameters
        du_h = delT*u_h;
		p_e = exp(-1*delT/t_e); //trapping parameters
        p_h = exp(-1*delT/t_h);
}
void directrans_(double *x, double *y, double *z, double *ee, double *ed, int *nehp_return, double *nehpdet_return) {
	init(delTL);

	//initialize random number generator (Twister)
    unsigned long int seed = 224;
    init_genrand( seed );

	//sample # of electron-hole pairs

    int nEHP = poissonSampling(*ed);
    limit_z = *z; //interaction depth

	//calculate burst size & thermalization distance
	double burst = calcBurstSizeFWHM(*ee);
	double ro = cubicFunction( h*vp*vp/q/D_h, -q*ufield_.EFZ, Eg - *ed/nEHP, -q*q/(4*pi*e_o*e_Se) );
	//initialize the particle array
	struct particle testE[nEHP];
    struct particle *ptrE = testE;
    struct particle testH[nEHP];
    struct particle *ptrH = testH;
	particle_init_All(ptrH, ptrE, nEHP, *x, *y, *z, burst, ro);

	//initialize while loop conditions for transport simulation (recombination and trapping)
	still_AliveH = nEHP; //if all EHPs are not removed. check "remove" in particle struct
    still_AliveE = nEHP; //if all EHPs are not removed. check "remove" in particle struct

	while ((still_AliveH > 0) && (still_AliveE > 0)) {
		still_AliveH = 0; //if all EHPs are not removed. check "remove" in particle struct
		still_AliveE = 0;

		//reset the Ecoul elements of particle struct for new iteration
		//calculate the Coulombic Force between all elements
		//update the locations of all elements
		clear_Ecoul_All(ptrH,nEHP);
        clear_Ecoul_All(ptrE,nEHP);
		calc_Ecoul_All(ptrH, ptrE, nEHP, nEHP);
		updateLocation_All(ptrH,nEHP);
        updateLocation_All(ptrE,nEHP);

	}

	//initialize while loop conditions for transport simulation (trapping only)
	trappingInit(ptrH,ptrE,nEHP);

		if ((still_AliveH > 0) || (still_AliveE > 0)) {
		init(delTH);
		clear_Ecoul_All(ptrH,nEHP);
		clear_Ecoul_All(ptrE,nEHP);

    	while ((still_AliveH > 0) || (still_AliveE > 0)) {
            	still_AliveH = 0; //if all EHPs are not removed. check "remove" in particle struct
           		still_AliveE = 0; //if all EHPs are not removed. check "remove" in particle struct

           		//update the locations of all elements
           	 	updateLocation_All(ptrH,nEHP);
            	updateLocation_All(ptrE,nEHP);
    	}
	}

	//update number of created and detected EHPs from simulations
	*nehp_return = nEHP;
	*nehpdet_return = (double)(particles_detected(ptrH,nEHP)+particles_detected(ptrE,nEHP))/2;
}

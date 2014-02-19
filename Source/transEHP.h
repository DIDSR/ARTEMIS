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
//	     @file    		transEHP.h
//       @Description	Library definition file
//       @author  		Yuan Fang (Yuan.Fang@fda.hhs.gov)
//       @date    		Jan 02, 2014
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "twister.h"

#define e_o 8.854e-14 //permittivity of free space (F/cm)
#define e_Se 6.3 //permittivity of a-Se
#define q 1.6e-19 //elementary electric charge (Coulomb)
#define m 9.11e-31 //electron mass (kg)
#define h 6.626e-34 //Planck's const (Js)
#define u_e 3e-3 //electron drift mobility (cm^2/Vs)
#define u_h 1.2e-1 //hole drift mobility (cm^2/Vs)
//#define u_e 6.3e-3 //electron drift mobility (cm^2/Vs)
//#define u_h 1.9e-1 //hole drift mobility (cm^2/Vs)
#define kb 1.38065e-23 //boltzmann constant (J/K)
#define T 300 //temperature (Kelvin)
#define pi 3.1415926
#define omega_pe 9.7e15 //plasma frequency
#define vp 3.6221e12 //phonon frequency: 1e12 - 1e13 (Hz)
#define Eg 2.3 //band gap (eV)
#define HOLE 1 //electron type
#define ELECTRON -1 //hole type
#define radius_recomb 1e-7 //distance between EHPs where recombination happen, and particle removed
#define limit_recomb 5e-6 //distance between EHPs where recombination is considered could no longer happen
#define delTH 1e-9
#define delTL 1e-13

#define W 5.00
//#define t_e 1e-5 //high trapping 10us
//#define t_h 1e-4 //high trapping 100us
#define t_e 5e-5 //low trapping 50us
#define t_h 1e-3 //low trapping 1000us
#define topElectrode 0
#define pixElectrode 0.015 //cm = 150um

static double delT;
static double Ldiff_e; //maximum diffusion distance for electron  !!DeBuG!!
static double Ldiff_h; //maximum diffusion distnace for hole  !!DeBuG!!
static double p_e; //electron trapping probability
static double p_h; //hole trapping probability
static double du_e; //constant used for diffusion of electron
static double du_h; //constant used for diffusion of hole

double const_Coulomb = q/(4*pi*e_o*e_Se); //constant used for calculation of Ecoul
double D_e = u_e*kb*T/q; //electron diffisivity of a-Se
double D_h = u_h*kb*T/q; //hole diffisivity of a-Se
double limit_z;
int still_AliveH; //if all holes are not removed. check "remove" in particle struct
int still_AliveE; //if all electrons are not remove.

//Electric field parameters from input file
extern struct {
	double EFX, EFY, EFZ, BFX, BFY, BFZ;
}ufield_;

//definition of structure for particles (either electron or hole)
struct particle {
	int type;
	int remove;
	double x_coord;
	double y_coord;
	double z_coord;
	double x_Ecoul;
	double y_Ecoul;
	double z_Ecoul;
};

//initializatio of location of an electron-hole pair with a separation distance of r_thermal randomly in space
void location(struct particle *hole, struct particle *elec, double x, double y, double z, double burst, double ro) {
	//sample randomly two angles, theta (between 0 and pi, and gamma (between 0 and 2*pi)
	double theta1 = randMT(0,1)*pi;
	double gamma1 = randMT(0,1)*2*pi;
	double theta2 = randMT(0,1)*pi;
	double gamma2 = randMT(0,1)*2*pi;

	//determine burst location
	double tempBurst = burst;
	double x_temp = x + tempBurst*sin(theta1)*cos(gamma1);
	double y_temp = y + tempBurst*sin(theta1)*sin(gamma1);
	double z_temp = z + tempBurst*cos(theta1);

	//particle isotropically initialized 0.5*r_thermal away from origin (0,0,0)
	hole->x_coord = x_temp + 0.5*ro*sin(theta2)*cos(gamma2);
	hole->y_coord = y_temp + 0.5*ro*sin(theta2)*sin(gamma2);
	hole->z_coord = z_temp + 0.5*ro*cos(theta2);
	elec->x_coord = x_temp - 0.5*ro*sin(theta2)*cos(gamma2);
	elec->y_coord = y_temp - 0.5*ro*sin(theta2)*sin(gamma2);
	elec->z_coord = z_temp - 0.5*ro*cos(theta2);
}

//initialization of a particle (set everything to zero, except the type)
void particle_init (struct particle *a, int type) {
	a->type = type;
	a->remove = 0;
	a->x_coord = 0;
	a->y_coord = 0;
	a->z_coord = 0;
	a->x_Ecoul = 0;
	a->y_Ecoul = 0;
	a->z_Ecoul = 0;
}

//initialization of all particles
void particle_init_All (struct particle *arrayH, struct particle *arrayE, int array_size, double x, double y, double z, double burst, double ro) {
	int i;
	for (i=0; i<array_size; i++){
		//create hole and electrons alternatively
		particle_init(&arrayH[i], HOLE);
    	particle_init(&arrayE[i], ELECTRON);

    	//initialize location of EHP with distance r_o apart
    	location(&arrayH[i], &arrayE[i], x, y, z, burst, ro);
	}
}

//clear the Coulomb Electric Field
void clear_Ecoul (struct particle *a) {
		a->x_Ecoul = 0;
    	a->y_Ecoul = 0;
    	a->z_Ecoul = 0;
}

//clear all Coulomb Electric Field entries
void clear_Ecoul_All (struct particle *array, int array_size) {
		int i;
    	for (i=0; i<array_size; i++){
		clear_Ecoul(&array[i]);
	}
}

//calculate the coulombic forces between all particles
void calc_Ecoul (struct particle *a, struct particle *b) {
	double dist = sqrt(pow((a->x_coord-b->x_coord),2)+pow((a->y_coord-b->y_coord),2)+pow((a->z_coord-b->z_coord),2));

	//check recombination, if satisfied, set remove for both particles
	if (dist < radius_recomb) {
		a->remove = 2;
		b->remove = 2;
	} else {

		//if not recombination, calculate the Coulombic Force between them:
			double dist_cubed = pow(dist,3);
	    	double x_temp, y_temp, z_temp;

	    	if ( ((a->type != b->type) && (a->type == HOLE)) || (a->type == ELECTRON) ) {
		    x_temp = const_Coulomb*((b->x_coord - a->x_coord)/dist_cubed);
		    y_temp = const_Coulomb*((b->y_coord - a->y_coord)/dist_cubed);
		    z_temp = const_Coulomb*((b->z_coord - a->z_coord)/dist_cubed);
		} else {
		    x_temp = const_Coulomb*((a->x_coord - b->x_coord)/dist_cubed);
		    y_temp = const_Coulomb*((a->y_coord - b->y_coord)/dist_cubed);
		    z_temp = const_Coulomb*((a->z_coord - b->z_coord)/dist_cubed);
		}
	    	a->x_Ecoul = a->x_Ecoul + x_temp;
	    	a->y_Ecoul = a->y_Ecoul + y_temp;
	    	a->z_Ecoul = a->z_Ecoul + z_temp;
			b->x_Ecoul = b->x_Ecoul - x_temp;
			b->y_Ecoul = b->y_Ecoul - y_temp;
			b->z_Ecoul = b->z_Ecoul - z_temp;
	}
}

//calculate the Coulombic Force between all elements
void calc_Ecoul_All (struct particle *arrayH, struct particle *arrayE, int array_sizeH, int array_sizeE) {
	int i, j;
    	for (i=0; i<array_sizeH; i++){
		for (j=0; j<array_sizeE; j++){
			if ((arrayH[i].remove == 0)&&(arrayE[j].remove == 0)) {
				calc_Ecoul(&arrayH[i], &arrayE[j]);
    			}
    		}
	}
}

int trapping (double prob) {
	if ( randMT(0,1) > prob ) {
		return 1;
	} else {
		return 0;
	}
}

void updateLocation (struct particle *a) {
	if ((a->remove == 1) || (a->remove == -1)) {

		//do nothing
    	} else if (a->remove == 2) {
	 	//if remove==2, (it means particle was removed from previous history
		//but not completely, because it had to be considered for the Ecoul calculation from prev. history
		//set remove=1, so that it is truely removed
		a->remove = 1;
    	} else {

	    //update the locations
	    if (a->type==ELECTRON){
	    	//Two components: drift due to Coulombic field and electric fied, and random diffusion
	    	a->x_coord = a->x_coord + du_e*a->x_Ecoul + Ldiff_e*randMT(-1,1);
	    	a->y_coord = a->y_coord + du_e*a->y_Ecoul + Ldiff_e*randMT(-1,1);
	    	a->z_coord = a->z_coord + du_e*(a->z_Ecoul-ufield_.EFZ) + Ldiff_e*randMT(-1,1);

			if ((delT==delTL ) && (a->z_coord < (limit_z-limit_recomb)) ) { //escape recombination
				a->remove = -1;
			}
			if ((delT==delTH) && (a->z_coord <= topElectrode) ){
				a->remove = -1; //should be another code for this
			}
			if ( trapping(p_e) ) {
				//printf("update Elec: trapped");
				a->remove = 1;
			}
			//check if the particle is removed, if so, return 0. if not, return 1.
			if (a->remove == 0){
				still_AliveE++;
			}
		} else if (a->type == HOLE) {
		    //There are two components: drift due to Coulombic field and electric fied, and random diffusion
			a->x_coord = a->x_coord + du_h*a->x_Ecoul + Ldiff_h*randMT(-1,1);
		        a->y_coord = a->y_coord + du_h*a->y_Ecoul + Ldiff_h*randMT(-1,1);
		        a->z_coord = a->z_coord + du_h*(a->z_Ecoul+ufield_.EFZ) + Ldiff_h*randMT(-1,1);

			if ((delT==delTL ) && (a->z_coord > (limit_z+limit_recomb))) { //escape recombination
				a->remove = -1;
			}
                        if ((delT==delTH) && (a->z_coord >= pixElectrode)) {
                                a->remove = -1; //should be another code for this
                        }
			if ( trapping(p_h) ) {
                                //printf("update Hole: trapped");
				a->remove = 1;
            }

            //check if the particle is removed.
			if (a->remove == 0){
				still_AliveH++;
			}
		 }
	 }
}

//update the locations of all elements
void updateLocation_All (struct particle *array, int array_size) {
	int i;
    	for (i=0; i<array_size; i++){
	    updateLocation(&array[i]);
	}
}

//loop through particle array, check if all the particles are removed or not
int particles_detected (struct particle *array, int array_size) {
	int size = array_size;
	int i;
	int result = 0;
	for (i=0; i<size; i++){
		if ( (array[i].remove == -1) || (array[i].remove == 0) ){
			result = result + 1;
		}
	}
	return result;
}

//initialize parameters for trapping
void trappingInit(struct particle *arrayH, struct particle *arrayE, int array_size) {
	still_AliveH = 0;
        still_AliveE = 0;
        int i;
        for (i=0; i<array_size; i++){
                if ((arrayH[i].remove == -1) || (arrayH[i].remove == 0)){
                        still_AliveH = still_AliveH + 1;
                        arrayH[i].remove = 0;
                }
                if ((arrayE[i].remove == -1) || (arrayE[i].remove == 0)) {
                        still_AliveE = still_AliveE + 1;
                        arrayE[i].remove = 0;
                }
        }
}

/* Functions needed for simulatio with PENELOPE output such as:
 * -Sampling of burst size (Gaussian FWHM)
 * -Poisson sampling
 * -cubic function used for calculation of thermalization distance
 */
double gaussianSampling (double sigma) {
	double result = box_muller(0, sigma);
	if (result < 0) {
		result = result * -1;
	}
	return result;
}
//Sampling of burst size (Gaussian FWHM)
double calcBurstSizeFWHM(double ke) {
	return gaussianSampling(sqrt(2*ke*q/m)*1e2/omega_pe);
}

int poissonSampling (double ed) {
	double edEhp = 0;
	int nEHP;
        while (edEhp < Eg) {
                nEHP = (int)poidev((float)ed/W);
                edEhp = ed/nEHP;

                if (nEHP == 1) {
                        break;
                }
        }
	return nEHP;
}

double cubicFunction(double A, double B, double C, double D) {
	double p, qq, DD, phi, temp1, temp2, y1, y2, y3, u, v, x;

	//Cubic equation with 3 roots
    	int nroot = 3;

	// Step 1: Calculate p and q
    	p  = C/A - B*B/A/A/3. ;
    	qq  = (2.*B*B*B/A/A/A - 9.*B*C/A/A + 27.*D/A) / 27. ;

	// Step 2: Calculate DD (discriminant)
    	DD = p*p*p/27. + qq*qq/4. ;

	// Step 3: Branch to different algorithms based on DD
    	if (DD < 0.) {
       		//Step 3b:
       		//3 real unequal roots -- use the trigonometric formulation
        	phi = acos(-qq/2./sqrt(fabs(p*p*p)/27.));
        	temp1=2.*sqrt(fabs(p)/3.);
        	y1 =  temp1*cos(phi/3.);
        	y2 = -temp1*cos((phi+pi)/3.);
        	y3 = -temp1*cos((phi-pi)/3.);
	} else {
		//Step 3a:
		//1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
        	temp1 = -qq/2. + sqrt(DD);
        	temp2 = -qq/2. - sqrt(DD);
        	//u = abs(temp1)^(1./3.);
        	//u = 1/(abs(temp1)*abs(temp1)*abs(temp1));
        	u = pow(fabs(temp1), 1./3.);

        	//v = abs(temp2)^(1./3.);
        	//v = 1/(abs(temp2)*abs(temp2)*abs(temp2));
        	v = pow(fabs(temp2), 1./3.);

        	if (temp1 < 0.) {
	        	u=-u;
        	}
        	if (temp2 < 0.) {
	        	v=-v;
        	}
        	y1  = u + v;
        	//y2r = -(u+v)/2.;
        	//y2i =  (u-v)*sqrt(3.)/2.;
	}

	// Step 4: Final transformation
    	temp1 = B/A/3.;
    	y1 = y1-temp1;
    	if (DD < 0.) {
        y2 = y2-temp1;
        y3 = y3-temp1;
  	} else {
        //y2r=y2r-temp1;
	}

	// Assign answer
    	if (DD < 0.) {
		//x(1) = y1;
		//x(2) = y2;
		//x(3) = y3;
        //three distinct real roots
        if ((y2 >= y1) && (y2 >= y3)) {
            return y2;
    	} else if ((y3 >= y2) && (y3 >= y1)) {
            return y3;
    	} else {
            return y1;
    	}
	} else if (DD == 0.) {
		//x(1) =  y1;
		//x(2) = y2r;
		//x(3) = y2r;
        //three real roots, at least two equal
        if ((y2 >= y1) && (y2 >= y3)) {
            return y2;
    	} else if ((y3 >= y2) && (y3 >= y1)) {
            return y3;
    	} else {
            return y1;
    	}
  	} else {
		//x(1) = y1;
		//x(2) = y2r + y2i*i;
		//x(3) = y2r - y2i*i;
        //only one real root, two imaginary
        return y1;
	}
}

 /* reverse:  reverse string s in place */
void reverse(char s[]) {
     int i, j;
     char c;

     for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
         c = s[i];
         s[i] = s[j];
         s[j] = c;
     }
}

/* itoa:  convert n to characters in s */
void itoa(int n, char s[]) {
     int i, sign;

     if ((sign = n) < 0)  /* record sign */
         n = -n;          /* make n positive */
     i = 0;
     do {       /* generate digits in reverse order */
         s[i++] = n % 10 + '0';   /* get next digit */
     } while ((n /= 10) > 0);     /* delete it */
     if (sign < 0)
         s[i++] = '-';
     s[i] = '\0';
     reverse(s);
}









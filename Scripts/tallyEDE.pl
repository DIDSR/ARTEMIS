##########################################################################################################################################################
##
## ****Disclaimer****
##  This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in
##  the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection
##  and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software
##  without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the
##  Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
##  parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality,
##  reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory
##  decisions. Although this software can be redistributed and#or modified freely, we ask that any derivative works bear some notice that they are
##  derived from it, and any modified versions bear some notice that they have been modified.
##
##	     @file    		tallyEDE.pl
##       @Description	Perl script file for pulse-height output and Swank noise calculations
##       @author  		Yuan Fang (Yuan.Fang@fda.hhs.gov)
##       @date    		Jan 02, 2014
##
##########################################################################################################################################################


#!/usr/bin/perl -w
use POSIX;

$a=0;
$b=0;
$c=0;
@efficiencyArray=();
@detEHPArray=();

$energy=100e3; #eV
$maxEHP=$energy*1.2/5;
$stepEHP=1;
$arraySize=$maxEHP/$stepEHP;
$totalCount=0;
for ($i = 0; $i < $arraySize; $i++) {
	$phsEHP[$i]=0;
	$phsEHPT[$i]=0
}

print ("Code Running..\n");

$fname = $ARGV[0];
open (OUTFILE2, ">$fname".".phs") || die "Can't.\n";
open (OPENFILE,$fname) || die "Can't Open File: $fname\n";
	while (<OPENFILE>)  # While still input lines in the file...
		{ 
			$totalCount = $totalCount + 1;
			($a,$b,$c)=split; 
			#print OUTFILE1 "$a $b $c \n";
			if($c == 0){
				push(@efficiencyArray,0);
				$phsEHP[0]=$phsEHP[0]+1;
			} else {
				push(@efficiencyArray,$c/$b);
				$phsEHP[ceil($c/$stepEHP)]=$phsEHP[ceil($c/$stepEHP)]+1;
			}
			$phsEHPT[ceil($b/$stepEHP)]=$phsEHPT[ceil($b/$stepEHP)]+1;
			push(@detEHPArray,$c);
		}
	close (OPENFILE);   # Close the file }
	
	#mean calculation for percentage of EHPs detected 
	$mean = 0;
	for ($i = 0; $i < @efficiencyArray; $i++) {
		$mean = $mean + $efficiencyArray[$i];
	}
	$mean = $mean/@efficiencyArray;

	#standard deviation calculation for percentage of EHPs detected
	$stdev = 0;
	for ($i = 0; $i < @efficiencyArray; $i++) {
                $stdev = $stdev + ($efficiencyArray[$i]-$mean)**2;
        }
	$stdev = sqrt($stdev/@efficiencyArray);
	#print OUTFILE2 "#mean = ", $mean, ", standard deviation = ", $stdev, ".\n";

        #mean calculation for # of EHPs detected
        $mean = 0;
        for ($i = 0; $i < @detEHPArray; $i++) {
                $mean = $mean + $detEHPArray[$i];
        }
        $mean = $mean/@detEHPArray;

        #standard deviation for # of EHPs detected
        $stdev = 0;
        for ($i = 0; $i < @detEHPArray; $i++) {
                $stdev = $stdev + ($detEHPArray[$i]-$mean)**2;
        }
        $stdev = sqrt($stdev/@efficiencyArray);
        print OUTFILE2 "#mean = ", $mean, ", standard deviation = ", $stdev, ".\n";
	
	$M0 = 0;
	$M1 = 0;
	$M2 = 0;
	#$M0T = 0;
	#$M1T = 0;
	#$M2T = 0;
	for ($i = 0; $i < $arraySize; $i++) {
		if ($phsEHP[$i]==0) {
			$temp = 0;
		} else {
			$temp = $phsEHP[$i]/$totalCount;
		}
		$phsEHP[$i]=$temp; 
		
		if ($phsEHPT[$i]==0) {
                        $temp = 0;
                } else {
                        $temp = $phsEHPT[$i]/$totalCount;
                }
                $phsEHPT[$i]=$temp;

		$M0 = $M0 + $phsEHP[$i];
		$M1 = $M1 + $phsEHP[$i]*($i*$stepEHP);
		$M2 = $M2 + $phsEHP[$i]*($i*$stepEHP)*($i*$stepEHP);
		#$M0T = $M0T + $phsEHPT[$i];
                #$M1T = $M1T + $phsEHPT[$i]*($i*$stepEHP);
                #$M2T = $M2T + $phsEHPT[$i]*($i*$stepEHP)*($i*$stepEHP);


	}
	$I = $M1*$M1/$M0/$M2;
	#$IT = $M1T*$M1T/$M0T/$M2T;
	print OUTFILE2 "#Swank results: I = ", $I, "M0 = ", $M0, " M1 = ", $M1, " M2 = ", $M2, ".\n";
	#print OUTFILE2 "#M0 = ", $M0T, " M1 = ", $M1T, " M2 = ", $M2T, " I = ", $IT, ".\n";

	for ($i = 0; $i < $arraySize; $i++) {
        print OUTFILE2 $i*$stepEHP, " ", $phsEHP[$i]*$totalCount, "\n";
    }
close (OUTFILE2);	

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
##	     @file    		PHS.gnu
##       @Description	GNUPLOT file for plotting pulse-height results
##       @author  		Yuan Fang (Yuan.Fang@fda.hhs.gov)
##       @date    		Jan 02, 2014
##
##########################################################################################################################################################

# Gnuplot MS-Windows 32 bit version 4.0
#! tables

#set zero 1.0e-35
set key right top

#set term postscript enhanced color lw 4 "Helvetica" 23
#set output "phs.eps"

set xlabel 'Electron-Hole Pairs'
set ylabel 'Count'

set xrange [0:1500]
set yrange [0:*]
plot 'tallyEDE-Mo4V1e4.dat.phs' u 1:2 title "Mo Spectra: 4V/{/Symbol m}m 1e4 Hist" w l lc 3, \
 	 'tallyEDE-Mo30V1e4.dat.phs' u 1:2 title "Mo Spectra: 30V/{/Symbol m}m 1e4 Hist" w l lc 4
pause -1

plot 'tallyEDE-Mo4V1e5.dat.phs' u 1:2 title "Mo Spectra: 4V/{/Symbol m}m 1e5 Hist" w l lc 3, \
 	 'tallyEDE-Mo30V1e5.dat.phs' u 1:2 title "Mo Spectra: 30V/{/Symbol m}m 1e5 Hist" w l lc 4
pause -1

plot 'tallyEDE-Mo4V1e6.dat.phs' u 1:2 title "Mo Spectra: 4V/{/Symbol m}m 1e6 Hist" w l lc 3, \
 	 'tallyEDE-Mo30V1e6.dat.phs' u 1:2 title "Mo Spectra: 30V/{/Symbol m}m 1e6 Hist" w l lc 4
pause -1

plot 'tallyEDE-W4V1e4.dat.phs' u 1:2 title "W Spectra: 4V/{/Symbol m}m 1e4 Hist" w l lc 3, \
 	 'tallyEDE-W30V1e4.dat.phs' u 1:2 title "W Spectra: 30V/{/Symbol m}m 1e4 Hist" w l lc 4
pause -1

plot 'tallyEDE-W4V1e5.dat.phs' u 1:2 title "W Spectra: 4V/{/Symbol m}m 1e5 Hist" w l lc 3, \
 	 'tallyEDE-W30V1e5.dat.phs' u 1:2 title "W Spectra: 30V/{/Symbol m}m 1e5 Hist" w l lc 4
pause -1

plot 'tallyEDE-W4V1e6.dat.phs' u 1:2 title "W Spectra: 4V/{/Symbol m}m 1e6 Hist" w l lc 3, \
 	 'tallyEDE-W30V1e6.dat.phs' u 1:2 title "W Spectra: 30V/{/Symbol m}m 1e6 Hist" w l lc 4
pause -1

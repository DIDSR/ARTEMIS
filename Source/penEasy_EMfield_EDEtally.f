CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
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
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!*******************************************************************
!*                         PENEASY                                 *
!*                                                                 *
!* Short description:                                              *
!*   General-purpose PENELOPE main program. Please refer to the    *
!*   README.txt file for detailed instructions.                    *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> routines CLEANS,START,SECPAR                               *
!*   from PENVARED:                                                *
!*   -> routine JUMPF                                              *
!*   from other penEasy files:                                     *
!*   -> routines in penaux.f, penpatch.f, penvox.f,                *
!*      sourceXX.f, tallyXX.f and timing.f                         *
!*   In particular, these routines supersede their PENELOPE        *
!*   equivalents:                                                  *
!*   -> KNOCKX (necessary for the fluence tally) replaces KNOCK;   *
!*      it is included in penpatch.f                               *
!*   -> STEPX (from penvox.f) tracks particle trajectories in      *
!*      quadric+voxel geometries, replacing STEP; it is included   *
!*      in penvox.f                                                *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2006                                                          *
!*                                                                 *
!* Josep Sempau, email: josep.sempau@upc.es                        *
!* Universitat Politecnica de Catalunya, Barcelona, Spain          *
!* SEE COPYRIGHT NOTICE IN FILE README.txt                         *
!*******************************************************************

!*******************************************************************
!*    Includes                                                     *
!*******************************************************************
      ! PENELOPE routines:     !!DeBuG!! PENELOPE files included in Makefile to avoid recompilation
!       include 'penelope.f'
!       include 'pengeom.f'
!       include 'penvared.f'

      ! Auxiliary routines:
      include 'penaux.f'
      include 'penpatch.f'
      include 'penvox.f'
      include 'timing.f'

      ! Source models (see documentation for a detailed description):
      include 'sourceBoxIsotropicGaussSpectrum.f'
      include 'sourcePhaseSpaceFile.f'
      ! <you may add your own here and in the source routines below>

      ! Tallies (see documentation for a detailed description):
      include 'tallyVoxelDoseDistrib.f'
      include 'tallySpatialDoseDistrib.f'
      include 'tallyCylindricalDoseDistrib.f'
      include 'tallySphericalDoseDistrib.f'
      include 'tallyEnergyDepositionPulseHeightSpectrum.f'
      include 'tallyFluenceTrackLength.f'
      include 'tallyPhaseSpaceFile.f'
      include 'tallyParticleCurrentSpectrum.f'
      include 'tallyParticleTrackStructure.f'
      ! <you may add your own here and in the tally routines below>
      

!*******************************************************************
!*    MAIN                                                         *
!*******************************************************************
      program main
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical endsim,absorb,forcing,isforcing
      integer*4 ncross,icol,left
      real*8 n,ds,dsef,de,dsmax

C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      real*8 ULDE, ULDV, ULEM, TMAX, emax_track  !!DeBuG!!  EM field
              
      real*8 Zmax, Zmin                     !!DeBuG!!  EM field: Z extend of the field, passed from main to GETEMF
      COMMON/UFIELD_Z/ Zmin, Zmax
      
C  ****  EM field.                          !!DeBuG!!  EM field
      real*8 EFX,EFY,EFZ,BFX,BFY,BFZ
      COMMON/UFIELD/EFX,EFY,EFZ,BFX,BFY,BFZ
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


      write(*,'(a)')' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)')
     & '>>>> This is penEasy v.2008-06-15 >>>>'
      write(*,'(a)')                                  !!DeBuG!!  EM field
     & '>>>>   !!EXTENDED TO EM FIELDS!!  >>>>'       !!DeBuG!!  EM field
      write(*,'(a)')                                  !!DeBuG!!  EM field
     & '>>>> !Voxel geometry disregarded! >>>>'       !!DeBuG!!  EM field
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)')' '

      call initime ! Write date on the screen
      
     
C>>>>>>>>> EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  !!DeBuG!!  EM field
C  (Code based on the penelope package example 'pm-field.f')
C  ************  Initializing the electromagnetic field routines.         !!DeBuG!!  EM field
C
      ! Init EM field to 0:
      EFX=0.0D0
      EFY=0.0D0
      EFZ=0.0D0
      BFX=0.0D0
      BFY=0.0D0
      BFZ=0.0D0
      ! Read the Z component of the E field from standard input:
      write(*,'(a)')" ** ELECTRIC FIELD INITIALIZATION:"
      read(*,*)  ! Skip title line in input
      write(*,'(a)')"    -- Uniform field along Z: EFZ [V/cm] = "
      read(*,*) EFZ
      write(*,'(5x,1pe12.5)') EFZ
      write(*,'(a)')"    -- Z interval: [Zmin,Zmax] = "
      read(*,*) Zmin, Zmax
      write(*,'(5x,1pe12.5,2x,1pe12.5)') Zmin, Zmax
      write(*,'(a)')"    -- Simulation accuracy (default 0.02): "    
      read(*,*) ULDV
      write(*,'(5x,1pe12.5)') ULDV  ! ULDV == Upper limit on the amount of deflection over the step due to the em-field
      ULDE=ULDV                   ! ULDE == Upper limit on the relative energy variation over the step in the em-field
      ULEM=ULDV                   ! ULEM == Upper limit on the amount of change of the em-field over the step

      emax_track = 0.0d0   !!DeBuG!!  EM field: tally the maximum electron energy in the simulation (may be higher than initial due to electric field acceleration).
C<<<<<<<<< EM field >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  !!DeBuG!!  EM field



      call init    ! Initialize the penEasy+PENELOPE+PENGEOM system
      call treset  ! Reset simulation timer
      n = 0.0d0    ! Reset history counter

      history: do              ! Each iteration simulates a new history
        n = n+1.0d0            ! Update history counter
        call noverflow(n)      ! Check that N does not overflow
        call cleans            ! Empty the stack, just in case
        call tally(1,n)        ! The simulation of this history begins
        call source(n)         ! Put primary particles (from the same history) in stack

        particle: do                       ! Each iteration simulates a new particle
          call secpar(left)                ! Retrieve a particle from the stack
          if (left.eq.0) exit particle     ! Stack was empty
          call tally(-99,-e)               ! The simulation of this particle begins

          if(kpar.eq.1 .and. e.gt.emax_track) emax_track=e   !!DeBuG!!  EM field: tally the maximum electron energy in the track

          if (absorb()) cycle particle     ! Check particle absorption
          call start                       ! Reset transport mechanics
          forcing = isforcing()            ! Set interaction forcing (variance reduct)

          interact: do                     ! Each iteration simulates an interaction
            if (absorb()) exit interact    ! Check particle absorption


            call TPEMF0(ULDV,ULDE,ULEM,TMAX)     !!DeBuG!!  EM field. Output: TMAX = maximum allowed step length.
            if (forcing) then
              call jumpf(min(TMAX,dsmax()),ds)   !!DeBuG!!  EM field.
                !! call jumpf(dsmax(),ds)       
            else
              call jump(min(TMAX,dsmax()),ds)    !!DeBuG!!  EM field.
                !! call jump(dsmax(),ds)        
            endif
            call TPEMF1(ds,dsef,ncross)          !!DeBuG!!  EM field. TPEMF1 will call step; voxels not considered
              !! call stepx(ds,dsef,ncross)

            if(kpar.eq.1 .and. e.gt.emax_track) emax_track=e   !!DeBuG!!  EM field: the electric field may increase the electron's energy!
              
            if (absorb()) exit interact          !!DeBuG!!  EM field: e- moving against the field will reduce energy, eventually below Eabs.

            if (ncross.eq.0) then
              call tally(3,ds)             ! Moved a distance DS, no interface crossed
            else
              call tally(4,dsef)           ! Moved a distance DSEF, interface found
              if (mat.eq.0) exit interact  ! New material is vacuum => gone
              call start                   ! New material => reset transport mechanics
              forcing = isforcing()        ! Set interaction forcing (variance reduct)
              cycle interact
            endif
            if (forcing) then
              call knockfx(de,icol)        ! Interaction forcing (see PENELOPE manual)
            else
              call knockx(de,icol)         ! Simulate an interaction
            endif
            call tally(-int(icol),de)      ! Tally kinetic energy released
          enddo interact
        enddo particle

        call tally(6,n)                    ! End-of-history bookkeeping

        if (endsim(n)) exit history        ! Simulation is finished
      enddo history

      call report(n)                       ! Write final report

      write(*,'(a)')"  "
      write(*,'(a)')
     & ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,'(a)')"report: Electromagnetic field"
      write(*,'(a)')"  Maximum electron energy in the simulation [eV]: "
      write(*,'(5x,1pe12.5)') emax_track
      write(*,'(a)')
     & ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"      
      write(*,'(a)')"  "

      end


      subroutine init
!*******************************************************************
!*    Initializes                                                  *
!*******************************************************************
      implicit none
      integer*4 nmat
      real*8 emax,realtime,cputime

      call iniconfig              ! Simulation config
      call inisource(emax)        ! Source models
      call inigeo(nmat)           ! Geometry: PENGEOM & penVox


      !!!emax = 100.0d0*emax  !!DeBuG!!  EM field:  increase emax: e- will get energy from the E field!
      
      call inipen(emax,nmat)      ! PENELOPE
      call initally               ! Tallies
      call iniforce(emax)         ! Interaction forcing
      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'init: INITIALIZATION ENDED'
      write(*,'(a,f9.2,a)') 'Elapsed real time:',realtime(),' s'
      write(*,'(a,f9.2,a)') 'Elapsed CPU time :',cputime(),' s'
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      end


      subroutine report(n)
!*******************************************************************
!*    Reports final results                                        *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*******************************************************************
      implicit none
      integer unc
      real*8 n,cputime,realtime,nowcpu,nowcpu2
      integer*4 seed1,seed2
      common/rseed/seed1,seed2

      nowcpu = cputime()
      call tallyreport(n,nowcpu,unc)

      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'report: SIMULATION ENDED'
      write(*,'(a)')
     & 'Results have been written to the corresponding DAT files.'
      select case(unc)
      case(0)
        write(*,'(a)')
     &   'The requested uncertainty has NOT been reached.'
      case(2)
        write(*,'(a)')
     &   'The requested uncertainty has been reached.'
      end select
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      ! Write generic report to the screen:
      write(*,*) ''
      write(*,'(a)') 'Last random seeds:'
      write(*,'(2(1x,i0))') seed1,seed2

      write(*,'(a)') 'Elapsed real time (s), excluding init:'
      write(*,'(1x,1pe12.5)') realtime()

      nowcpu2 = cputime()
      write(*,'(a)') 'Elapsed CPU time (s), excluding init:'
      write(*,'(1x,1pe12.5)') nowcpu2

      write(*,'(a)') 'Each report update took (in CPU s):'
      write(*,'(1x,1pe12.5)') nowcpu2-nowcpu

      write(*,'(a)') 'No. of histories simulated:'
      write(*,'(1x,f18.0)') n

      if (nowcpu.gt.0.0d0) then
        write(*,'(a)') 'CPU Speed (histories/s):'
        write(*,'(1x,1pe12.5)') n/nowcpu
      endif

      call endtime  ! Report timings
      end


!*******************************************************************
!*******************************************************************
!*    Source routines start here.                                  *
!*    Add your own models or delete the unwated ones.              *
!*    Source models usually require:                               *
!*     i) an initialization routine that must be called by         *
!*        INISOURCE                                                *
!*    ii) a particle generation routine that must be called        *
!*        by SOURCE                                                *
!*******************************************************************
!*******************************************************************

      subroutine inisource(emax)
!*******************************************************************
!*    Init routines for source models                              *
!*                                                                 *
!*    Output:                                                      *
!*      emax -> max source energy (eV)                             *
!*******************************************************************
      implicit none
      real*8 emax
      logical active
      integer nsrc

      nsrc = 0
      call BIGSinisrc(active,emax)
      if (active) nsrc = nsrc+1
      call PSFinisrc(active,emax)
      if (active) nsrc = nsrc+1

      if (nsrc.ne.1) then
        write(*,*) ''
        write(*,'(a)')
     &    'inisource:ERROR: There must be ONE active source'
        stop
      endif
      end


      subroutine source(n)
!*******************************************************************
!*    Source models                                                *
!*                                                                 *
!*    Input:                                                       *
!*      n -> top history counter                                   *
!*******************************************************************
      implicit none
      real*8 n

      call BIGSsource(n)
      call PSFsource(n)
      end


!*******************************************************************
!*******************************************************************
!*    Tally routines start here.                                   *
!*    Add your own tallies or delete the unwated ones              *
!*    Tallies usually require:                                     *
!*    i) an initialization routine that must be called by INITALLY *
!*    ii) a tally routine that must be called by TALLY             *
!*    iii) a reporting routine that must be called by TALLYREPORT  *
!*                                                                 *
!*    Notice that the ordering of the tally initialization routines*
!*    must coincide with the ordering of the corresponding sections*
!*    in the input file.                                           *
!*******************************************************************
!*******************************************************************

      subroutine initally
!*******************************************************************
!*    Init tallying routines.                                      *
!*                                                                 *
!*    Comments:                                                    *
!*      -> VDDinitally sets variables that are needed for          *
!*         proper particle transport in voxelized geometries.      *
!*         Therefore, it should not be deactivated even if the     *
!*         tally is not used.                                      *
!*******************************************************************
      implicit none
      call VDDinitally
      call SDDinitally
      call CDDinitally
      call SPDinitally
      call EPSinitally
      call FTLinitally
      call PSFinitally
      call PCSinitally
      call PTSinitally

      call EDEinitally  !!DeBuG!! Calling the new EDE tally!

      end


      subroutine tally(mode,arg)
!*******************************************************************
!*    Tallying routines.                                           *
!*                                                                 *
!*    Comments:                                                    *
!*      -> VDDtally sets state variables that are needed for       *
!*         proper particle transport in voxelized geometries.      *
!*         Therefore, it should not be deactivated even if the     *
!*         tally is not used.                                      *
!*      -> Furthermore, these variables could be used by other     *
!*         tallies and, in consequence, VDDtally should be the     *
!*         first tally to be called.                               *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      call VDDtally(mode,arg)  ! Must be 1st tally to be called
      call SDDtally(mode,arg)
      call CDDtally(mode,arg)
      call SPDtally(mode,arg)
      call EPStally(mode,arg)
      call FTLtally(mode,arg)
      call PSFtally(mode,arg)
      call PCStally(mode)
      call PTStally(mode,arg)

      call EDEtally(mode,arg)  !!DeBuG!! Calling the new EDE tally!
      
      end


      subroutine tallyreport(n,cputim,unc)
!*******************************************************************
!*    Calls report routines for all tallies                        *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      simtim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      unc -> larger than 1 if requested uncert has been reached  *
!*******************************************************************
      integer unc,uncdone
      real*8 n,cputim

      ! Write partial reports to corresponding data files:
      unc = 1
      call VDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call SDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call CDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call SPDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call EPSreport(n,cputim,uncdone)
      unc = unc*uncdone
      call FTLreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PSFreport(n,cputim)           ! No uncertainty for this tally
      call PCSreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PTSreport                     ! No arguments for this tally

      call EDEreport(n)  !!DeBuG!! Calling the new EDE tally! 
      
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
!*                          TALLY                                  *
!*                 ENERGY DEPOSITION EVENTS                        *
!*                                                                 *
!* Short description:                                              *


!*                                                                 *
!*   Writes to a file the position and value of energy loss events *
!*   so that particle tracks can be visualized.                    *
!*                                                                 *
!* Dependencies:                                                   *
!*   from PENELOPE:                                                *
!*   -> common /TRACK/                                             *
!*   -> common /RSEED/                                             *
!*   from other penEasy files:                                     *
!*   -> routine GETLINE,FINDUF                                     *
!*                                                                 *
!* Compatible with PENELOPE versions:                              *
!*   2005,2006                                                     *
!*******************************************************************


      subroutine EDEtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      implicit none
      integer*4 mode, i, nehp
      real*8 arg, de, nehpdet

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)

C  ****  Secondary stack (PENELOPE).
      integer*4 KS,IBODYS,MS,ILBS,NSEC,NMS
      real*8 ES,XS,YS,ZS,US,VS,WS,WGHTS
      PARAMETER (NMS=2000)
      COMMON/SECST/ES(NMS),XS(NMS),YS(NMS),ZS(NMS),US(NMS),
     1   VS(NMS),WS(NMS),WGHTS(NMS),KS(NMS),IBODYS(NMS),MS(NMS),
     2   ILBS(5,NMS),NSEC

      logical active
      integer edeunit
      integer*4 num_ede, num_real, mat_detector, secondaries_bef_knock,
     &          totnehp
      real*8 energy_before_knock, total_edep_history, distance_step,
     &       totnehpdet
      common /comEDE/ energy_before_knock, total_edep_history,
     &                distance_step, totnehpdet, totnehp,
     &                secondaries_bef_knock,
     &                num_ede, num_real, mat_detector, edeunit, active

      if (.not.active) return

      if (mode.eq.0) then
        ! ** New primary stored by the SOURCE at the secondary stack:
        total_edep_history = 0.0d0    ! Reset the history edep counter
        totnehp = 0
        totnehpdet = 0.0d0
      endif
          
      if (mode.eq.-99) then
        ! ** Particle track begins for next primary or secondary:
        if (ilb(3).eq.5) then
          ! -- This particle was created in an inner shell ionization event: do not simulate.   !!DeBuG!!
          e = 0.0
!        else
!          write(edeunit,'(a)')" "     ! Separate particles
        endif  

      else if (mode.eq.3) then
        ! ** An interaction will take place:
        secondaries_bef_knock = NSEC   ! Store the number of particles in secondary stack
        energy_before_knock = e
        distance_step = arg
        
      else if (mode.lt.0     .and.
     &          arg.gt.0.0d0 .and.           ! (Neglect elastic and inner shell events)
     &          mat.eq.mat_detector) then
        ! ** An interaction happened with energy deposition inside the detector:
        num_ede = num_ede + 1
        
        ! Subtracting the energy of the new secondaries from ARG:
        de = arg
        do i = (secondaries_bef_knock+1), NSEC
          de = de - ES(i)   ! Subtract the energy of the new secondaries            
        enddo

        total_edep_history = total_edep_history + de
  
      !- if the interaction created secondary particles but did not deposit energy !!DeBuG!!
      if (de.gt.0.0d0) then !

	call directrans(x, y, z, energy_before_knock, de, nehp, nehpdet)
       
        totnehp = totnehp + nehp
        totnehpdet = totnehpdet + nehpdet

      endif
     
          
      else if (mode.eq.6 .and. total_edep_history.gt.0.0d0) then
   		num_real = num_real + 1
   		
        ! ** End-of-history bookkeeping
        write(edeunit,'(1pe16.4,1x,i8,1x,1pe16.4)')
     &        total_edep_history,totnehp,totnehpdet

!         call report_globals_debug   !!DeBuG!!

      endif
      
      end


      subroutine EDEreport(n)
!*******************************************************************
!*    Report number of deposition events per primary               *
!*******************************************************************
      implicit none

      real*8 n

      logical active
      integer edeunit
      integer*4 num_ede, num_real, mat_detector, secondaries_bef_knock,
     &          totnehp
      real*8 energy_before_knock, total_edep_history, distance_step,
     &       totnehpdet
      common /comEDE/ energy_before_knock, total_edep_history,
     &                distance_step, totnehpdet, totnehp,
     &                secondaries_bef_knock,
     &                num_ede, num_real, mat_detector, edeunit, active

      write(*,'(a)') " "
      write(*,'(a)')   ">>>>>>>>> EDEreport >>>>>>>>>"
      write(*,'(a,i8)')"    Number of energy deposition"//
     &                 " events written to the output file: ", num_ede
      write(*,'(a,1pe12.5)')"    Depositions per primary = ",
     &                           dble(num_ede)/n
     
      end


      subroutine EDEinitally
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      implicit none
      
      logical active
      integer edeunit
      integer*4 num_ede, num_real, mat_detector, secondaries_bef_knock,
     &          totnehp
      real*8 energy_before_knock, total_edep_history, distance_step,
     &       totnehpdet
      common /comEDE/ energy_before_knock, total_edep_history,
     &                distance_step, totnehpdet, totnehp,
     &                secondaries_bef_knock,
     &                num_ede, num_real, mat_detector, edeunit, active

      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY ENERGY DEPOSITION EVENTS v.2010-06-03]')
      parameter (eos='[END OF EDE SECTION]')
      character*80 buffer
      integer finduf,error

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'EDEinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(1x,a3)') buffer
      if (adjustl(buffer(1:3)).eq.'ON') then
        active = .true.
      else if (buffer(1:3).eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Particle Track Structure is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'EDEinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'EDEinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer(1:3)
        stop
      endif

      ! Init counters:
      total_edep_history  = 0.0d0
      energy_before_knock = 0.0d0
      secondaries_bef_knock = 0
      num_ede = 0
      num_real = 0

      write(*,'(a)') ' Detector sensitive material number:'
      read(*,*) mat_detector
      write(*,'(2x,i3)') mat_detector

      write(*,'(a)') ' EDE output file name:'
      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      write(*,'(1x,a)') buffer

      edeunit = finduf()
      open(edeunit,file=buffer,iostat=error)
      if (error.ne.0) then
        write(*,'(a)')
     &    'EDEinitally:ERROR: cannot open output data file:'
        write(*,'(a80)') buffer
        stop
      endif
      
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'EDEinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> EDE tally initialization finished >>>>'
      
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


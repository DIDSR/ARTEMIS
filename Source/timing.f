!*******************************************************************
!*                         TIMING                                  *
!*                                                                 *
!* Short description:                                              *
!*   Time routines complying with the F95 standard. They have been *
!*   tested on the fortran compilers g77, g95, Compaq Visual       *
!*   Fortran Pro 6.1.0 and Intel Fortran.                          *
!*                                                                 *
!* Dependencies:                                                   *
!*   The following implicits are used:                             *
!*   - date_and_time(char*8 sdate,char*10 stime,char*5 zone,       *
!*                   int values(8))                                *
!*     where                                                       *
!*       sdate = YYYYMMDD                                          *
!*       stime = hhmmss.sss                                        *
!*       zone = hhmm , difference with respect to UTC              *
!*       values(1) = year,                                         *
!*       values(2) = month (1..12)                                 *
!*       values(3) = day (1..31)                                   *
!*       values(4) = difference with respect to UTC in minutes     *
!*       values(5) = hour of the day (0..23)                       *
!*       values(6) = minutes of the hour                           *
!*       values(7) = seconds of the minute                         *
!*       values(8) = milliseconds of the second                    *
!*   - cpu_time(real time) (in s)                                  *
!*                                                                 *
!* Last update:                                                    *
!*   2008-03-01 JS                                                 *
!*******************************************************************


      subroutine initime
!*******************************************************************
!*    Initializes the timers (real and CPU) and writes the current *
!*    date to stdout in text format. It should be called at the    *
!*    beginning of the main program.                               *
!*******************************************************************
      implicit none
      character*100 fdate
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call treset
      call date_and_time(sdate,stime,zone,values)
      call dateString(values,fdate)
      write(*,*)
      write(*,'(a,a)') 'Program timer started on ',fdate
      end


      subroutine endtime
!*******************************************************************
!*    Writes current date to stdout in text format. Should be      *
!*    called just before your program ends.                        *
!*******************************************************************
      implicit none
      character*100 fdate
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)
      ! real*8 realtime,cputime

      call date_and_time(sdate,stime,zone,values)
      call dateString(values,fdate)
      write(*,*)
      write(*,'(a,a)') 'Program ended on ',fdate
      ! write(*,*)
      ! write(*,'(a,f12.2,a)') 'Real time:',realtime(),' s'
      ! write(*,'(a,f12.2,a)') 'CPU time :',cputime(),' s'
      write(*,*)
      write(*,'(a)') 'Have a nice day.'
      end


      subroutine treset
!*******************************************************************
!*    Resets the real and CPU timers                               *
!*******************************************************************
      implicit none
      real utime0
      real*8 rtime0,lasthour,tshift
      common /tim001/ lasthour,tshift,rtime0,utime0
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call date_and_time(sdate,stime,zone,values)
      lasthour = values(5)
      rtime0 = (((values(5)*60.0d0)+values(6))*60.0d0)+values(7)+
     &            1.0d-3*values(8)
      tshift = 0.0d0
      call cpu_time(utime0)
      if (utime0.lt.0.0) then
        write(*,*) 'treset: WARNING: intrinsic time routine failed'
        write(*,*) '        Timings will not be reliable.'
      endif
      end


      real*8 function cputime()
!*******************************************************************
!*    Returns CPU (user) time elapsed since last call to TRESET    *
!*                                                                 *
!*    Output:                                                      *
!*      CPU time in s                                              *
!*******************************************************************
      implicit none
      real utime0
      real*8 rtime0,lasthour,tshift
      common /tim001/ lasthour,tshift,rtime0,utime0
      real time

      call cpu_time(time)
      if (time.lt.0.0) then
        write(*,*) 'cputime: WARNING: intrinsic time routine failed'
        write(*,*) '         Timings will not be reliable.'
      endif
      cputime = time-utime0
      end


      real*8 function realtime()
!*******************************************************************
!*    Returns the real time elapsed since last call to TRESET      *
!*                                                                 *
!*    Output:                                                      *
!*      Real time in s                                             *
!*    Comments:                                                    *
!*      -> MUST be called at least once per day or it will not     *
!*         return reliable data                                    *
!*******************************************************************
      implicit none
      real utime0
      real*8 rtime0,lasthour,tshift
      common /tim001/ lasthour,tshift,rtime0,utime0
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call date_and_time(sdate,stime,zone,values)
      if (dble(values(5)).lt.lasthour) tshift = tshift+86400.0d0  ! The day changed
      lasthour = values(5)
      realtime = (((values(5)*60.0d0)+values(6))*60.0d0)+values(7)+
     &            1.0d-3*values(8)+tshift-rtime0
      end


      subroutine dateString(values,fdate)
!*******************************************************************
!*    Returns a string with the date and time                      *
!*                                                                 *
!*    Input:                                                       *
!*      values(8) -> date array in standard Fortran (see header)   *
!*    Output:                                                      *
!*      fdate -> formatted date (and time) in a char string        *
!*******************************************************************
      implicit none
      integer values(8)
      character*100 fdate
      character*3 month(12)

      month(1)  = 'Jan'
      month(2)  = 'Feb'
      month(3)  = 'Mar'
      month(4)  = 'Apr'
      month(5)  = 'May'
      month(6)  = 'Jun'
      month(7)  = 'Jul'
      month(8)  = 'Aug'
      month(9)  = 'Sep'
      month(10) = 'Oct'
      month(11) = 'Nov'
      month(12) = 'Dec'
      write(fdate,'(i2,1x,a3,1x,i4,2x,i2.2,a1,i2.2,a1,i2.2)')
     &  values(3),month(values(2)),values(1),values(5),':',
     &  values(6),':',values(7)
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

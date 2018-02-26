*     +NUSTAR_LIMIT
      subroutine NUSTAR_LIMIT( results,n_res)
      implicit none

      integer n_res
      real results(0:n_res)


*     Description: 
*     Calculates exposure time for a 3 and 5 sigma detection
*     assuming a background rate 0.3
*     
*     Arguments:
*     results (i)  : Predicted count rates
*     n_res   (i)  : Number of standard bands
*     
*     Dependencies:
*     PWRITE
*
*     Major re-write following discussion with Kristin Madsen et al.
*     of NuSTAR team, using 2013 August simulation response & background
*
*     Modified 2014 July 24 - changed upper limit from 80 to 79 keV
*     and emphasized the two-module nature of the count rates, for v4.7a
*
*     Modified 2014 July 25 - cleaned up dead time fraction display for >99%
*     
*     -NuSTAR_LIMIT

      integer r_index,n_bands
      parameter (r_index=0)

      parameter (n_bands=5)

      real DT
      parameter (DT=0.0025)

      real frac_limit
      parameter (frac_limit=0.15)
      
      real dtadjust
      parameter( dtadjust = 1.0 )
c     To calculate the dead time, we want the total count rate
c     per telescope.  PIMMS estimate counts in the standard
c     extraction region (containing 50% of the PSF) but for both
c     telescopes, so these cancel out.  This parameter allows future
c     changes in this.  For the purpose of dead-time calculation,
c     we're currently ignoring the background, which is probably a
c     good approximation.


c      real countrate_lim
c      parameter (countrate_lim=60)
      
      integer highlimit_flag(0:n_bands-1)
      integer highlimit_flag_dt(0:n_bands-1)

      character*80 nu_strng
      character*200 info_strng

      character*8 bounds(0:n_bands-1)

      integer counter

      real bgd_rate(0:n_bands-1)
      real cxb_rate(0:n_bands-1)
      real results_dt(0:n_bands-1)

      real tot_counts, dt_frac
c      real dt_percent(0:n_bands-1)
      real dt_percent

      integer sigs(2)
      
c      integer printmsg
c      integer highlimit_flag_all
c      integer dt_percent_all

      real temp,work
      
      data bounds /'3-79','3-10','10-20','20-40','40-79'/
      
c      data bgd_rate /1.876e-03,6.539e-04,3.676e-04,4.323e-04,
c     &     4.133e-04/
c      data bgd_rate / 4.739e-3, 1.393e-3, 5.960e-4, 1.430e-3, 1.320e-3 /
c                  the above for 80 keV upper bound; the below, 79 keV
      data bgd_rate / 4.715e-3, 1.393e-3, 5.960e-4, 1.430e-3, 1.296e-3 /

ccc Add a statment if the count rate across the entire band is > 60
c      if (results(0) .gt. 60) then
c         write(info_strng(1:200),'(a)')'* Count rate > 60 c/s'//
c     &        ' adjusted  by 7/5'
c         call PWRITE(info_strng)
c         
c         
c         write(info_strng(1:200),'(a,f8.3,a)')'* PIMMS '//
c     & 'predicts ', results(0)*7/5,' cps with NUSTAR'
c         call PWRITE(info_strng)
c      end if

      
      write(info_strng(1:200),'(a)')' '  
      call PWRITE(info_strng)


      write(info_strng(1:200),'(a)')'* Count rates in different '//
     & 'energy bands'
      call PWRITE(info_strng)





ccc Modify c/s if the number is greater than 60 c/s
c      do counter=0,n_bands-1
c         if (results(counter) .gt. 60) then
c            results(counter) = results(counter)*7/5
c         end if
c      end do   



cc    Calculate the dead time count rate fraction for the whole band 3-79 keV
cc    This resulting ratio will be applied to all sub bands

c     The total count rate per telescope, ignoring background, is exactly
c     as calculated by PIMMS - notionally, multiply by 2 to get the
c     total source count rate (inc. photons outside the 50% PSF extraction
c     region), then divide by 2 to get the per telescope number.
c     It would be more accurate to add the total background rate per telescope
c     but that's probably a small effect.
      tot_counts = results( 0 )
      dt_frac = 1.0 / ( 1.0 + tot_counts * dtadjust * DT )
      dt_percent = ( 1.0 - dt_frac ) * 100.0
      
c      results_dt( 0 ) = results( 0 ) / ( 1 + ( results (
c     &     0 ) * DT) )
c
c      dt_frac = results_dt( 0 ) / results( 0 )
c
c      dt_percent(0) = (results(0) - results_dt(0))/results(0)




cc    Whether or not to print the warning that the countrate
cc    is greater than 60 c/s
c      printmsg = 0
c      highlimit_flag_all = 0
c      dt_percent_all = 0

      do counter=0,n_bands-1

cc    Must calculate the actual recorded rate corrected for dead time
cc    Rr = Rt / (1 + Rt*DT)    ==> DT is dead time, Rt is true rate
cc    Rr is recorded rate
cc    c/s with dead time most likely calculated for 3-79 keV band
cc    Take this ratio and apply it to all of the sub-bands rather than
cc    calculating the dead time for each band individually


         results_dt( counter ) = results( counter ) * dt_frac

         
         
c         dt_percent(counter) = (results(counter) - 
c     &        results_dt(counter)) /results(counter)
c
c         
c         if (dt_percent(counter) .gt. 0.15) then
c            dt_percent(counter) = 1
c            dt_percent_all = 1
c            printmsg = 1
c         else
c            dt_percent(counter) = 0
c         end if

cc    If the count rate is higher than 60 c/s, then set the high limit flag
cc    This means that the extraction radius needs to be adjusted by 7/5         

ccc       Without dead time
c         if (results(counter) .gt. countrate_lim) then
c            highlimit_flag(counter) = 1
c            highlimit_flag_all = 1
c         else
c            highlimit_flag(counter) = 0
c         end if

cc       *With* dead time
c         if (results_dt(counter) .gt. countrate_lim) then
c            highlimit_flag_dt(counter) = 1
c            highlimit_flag_all = 1 
c         else
c            highlimit_flag_dt(counter) = 0
c         end if


      end do
ccccccccccccccccccccccc


c      if (highlimit_flag_all .eq. 1) then
c
c         write(info_strng(1:200),'(a)')'* Count rate > 60 c/s '//
c     &        'adjusted by 7/5' 
c         call PWRITE(info_strng)
c      end if


         
c      if (printmsg .eq. 1) then
         write(info_strng(1:200),'(a)')'* The dead time value '//
     &        ' assumed is 0.0025 sec'
         call PWRITE(info_strng)

         info_strng =
     &               '* The dead time fraction is estimated to be 12.3%'
         if( dt_percent .ge. 0.99 ) then
           info_strng( 45: 48 ) = '> 99'
         else
           write( info_strng( 45: 48 ), '(f4.1)' ) dt_percent
         end if
         call PWRITE( info_strng )

         info_strng = ' '
c         write(info_strng(1:200),'(a,a)')'                                  '//
c     &        '                                                                  '

         call PWRITE(info_strng)
         
         call writeheader(nu_strng)
c      else
c         call writeheader_wo(nu_strng)
c      end if
         
      do counter=0,n_bands-1
cc Write count rates (without dead time) to the screen
         call writecounts(bounds(counter),results(counter),
     &        bgd_rate(counter),results_dt(counter),nu_strng )
c     &        highlimit_flag_dt(counter),highlimit_flag(counter),
c     &        dt_percent(counter))

      end do


cc    Choosing the entire band to determine whether or not the signal
cc    is detectable at 3, 5 sigma
      counter = 0

cc Start the calculation for time to reach an 'n'-sigma detection

cc   results(counter-1) goes with bgd_rate(counter) and bounds(counter)
cc   square number of counts
         temp = results_dt( counter ) * results_dt( counter )

         
cc     If the square of the counts is high enough, the calculate the time
cc     to reach 3, 5 sigma detection         
         if (temp .gt. 1.0e-30) then
            call calctime(temp,results_dt( counter ),
     &           bgd_rate( counter ) * dt_frac, bounds(counter),
     &           DT,results(counter))
            
         else
            nu_strng = '  The source is undetectable at 3-sigma'
            call PWRITE(nu_strng)
         end if 



c      if (dt_percent_all .eq. 1) then
c         write(info_strng(1:200),'(a)')'* WARNING! The dead time'//
c     &' fraction is greater than 15%'
c         call PWRITE(info_strng)
c      end if

      end
      




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeheader(nu_strng)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      character*80 nu_strng

      write( nu_strng(1:4),'(a)')''
      write( nu_strng(5:17),'(a)')'band'
      write( nu_strng(18:33), '(a)') 'counts/s'
      write( nu_strng(34:45),'(a)') 'counts/s'
      write( nu_strng(46:60), '(a)') 'counts/s'
      write( nu_strng(61:80),'(a)')' '
      
      
      call PWRITE(nu_strng)

      nu_strng = '                 ---------    both modules   --------'
      call PWRITE( nu_strng )
      
      write( nu_strng(1:4),'(a)')''
      write( nu_strng(5:17),'(a)')'(keV)'
      write( nu_strng(18:33), '(a)') 'no dead time'
      write( nu_strng(34:45),'(a)') 'background'
      write( nu_strng(46:60), '(a)') 'with dead time'
      write( nu_strng(61:80),'(a)')' '
      call PWRITE(nu_strng)
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine writeheader_wo(nu_strng)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c      character*80 nu_strng
c
c      write( nu_strng(1:4),'(a)')''
c      write( nu_strng(5:17),'(a)')'band'
c      write( nu_strng(18:33), '(a)') 'counts/s'
c      write( nu_strng(34:45),'(a)') 'counts/s'
c      write( nu_strng(46:60), '(a)') '      '
c      write( nu_strng(61:80),'(a)')' '
c      
c      
c      call PWRITE(nu_strng)
c      
c      
c      
c      write( nu_strng(1:4),'(a)')''
c      write( nu_strng(5:17),'(a)')'(keV)'
c      write( nu_strng(18:33), '(a)') 'no dead time'
c      write( nu_strng(34:45),'(a)') 'background'
c      write( nu_strng(46:60), '(a)') '  '
c      write( nu_strng(61:80),'(a)')' '
c      call PWRITE(nu_strng)
c      
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writecounts(bounds,cts,bgd_rate,cts_dt,nu_strng)
c      subroutine writecounts(bounds,cts,bgd_rate,cts_dt,nu_strng,
c     & highlimit_flag_dt,highlimit_flag,dt_percent_flag)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      character*8 bounds
      character*80 nu_strng
      real cts,cts_dt,bgd_rate
c      integer highlimit_flag
c      integer highlimit_flag_dt
c      real dt_percent_flag

      
      write( nu_strng(1:4),'(a)')''
      write( nu_strng(5:17),'(a)')bounds
      write( nu_strng(18:33), '(1p,e9.3)') cts



      write( nu_strng(34:45),'(1p,e9.3)') bgd_rate


C  Current recommendation: always report dead time 
cc Only print out dead time if the dead time is greater than 15%
c      if (dt_percent_flag .eq. 1) then
         write( nu_strng(46:62), '(1p,e11.5)') cts_dt
c      else if (dt_percent_flag .eq. 0) then
c         write( nu_strng(46:62), '(a)') ' '
c      end if



      write( nu_strng(63:80),'(a)')' '

      call PWRITE(nu_strng)


      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calctime(temp,cts_dt,bgd_rate,bounds,DT,cts)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      character*8 bounds
c                   unused, (3-80 keV) hardwired
      character*80 nu_strng
      real temp,cts_dt, bgd_rate,t_3,t_5,work3,work5,DT,cts
      integer sig


cc    Three sigma detection limit
      sig = 3
      
      work3 = sig * sig * (cts_dt + bgd_rate)

cc    Five sigma detection limit
      sig = 5
      work5 = sig * sig * (cts_dt + bgd_rate)



cc  Is the source detectable at 5-sigma?


c     NOTE: work5, in the absence of systematic error term, is always positive
      if (work5 .lt. 0.0) then

cc If not, print B  that to the screen and check to see if there
cc could be a three sigma detection
         nu_strng = '  The source is undetectable at 5-sigma'
      else

cc If it is detectable at 5- sigma, print B  how long it will take
cc to reach such a detection and also how long to reach a three
cc detection

         t_5 = work5 / temp
         write(nu_strng(1:41),'(a)')'The source is detectable at  '//
     &  '5 sigma in'
         write(nu_strng(42:54),'(1p,e11.5)') t_5
         write(nu_strng(55:80),'(a)')'seconds (3-79 keV)'
         call PWRITE(nu_strng)
      end if


c     NOTE: work3, in the absence of systematic error term, is always positive
cc Is the source detectable at 3 sigma?         
      if (work3 .lt. 0.0) then
         nu_strng = '  The source is undetectable at 3-sigma'
         call PWRITE(nu_strng)
      else
cc If the source is detectable at 3 sigma, how long until this 
cc detection can be reached?
         t_3 = work3 / temp
c         t_5 = work5 / temp

         write(nu_strng(1:41),'(a)')'The source is detectable at  '//
     &  '3 sigma in'
         write(nu_strng(42:54),'(1p,e11.5)') t_3
         write(nu_strng(55:80),'(a)')'seconds (3-79 keV)'

         call PWRITE(nu_strng)


      endif




      
      end



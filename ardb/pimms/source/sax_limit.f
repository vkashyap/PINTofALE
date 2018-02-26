*+SAX_LIMIT
        subroutine SAX_LIMIT( special, results, n_res)

*       Description:
*         Given a point source count rate calculates the exposure time
*         for 5 sigma detection assuming a background rate.
*
*       Arguments:
*         special      (i)   : flag for instrument 
*                               401 LECS
*                               402 MECS
*                               403 PDS
*                               404 HPGSPC
*                               405 WFC
*         results      (i)   : count rate of the point source 
*
*       Origin:
*
*
*-SAX_LIMIT
        implicit none
*Import
        integer special, n_res
        real results( 0: n_res )
*Local
        real time, f, sigma, back, back_low, back_total, back_high
        real cosmic, back_med
        character*80 xw_strng
        character*7 cband
*
        sigma=5.
*
        if (special.eq.401) then
* LECS
*------
cccc version 1996
* 0.1-10 keV internal background is obtained from 
* (1.6e-4 count/s/arcmin**2) X  (pi*6*6arcmin)
* where 6 arcmin is the radius of the region containing 
* 80 % of the counts at 0.284 keV
*
c           back_total=0.0181
c           cosmic=8.1e-3
c           back_total=back_total+cosmic
*
* 0.1-1 keV internal background is obtained from 
* (1.6e-5 count/s/arcmin**2) X  (pi*6*6arcmin)
* where 6 arcmin is the radius of the region containing 
* 80 % of the counts at 0.284 keV
*
           back_low=1.81e-3
           cosmic=2.e-3
           back_low=back_low+cosmic
c
           f=0.80
cccc
c
c version 31 Aug 1991
c
c region of 8 arcmin radius.
c from blanck sky observations:
c total background 0.1-10 keV = 2.5E-02           
c 0.1-1  keV background = 0.005 

           back_total=0.025
           back_low=0.005

           f=0.99 
           results(1)=results(1)*f
c
           call PWRITE( ' ')
           call PWRITE( '* Results in 2 SAX LECS bands are:')
           call PWRITE( '* Band    Source     BGD      5-sigma')
           call PWRITE( '* keV      cps       cps     detection (s)')
           call PWRITE( ' ')
*
           cband=' 0.1-10'
           time= (sigma**2)*(results(0)+back_total)/results(0)**2
           write( xw_strng, 401) cband,results(0), back_total, time
           call PWRITE( xw_strng( : 60 ))
c
           cband=' 0.1-1'
           time= (sigma**2)*(results(1)+2.*back_low)/results(1)**2
           write( xw_strng, 401) cband, results(1), back_low, time
           call PWRITE( xw_strng( : 60 ))
*
           call PWRITE( ' ')
           call PWRITE( '* Assumed 8 arcmin radius '//
     &                  '(95 % source region at 0.284 keV)')
           call PWRITE( '* 0.1-10 keV bgd rate of 0.025 cps '//
     &                  'for 1 LECS')
*
        elseif (special.eq.402) then
*
* MECS
*------
cccc version 1996
c
* 2-10 keV internal background is obtained from 
* (1.4e-4 count/s/arcmin**2) X  (pi*2.5*2.5arcmin)
* where 2.5 arcmin is the radius of the region containing 
* 80 % of the counts at 6.4 keV
*
           back=2.83e-3
           cosmic=1.6e-3
           back=back+cosmic
*
           f=0.8
ccc end version 1996

c  version 31 Aug 1997  !!!!2MECS ONLY, mecs2+mecs3!!!
c
c region of 4 arcmin radius
c from blanck sky observations:
c total background 2-10 keV = 0.0058          
c total background 2-10 keV 3 mecs units = 0.0087

           back=0.0058

 
           f=0.95
           results(0)=results(0)*f
*
           time= (sigma**2)*(results(0) + back)/results(0)**2
*
*     Format 201 tweaked eliminate a warning message - KM, 2001 Oct 12
           write( xw_strng, 201 ) time
201        format( '* An exposure of ', f10.2,
     &                         's is required for a 5-sigma detection' )
           call PWRITE( xw_strng( : 66 ) )
           call PWRITE( ' ')
           call PWRITE( '* Assumed 4 arcmin radius '//
     &                  '(90 % source region at 5 keV)')
           call PWRITE( '* 2-10 keV bgd rate of 0.0058 cps '//
     &                  'for 2 MECS')
* 
       elseif ( special.eq.403) then
*
*PDS
*------
* The Effective area is for four photswich. 
* Total geometric area of 795 cm2
* NOTE the background value is for four detectors.
*
ccc version 1996
* 15-300 keV total background of 32.8 counts s-1
           back_total=32.8
* 15-30 keV background of 5.4 counts s-1
           back_low=5.4
* 30-80 keV background of 5.0 counts s-1
           back_med=5.0
* 80-300 keV background of 22.4 counts s-1
           back_high=22.4
ccc end version 1996

c
c version 31 Aug 1997
c
* 13-200 keV total background of 31.5 counts s-1
           back_total=31.5
* 13-30 keV background of 7.4 counts s-1
           back_low=7.4
* 30-80 keV background of 12 counts s-1
           back_med=12.0
* 80-200 keV background of 11.2 counts s-1
           back_high=11.2

*
c           time= (sigma**2)*(results(0)+2.*back_total)/results(0)**2
c           write( xw_strng, 201 ) time
c           call PWRITE( xw_strng( : 50 ) ) 

           call PWRITE( ' ')
           call PWRITE( '* Results in 4 SAX PDS bands are:')
           call PWRITE( '* Band    Source     BGD      5-sigma')
           call PWRITE( '* keV      cps       cps     detection (s)')
           call PWRITE( ' ')
c
           cband=' 13-200'
           time= (sigma**2)*(results(0)+2.*back_total)/results(0)**2
           write( xw_strng, 401) cband,results(0), back_total, time
           call PWRITE( xw_strng( : 60 ))
c
           cband=' 13-30'
           time= (sigma**2)*(results(1)+2.*back_low)/results(1)**2
           write( xw_strng, 401) cband, results(1), back_low, time
           call PWRITE( xw_strng( : 60 ))
c
           cband=' 30-80'
           time= (sigma**2)*(results(2)+2.*back_med)/results(2)**2
           write( xw_strng, 401) cband, results(2), back_med, time
           call PWRITE( xw_strng( : 60 ))
c
           cband=' 80-200'
           time= (sigma**2)*(results(3)+2.*back_high)/results(3)**2
           write( xw_strng, 401) cband, results(3), back_high, time
           call PWRITE( xw_strng( : 60 ))
c
           call PWRITE( ' ')
           call PWRITE( '* !!!NOTE!!! that using the PDS in the '//
     &                  'standard rocking mode half of the')
           call PWRITE( '  observing time is allocated to monitor '//
     &                  'the bgd')
        elseif ( special.eq.404) then
*
*HPGSPC
*------
* Estimated internal background of 0.002 counts cm-2 s-1 keV-1 
* Total geometric area of 450 cm2
* 4-120 keV total background of 120 counts s-1
           back_total=120.
* 4-34 keV background of 50 counts s-1
           back_low=50.
* 34-120 keV background of 70 counts s-1
           back_high=70.
*
c           time= (sigma**2)*(results(0)+2.*back_total)/results(0)**2
c           write( xw_strng, 201 ) time
c           call PWRITE( xw_strng( : 50 ) ) 

           call PWRITE( ' ')
           call PWRITE( '* Results in 3 SAX HPGSPC bands are:')
           call PWRITE( '* Band    Source     BGD      5-sigma')
           call PWRITE( '* keV      cps       cps     detection (s)')
           call PWRITE( ' ')
c
           cband=' 4-120'
           time= (sigma**2)*(results(0)+2.*back_total)/results(0)**2
           write( xw_strng, 401) cband,results(0), back_total, time
           call PWRITE( xw_strng( : 60 ))
c
           cband=' 4-34'
           time= (sigma**2)*(results(1)+2.*back_low)/results(1)**2
           write( xw_strng, 401) cband, results(1), back_low, time
           call PWRITE( xw_strng( : 60 ))
c
           cband=' 34-120'
           time= (sigma**2)*(results(2)+2.*back_high)/results(2)**2
           write( xw_strng, 401) cband, results(2), back_high, time
           call PWRITE( xw_strng( : 60 ))
c
           call PWRITE( ' ')
           call PWRITE( '* !!!NOTE!!! that using the HPGSPC in the '//
     &                  'standard rocking mode half of the')
           call PWRITE('  observing time is allocated to monitor '//
     &                 'the bgd')

        elseif ( special.eq.405) then
*
*WFC
*------
*
           call PWRITE( 'no Special for now')
        endif
c
 401       format(a7,f7.3,3x,f7.3,3x,f10.2)
c
        end

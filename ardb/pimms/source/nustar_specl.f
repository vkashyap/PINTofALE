        else if ( special .eq. 1001 ) then
           call PWRITE('* NOTE: The count rates are estimated for two'//
     &                ' modules with a 50%')
           call PWRITE('        extraction region')
           call PWRITE('* The count rate does not account for '//
     &                  'deadtime')


          if( restr ) then
             call PWRITE('* PIMMS provides NuSTAR background and '//
     &           'deadtime information if run without ')
             call PWRITE('  specifying a limited output energy range')

c            call PWRITE( '% Further instrument-specific information ' //
c     &                      'can be obtained by re-running PIMMS' )
c            call PWRITE( '  without specifying a limited ' //
c     &                                    '(non-default) energy range' )

          else
           call NUSTAR_LIMIT( results,n_res)
          end if

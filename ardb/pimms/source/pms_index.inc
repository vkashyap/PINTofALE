*       This is the include file defining the Mission/Detector/Filter index
*       system used within PIMMS.
*
*       Modified - 1997 Dec, added the parameter ang_kev
*       Last modified - 2012 Oct, increased the number of missions allowed

        integer max_len, max_index
        integer max_mssn, max_dtct, max_fltr
        real ang_kev
        parameter( max_len = 32 )
        parameter( max_index = 4096 )
        parameter( max_mssn = 128 )
        parameter( max_dtct = 32 )
        parameter( max_fltr = 64 )
*            Maximum number of mission allowed is 32, maximum number of
*            detectors for any single mission is 8, maximum number of
*            filters per detector is 16 --- that's what it says, but
*            looks like all it matters is the combined total not to exceed
*            512.  Each name can be up to 32 characters long.
        parameter( ang_kev = 12.39854 )
*            If lambda = l A, E = ( ang_kev / l ) keV

        character*( max_len ) nm_array( max_index )
        common / pms_mndx0 / nm_array
*            This common stores all the names, mission, detectors and filters

        integer ln_array( max_index )
        integer st_array( max_index )
        common / pms_mndx1 / ln_array, st_array
*            This common stores length of the names and "status" associated
*            with each --- for mission/detector, it's the actual number of
*            detectors/filters associated with it, for filters it's the
*            encoded status for things like if it has an effecitve area curve.

        integer ind_fltr( max_index )
        integer ind_dtct( max_index )
        integer ind_mssn( max_mssn )
        integer n_mssn
        common / pms_mndx2 / ind_fltr, ind_dtct, ind_mssn, n_mssn
*            A three-tiered index array structure so that you can recover
*            these information.

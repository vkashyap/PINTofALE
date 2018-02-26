function mixie,lstr,sstr,lwdt=lwdt,nsig=nsig,lthr=lthr,rmfstr=rmfstr,mixstr=mixstr,$
 contamfrac=contamfrac,victim=victim,culprit=culprit,contamint=contamint,idtag=idtag,$
 vidtag=vidtag,vint=vint,mproxy=mproxy,mplet=mplet,verbose=verbose,$
 nhne=nhne,abund=abund,dem=dem,dlogt=dlogt, obsflx=obsflx,ofsig=ofsig, brnstr=brnstr,$
 abndupdt=abndupdt,multproxy=multproxy,multplet=multplet,multiproxy=multiproxy,$
 multiplet=multiplet,_extra=e

;+
;function mixie
;
;      MIXIE() computes correction factors for line blends given a set of
;      intersting lines (LSTR) and a set of contaminants (SSTR). 
;      LSTR and SSTR should both be RD_LIST() style structures, which 
;      contain necessary line information (emissivity, wvl, z, etc.) 
; 
;      MIXIE()  will identify which lines in SSTR are possible contaminants
;      for each line in LSTR.  SSTR lines must meet two criteria to be
;      included in the correction fraction computation of each LSTR line: 
;
;      1) The SSTR line lies within 2*NSIG*LWDT on either side of the LSTR line 
;      2) The total contribution of the line is  greather than a 
;         user-defined threshold: LTHR * INTENSITY OF LINE (LTHR is 
;         keyword input while intensity of the interesting line is
;         computed internally.)
; 
;      Correction factors are estimated by computing complimentary 
;      error functions with an IDL native routine ERRORF() for each 
;      contaminant (i.e. gaussian  profiles are assumed). If an
;      OGIP-compliant rmf is fed via the RMFSTR keyword, the rmf 
;      will be used to estimate correction factors rather than
;      ERRORF() <<RMFSTR NOT IMPLEMENTED YET>>.
;      The correction factor is defined as follows: 
; 
;        CORRECTION FACTOR = INT / (INT + total contributions from blends)
;        where INT = intensity of interesting line
;
;      If set, the MIXSTR keyword can return a structure containing the
;      identification, emissivities, and wavelengths of each contaminant
;      line, together with LSTR and SSTR indices which match contaminated line 
;      to contaminant (VICTIM and CULPRIT respectively). Also to be contained in
;      the structure are the fractional contributions CONTAMFRAC and intensity of
;      of the contaminant CONTAMINT (i.e. the total contribution of the contaminant is
;      CONTAMINT*CONTAMFRAC). The fields in the output structure will be the
;      the following (see KEYWORDS for individual descriptions):
;
;      {CONTAMFRAC,CONTAMINT,VICTIM,CULPRIT,WVL,IDTAG}
;
; 
;syntax 
;      x=mixie(lstr,sstr,lwdt=lwdt,nsig=nsig,lthr=lthr,rmfstr=rmfstr,mixstr=mixstr,$
;           contamfrac=contamfrac,victim=victim,culprit=culprit,contamint=contamint,$
;           idtag=idtag,vint=vint,mproxy=mproxy,mplet=mplet,_extra=e)
;
;parameters
;      lstr   [INPUT;required] Line (and/or conitnuum) structure of
;                      RD_LIST() format containing info (wvl,emis,z) on 
;                      lines (and/or continuum) for which possible contaminants will be identified 
;                      and for which correction factors  will be
;                      calculated. Assumes that ion equilibria are
;                      included. (see RD_LIST() keyword INCIEQ)
;      sstr   [I/O]    Line structure of RD_LIST() format 
;                      containing user specified line contaminant suspect
;                      info (wvl, emis, z).  Assumes that ion equilibria are
;                      included. (see RD_LIST() keyword INCIEQ). If SSTR is
;                      ommitted or undefined on input, the emissivity database will
;                      read with RD_LIST() in the neccessary wavelength range
;                      and including all species to create a list of  suspects.
;                      If undefined on input, will contain RD_LIST() structure
;                      on output.                
;                     
;keywords 
;      lwdt   [INPUT]  Std. deviation to be used as an estimate for 
;                      line broadening. Only contamination candidate profiles (sstr) 
;                      within 2*lwdt*nsig of the line in question will
;                      be considered. May be single value or an array
;                      with values for each element in LSTR. (must
;                      match number of lines+continuum specified in LSTR)
;                      default = 0.2 pixels
;      nsig   [INPUT]  number of standard deviations from line of interest that
;                      defines wavelength range within which to consider contamination.
;                      If not set, the default is determined by LTHR.  
;                      May be single value or an array of values for
;                      each element in lstr. If array, must match the
;                      number of lines+continuum specified in LSTR.
;                      
;                      So the range within which to consider
;                      contamination for each line is defined as: 
;
;                         WRANGE = [LSTR.WVL(x)-LWDT*NSIG,LSTR.WVL(x)+LWDT*NSIG]                     
;
;                      or if LWDT and NSIG are arrays: 
;
;                         WRANGE = [LSTR.WVL(x)-LWDT(x)*NSIG(x),LSTR.WVL(x)+LWDT(x)*NSIG(x)]
; 
;                      NOTE: To handle continuum measurements, 
;                      choose an NSIG which adequatly covers your
;                      continuum measurement wavelength range.        
;         
;            
;      lthr   [INPUT]  fraction of intensity of the line of interest to
;                      use as cut-off intensity value for contaminant
;                      candidate fractional flux  contributions. default = 1e-4  
;                      i.e. cutoff = (intnsity of line of interest)/100
;      rmfstr [INPUT]  <<NOT IMPLEMENTED YET>> structure containing OGIP compliant RMF. see
;                      RD_OGIP_RMF(). If set, contamination will be
;                      estimated by analyzing RMF instead of using ERRORF()
;      mixstr [OUTPUT] Structure containing detailed information for
;                      each pair of blends. The structure fields are
;                      the following: 
;
;                      {CONTAMFRAC,CONTAMINT,VINT,VICTIM,CULPRIT,WVL,IDTAG}
;
;                      Each field can also be returned individually:
;
;      contamfrac [OUTPUT] array containting contribution fractions for
;                      each culpable line
;      contamint [OUTPUT] the contaminant line intensities         
;      vint [OUTPUT] the victim line intensities 
;      victim [OUTPUT] lstr index of intersting line being contaminated 
;      culprit[OUTPUT] sstr index of indices of contaminant lines
;      idtag  [OUTPUT] string containing derscription of each culprit.
;                      the description is as follows:
;
;                      SPECIES, WVL, DATABASE, DESCRIPTION
;                      
;                      SPECIES gives atomic symbol and ionic
;                      state. WVL gives Wavelength. DATABASE gives 
;                      the name of the DATABASE used. And DESCRIPTION 
;                      gives the level designations and electron
;                      configurations of the transition. e.g.
; 
;                      MgXII 8.41920 $CHIANTI (1s) 2S_1/2 (2p) 2P_3/2 
;                      OVIII 18.9726 $CHIANTI (1s) 2S_1/2 (2p) 2P_1/2 
;
;      vidtag [OUPUT]  string containing description of each victim
;                      and of same format as idtag 
;      mproxy [INPUT] If a correction factor is needed for  a
;                       multiplet measured using a single profile,
;                       then MPROXY will specify the
;                       multiplet components on input, and contain
;                       the multiplet correction factor on output.
; 
;                       One component of the multiplet is to stand
;                       proxy for the group, and must be specified in 
;                       LSTR. 
;                       
;                         MPROXY is defined  as follows: 
;                     
;                         mproxy = strarr(2,N) or 
;                         mproxy = intarr(2,N)
;                      
;                         N is the number of component-proxy
;                         relationships. Each component-proxy
;                         relationship  is expressed as a two
;                         element string array. The first element should
;                         contain the LSTR index of the proxy. The
;                         second element should contain either the LSTR index of the
;                         component if it is included in LSTR, or a
;                         string description of the component in IDTAG format. 
;                         One may optionally leave out:
;
;                            o the DESCRIPTION field (e.g. just use SPECIES, WVL,DATABASE). 
;                            o the DATABASE and WVL fields (e.g. just use
;                              SPECIES, DATABASE)
; 
;                       On output, one correction factor will be given 
;                       for each multiplet. The index of this
;                       correction factor in the output correction
;                       factor array is determined by the position of
;                       the proxy in LSTR. (See examples below)
;                                                 
;      mplet  [INPUT] If blending between one or
;                       more multiplet components is already accounted
;                       for (e.g. a multiplet is modelled with
;                       serveral profiles)  then MPLET can
;                       specify which multiplet component to leave 
;                       out of the correction factor calculation of 
;                       another component. 
; 
;                         MPLET is defined as follows: 
; 
;                         mplet = strarr(2,N)
;                         mplet = intarr(2,N)
;          
;                         N is the number of component-component 
;                         relationships. Each component-component 
;                         relationship is expressed as a two element 
;                         string array. The first element should
;                         contain  the LSTR index of the component 
;                         whose correction factor is needed. 
;                         The second element should contain either the
;                         LSTR index of the component to leave out, or
;                         a string description of the component to
;                         leave out in IDTAG format. 
;                         One may optionally leave out:
;
;                            o the DESCRIPTION field (e.g. just use SPECIES, WVL,DATABASE). 
;                            o the DATABASE or WVL fields (e.g. just use
;                              SPECIES, DATABASE)
;
;                      Correction factors will not be calculated for 
;                      lines specified to be left out of another
;                      lines' correction factor calculation. See
;                      examples below.
;      brnstr [I/O]    if mixie is to be run multiple times
;                      (e.g. using MCMC_DEM() ) the brnstr keyword can 
;                      be set so that MIXIE() creates a sturcture
;                      containing internally calculated and recycle-able information, to be fed
;                      to subsequent MIXIE() calls. The structure is
;                      of the form: 
;
;                      {FRACS,NDX,INTNDX,NSS,NNSS,PRXBSE,MULTSTR} 
;
;                      where: 
;
;                      FRACS   = float array contribution fractions of all non-exempt lines
;                      NDX     = long array  sstr index of each of thes lines 
;                      INTNDX  = float array ndx and fracs index intervals corresponding to each victim
;                      NSS     = float array lstr index of each victim  
;                      NNSS    = long array  number of actual victims i.e. n_elements(nss) 
;                      PRXBSE  = long array  lstr index of the proxy in each proxy relation. 
;                      MULTSTR = STRUCT      rd_line structure of proxee in each relation 
;
;                      When fed as input, BRNSTR allows MIXIE to avoid
;                      explicit contribution factor calculations for
;                      each suspect and cumbersome MPROXY and
;                      MPLET sting parsing. 
;      DEM    [INPUT]  differential emission measure  at each T [cm^-5/logK]
;                      * default is a constant=1e14!!
;                      * if defined, DLOGT must also be specified
;                      * if LSTR and SSTR emissivity temperature grids do not
;                        match  DLOGT then they are rebinned and
;                        resampled so that they do 
;      DLOGT  [INPUT   logarithmic temperature grid at which DEM is defined
;                      * default is 4.0 to 8.0 in steps of 0.05
;      ABUND  [I/O]  abudnanaces relative to H abund(0) = 1
;      	               * abund(Z-1) contains the abundance for element Z
;                      * default: Grevesse & Sauval(1998)
;      ABNDUPDT [INPUT] If set, then abundnaces will be updated on
;                      output. Abundances are calculated after taking 
;                      into account the correction factors 
;      OBSFLX [INPUT]  Obseved fluxes used to calculate abundnaces,
;                      must match the number of correction factors
;                      calculated                         
;      OFSIG  [INPUT]  errors for OBSFLX
;      MULTIPLET [INPUT] temporary, same as MPLET
;      MULTPLET [INPUT]  temporary, same as MPLET 
;      MULTPROXY [INPUT] temporary, same as MPROXY
;      MULTIPROXY[INPUT] temporary, same as MPROXY
;      verbose[INPUT]  set verbosity 
;      _extra [INPUT]  pass defined keywords to
;	               LINEFLX (DEM, ABUND, EFFAR, WVLAR, TEMP, NOPH, IKEV, REGRID)
;                      RD_LIST (DBDIR,LOGP,PRES,N_E,EQFILE)
;      
;subroutines
;      LINEFLX
;      RD_LIST
;      INICON
;      CAT_LN
;      REBINX
;      MCMC_ABUND
;      MCMC_EROR 
;      IS_KEYWORD_SET
;
;restrictoins
;      Is is suggested that LSTR and SSTR be read from a common
;      database to avoid accidentally identifying a line as 'self blending' 
;      MIXIE() will NOT compute correct correction factors for DR lines. 
;
;usage examples
;
;      * !LDBDIR = '$CHIANTI'
;        Ne9   = rd_list('Ne9|[13.44,13.45]',/incieq,sep='|')
;        cand  = rd_list('ALL|[12.44,15.45]',/incieq,sep='|')
;        factor = mixie(Ne9,Cand,lwdt=0.008,nsig=5,lthr=0.02,contamfrac=contamfrac,victim=victim,$
;                      culprit=culprit,contamint=contamint,idtag=idtag, mixst=mixstr)
;        for j = 0,n_elements(culprit)-1 do print,mixstr.idtag(j), mixstr.contamint(j)* mixstr.contamfrac(j)
;
;      * !LDBDIR = '$CHIANTI'
;        Ne9    = rd_list('Ne9|[13.44,13.45]',/incieq,sep='|')
;        factor = mixie(Ne9,lwdt=0.008,nsig=5,lthr=0.02,contamfrac=contamfrac,victim=victim,$
;                 culprit=culprit,contamint=contamint,idtag=idtag, mixst=mixstr)
;        for j = 0,n_elements(culprit)-1 do print,mixstr.idtag(j), mixstr.contamint(j)* mixstr.contamfrac(j)
;
;      * !LDBDIR = '$CHIANTI'
;        T_components  = [6.6,6.9,7.4]           ; log(T[K]) components in the EM
;        EM_components = [3.05,3.05,3.55]*1d11   ; Emission Measure [cm^-3]
;        !DEM=mk_dem('delta', logT = !LOGT, pardem=T_components, indem=EM_components)
;        Ne9   = rd_list('Ne9|[13.44,13.45]',/incieq,sep='|')
;        factor = mixie(Ne9,lwdt=0.008,nsig=5,lthr=0.02,contamfrac=contamfrac,victim=victim,$
;                      culprit=culprit,contamint=contamint,idtag=idtag, mixst=mixstr,DEM=!DEM, DLOGT=!LOGT)
;        for j = 0,n_elements(culprit)-1 do print,mixstr.idtag(j),mixstr.contamint(j)* mixstr.contamfrac(j)
; 
;      * !LDBDIR = '$CHIANTI' 
;        N7  = rd_list('N7|[24.77,24.79]',/incieq,sep='|'); this list specifies a triplet
;        mproxy = [ [2,1], [2,0] ]  ; We use one triplet component as proxy       
;        factor = mixie(N7,lwdt=0.008,nsig=5,lthr=0.008,contamfrac=contamfrac,victim=victim,$
;                 culprit=culprit,contamint=contamint,idtag=idtag,mproxy=mproxy,$
;                 mixstr=mixstr)
;
;      * !LDBDIR = '$CHIANTI' 
;        N7 = rd_list('N7|[24.7846,24.7847]',/incieq,sep='|'); this list contains one N7 triplet component 
;        ; use MPROXY to specify N7 triplet components the one component are to stand proxy for
;        mproxy = [['0','NVII 24.7793 $CHIANTI (1s) 2S_1/2 (2p) 2P_3/2'], $ 
;                  ['0','NVII 24.7844 $CHIANTI (1s) 2S_1/2 (2s) 2S_1/2']]  
;        factor = mixie(N7,njk,lwdt=0.008,nsig=5,lthr=0.008,contamfrac=contamfrac,victim=victim,$
;                       culprit=culprit,contamint=contamint,idtag=idtag,$
;                       mixst=mixstr,mproxy=mproxy )
;
;      * !LDBDIR = '$CHIANTI' 
;        N7 = rd_list('NVII|24.7847|$CHIANTI|(1s) 2S_1/2 (2p) 2P_1/2',sep='|',/incieq)
;        ; use MPROXY to specify multiplet components the N7 at 24.7847 is to stand proxy for
;        mproxy = [['0','NVII 24.7793 $CHIANTI (1s) 2S_1/2 (2p) 2P_3/2'], $ 
;                  ['0','NVII 24.7844 $CHIANTI (1s) 2S_1/2 (2s) 2S_1/2']]  
;        factor = mixie(N7,j,lwdt=0.008,nsig=5,lthr=0.008,contamfrac=contamfrac,victim=victim,$
;                       culprit=culprit,contamint=contamint,idtag=idtag,$
;                       mixst=mixstr,mproxy=mproxy)
;                 
;history
;                INTENDED as an improvement on bland() (LL:Jul'03)
;                NEW documentation new default nsig 
;                  output now correction factor and not mixstr (LL:Aug'03)
;                BUG FIX: 1/SQRT(2) factor needed for arguement of ERRORF so
;                  lwdt is actually standard deviation
;                now prints message when no blends are found
;                now robust to theoretical lines (designated by
;                negative wavelengths)                      
;                BUG FIX correction factor calculations for multiple
;                victim list was incorrect (LL/DG:Sep'03) 
;                ADDED keyword COMPOSITE to handle doublets (LL:Sep'03) 
;                REPLACED keyword COMPOSITE with MUTLPROXY which now handles
;                  multiplets. ADDED related keyword MULTIPLET. (LL:Oct'03)
;                ADDED correction for 'over correction' i.e. a CULPRIT now
;                  cannot 'overdistribute' its flux when contaminating
;                  several VICTIMS at once. (LL:Oct'03)  
;                CHANGED LTHR and NSIG so that they can handle arrays
;                  so correction factors for continuum measurements
;                  are possible (LL:Oct'03) 
;                CHANGED default RD_LIST() call so that makes
;                  intelligent efficient wrange decision (LL:Oct'03) 
;                ADDED BRNSTR keyword so that doesn't bog down
;                  MCMC_DEM(). ALSO included explicit intensity
;                  calculations so to avoid being slowed by LINEFLX(). (LL:Oct'03) 
;                ADD ABNDUPDT keyword to toggle ABUND updating on
;                  output (LL:Jan'04) 
;                CHANGED location of ABUNDANCE calculation so that it
;                  abundances get updated after contamination
;                  correction. ADDED keyword OFSIG (LL:Jan'04) 
;                CHANGED MULTIPLET KEYWORD to MPLET and MULTPROXY to
;                  MPROXY (LL/DG:Jan'04) 
;		 changed default of LTHR to 1e-4 (VK;Jan'04)
;                ADDED new MPLET and MPROXY format (LL:Feb'04) 
;                BUG FIX over correction mechanism (i.e. considering
;                  situations of one culprit and two victims) was not
;                  working properly, leading to correction factors >
;                  1. (LL/DG:Feb'04) 
;                BUG FIX total(contamfrac) gt 1 changed to
;                  n_elements(culprit) gt 1  (LL:Feb'04) 
;                BUG FIX keywords idtag and mixstr were failing 
;                  because need to handle with keyword_set()  not
;                  arg_present() (LL:Mar'04) 
;                VARIOUS BUG FIXES 
;                  sstr enteed as dummy argument = crash 
;                  mixstr and idtag fix (LL:Mar'05)  
;		 replaced direct reference to !LOGT with calls to
;		   defsysv and setsysval (VK; Apr'05)
;                BUG FIX if sstr is one component only it was 
;                  ignored. (LL:Apr'05) 
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;-               

ok='ok' & np=n_params() & nl=n_elements(lstr)
ml=n_tags(lstr) & mss=size(sstr)
if np lt 1 then ok='Insufficient parameters' else $
 if nl eq 0 then ok='LSTR is not defined' else $
  if ml eq 0 then ok='LSTR must be a structure' else begin
      lname=tag_names(lstr)
      o1=(where(strpos(lname,'LINE_INT',0) ge 0))[0]
      o2=(where(strpos(lname,'LOGT',0) ge 0))[0]
      o3=(where(strpos(lname,'WVL',0) ge 0))[0]
      o4=(where(strpos(lname,'Z',0) ge 0))[0]
      o5=(where(strpos(lname,'ION',0) ge 0))[0]
      o6=(where(strpos(lname,'DESIG',0) ge 0))[0]
      o7=(where(strpos(lname,'CONFIG',0) ge 0))[0]
      o8=(where(strpos(lname,'SRC',0) ge 0))[0]
      o9=(where(strpos(lname,'JON',0) ge 0))[0]
    if o1 lt 0 or o2 lt 0 or o3 lt 0 or o4 lt 0 or $
      o5 lt 0 or o6 lt 0 or o7 lt 0 or o8 lt 0 or o9 lt 0 then $
      ok='LSTR not in standard format'
endelse 
;    if (mss(mss(0)+1) ge 1) and (mss(mss(0)+1) ne 8) then ok='SSTR must be a structure' 
       if mss(mss(0)+1) eq 8 then begin 
       sname=tag_names(sstr)
       o1=(where(strpos(sname,'LINE_INT',0) ge 0))[0]
       o2=(where(strpos(sname,'LOGT',0) ge 0))[0]
       o3=(where(strpos(sname,'WVL',0) ge 0))[0]
       o4=(where(strpos(sname,'Z',0) ge 0))[0]
       o5=(where(strpos(sname,'ION',0) ge 0))[0]
       o6=(where(strpos(sname,'DESIG',0) ge 0))[0]
       o7=(where(strpos(sname,'CONFIG',0) ge 0))[0]
       o8=(where(strpos(sname,'SRC',0) ge 0))[0]
       o9=(where(strpos(sname,'JON',0) ge 0))[0]
    if o1 lt 0 or o2 lt 0 or o3 lt 0 or o4 lt 0 or $
       o5 lt 0 or o6 lt 0 or o7 lt 0 or o8 lt 0 or o9 lt 0 then $
       ok='SSTR not in standard format'
    endif

if keyword_set(rmfstr) then begin 
    if n_tags(rmfstr) eq 0 then ok = 'RMFSTR must be a structure'
       rname=tag_names(rmfstr)
       o1=(where(strpos(rname,'NNRG',0) ge 0))[0]
       o2=(where(strpos(rname,'ELO',0) ge 0))[0]
       o3=(where(strpos(rname,'EHI',0) ge 0))[0]
       o4=(where(strpos(rname,'NCHAN',0) ge 0))[0]
       o5=(where(strpos(rname,'EMN',0) ge 0))[0]
       o6=(where(strpos(rname,'EMX',0) ge 0))[0]
       o7=(where(strpos(rname,'N_GRP',0) ge 0))[0]
       o8=(where(strpos(rname,'F_CHAN',0) ge 0))[0]
       o9=(where(strpos(rname,'N_CHAN',0) ge 0))[0]
       o10=(where(strpos(rname,'MATRIX',0) ge 0))[0]
       o11=(where(strpos(rname,'FIRSTCHAN',0) ge 0))[0]
   if o1 lt 0 or o2 lt 0 or o3 lt 0 or o4 lt 0 or o5 lt 0 or o6 lt 0 $
       or o7 lt 0 or o8 lt 0 or o9 lt 0 or o10 lt 0 or o11 lt 0 then $
       ok='RMFSTR not in standard format'
    endif

if ok ne 'ok' then begin 
print, 'Usage: mixie,lstr,sstr,lwdt=lwdt,nsig=nsig,lthr=lthr,rmfstr=rmfstr,mixstr=mixstr,$'
print, 'contamfrac=contamfrac,victim=victim,culprit=culprit,contamint=contamint,idtag=idtag,$'
print, 'vidtag=vidtag,vint=vint,mproxy=mproxy,mplet = mplet,verbose=verbose,$'
print, 'nhne=nhne,abund=abund,dem=dem,dlogt=dlogt, obsflx=obsflx,ofsig=ofsig, brnstr=brnstr,$'
print, 'abndupdt=abndupdt,_extra=e'
   if np ne 0 then message, ok, /info 
   return,-1
endif

vv=0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1
if not keyword_set(lwdt) then lwdt =0.2
if not keyword_set(lthr) then lthr = 1e-4 else lthr = (lthr[0]>1e-10)<1
if not keyword_set(nsig) then nsig = interpol(findgen(15),errorf(findgen(15)/sqrt(2)),1-lthr)
if not keyword_set(abund) then abund = getabund('grevesse & sauval') 
if not keyword_set(nhne) then nhne= 0.83 & nhne = nhne(0) 
;initiate element symbols and roman numerals
inicon, atom=atom, roman=roman
sorc = [' ?', '$SPEX','$CHIANTI','$APED'] 

;parse LSTR 
   i_wvl  = lstr.wvl  & i_emis = lstr.line_int
   i_logt = lstr.logt & i_z    = lstr.z  
   i_nw   = n_elements(i_wvl) & i_ion = lstr.ion & i_jon = lstr.jon
   i_desig = lstr.desig & i_config = lstr.config 
   i_nwo = i_nw & i_src = lstr.src
   i_descrip =  '('+lstr.config[0,*]+') '+lstr.desig[0,*]+' ('+lstr.config[1,*]+') '+lstr.desig[1,*]
if keyword_set(dem) and not keyword_set(dlogt) then message,'If DEM is set then DLOGT must be set' 
if keyword_set(dem) and keyword_set(dlogt) then begin 
   ndt = n_elements(dlogt) & nd  = n_elements(dem) 
   cdt = n_elements(i_logt)
   if ndt ne nd then message,'DEM and DLOGT arrays do not match' 
   if (ndt ne cdt) or total(dlogt-i_logt) ne 0 then begin 
     ni_emis = fltarr(ndt,i_nw)  
     for i=0,i_nw-1 do ni_emis(*,i) = rebinx(i_emis(*,i),i_logt,dlogt)
     i_emis = ni_emis & i_logt = dlogt
     lstr=create_struct('LINE_INT',ni_emis,'LOGT',dlogt,'WVL',i_wvl,'z',i_z,'ION',i_ion, $
     'DESIG', i_desig,'CONFIG',i_config,'SRC',i_src, 'JON',i_jon)
   endif
endif else begin 
   defsysv,'!LOGT',exists=ivar
   if ivar eq 0 then dlogt=findgen(81)*0.05+4. else $
	setsysval,'LOGT',dlogt,/getval
   dem = dlogt*0+1d14
   ndt = n_elements(dlogt) & nd  = n_elements(dem) 
   cdt = n_elements(i_logt)
   if (ndt ne cdt) or total(dlogt-i_logt) ne 0 then begin 
     dem = rebinx(dem, dlogt, i_logt)
     dlogt = i_logt
   endif
endelse

       ndt = n_elements(dlogt) 
       h=6.626176e-27 & c= 2.9979e10
       nt = n_elements(i_logt) 
       enrg = (h*c*1e8/(abs(i_wvl)))
       mdt=alog(10)*total(i_logt(1:*)-i_logt)/(nt-1)
       tmp = fltarr(nt)+1 & tmp(0)= 0.5 & tmp(nt-1)= 0.5
       if i_nw eq 1 then i_int = $ 
       mdt*nhne*(abund(i_z-1)*((1d-23)*total(i_emis*(tmp*dem))))/enrg else $ 
       i_int = mdt*nhne*(abund(i_z-1)*((1d-23)*i_emis##(tmp*dem)))/enrg
  ;    i_int2  = lineflx(i_emis, i_logt, i_wvl, i_z,dem=dem,abund=abund, _extra = e)   

;check if lwdt and nsig are arrays and if they match
nlwdt=n_elements(lwdt) & nnsig=n_elements(nsig) 
if nlwdt eq 1 then lwdt = fltarr(i_nw)+lwdt
if nnsig eq 1 then nsig = fltarr(i_nw)+nsig
nlwdt=n_elements(lwdt) & nnsig=n_elements(nsig) 
if nlwdt ne i_nw then ok = 'LWDT array does not match LSTR' else $ 
if nnsig ne i_nw then ok = 'NSIG array does not match LSTR' else $ 
if nnsig ne nlwdt then ok = 'LWDT size does not match NSIG size' 
if ok ne 'ok' then begin 
   message, ok, /info
   return, -1 
endif 
if vv gt 7 then message, 'deblending....',/info 

;if SSTR not present or dummy variable, generate it according to necessary wavelength range
    if not (arg_present((sstr))) or (mss(2) ne 8) then begin 
       if vv gt 5 then message, 'SSTR not present or undefined in call.',/info
       if vv gt 5 then message, 'Will now read data base to create list of possible contaminants.',/info
    ;find wvl intervals with which to call rd_list and deal with overlaps
       ints = fltarr(2,i_nw)
       ordr = sort(i_wvl) & tmp0 = i_wvl-2*nsig*lwdt & tmp1 =  i_wvl+2*nsig*lwdt 
       ints(0,*) = tmp0(ordr)  &  ints(1,*) = tmp1(ordr) ;downs & ps  
       qup = intarr(i_nw)
       for qq = 0, i_nw-1 do qup(qq) = total(where(ints(0,qq:*) lt ints(1,qq)))  
    ; 0s in qup indicate borders between islands 
       qupb = [-1,where(qup eq 0)] & nqupb = n_elements(qupb)-1 
       nuints = fltarr(2,nqupb)  
       for ww = 0, nqupb-1 do nuints(*,ww)=[ints(0,qupb(ww)+1),ints(1,qupb(ww+1))]
       for qq = 0, nqupb-1 do begin 
           if qq eq 0 then begin 
              lnlist = 'All|'+'['+string(nuints(0,0))+','+string(nuints(1,0))+']'
              sstr =  rd_list(lnlist,/incieq,sep='|', _extra=e) 
           endif else begin 
              lnlist = 'All|'+'['+string(nuints(0,qq))+','+string(nuints(1,qq))+']'
              sstrx   = rd_list(lnlist,/incieq,sep='|', _extra=e) 
              sstr = cat_ln(sstr,sstrx)
           endelse
       endfor 
   endif
;parse SSTR 
   c_wvl  = sstr.wvl          &   c_emis = sstr.line_int
   c_logt = sstr.logt         &   c_z    = sstr.z
   c_nw   = n_elements(c_wvl) &   c_src  = sstr.src
   c_ion  = sstr.ion          &   c_jon  = sstr.jon
   c_desig = sstr.desig  & c_config = sstr.config 
if keyword_set(dem) and not keyword_set(dlogt) then message,'If DEM is set then DLOGT must be set' 
if keyword_set(dem) and keyword_set(dlogt) then begin 
   ndt = n_elements(dlogt) & nd  = n_elements(dem) 
   cdt = n_elements(c_logt)
   if ndt ne nd then message,'DEM and DLOGT arrays do not match' 
   if (ndt ne cdt) or total(dlogt-c_logt) ne 0 then begin 
     nc_emis = fltarr(ndt,c_nw)  
     for i=0,c_nw-1 do nc_emis(*,i) = rebinx(c_emis(*,i),c_logt,dlogt)
     c_emis = nc_emis & c_logt = dlogt
     sstr=create_struct('LINE_INT',nc_emis,'LOGT',dlogt,'WVL',c_wvl,'z',c_z,'ION',c_ion, $
     'DESIG', c_desig,'CONFIG',c_config,'SRC',c_src, 'JON',c_jon)
   endif
endif else begin 
   defsysv,'!LOGT',exists=ivar
   if ivar eq 0 then dlogt=findgen(81)*0.05+4. else $
	setsysval,'LOGT',dlogt,/getval
   dem = dlogt*0+1d14
endelse
ndt = n_elements(dlogt) 

       h=6.626176e-27 & c= 2.9979e10
       nt = n_elements(c_logt) 
       enrg = (h*c*1e8/(abs(c_wvl)))
       mdt=alog(10)*total(c_logt(1:*)-c_logt)/(nt-1)
       tmp = fltarr(nt)+1 & tmp(0)= 0.5 & tmp(nt-1)= 0.5
       if c_nw eq 1 then $ 
     c_int =mdt*nhne*(abund(c_z-1)*((1d-23)*total(c_emis*(tmp*dem))))/enrg else $
     c_int =mdt*nhne*(abund(c_z-1)*((1d-23)*c_emis##(tmp*dem)))/enrg

exempt_list = -1 ; initiate exempt list.. these guys are above the law for select cases 
nmlt = 0 
nprx = 0 
bse  = -1

if keyword_set(multproxy) then begin 
  message, 'MULTPROXY is now obsolete, use MPROXY instead.',/info 
  return,-1L
endif 
if keyword_set(multiproxy) then begin 
  message, 'MULTIPROXY is obsolete, use MPROXY instead.',/info 
  return,-1L
endif 
if keyword_set(multiplet) then begin 
  message, 'MULTIPLET is now obsolete, use MPLET instead.',/info 
  return,-1L
endif 
if keyword_set(multplet) then begin 
  message, 'MULTPLET is obsolete, use MPLET instead.',/info 
  return,-1L
endif 

if keyword_set(mproxy) then begin
   t1 = size(mproxy) 
   if t1(0) eq 1 then mproxy = reform(mproxy,2,1)
   t1 = size(mproxy)
   if t1(3) ne 7 then mproxy = strcompress(string(long(mproxy)),/remove_all)  
   if t1(0) ne 2 then begin 
      message, 'MPROXY must be two-dimensional string array', /info 
      return, -1 
   endif 
   mult = mproxy & nprx = t1(2) 
endif
if keyword_set(mplet) then begin 
   t1 = size(mplet)  
   if t1(0) eq 1 then mplet = reform(mplet,2,1)
   t1 = size(mplet) 
   if t1(3) ne 7 then multplet = strcompress(string(long(mplet)),/remove_all)  
   if t1(0) ne 2 then begin 
       message, 'MPROXY must be two-dimensional string array', /info 
       return, -1 
   endif 
   mult = mplet & nmlt = t1(2) 
endif
if keyword_set(mproxy) and keyword_set(mplet) then begin 
   mult = transpose([transpose(mproxy),transpose(mplet)])
endif

nss = findgen(i_nw) & nnss=n_elements(nss)
witness_list = [-1] ; initiate list of lstr members who are not direct victims 

if not keyword_set(brnstr) then begin
multstr = cat_ln(sstr, pick = [0]) 
c_descrip =  '('+sstr.config[0,*]+') '+sstr.desig[0,*]+' ('+sstr.config[1,*]+') '+sstr.desig[1,*]
if keyword_set(mult)  then begin
      nR = nprx+nmlt  & numlst = ['1','2','3','4','5','6','7','8','9','0']
      tst2 = fltarr(nR) & tst3 = tst2 & mee_wvl = fltarr(nR) & mee_flx = fltarr(nR)  
      mee_z = fltarr(nR) & mee_ion = fltarr(nR) & mee_descrip = strarr(nR)
      for r = 0,nR-1 do begin ; loop through mult relationships
         ;test 1st elements 
         bse =  strtrim(mult(0,r),2) 
         lenp = strlen(bse) & tmpA = intarr(lenp) & tmpB = intarr(10)  
         for ww = 0,lenp-1 do begin 
            for qq = 0,9 do begin 
               tmpB(qq) = strmatch(strmid(bse,ww,1),numlst(qq))  
            endfor  
            tmpA(ww) = total(tmpB)
        endfor
        if total(tmpA) ne lenp then begin  
           message, 'MPROXY or MPLET format incorrect.', /info 
           return, -1 
         endif 
         bse  = long(bse)
         if bse gt i_nw-1 or bse lt 0 then begin 
           message, "MPROXY or MPLET component's LSTR index out of range"
           return, -1 
         endif 
         ;test 2nd elements 
          cmp = strtrim(mult(1,r),2) ;remove leading and trailing blanks
          tagtst = strsplit(cmp, /extract) 
          if n_elements(tagtst) eq 1 then begin ; is spcified string an index?
            lenc = strlen(cmp) & tmpA = intarr(lenc) & tmpB = intarr(10) 
            for ww = 0, lenc-1 do begin 
               for qq = 0,9 do begin 
                  tmpB(qq) = strmatch(strmid(cmp,ww,1),numlst(qq)) 
               endfor
            tmpA(ww) = total(tmpB) 
            endfor
            if total(tmpA) ne lenc then begin 
              message, 'MPROXY or MPLET format is incorrect.' , /info 
              return, -1 
            endif
            lngc = long(cmp)
            if lngc gt i_nw-1 or lngc lt 0 then begin 
              message, "MPROXY or MPLET component's LSTR index out of range.",/info
              return, -1 
            endif 
            witness_list=[witness_list,lngc]
            mee_flx(r) = lineflx(i_emis(*,lngc), i_logt, i_wvl(lngc), i_z(lngc),dem=dem,abund=abund, _extra=e) 
            mee_wvl(r) = i_wvl(lngc) & mee_z(r) = i_z(lngc) & mee_ion(r) = i_ion(lngc)  
            mee_descrip(r) = i_descrip(lngc)
            tmpstr = cat_ln(lstr, pick = [r]) 
        endif else begin ; else assume IDTAG format
            ; check idtag fields against LSTR and SSTR to see if we can avoid RD_LIST() 
            avoid = -1
            if n_elements(tagtst) eq 7 then begin $ ; if tagtst adeqaute, examine LSTR 
               probono = where(i_wvl eq tagtst(1)) & scndcheck  = 0 
               if total(probono) ge 0 then begin  
                  descm=where(strcompress(i_descrip(probono),/remove_all) $
                                eq strcompress(strjoin(tagtst(3:*)),/remove_all)) 
                  if total(descm) ge 0 then begin          
                     symb2zion, tagtst(0), xz, zion, _extra = e  
                     speciesm = where((i_z(probono(descm))) eq xz) and (i_ion(probono(descm)) eq zion) 
                     if total(speciesm) ge 0 then begin  
                        fc = long(probono(descm(speciesm(0))))
                        if strmatch(tagtst(2),sorc(i_src(fc))) then begin 
                          avoid = fc ; we can
                          witness_list=[witness_list,fc]
                          tmpstr = cat_ln(lstr, pick = [avoid])
                        endif 
                    endif 
                endif
            endif 
                if avoid eq -1 then begin ; now check SSTR 
                    probono = where(c_wvl eq tagtst(1)) & scndcheck = 1 
                    if total(probono) ge 0 then begin  
                        descm=where(strcompress(c_descrip(probono),/remove_all) $
                                   eq strcompress(strjoin(tagtst(3:*)),/remove_all))          
                        if total(descm) ge 0 then begin 
                            symb2zion, tagtst(0), xz, zion, _extra = e  
                            speciesm = where((c_z(probono(descm)) eq xz) and (c_ion(probono(descm)) eq zion))          
                            if total(speciesm) ge 0 then begin 
                                fc = (probono(descm(speciesm(0))))
                                if strmatch(tagtst(2),sorc(c_src(fc))) then begin 
                                  avoid = fc ; we can 
                                  tmpstr = cat_ln(sstr, pick = [avoid])
                                endif
                            endif 
                        endif
                    endif 
                endif 
            endif 
            if n_elements(tagtst) eq 5 then begin ; first check LSTR
                 symb2zion, tagtst(0), xz,zion, _extra = e & scndcheck = 0
                 speciesm = where((i_z eq xz) and (i_ion eq zion))
                 if total(speciesm) ge 0 then begin 
                    descm = where(strcompress(i_descrip(speciesm),/remove_all) $
                              eq strcompress(strjoin(tagtst(1:*)),/remove_all))
                    if total(descm) ge 0 then begin 
                       avoid = speciesm(descm(0)) & tmpstr = cat_ln(lstr, pick = [avoid])
                    endif
                 endif
                 if avoid eq -1 then begin ; now check SSTR
                    symb2zion, tagtst(0), xz,zion, _extra = e & scndcheck = 1
                    speciesm = where((c_z eq xz) and (c_ion eq zion))
                    if total(speciesm) ge 0 then begin 
                       descm = where(strcompress(c_descrip(speciesm),/remove_all) $
                              eq strcompress(strjoin(tagtst(1:*)),/remove_all))
                    if total(descm) ge 0 then begin 
                       avoid = speciesm(descm(0)) & tmpstr = cat_ln(sstr, pick = [avoid])
                    endif
                 endif                    
                 endif
            endif 
            if avoid ne -1 then begin ; line in SSTR or LSTR 
               if scndcheck then begin 
                  mee_flx(r) = lineflx(c_emis(*,avoid), c_logt, c_wvl(avoid), c_z(avoid), $
                                       dem=dem,abund=abund,_extra=e) 
                  mee_wvl(r) = c_wvl(avoid) & mee_z(r) = c_z(avoid) & mee_ion(r) = c_ion(avoid)
                  mee_descrip(r) = c_descrip(avoid)
               endif else begin 
                  mee_flx(r) = lineflx(i_emis(*,avoid), i_logt, i_wvl(avoid), i_z(avoid), $
                                       dem=dem,abund=abund,_extra=e) 
                  mee_wvl(r) = i_wvl(avoid) & mee_z(r) = i_z(avoid) & mee_ion(r) = i_ion(avoid) 
                  mee_descrip(r) = c_descrip(avoid)
               endelse 
            endif else begin ; RD_LIST() can't be avoided
               cmpstr =  rd_list(cmp, /incieq, sep=' ') 
               if n_elements(cmpstr.wvl) gt 1 then begin ; too many lines?
                  message, 'MPROXY or MPLET format incorrect. Each IDTAG must specify only one line.' 
                  return, -1         
               endif else begin  ; if too many lines  
                  if total(cmpstr.line_int) eq -1 then begin 
                      message, 'MPROXY or MPLET specification incorrect', /info 
                      return, -1 
                  endif else begin 
                      tmpstr = cmpstr 
                      tmp_emis = rebinx(cmpstr.line_int[*,0], cmpstr.logt, dlogt)
                      mee_flx(r) = lineflx(tmp_emis, dlogt,cmpstr.wvl,cmpstr.z,dem=dem,abund=abund, _extra=e)
                      mee_wvl(r) = cmpstr.wvl & mee_z(r) = cmpstr.z & mee_ion(r) = cmpstr.ion 
                      mee_descrip(r) = '('+cmpstr.config[0,*]+') '+cmpstr.desig[0,*]+ $ 
                                   ' ('+cmpstr.config[1,*]+') '+cmpstr.desig[1,*]
                  endelse
               endelse; else compute flx and populate variables
            endelse;  else use RD_LIST()
      endelse; assume IDTAG format 
  multstr = cat_ln(multstr, tmpstr)
   endfor ;loop through mult relationships       
  nss=where(histogram(witness_list,min=0,max=i_nw-1) eq 0, nnss) ;indeces of actual victims (exclude witnesses)
  bse = long(mult(0,*)) 
  multstr = cat_ln(multstr, pick = 0)
  if nprx gt 0 then  multstr = cat_ln(multstr, pick = [findgen(nprx)])
endif ; if mult is set 
  brnin = 1 
endif else begin ; if brnstr not set 
  brnin  = -1 
  fracs   = brnstr.fracs 
  ndx     = brnstr.ndx 
  intndx  = brnstr.intndx
  nss     = brnstr.nss 
  nnss    = brnstr.nnss
  prxbse  = brnstr.prxbse
  multstr = brnstr.multstr
  enrg    = h*c*1e8/(abs(multstr.wvl))
  if total(prxbse) ne -1 then begin 
    if nprx  eq 1 then multflx = $
    mdt*nhne*(abund(multstr.z-1))*((1d-23)*total(multstr.line_int*(tmp*dem)))/enrg
    if nprx gt 1 then multflx = $
    mdt*nhne*(abund(multstr.z-1))*((1d-23)*multstr.line_int##(tmp*dem))/enrg
  endif
endelse 


contamfrac = [-1] & victim = [-1] & culprit= [-1] & contamint = [-1] & jon = [-1] 
z    = [-1] & wvl    = [-1] & src    = ' ' & ion = -1    & tmax = [-1] 
 desig = ' '& config = ' ' & descrip=' '   & allfrac = [-1] & allsub = [-1] & nsub= [0]
crctn_factor = fltarr(nnss)

;begin trial '...nothing but the truth so help me God'.

;Use RMF if one is given 
if keyword_set(rmfstr) then begin 
  elo=rmfstr.elo & ehi=rmfstr.ehi; photon nrgis
  emn=rmfstr.emn & emx=rmfstr.emx; channel bin nrgis
  matrix=rmfstr.matrix 
  if rmfstr.firstchan eq 1 then f_chan = rmfstr.f_chan -1 > 0 else f_chan = rmfstr.firstchan
  n_grp=rmfstr.n_grp & n_chan = rmfstr.n_chan & nnrg=rmfstr.nnrg    ;
  for tt = 0, i_nw-1 do begin ;loop through victims
     i_w = abs(i_wvl(tt))
     sub = where((abs(c_wvl) gt i_w-2*nsig(tt)*lwdt(tt)>0) and (abs(c_wvl) lt i_w+2*nsig(tt)*lwdt(tt)))
    if total(sub) ge 1 then begin ; if there are suspects within nsig*lwdt then start 
     uu  = (where(hastogram(!fundae.kevang/i_w,[elo[0],ehi[1:*]]) eq 1))[0] ;phtn nrg grid index for victim    
     profv = fltarr(nnrg)                                                 
        if n_grp(uu) gt 1 then begin ; if more than one grp                  
           A = min(f_chan[*,uu])       ;
           B = max(f_chan[*,uu],mndx)  ;create jury(min max indices of profv support)
           B = B+n_chan[mndx,uu]-1     ;
        endif else begin                                                     
           A = f_chan[uu]              ;                                    
           B = A+n_chan[uu]-1          ; 
        endelse                                                              
     for  nn = 0, n_elements(sub)-1 do begin ;loop through suspect line subset
        cndx = sub(nn) 
        c_w  = abs(c_wvl(cndx))
        cc   = (where(hastogram(!fundae.kevang/c_w,[elo[0],ehi[1:*]]) eq 1))[0] ;phtn nrg grid index for suspect
        profc = fltarr(nnrg)                                                                    ;
        if n_grp(cc) gt 1 then begin ; if more than one grp                                     ; 
           profc[f_chan(0,cc):f_chan(0,cc)+n_chan(0,cc)-1] = matrix[0:n_chan(0,cc)-1,cc]        ;create suspect 
           for pp = 1, n_grp(cc)-1 do profc[f_chan(pp,cc):f_chan(pp,cc)+n_chan(pp,cc)-1] = $    ;profile (profc)
                                      matrix[n_chan(pp-1,cc):n_chan(pp-1,cc)+n_chan(pp,cc)-1,cc]; 
           C = min(f_chan[*,cc])        ;
           D = max(f_chan[*,cc],mndx)   ;create jury(min max indices of profc support)
           D = D+n_chan[mndx,cc]-1      ;
        endif else begin                                                                        ;
           C = f_chan[cc]               ;
           D = C+n_chan[cc]-1           ;
           profc[f_chan(cc):f_chan(cc)+n_chan(cc)-1] = matrix[0:n_chan(cc)-1,cc]                ;
        endelse                                                                                 ;
         ;Cutoff according to lthrold
        if (A>C) lt B and (B<D) gt A then cntrbtn_frc = total(profc[(A>C):(B<D)]) else cntrbtn_frc=0 ; verdict 
        if (cntrbtn_frc*c_int(cndx) gt lthr*i_int(tt)) and (c_w ne i_w) then begin ;if suspect is guilty 
           contamfrac   = [contamfrac    , cntrbtn_frc] ; contain contrubution fractions     
           victim = [victim  , tt         ] ; lstr indices of interesting line being contaminated
           culprit= [culprit , sub(nn)    ] ; sstr indices of contaminant lines 
           contamint = [contamint  , c_int(cndx)] ; the contaminant line intensities
           jnk    = max(c_emis(*,cndx),ji)  ; find TMAX for contaminant
           tmax   = [tmax, c_logt(ji)]      ; temperature of max formation 
           z      = [z       , c_z(cndx)  ] ; relevant atomic number
           wvl    = [wvl     , c_wvl(cndx)] ; relevant wavelength
           src    = [src     , sorc(c_src(cndx))] ; source (which database?)
           ion    = [ion     , c_ion(cndx)] ; relevant ionic state 
           jon    = [jon     , c_jon(cndx)] ; relevant significant (ionic-equilibria wise) ionic state  
        endif ;end guilty verdict if 
     endfor ;end loop through suspects
    endif; if there are suspects within nsig*lwdt 
   if total(z) eq -1 then begin 
       crctn_factor = 1D 
       if vv gt 5 then print, 'No blends found for '+ atom[i_z(tt)-1]+roman[i_ion(tt)-1]$
              +' at '+strcompress(string(i_wvl(tt)),/remove_all)
    ;calculate correction factor for each victim 
    endif else begin 
       crctn_factor(tt) = i_int(tt)/(i_int(tt)+total(contamfrac[1:*]*contamint[1:*]))
    endelse
  endfor ;end loop through victims
  
 endif else begin ;end rmfstr set if  

 ;Assume Gaussians otherwise 
 for jjp = 0, nnss-1 do begin ;loop through intersting lines 
    jj      = nss(jjp)
    i_w     = abs(i_wvl(jj)) 
    subfrac = [-1]
    nusub   = [-1]
    ; if burn is necessary 
    if brnin eq 1 then begin 
         ee = where( bse eq (jj) ) ; are there any other exemptions from mproxy or mplet? 
         exempt_w       = i_wvl(jj)    
         exempt_z       = i_z(jj) 
         exempt_ion     = i_ion(jj)  
         exempt_descrip = i_descrip(jj)
         ee = where( bse eq (jj) ) ; are there any other exemptions from mproxy or mplet? 
         if total(ee) ge 0 then begin 
             exempt_w       = [exempt_w,mee_wvl(ee)]
             exempt_z       = [exempt_z,mee_z(ee)] 
             exempt_ion     = [exempt_ion, mee_ion(ee)] 
             exempt_descrip = [exempt_descrip, mee_descrip(ee)]
         endif  
      sub = where((abs(c_wvl) gt i_w-2*nsig(jj)*lwdt(jj)>0) and (abs(c_wvl) lt i_w+2*nsig(jj)*lwdt(jj)) $
            and (c_int gt lthr*i_int(jj)))   
      if total(sub) ge 0 then begin ; if there are suspects within nsig*lwdt then start
        for ns = 0, n_elements(sub)-1 do begin 
           cndx = sub(ns) 
           c_w  = (abs(c_wvl(cndx))) 
           ;exclude exempt suspects from cross-exmination (after which mult
           ;parsing becomes obselete
           exempt = 1 
           if (total(where(exempt_w eq  c_wvl(cndx))) lt 0) then  exempt = -1 else $ 
           if (total(where(exempt_z eq  c_z(cndx)))  lt 0)  then  exempt = -1 else $ 
           if (total(where(exempt_ion eq c_ion(cndx))) lt 0) then exempt = -1 else $ 
           if (total(where(exempt_descrip  eq  c_descrip(cndx))) lt 0)  then   exempt = -1 
           if exempt eq -1 then begin 
               S = ( abs(c_w-i_w)-nsig(jj)*lwdt(jj) )/(sqrt(2)*lwdt(jj))
               if S lt 0 then frc = 0.5+errorf(abs(S))/2 else frc= (1-errorf(S))/2
               subfrac = [subfrac, frc]
               nusub = [nusub,cndx]
           endif
        endfor        
        if n_elements(nusub) eq 1 then nusub = [nusub, -1] 
        if n_elements(subfrac) eq 1 then subfrac = [subfrac, 0] 
        sub     = nusub[1:*]   & nsub    = [nsub, n_elements(sub)] 
        subfrac = subfrac[1:*] & allfrac = [allfrac,subfrac] 
        allsub  = [allsub,sub]             
      endif else begin 
        sub     = nusub        & nsub = [nsub, 1]  
        subfrac = 0  &  allfrac = [allfrac, subfrac] 
        allsub  = [allsub, sub]   
     endelse
           ;is it a proxy? 
           if keyword_set(mproxy) then   px = where( bse eq (jj) ) else px = -1 
           if total(px) ge 0 then  px = px(where(px le nprx -1)) 
           if total(px) ge 0 then i_int(jj) = i_int(jj) + total(mee_flx(px))
                                  ;create exempt list        
    endif else begin ; burn uncessecary (parse brnin keyword)
       sub     = ndx[intndx(jjp):intndx(jjp+1)-1]
       subfrac = fracs[intndx(jjp):intndx(jjp+1)-1] 
       px = where(prxbse eq (jj)) 
       if total(px) ge 0 then i_int(jj) = i_int(jj) + total(multflx(px))
    endelse

   if total(sub) ge 0 then begin ;if there are suspects within nsig*lwdt 
     for nn = 0, n_elements(sub)-1 do begin ;loop through suspect line subset
        cndx = sub(nn) 
        c_w  = abs(c_wvl(cndx))   
        ;compute contribution fractions via errf:
        cntrbtn_frc = subfrac[nn] 
        ;if intensity is adequate and candidate not on exempt list then encarcerate
        if (c_int(cndx)*cntrbtn_frc gt i_int(jj)*lthr) then begin 
          ;append contaminant line info to output arrays
          contamfrac   = [contamfrac    , cntrbtn_frc] ; contain contrubution fractions     
          victim = [victim  , jj         ] ; lstr indices of interesting line being contaminated
          culprit= [culprit , sub(nn)    ] ; sstr indices of contaminant lines 
          contamint = [contamint  , c_int(cndx)] ; the contaminant line intensities
          wvl    = [wvl     , c_wvl(cndx)] ; relevant wavelength
          if arg_present(mixstr) or arg_present(idtag) then begin 
            z      = [z       , c_z(cndx)  ] ; relevant atomic number
            src    = [src     , sorc(c_src(cndx))] ; source (which database?)
            ion    = [ion     , c_ion(cndx)] ; relevant ionic state 
            jon    = [jon     , c_jon(cndx)] ; relevant significant (ionic-equilibria wise) ionic state 
            descrip= [descrip , c_desig[1,cndx]+' '+c_config[1,cndx]+ $
                     ' - '+c_desig[0,cndx]+' '+c_config[0,cndx] ]; description 
          endif
        endif
      endfor ;end loop through suspects
   endif ; if ther are suspects to begin with 
   if total(where(victim eq jj)) lt 0 then begin 
      crctn_factor(jjp) = 1D 
      if vv gt 5 then print, 'No blends found for '+ atom[i_z(jj)-1]+roman[i_ion(jj)-1]$
             +' at '+strcompress(string(i_wvl(jj)),/remove_all)
   ;calculate correction factor for each victim 
   endif else begin 
      qq = where(victim eq jj)
      crctn_factor(jjp) = i_int(jj)/(i_int(jj)+total(contamfrac(qq)*contamint(qq)))
   endelse
 endfor ;end loop through victims
endelse ;end rmfstr present else  

;create brnstr 
if brnin eq 1 then begin   
  intndx = fltarr(nnss+1) 
  for j = 1, nnss do intndx(j) = total(nsub[0:j])
  brnstr = create_struct('FRACS'  , allfrac[1:*], $ ; contribution fractions of all non-exempt lines
                         'NDX'    , allsub[1:*] , $ ; sstr index of each of thes lines 
                         'INTNDX' , intndx , $ ; ndx and fracs index intervals corresponding to each victim
                         'NSS'    , nss    , $ ; lstr index of each victim  
                         'NNSS'   , nnss   , $ ; number of actual victims i.e. n_elements(nss) 
                         'PRXBSE' , bse[0:(nprx-1)>0], $; lstr index of the proxy in each proxy relation. 
                         'MULTSTR', multstr )  ; rd_line structure of proxee in each relation 
endif


if n_elements(culprit) gt 1 then begin 
  ;get rid of junk element -1 
  contamfrac = contamfrac[1:*] & victim = victim[1:*] & culprit = culprit[1:*] 
  contamint = contamint[1:*] & wvl = wvl[1:*] & vint = i_int
  ;correct for over-correction (e.g. one culprit contaminates two victims)
  ochist = histogram(culprit)  ; how many victims per culprit
  oo = where(ochist gt 1)  ; cull the culprits with multiple offenses
  culiq = oo + min(culprit) 
if total(oo) ge 0 then begin ; any overcorrections?
  for qq  = 0, n_elements(culiq)-1 do begin  ; loop throu the culprits with multiple offences
      zoo = where((culprit) eq culiq(qq)) 
      viciq = victim(zoo)    
      viciqw= (i_wvl(viciq)) & culiqw=(c_wvl(culiq(qq)))     
      octst = abs(culiqw)-abs(viciqw) 
      zertst = where(octst eq 0) 
      if total(zertst) ge 0 then begin  ; if exact wavlength match between culprit and victim(s) 
          for ww = 0,n_elements(zoo)-1 do if total(where(zertst eq zoo(ww))) ge 0 then $
          contamfrac(zoo(ww)) = 0 else contamfrac(zoo(ww))= 1/n_elements(zertst) ;???
      endif else begin 
          rblndrs = where(octst lt 0)  &  lblndrs = where(octst gt 0)  
          if total(rblndrs) ge 0 and total(lblndrs ge 0) then begin ; if there exist victims on both sides
              jnk = min(octst(rblndrs),rmaxbndx) & jnk = max(octst(lblndrs),lmaxbndx)
              rmaxb = viciqw(rblndrs(rmaxbndx))   &   lmaxb = viciqw(lblndrs(lmaxbndx))    
              rmaxvi= viciq(rblndrs(rmaxbndx)) & lmaxvi = viciq(lblndrs(lmaxbndx))   
              ; if only one group of victims (of two possible) 
              if lmaxb + nsig(lmaxvi)*lwdt(lmaxvi) gt rmaxb - nsig(rmaxvi)*lwdt(lmaxvi) then begin 
                  ;noctst = total(abs(octst)) 
                  nuwt = total(abs(octst))/abs(octst) 
                  nufact = nuwt/total(nuwt)  
                  for ww = 0,n_elements(zoo)-1 do contamfrac(zoo(ww))=nufact(ww)            
              endif else begin ; two groups of victims 
                  ;noctst = total(abs(octst(rblndrs))) 
                  nuwt = total(abs(octst(rblndrs)))/abs(octst) 
                  nufact = nuwt/total(nuwt)  
                  roo = zoo(rblndrs) & rmxfrac = contamfrac(zoo(rmaxbndx))
                  for ww = 0,n_elements(roo)-1 do contamfrac(roo(ww))=rmxfrac*nufact(ww)
                  ;noctst = total(abs(octst(lblndrs))) 
                  nuwt = total(abs(octst(lblndrs)))/abs(octst) 
                  nufact = nuwt/total(nuwt)  
                  loo = zoo(lblndrs) & lmxfrac = contamfrac(zoo(lmaxbndx))
                  for ww = 0,n_elements(loo)-1 do contamfrac(loo(ww))=lmxfrac*nufact(ww)  
              endelse
        endif else begin ; if there are victims from only one side we can assume one group
              ;noctst = total(abs(octst)) 
              nuwt = total(abs(octst))/abs(octst) 
              nufact = nuwt/total(nuwt)  
              for ww = 0,n_elements(zoo)-1 do contamfrac(zoo(ww))=nufact(ww)   
          endelse
      endelse
  endfor; any overcorrections?
  for jj = 0, n_elements(i_nw)-1 do begin 
      qq = where(victim eq jj)
      if total(qq) gt 0 then crctn_factor(jj-min(nss)) = $
        i_int(jj)/(i_int(jj)+total(contamfrac(qq)*contamint(qq))) else $ 
        crctn_factor(jj-min(nss)) = 1
  endfor
endif

;prepare outputs
 ;prepare tag string 

  vidtag = atom[i_z(nss)-1]+roman[i_ion(nss)-1]+' '+string(i_wvl(nss))+' '+sorc(i_src(nss)) +' '+ $
  i_descrip(nss) 

 ; idtag = atom(z[1:*]-1)+roman(ion[1:*]-1)+' '+ string(wvl)+' '+ src +' '+ descrip 
  if arg_present(mixstr) or arg_present(idtag) then begin 
 ;construct structure  
  idtag = atom[z[1:*]-1]+roman[ion[1:*]-1]+' '+ string(wvl)+' '+ src[1:*] +' '+ descrip[1:*] 
  mixstr = create_struct('CONTAMFRAC',contamfrac,'CONTAMINT',contamint,'VINT', vint, $
                         'VICTIM',victim,'CULPRIT',culprit,'WVL',wvl,'IDTAG',idtag,'VIDTAG',vidtag)
  endif
endif

if is_keyword_set(obsflx) and keyword_set(abndupdt) then begin 
   if n_elements(obsflx) eq i_nw then begin  
      if not keyword_set(ofsig) then ofsig = 1+sqrt(0.75+obsflx) else begin 
         if n_elements(obsflx) ne n_elements(ofsig) then ofsig = 1+sqrt(0.75+obsflx)  
     endelse 

      nuabund = mcmc_abund(dem, dlogt,i_wvl, obsflx, ofsig, i_emis, i_z, asig, zout, $ 
                           crctn_factor=crctn_factor,/nosol) 
      abund(zout-1) = nuabund 
   endif else begin 
      message, "Number of lines specified in OBSFLX does't match number of lines in LSTR, sorry", /info 
      return, -1L 
  endelse 
endif 

return, crctn_factor
end

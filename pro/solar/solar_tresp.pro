function solar_tresp,instrument,bandpass,tgrid=tgrid,ldbdir=ldbdir,date=date,$
        cdbdir=cdbdir,filter=filter,rspstr=rspstr,abref=abref,wrange=wrange,$ 
        abund=abund,eqfile=eqfile,perpix=perpix,lstr=lstr,cstr=cstr,units=units,$ 
        det=det,rhessirmf=rhessirmf,nocont=nocont,zresp=zresp,_extra=e

;+ 
;function	solar_tresp
;	Computes and returns temperature response for a selection of
;       bandpasses from a variety of solar X-Ray and EUV diagnostic
;       instruments. The explicitly supported instruments and bandpasses and their
;       corresponding ouput units (at the sun) are as follows:
;
;       TRACE 171,195,284 in [DN s-1 cm-2] per [cm-5 ln(K)-1] 
;       EIT   171,195,284 in [DN s-1 cm-2] per [cm-5 ln(K)-1] 
;       SXT   Al1, AlMg  in [DN s-1 cm-2] per [cm-5 ln(K)-1] 
;       RHESSI GD, in [Counts s-1] per [cm-5 ln(K)-1] 
;       CDS all GIS and NIS lines in [ergs s-1 cm-2] per [cm-5 ln(K)-1] 
;
;       The user may also specify an arbitrary instrument bandpass combination 
;       (one other than the explicitly supported ones listed above) with the
;       limitation that the WRANGE keyword must be set. The units in this case 
;       would be [Counts s-1] per [cm-5 ln(K)-1] if and effective area is included and 
;       [ph s-1 cm+2] per [cm-5 ln(K)-1]  if it is not. 
;
;       A choice of atomic databases, abundances, and ion equilibria
;       may be set via keywords LDBDIR,ABUND,and EQFILE respectively. 
;       the temperature response calculations are made like this: 
;
;  
;
;syntax 
;       resp=solar_tresp(instrument,filter,date,tgrid=tgrid,ldbdir=ldbdir,$
;       abund=abund,eqfile=eqfile) 
;
;parameters        
;       instrument  [INPUT;required] name of the relevant instrument
;                   * can be 'TRACE','EIT','SXT','CDS','RHESSI'
;                   * may also be set to a string descrption of 
;                     arbitrary instruments e.g. 'RESIK'. keyword 
;                     wrange must be set.
;       bandpass    [INPUT;required] select the relevant bandpass
;                   * can be '171','195','284' for 'TRACE' and 'EIT'
;                   * can be 'Open','Al1','AlMg','Be119','AL12','MG3' for 'SXT'
;                   * can be 'GD' for 'RHESSI' 
;                   * for 'CDS' must be a string specifiying an ionic 
;                     state (element abbreviation and roman
;                     numeral) and a transition wavelength.                     
;                     e.g. 
;                          'Mg VII 319.03' 
;                          'Ne VI 401.14' 
;                          'Fe IX 171.07' 
;                      
;                     CDS lines are treated as very narrow
;                     bandpasses, so that possible contributions 
;                     from weak lines are accounted for. 
;
;                     The transition wavelength must fall into one of
;                     the wavelength regions covered by GIS or NIS. 
;                     
;                     If the transition wavelength is designated WVL
;                     and the detector resolution as RES, then the  
;                     bandpass wavelength range is determined
;                     to be [WVL-3*RES,WVL+3*RES]. RES is 
;                     hardcoded as 0.18 Ang for NIS1, 0.32A for NIS12
;                     and 0.21 for GIS. This range may be overridden by 
;                     using keyword WRANGE.
;
;                  * can be set to and arbitrary string descrption 
;                    of bandpass if instrument is set to something other 
;                    than the explicitley supported instruments. e.g. 
;                    '1.123 - 1.293 AA' for instrument 'RESIK'
;                    keyword WRANGE must be set. 
;
;keywords 
;       filter      [INPUT] Filter wheel selection for the EIT and 
;                   TRACE passbands 
;                   * default is 'clr' for EIT can also be 'al1', 'al2'
;                   * default is 'AO' for TRACE can also be 'OO',
;                     'OA' or 'AA'                            
;                   * ignored if instrument selection is SXT or CDS
;       date        [INPUT] time in any Yohkoh format or a Yohkoh 
;                     index structure. 
;                   * see e.g. ANYTIME2INTS() SolarSoft function for 
;                     standard Yohkoh formats
;                   * SXT temperature responses time variable due to an enterance 
;                     filter failiure (Oct 1992). The value of this keyword will 
;                     be passed to GET_YO_DATES to determine appropriate correction.
;                   * if unset no correction is made. 
;       wrange      [INPUT] a two elements floating point vector which 
;                     denotes a wavelength range. overrides the internally
;                     determined bandpass range.           
;       tgrid       [INPUT] Log(T) grid on which the ouput response 
;                     is to be defined. 
;                   * default is PoA system variable !LOGT 
;                   * if !LOGT does not exist then tgrid=4.0+findgen(41)*0.1
;       lstr        [I/O] line emissivity structure  
;       cstr        [I/O] continuum emissivity structure
;       ldbdir      [INPUT] line emissivity database directory 
;                   * default is '$CHIANTI' can also be '$SPEX', '$APED'
;                     or an explicit directory specification
;                     e.g. '/data/foo/bar/spex/' or '../foo/bar/aped' etc. 
;       cdbdir      [INPUT] continuum emissivity database directory
;                   * default is '$CONT'
;       abref       [INPUT] string reference to the set of abundances to use 
;                   * default is 'feldman' which will call GETABUND() to 
;                   retrieve Feldman et al. (1992) abundances 
;                   * see GETABUND() documentation for other options.
;       eqfile      [INPUT] string reference to the ion equilibria file to use 
;                   * default is 'ioneq/mazzotta_etal.ioneq' 
;                   * see local CHIANTI database for other options
;       perpix      [INPUT] if set, responses are returned as 
;                   [DN s-1 pixel-1] per [cm-5 ln(K)-1]. 
;                   * ignored for SOHO/CDS 
;       det         [INPUT] if set then responses are computed at detector 
;                     and not at the SUN
;       rhessi_rmf  [INPUT] RHESSI spectrometer OGIP-compliant rmf(srm) fits file from which 
;                     and effective area will be extracted via rd_ogip_rmf() 
;       rspstr      [OUTPUT] output structure containing all
;                   ouputs and inputs as well as general info
;                   concerning what went into response calculation. 
;                   * The structure is defined with these tags: 
;                     INSTR           STRING    instrument parameter
;                     BANDPASS        STRING    bandpass parameter 
;                     FILTER          STRING    filter keyword
;                     RESP            DOUBLE    computed response
;                     TGRID           FLOAT     response Log(T) grid
;                     EFFAR           DOUBLE    effective area
;                     WVLAR           DOUBLE    effective area wavelength grid
;                     LDBDIR          STRING    ldbdir keyword;                     
;                     ABREF           STRING    abref keyword
;                     EQFILE          STRING    eqfile keyword
;                     ABUND           FLOAT     abund keyword
;       nocont      [INPUT] if set, continuum will not be included in response
;       zresp       [INPUT] if set, 
;       _extra	    [INPUT] use to transmit defined keywords to called functions
;
;restrictions 
;       * This function uses various SolarSoft routines 
;         and PINTofALE routines. It also employs files 
;         from Solar Soft databases. Note also that SolarSoft 
;         installations can vary in what missions/instruments 
;         are supported. 
;       * employs PINTofALE subroutines
;         LICOSPEC 
;         LINEFLX 
;         GETABUND
;         RD_LINE
;       * employs SolarSoft subroutines  
;         SSW_INSTRUMENTS 
;         CONCAT_DIR 
;         RESTGEN 
;         GET_YO_DATES
;         TRACE_EUV_RESP 
;         EIT_PARMS
;       * If the DATE keyword is set, user cannot set the verbose keyword 
;         due to a conflict in GET_YO_DATES
;history
;      liwei lin (Aug 2004)
;       BUGFIX (Jan 2005) had assumed that licospec output emissivities were in 
;        ergs cm-2 s-1. now /noph is fed through _extra to lineflx()
;       BUGFIX (LL Sep 2005)  Thanks to Harry Warren who points out that upon reading the SXT
;        effective area file, the FILTER entry already includes ENTRANCE, CCD,
;        and MIRROR effective areas. 
;       BUGFIX (LL Sep 2005) the reported units, both in documentation and
;        output structures were incorrect. Missing were terms 1d44 and cm^-5 
;       EDIT (LL Sep 2005) so that there is no hard reliance on PoA system
;        variables. RHESSI default bandpass no longer 'G4'; user is asked 
;        if RHESSI SRM file not fed via RHESSI_RMF keyword.  
;       ADDED (LL Sep 2005) ZRESP keyword 
;       BUGFIX (LL Oct 2005) wvlar returned wmid rather than wvlar
;       EDIT (LL Oct 2006) corrected reported UNITS
;- 
;     check and parse input
np = n_params(0)
if np eq 0 then begin
  print,'Usage:resp=solar_tresp(instrument,bandpass,tgrid=tgrid,ldbdir=ldbdir,$'
  print,'cdbdir=cdbdir,filter=filter,rspstr=rspstr,_extra=e)'
  return, -1L
endif
instr = strcompress(instrument,/remove_all)
if instr ne 'CDS' then bdps  = strcompress(bandpass,/remove_all) else $ 
bdps = strcompress(bandpass)
mshn  = 'IMPSSBL' & xplct = 1
if instr eq 'EIT' or instr eq 'CDS' then mshn ='SOHO' 
if instr eq 'SXT' then mshn = 'YOHKOH' 
if instr eq 'RHESSI' then mshn = 'RHESSI' 
if instr eq 'TRACE' then mshn = 'TRACE' 
if mshn eq 'IMPSSBL' then begin 
   if not keyword_set(wrange) then begin
      print, 'phooey' 
      message, 'Sorry, WRANGE keyword must be set when INSTRUMENT'+$
      ' is not explicitly supported.',/info  
      return, -1L 
   endif 
   mshn = instr  & xplct = 0
endif 
;if mshn  eq 'IMPSSBL' then begin 
;  message,'Sorry, '+instr+' is not supported by this funcion',/info
;  return, -1L 
;endif 
if instr eq 'EIT' or instr eq 'TRACE' then begin 
  bparr = ['171','195','284'] 
  bptst  = strmatch(bparr,bdps) 
  if total(bptst) eq 0 then begin 
    message,'Sorry, bandpass '+bdps+' not understood for '+instr+'.',/info
    return, -1L
  endif
endif 
if instr eq 'SXT' then begin 
  bparr = ['Open','Al1','AlMg','Mg3','AL12','Be119']
  bptst  = strmatch(bparr,bdps)
  if total(bptst) eq 0 then begin 
    message,'Sorry, bandpass '+bdps+' not understood for '+instr+'.',/info
    return, -1L 
  endif
endif 
if instr eq 'RHESSI' then begin 
  bparr = ['GD']
  bptst  = strmatch(bparr,bdps)
  if total(bptst) eq 0 then begin 
    message,'Sorry, bandpass '+bdps+' not understood for '+instr+'.',/info
    return, -1L 
  endif
endif 
if instr eq 'CDS' then begin 
  carr = strsplit(bdps,/extract) 
  ;tst = execute('jnk=float(carr(3))')  <-- IDL BUG !?
  wvl  = float(carr(2))
  smwr = [[151,221],$;GIS1   (these are launch values and may have changed 
          [256,338],$;GIS2    apologies to the CDS specialist)
          [393,493],$;GIS3
          [656,785],$;GIS4
          [308,381],$;NIS1 
          [513,633]];NIS2 
  smtr = where(smwr[0,*]-wvl le 0 and smwr[1,*]-wvl ge 0) & smtr = smtr(0)
  if smtr eq -1L then begin 
    message,'WARNING: the line you specified is out of CDS range.',/info
  endif 
endif 

;     check if required SolarSoft packages are available 
if xplct then rqib = mshn+'/'+instr else begin 
rqib = mshn & bdps = bandpass
endelse
aqib = ssw_instruments()

if total(strmatch(aqib,rqib)) eq 0 and rqib ne 'RHESSI/RHESSI' and $
xplct eq 1 then begin 
  message,'Sorry, the required SolarSoft package is not available for '$
  +rqib+'.',/info
  return, -1L
endif

;     check and parse keywords 
if not keyword_set(filter) and rqib eq 'SOHO/EIT' then filter = 'clr'
if rqib eq 'SOHO/EIT' then begin 
  if total(strmatch(['clr','al1','al2'],filter)) eq 0 then begin 
    message,'Filter instrument mismatch',/info
    return,-1L
  endif
endif
if not keyword_set(filter) and rqib eq 'TRACE/TRACE' then filter = 'AO' 
if rqib eq 'TRACE/TRACE' then begin 
  if total(strmatch(['AO','OO','OA','AA'],filter)) eq 0 then begin 
    message,'Filter instrument mismatch',/info 
    return,-1L 
  endif 
endif 
pi = 3.14159 
au = 1.49599e+13 ;[cm]
rsun = 6.96900e+10 ;[au] 
kevang = 12.398521 ;[cm] 
if not keyword_set(tgrid) then begin 
  defsysv,'!LOGT',exists=logte 
  if logte eq 0 then tgrid=4.0+findgen(41)*0.1 else tgrid=!LOGT 
endif  
if not keyword_set(filter) then filter = '-' 
if not keyword_set(abref)  then abref  = 'feldman' 
if not keyword_set(abund)  then abund  = getabund(abref,_extra=e) 
if not keyword_set(eqfile) then eqfile = 'ioneq/mazzotta_etal.ioneq' 
if not keyword_set(nh_ne) then nh_ne = 0.83

;     read in effective area and establish wrange and wgrid 
if rqib eq 'SOHO/EIT' then begin 
  wvlar  = 170 + 0.1d*findgen(1300)
  effar  = eit_parms(wvlar,bdps,filter)
  ow     = where(effar gt 0,mow) 
  wgrid  = wvlar(ow)
  wrange = minmax(wgrid)
endif
if rqib eq 'SOHO/CDS' then begin
  if smtr lt 3 then res = 0.21 ; GIS 1-4
  if smtr eq 4 then res = 0.18 ; NIS 1 
  if smtr eq 5 then res = 0.32 ; NIS 2
  wrange = [wvl-3*res,wvl+3*res] 
  wgrid  = wrange(1)+findgen(30)*((wrange(1)-wrange(0))/30)
  ; create dummy cstr structure so licospec won't choke if wvl out of cont db range
  crange = [0.8*wrange(0),wgrid,1.2*wrange(1)] & crange1= 0.5*(crange[1:*]+crange)
  cstr = {cont_int:fltarr(81,n_elements(crange)-1),logt:tgrid,wvl:crange,$
          tkev:crange/1000d,ekev:crange/1000d,midwvl:crange1}
endif
if rqib eq 'TRACE/TRACE' then begin 
  effar  = trace_euv_resp(bdps,filter,lambda=wvlar) 
  ow     = where(effar gt 0,mow) 
  wgrid  = wvlar(ow)
endif 
if rqib eq 'YOHKOH/SXT' then begin 
  afile  = concat_dir('$DIR_SXT_SENSITIVE','sra*genx')
  afile  = findfile(afile) 
  restgen,file=afile, str=area, head = head,text=text;CHIANTI
  filt  = where(bptst eq 1) 
  effar = area.filter[filt,*];*area.mirror*area.ccd*area.entrance
  effarO= area.filter[0,*];*area.mirror*area.ccd 
  wvlar=area.wave  
  ow    = where(effar gt 0)
  wgrid = area.wave(ow)
endif 
if rqib eq 'RHESSI/RHESSI' then begin 
  if keyword_set(rhessi_rmf) then afile = rhessi_rmf else begin 
    ;afile = findfile(!ardb+'*srm*')  
    afile = ''
    message, $
    'Please enter RHESSI SRM file (to avoid this message use keyword RHESSI_RMF):',/info
    read, afile
    afile = string(afile)
  endelse
  rrmf = rd_ogip_rmf(afile(0), effar = effar) 
  egrid = [rrmf.elo(0),rrmf.ehi] & egrid = (egrid+egrid[1:*])/2
  wgrid = kevang/egrid 
  wvlar = wgrid
endif
if keyword_set(wrange) then begin 
   wgrid  = wrange(0)+findgen(30)*((wrange(1)-wrange(0))/30)
 ; create dummy cstr structure so licospec won't choke if wvl out of cont db range
  crange = [0.8*wrange(0),wgrid,1.2*wrange(1)] & crange1= 0.5*(crange[1:*]+crange)
  if keyword_set(nocont) then $
  cstr = {cont_int:fltarr(81,n_elements(crange)-1),logt:tgrid,wvl:crange,$
          tkev:crange/1000d,ekev:crange/1000d,midwvl:crange1}
;  dw=0.1 & nw = long((wrange(1)-wrange(0))/dw+0.1) 
;  wgrid = findgen(nw+1L)*dw+wrange(0)
endif
wrange = minmax(wgrid)
if not keyword_set(ldbdir) then  ldbdir = '$CHIANTI'
if ldbdir eq '$APED' then fold_ioneq=0 else fold_ioneq=1 
if ldbdir eq '$APED' then abund=abund/getabund('anders & grevesse')
licospec,wgrid,lspec,cspec,lstr=lstr,lemis=lemis,cemis=cemis,$
cstr=cstr,abund=abund,ldbdir=ldbdir,eqfile=eqfile,_extra=e

;     compute flux for each Temperature  
wmid=0.5*(wgrid[1:*]+wgrid)
au2=au^2
;if keyword_set(det) then au2=au^2 else au2 = 1
unit_conv0=1d/(4*pi*au2)/alog(10)/0.05d ; for SXT 
unit_conv1=1d/(4*pi)/alog(10)/0.05d     ; for EIT, CDS,and TRACE
unit_conv2=1d/alog(10)/0.05d            ; for RHESSI and misc.
unit_conv3=(rsun^2)/(au^2)/alog(10)/0.05d
resp = tgrid
resp0 = resp
if rqib eq 'SOHO/CDS' then temis=lemis else temis=lemis+cemis
if keyword_set(nocont) then temis = lemis

for i=0L,n_elements(resp)-1L do begin
  emd=0.D*tgrid & emd[i]=1d44;*unit_conv
  if rqib eq 'YOHKOH/SXT' then begin 
    resp[i]=unit_conv0*total(lineflx(temis,$
    tgrid,wmid,dem=emd,effar=effar,wvlar=wvlar,/noph,_extra=e))
    resp0[i]=unit_conv0*total(lineflx(temis,$ 
    tgrid,wmid,dem=emd,effar=effarO,wvlar=wvlar,/noph,_extra=e))
    ; resp units are [ ergs s-1 cm-2 ] 
  endif 
  if rqib eq 'SOHO/CDS' then begin 
    resp[i]=unit_conv1*total(lineflx(temis,$
    tgrid,wmid,dem=emd,effar=effar,wvlar=wvlar,/noph,_extra=e))
    ; resp units are [ ergs s-1 cm-2 ] for CDS 
  endif 
  if rqib eq 'RHESSI/RHESSI' then begin 
    resp[i]=unit_conv2*total(lineflx(temis,$
    tgrid,wmid,dem=emd,effar=effar,wvlar=wvlar,/noph,_extra=e))
    ; resp units are [ Counts s-1 ] for RHESSI 
  endif
  if rqib eq 'TRACE/TRACE' or rqib eq 'SOHO/EIT' then begin 
    resp[i]=unit_conv1*total(lineflx(temis,$
    tgrid,wmid,dem=emd,effar=effar,wvlar=wvlar,_extra=e))
    ; resp units are [ DN s-1 pixel-1 ] for EIT and TRACE 
  endif  
  if xplct eq 0 then begin 

    resp[i]=unit_conv3*total(lineflx(temis,$
    tgrid,wmid,dem=emd,effar=effar,wvlar=wvlar,_extra=e))
 ;   resp[i]=unit_conv2*total(lineflx(temis,$
 ;   tgrid,wmid,dem=emd,effar=effar,wvlar=wvlar,/noph,_extra=e))    
    ; resp units are [PH s-1 cm-2] if effar is not set 
    ; resp units are [Counts s-1] if effar is set 
  endif 
endfor
;     convert to DN units and make entrance filter correction if necessary
if rqib eq 'YOHKOH/SXT' then begin  
  if keyword_set(date) then $ 
    fracopn=get_yo_dates(date,/entr,/value,_extra=e) else fracopn = 0d  
    resp = resp*(1-fracopn)+resp0*fracopn
  for i=0L,n_elements(resp)-1L do begin 
    erg2ev =  (1d-7)/(1.602d-19) ; 1e-7 [J/erg]  1.602e-19 [J/eV]->[eV/erg]  
    resp(i) = [(1/3.65d)*resp(i)*erg2ev]/100D ; 100 [DN/eV] SXT GAIN
  endfor
  ; resp units are here [ DN s-1 cm-2]
endif

;     convert to units of [DN s-1 cm-2] from [DN s-1 pixel-1] 
if rqib eq 'YOHKOH/SXT'  then arc_pix = 4.9218750 ; [arcsecs pixel-1]
if rqib eq 'TRACE/TRACE' then arc_pix = 0.501     ; [arcsecs pixel-1]
if rqib eq 'SOHO/EIT'    then arc_pix = 2.6       ; [arcsecs pixel-1]
if rqib ne 'SOHO/CDS' and rqib ne 'RHESSI/RHESSI' and xplct eq 1 then begin 
  ds  = arc_pix*pi/(60d*60d*180d) ; [radians pixel -1] 
  cm2 = (au*ds)^2  ; [cm+2 pixel-1]
endif else cm2 = 1

if keyword_set(perpix) then begin 
  if rqib eq 'YOHKOH/SXT' then resp = resp*cm2 
  ; resp units are here [ DN s-1 pixel-1]
  units = '[DN s-1 pixel-1] per [cm-5 ln(K)-1]' 
endif else begin 
  if rqib eq 'TRACE/TRACE' then resp = resp/cm2 
  if rqib eq 'SOHO/EIT' then resp = resp/cm2 
  ; resp units are here [ DN s-1 cm-2] 
  units = '[DN s-1 cm-3] per [cm-5 ln(K)-1]' 
endelse 
if rqib eq 'SOHO/CDS' then units = '[ergs s-1 cm-2] per [cm-5 ln(K)-1]' 
if rqib eq 'RHESSI/RHESSI' then units = '[Counts s-1] per [cm-5 ln(K)-1]'
if keyword_set(e) then effset= total(strmatch(tag_names(e),'EFFAR')) $ 
else effset = 0 
if keyword_set(e) then wvlset= total(strmatch(tag_names(e),'WVLAR')) $ 
else wvlset = 0 
if effset then effar = e.EFFAR
if wvlset then wvlar = e.WVLAR ;else wvlar = wmid 
if xplct eq 0 then if effset then $
units = '[Counts s-1] per [cm-5 ln(K)-1]' else units = '[ph s-1 cm-2] per [cm^-5 ln(K)-1]'

;    convert to units of [DN s-1 cm+3] from [DN s-1 cm-2] 
dem  = 0.D*tgrid+1d44 
nT = n_elements(tgrid) 
dem[0]=0.5*dem[0] & dem[nT-1]=0.5*dem[nT-1] ;accnt for trpzdl integration
resp = resp/dem 

;    create zresp array 
zrespi=1 
if n_elements(zresp) gt 0 then begin 
  if zresp[0] eq 0 then zrespi = 0 
endif
if keyword_set(zrespi) then begin
  nabn = n_elements(abund)
  zresp = fltarr(nT,nabn)
  for j = 0, nabn-1 do begin      
     oo = where(lstr.z eq j) 
     if oo[0] ne -1 then begin 
         zlstr=cat_ln(lstr,pick=[oo]) 
         zresp[*,j]=solar_tresp(instrument,bandpass,tgrid=tgrid,lstr=zlstr,$
         /nocont,zresp=0,/allah,wrange=wrange,_extra=e)
     endif
  endfor
endif  

;     create output structure 
if arg_present(rspstr) then begin 
   if rqib eq 'SOHO/CDS' or xplct eq 0 and not effset then begin
     rspstr = {instr:rqib,$               
               bandpass:bdps,$
               filter:filter,$
               resp:resp,$  
               zresp:zresp,$
               units:units,$ 
               tgrid:tgrid,$ 
               ldbdir:ldbdir,$ 
               abref:abref,$ 
               eqfile:eqfile,$ 
               abund:abund} 
   endif else begin 
     rspstr = {instr:rqib,$               
               bandpass:bdps,$
               filter:filter,$
               resp:resp,$
               zresp:zresp,$
               units:units,$   
               tgrid:tgrid,$ 
               effar:effar,$ 
               wvlar:wvlar,$ 
               ldbdir:ldbdir,$ 
               abref:abref,$ 
               eqfile:eqfile,$ 
               abund:abund}
   endelse
endif
return, resp 
end

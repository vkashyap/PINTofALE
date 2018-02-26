function line_chianti,elem,pres,temp,ang,trans,kev=kev,ikey=ikey,jkey=jkey,$
	econf=econf,wmn=wmn,wmx=wmx,n_e=n_e,istate=istate, _extra=e
;+
;FUNCTION	line_chianti
;	returns a 2D array containing line cooling emissivities at a set of
;	temperatures for all available spectral lines (in a CHIANTI style
;	database [http://wwwsolar.nrl.navy.mil/chianti.html]) for the
;	specified element at the specified pressure or density
;	[1e-23 erg cm^3 s^-1] (NTEMP,NANG)
;	NOTE: Ion Balance is NOT included!!!
;
;SYNTAX
;	fx=line_chianti(elem,pres,logT,ang,trans,kev=kev,ikey=ikey,jkey=jkey,$
;	econf=econf,wmn=wmn,wmx=wmx,n_e=n_e,istate=istate,/tex,chidir=chidir)
;
;PARAMETERS
;	elem	[INPUT; required] name of element e.g., "He", "C", etc.
;		* may specify ionic state (e.g., 'Fe XXII')
;		* if an array, only the first element is looked at
;	pres	[INPUT; default: 1e15 cm^-3 K] pressure at which to compute
;		level populations
;	temp	[INPUT/OUTPUT] log(Temperature[K]) at which to compute the
;		line fluxes.  if not given, then TEMP=4.0+FINDGEN(81)*0.05
;		* best results if evenly spaced in log(T)
;	ang	[OUTPUT] all the wavelengths [Angstrom]
;	trans	[OUTPUT] 2D array of transitions --
;		trans[0,*] lower levels; trans[1,*] upper levels
;
;KEYWORDS
;	kev	[OUTPUT] all the line energies (same as ANG, but in [keV])
;	ikey	[OUTPUT] ionic state of element producing the line at ANG
;	jkey	[OUTPUT] ionic state that matters to ion balance
;	econf	[OUTPUT] 2D string array of e-configurations --
;		econf[0,*] lower levels; econf[1,*] upper levels
;	wmn	[INPUT; default: 0 Angstrom] minimum value of wavelength
;		to include in output
;	wmx	[INPUT; default: 900 Angstrom] maximum value of wavelength
;		to include in output
;	n_e	[INPUT] electron density [/cm^3]
;		* OVERRIDES values determined using PRES and TEMP
;	istate	[INPUT; default: ALL; overridden by ELEM] limit calculations
;		to specified ionic states (1=neutral, 2=singly ionized, etc.;
;		can be an array)
;
;	_extra	[INPUT ONLY] used to specify
;		* /TEX -- TRANS in TeX format (1) or not (0)
;		* CHIDIR -- path name to the CHIANTI installation
;
;RESTRICTIONS
;	* uses the CHIANTI database
;	* requires subroutines:
;	  -- GET_IONLIST
;	  -- SYMB2ZION [LAT2ARAB]
;	  -- RD_CHIANTI [READ_WGFA2(%), READ_ELVLC(%), READ_SPLUPS(%), READ_IP(%)]
;	  -- POPSOL [DESCALE_UPS(%), NENH]
;	  -- TRANSLABEL
;	  -- KILROY
;	  -- IS_KEYWORD_SET
;	  (%): CHIANTI subroutine, used as is
;
;HISTORY
;	vinay kashyap (Nov96)
;	removed ion balance, changed input keyword CONF to output ECONF,
;	  changed mult. factor from 1e30 to 1e23 (VK; Dec96)
;	added level excitation check (VK; Jun97)
;	corrected bug with TEMP specification; forced ELEM to be scalar
;	  (VK; Nov98)
;	changed call from POPLEV to POPSOL (VK; SepMM)
;	added keyword JKEY; corrected calls to modifications of GET_IONLIST,
;	  RD_CHIANTI, and POPSOL to account for dielectronic recombination
;	  lines; converted syntax to IDL5 (VK; OctMM)
;	bug correction: CHIDIR was not being passed on;	caught and
;	  filtered -ve population levels (VK; DecMM)
;	bug correction: transition labels were being filtered on wavelength
;	  twice, once within TRANSLABEL and once here.  now forced the range
;	  in TRANSLABEL to be ignored (VK; Jun01)
;	bug correction: was failing if ELEM contained ION info (VK; May02)
;	modified call to POPSOL, per CHIANTI 4.2 changes (VK; Apr04)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	made loops go over long integers (VK; 11Jun2010)
;-

np=n_params(0)
if np eq 0 then begin
  print,'Usage: fx=line_chianti(element,pressure,logT,wavelengths,transitions,$'
  print,'       kev=energies,ikey=ion_wvl_key,jkey=ion_balance_key,$'
  print,'       econf=e_configurations,wmn=wmn,wmx=wmx,n_e=n_e,istate=istate,$'
  print,'       /tex,chidir=chidir)'
  print,'  returns 2D array of line fluxes for given element at'
  print,'  all wavelengths and given temperatures'
  return,-1L
endif

;pressure
if n_elements(pres) eq 0 then pres=1e15

;temperature array
nt=n_elements(temp)
if nt eq 1 then t=[float(temp)]
if nt lt 1 then begin			;undefined
  nt=81 & t=4.0+findgen(nt)*0.05 & temp=t
endif
if nt gt 1 then t=float(temp)

;initialize
fx=fltarr(nt) & ang=[-1.] & kev=ang & ikey=[0L] & jkey=[0L]	;output arrays
ktrans=-1L
if not keyword_set(wmn) then wmn=0.	;[Angstroms]
if not keyword_set(wmx) then wmx=900.	;[Angstroms]

;figure out which element and in which ionic state
symb2zion,elem[0],z,ist
if is_keyword_set(istate) and not keyword_set(ist) then ist=fix(istate)

;make sure that data exists for (Z,ISTATE)
availon=get_ionlist(z,dielec=dielec,_extra=e)	;default list of ions
if ist[0] eq 0 then begin
  ist=availon	;"ALL"
  ide=dielec	;dielectronic recombination flag
endif else begin
  ide=0 & oldist=ist
  ;for ii=0L,n_elements(availon)-1L do begin
  for ii=0L,n_elements(oldist)-1L do begin
    oo=where(availon eq oldist[ii],moo)
    if moo gt 0 then begin
      if not is_keyword_set(ide) then begin
	ide=dielec[oo] & ist=availon[oo]
      endif else begin
	ide=[ide,dielec[oo]] & ist=[ist,availon[oo]]
      endelse
    endif
  endfor
  if n_elements(ide) ne n_elements(ist) then ide=intarr(n_elements(ist))
endelse
nion=n_elements(ist) & kion=intarr(nion)+1	;any "bad" ions?
for i=0L,nion-1L do begin
  oo=where(ist[i] eq availon)
  if oo[0] eq -1 then begin			;bad ion.  baad, baad ion.
    message,'ION state '+strtrim(ist[i],2)+' has no data; eliminating!',/info
    kion[i]=-1
  endif
endfor
oo=where(kion gt 0)				;any good ions?
if oo[0] ne -1 then ist=ist[oo] else ist=-1L
if ist[0] eq -1 then begin
  message,'NO IONS FOUND',/info & return,fx
endif
nion=n_elements(ist)

for i=0L,nion-1L do begin			;for each ion...
  jon=ist[i]
  delec=ide[i]	;dielectric recombination flag (1=yes, 0=no)
  if is_keyword_set(delec) then print,'Ionic state (DR)',jon else $
  	print,'Ionic state ',jon

  ;	read in atomic data from CHIANTI database
  ;
  ;chline=rd_chianti(z,jon,delec=delec,chidir=chidir)
  chline=rd_chianti(z,jon,delec=delec, _extra=e)
  ;-----------------------------------------------------;
  ;	CHLINE is a structure which is described by	;
  ;	CHELP=RD_CHIANTI() & HELP,CHELP,/STRUCTURE	;
  ;-----------------------------------------------------;

  ;	are there any lines in the specified range??
  ww=chline.wvl & oo=where(abs(ww) gt wmn and abs(ww) le wmx,nw)
  print,'# of lines in wavelength range = '+strtrim(nw,2)

  if nw gt 0 then begin			;(there are some lines
    ;	get level populations for the various temperatures for this ion
    ;
    nlev=chline.nlev & pop=dblarr(nt,nlev) & popt=dblarr(nlev)
    for j=0L,nt-1L do begin
      tt=10.^(t[j]) & eden=pres/tt
      ;
      ;tmp=poplev(tt,pres,chline,n_e=n_e)	;original
      ;tmp=popsol(tt,pres,chline,n_e=n_e, _extra=e)	;chianti 4
      tmp=popsol(tt,pres,chline,n_e=n_e,delec=delec,z=z,jon=jon,$
	_extra=e)	;chianti 4.2
      ;
      ;	the following 2 lines added in Dec 2000 to avoid the problem
      ;	of -ve emissivities
      mtmp=max(tmp) & obad=where(tmp lt 1d-30*mtmp,mobad)
      if mobad gt 0 then tmp[obad]=0.
      ;
      popt[0]=tmp[*] & npop=n_elements(tmp)
      if keyword_set(n_e) then popt=popt/n_e else popt=popt/eden
      ;
      pop[j,0:npop-1]=popt[0:npop-1]*1e23 > 0
      		;multiply by 1e23 just to make the numbers large and get of
		;any -ve values which are surely just numerical errors
      kilroy; was here.
    endfor	;j=0,nt-1
    print,''

    ;	for each transition, get intensities for the set of temperatures
    ;
    ntrans=chline.ntrans			;how many transitions
    splups=chline.spups
    for j=0L,ntrans-1L do begin			;{for each line
      wvl=ww[j] & l2=(chline.lev2)[j] & a_val=(chline.a)[j]

      ;;	can this level be excited?
      ;upssiz=size(splups) & upssiz2=upssiz[2]-1
      ;if l2 le upssiz2 then begin
	;maxups=max(splups(*,l2-1,0))
	;if maxups gt 0. then upstst=1 else upstst=0
      ;endif else upstst=0
      ;if abs(wvl) gt wmn and abs(wvl) le wmx and upstst gt 0 then begin

      ;	include this line in calculation?
      ;	NOTE: CHIANTI's CH_SYNTHETIC selects lines based on the
      ;	criteria evaluated below.  We however ignore the one that
      ;	rules out lines of 0 intensity based on the population
      ;	fraction, because of the manner in which our emissivities
      ;	are stored -- can't have the wavelength grid depend on the
      ;	density at which the calculation is carried out!
      ok='ok'
      siz_chianti=-1 & siz_true=max([chline.splstr.lvl1,chline.splstr.lvl2])
      if total(pop[*,l2-1]) ne 0 then siz_chianti=siz_true
      i_spl=where(chline.splstr.lvl2 eq l2,mi_spl)

      ;THIS LINE COMMENTED OUT DELIBERATELY; if siz_chianti lt 0 then ok='no lines of 0 intensity gets into line list' else $
       if l2 gt siz_true then ok='not possible to excite this level?' else $
	if mi_spl eq 0 and delec eq 1 then ok='exclude lines already in '+$
		'standard list from turning up in DR spectrum as well' else $
	 if abs(wvl) lt wmn then ok='line falls below range' else $
	  if abs(wvl) gt wmx then ok='line falls above range' else $
	   if a_val eq 0 then ok='A=0'
	

      if ok eq 'ok' then begin	;(go ahead and compute emissivity

	;------------------------------------------------------------;
        ;	you are here.  this is why you are here.             ;
	;                                                            ;
	;	emissivity = (h*c/wvl)*A*pop_fraction                ;
	;	[1e-23 erg cm^3/s] =                                 ;
	;		[(erg s * cm/s)/cm * (cm^3/s) * (1e-23)]     ;
	;                                                            ;
	;	but please do keep in mind that pop_fraction does    ;
	;	NOT include ion balance!                             ;
	;------------------------------------------------------------;

        ;	A -> cm -> hz -> erg -> keV
        if wvl eq 0 then nrg=0. else $
		nrg=((6.626d-27*(2.998d10/(abs(wvl)*1d-8)))/1.602d-9)

        hc = 6.626d-27 * 2.998d+10 * 1.d+8 / abs(wvl)
        intens=dblarr(nt)
	if siz_chianti gt 0 then for k=0L,nt-1L do intens[k]=hc*a_val*pop[k,l2-1]

        ;	append each transition to output(s)
        fx=[fx,intens] & ang=[ang,wvl] & ikey=[ikey,jon]
	jkey=[jkey,jon+delec] & kev=[kev,nrg]
	if not is_keyword_set(ow) then ow=[long(j)] else ow=[ow,long(j)]

      endif					;OK)

    endfor					;j=0,ntrans-1}

    ;level designations and configurations
    ;tmp=translabel(chline,config=config,wmn=wmn,wmx=wmx,_extra=e)
    ;	NOTE: ABOVE SHOULD NOT HAVE WMN AND WMX SAME AS IN THIS PROGRAM
    ;	BECAUSE OF OW FILTER BELOW.
    tmp=translabel(chline,config=config,wmn=-1e10,wmx=1e10,_extra=e)
    if is_keyword_set(ow) then tmp=tmp[*,ow] else tmp=''
    if is_keyword_set(ow) then config=config[*,ow]
    if is_keyword_set(ow) then print,'# included = '+strtrim(n_elements(ow),2)
    ow=0L
    if not keyword_set(itrans) then begin	;first
      if tmp[0] ne '' then begin
        trans=tmp[*] & econf=config[*] & itrans=1
      endif
    endif else begin				;subsequent
      if tmp[0] ne '' then begin
        trans=[trans,tmp[*]] & econf=[econf,config[*]]
      endif
    endelse

  endif					;wmn<ww<=wmx)
endfor					;i=0,nion-1

;output
nw=n_elements(ang)-1
if nw gt 0 then begin
  fx=fx[nt:*] & ang=ang[1:*] & kev=kev[1:*] & ikey=ikey[1:*] & jkey=jkey[1:*]
  fx=reform(fx,nt,nw) & trans=reform(trans,2,nw) & econf=reform(econf,2,nw)
endif

return,fx
end

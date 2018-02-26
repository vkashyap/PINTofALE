function cecf,wvlar,effar,wrange=wrange,chnrng=chnrng,NH=NH,$
	rmfile=rmfile,rmstr=rmstr,lstr=lstr,cstr=cstr,abund=abund,$
	ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,verbose=verbose,$
	logT=logT,lspec=lspec,cspec=cspec,ctrate=ctrate,ergfx=ergfx,$
	eps=eps,metal=metal, _extra=e
;+
;function	cecf
;	compute and return a counts-to-energy conversion factor
;	as [ergs/cm^2/ct] for a variety of temperatures
;
;syntax
;	f=cecf(wvlar,effar,wrange=wrange,chnrng=chnrng,NH=NH,$
;	rmfile=rmfile,rmstr=rmstr,lstr=lstr,cstr=cstr,abund=abund,$
;	ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,verbose=verbose,$
;	logT=logT,lspec=lspec,cspec=cspec,ctrate=ctrate,ergfx=ergfx,$
;	eps=eps,logp=logp,n_e=n_e,eqfile=eqfile,chifil=chifil,$
;	fH2=fH2,He1=He1,HeII=HeII,/fano,/wam,/bam,/mam,noHeH=noHeH)
;
;parameters
;	wvlar	[INPUT; required] array of wavelengths [Ang]
;		at which effective area is defined
;	effar	[INPUT] effective area [cm^2]
;		* must be >0; all -ve values will be set to 0
;		* size must match WVLAR.  if it does not:
;		-- if N(EFFAR)=0, then EFFAR=fltarr(N(WVLAR))+1.
;		-- if N(EFFAR)=1, then EFFAR=fltarr(N(WVLAR))+abs(EFFAR[0])
;		-- if N(EFFAR)>1 and .ne.N(WVLAR), interpolated onto
;		   minmax(WVLAR) assuming a regular grid
;
;keywords
;	wrange	[INPUT] wavelength range of spectral energy distribution
;		to include in the calculation
;		* ignored if not 2-element array
;		* if smaller than minmax(WVLAR), the overhanging parts
;		  of the effective area are thrown away
;	chnrng	[INPUT] if given, totals counts within the specified
;		channel range
;		* ignored if not 2-element array
;		* if RMF is available (see RMFILE and RMSTR below), then
;		  -- if integers, assumed to be the range in PI or PHA
;		     (whatever is defined in the RMF), else
;		  -- assumed to be the range in [keV]
;		* if RMF is not available, then assumed to be the range
;		  in [keV] if given, and corresponds to WRANGE if not given
;	NH	[INPUT] H column density
;		* default is 1e18 cm^-2; explicitly set to 0 to not use NH
;		* if value is between 10 and 30, then assumed to have been
;		  given in log10(NH)
;	rmfile	[INPUT] name of response matrix file
;		* ignored if RMSTR is valid
;	rmstr	[I/O] a response matrix structure of the type
;		output by RD_OGIP_RMF()
;		* if it exists on input, this is used instead of
;		  trying to read in from RMFILE; once read in, is
;		  returned as output via this keyword 
;	lstr	[I/O] line emissivity structure of the type read
;		in from RD_LINE()
;		* if it exists on input, this is used instead of
;		  trying to read it in from the line emissivity
;		  database; once read in, is returned as output
;		  via this keyword
;		* WARNING: the output emissivity matrix will include
;		  ion balance (but not abundances; do not use output
;		  elsewhere unless you know what you are doing)
;	cstr	[I/O] as with LSTR, for continuum emissivities
;	abund	[INPUT] abundances (see GETABUND() for details)
;		* default is Grevesse et al. (1992)
;	ldbdir	[INPUT] line emissivity database location
;		* default is '$CHIANTI'
;		* used only if LSTR is not defined on input
;	cdbdir	[INPUT] continuum emissivity database location
;		* default is '$CONT'
;		* used only if CSTR is not defined on input
;	ceroot	[INPUT] continuum emissivity files rootname
;		* default is 'cie'
;		* used only if CSTR is not defined on input
;	verbose	[INPUT] controls chatter
;	logT	[OUTPUT] the temperature grid for which the cecf is
;		calculated
;	lspec	[OUTPUT] the contribution from line emission to the [erg/s/cm^2]
;	cspec	[OUTPUT] like LSPEC, for continuum emission [erg/s/cm^2]
;	ctrate	[OUTPUT] the countrate that corresponds to LSPEC+CSPEC [ct/s]
;	ergfx	[OUTPUT] the absorbed flux that corresponds to LSPEC+CSPEC [ergs/s/cm^2]
;	eps	[INPUT] a small number; set all values below
;		EPS*max(LSPEC+CSPEC) or EPS*max(CTRATE) to !VALUES.NAN
;		* explicitly set to 0 to avoid doing this
;		* default is 1e-6
;	metal	[JUNK] here to trap and discard possible keyword to RD_CONT()
;	_extra	[INPUT ONLY] pass defined inputs to subroutines:
;		RD_LINE: LOGP, N_E
;		FOLD_IONEQ: EQFILE, CHIFILE
;		RD_CONT: LOGP, N_E
;		ISMTAU: FH2, HE1, HEII, FANO, WAM, BAM, MAM, NOHEH
;
;description
;
;requires subroutines
;	RD_LINE, SYMB2ZION, KILROY, WHEE, FOLD_IONEQ, READ_IONEQ, RD_IONEQ,
;	RD_CONT, RD_OGIP_RMF, ISMTAU, BAMABS, MID2BOUND, GETABUND, INICON,
;	IDL-Astro
;
;history
;	vinay kashyap (MMV.III)
;	modified behaviour of CHNRNG to work even in absence
;	  of RMF (VK; MMV.IX)
;	added keyword ERGFX (VK; MMXI.III)
;	dwvl for continuum forced to be +ve; APEC abundances are now
;	  removed before being added back on via ABUND (VK; MMXIII.I)
;	added more tests to catch ATOMDB (VK; MMXIV.XI)
;-

;	usage
ok='ok' & np=n_params() & nw=n_elements(wvlar) & na=n_elements(effar)
if np lt 1 then ok='Insufficient parameters' else $
 if nw eq 0 then ok='WVLAR is undefined' else $
  if nw lt 2 then ok='WVLAR is insufficient'
if ok ne 'ok' then begin
  print,'Usage: f=cecf(wvlar,effar,wrange=wrange,chnrng=chnrng,NH=NH,$'
  print,'	rmfile=rmfile,rmstr=rmstr,lstr=lstr,cstr=cstr,abund=abund,$'
  print,'	ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,verbose=verbose,$'
  print,'	logT=logT,lspec=lspec,cspec=cspec,ctrate=ctrate,ergfx=ergfx,$
  print,'	eps=eps,logp=logp,n_e=n_e,eqfile=eqfile,chifil=chifil$'
  print,'	fH2=fH2,He1=He1,HeII=HeII,/fano,/wam,/bam,/mam,noHeH=noHeH)'
  print,'  compute a counts-to-energy conversion factor [ergs/cm^2/ct]'
  print,'  for a variety of temperatures'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
arwvl=[wvlar[*]] & awmin=min(arwvl,max=awmax)
;
if na eq nw then areff=[effar[*]>0] else begin
  if na eq 0 then areff=fltarr(nw)+1. else $
   if na eq 1 then areff=fltarr(nw)+abs(effar[0]) else begin
     if vv gt 1 then $
	message,'Interpolating EFFAR to match size of WVLAR',/informational
     dx=(wmax-wmin)/na & x=findgen(na)*dx+wmin
     areff=interpol(effar,x,arwvl)
   endelse
endelse
na=n_elements(areff)
;
wmin=awmin & wmax=awmax
if n_elements(wrange) eq 2 then wmin=min(wrange,max=wmax)
;
;if n_elements(wrange) eq 2 then begin
;  wmin=min(wrange,max=wmax)
;  if wmin lt awmin then begin
;    if vv ge 5 then $
;      message,'min(WRANGE) < min(WVLAR); ignoring min(WRANGE)',/informational
;    wmin=awmin
;  endif
;  if wmax gt awmax then begin
;    if vv ge 5 then $
;      message,'max(WRANGE) > max(WVLAR); ignoring max(WRANGE)',/informational
;    wmax=awmax
;  endif
;endif

oo=where(arwvl ge wmin and arwvl le wmax,moo)
if moo eq 0 then begin
  message,'WRANGE selected out entire WVLAR; returning',/informational
  return,0.
endif
arwvl=arwvl[oo] & areff=areff[oo]
wgrid=mid2bound(arwvl)
;
nrf=n_elements(rmfile) & szr=size(rmfile)
nrs=n_elements(rmstr) & mrs=n_tags(rmstr)
if mrs ne 0 then begin
  tnamrm=tag_names(rmstr) & k=0
  for i=0L,n_elements(tnamrm)-1L do begin
    if tnamrm[i] eq 'NNRG' then k=k+1
    if tnamrm[i] eq 'ELO' then k=k+1
    if tnamrm[i] eq 'EHI' then k=k+1
    if tnamrm[i] eq 'NCHAN' then k=k+1
    if tnamrm[i] eq 'EMN' then k=k+1
    if tnamrm[i] eq 'EMX' then k=k+1
    if tnamrm[i] eq 'N_GRP' then k=k+1
    if tnamrm[i] eq 'F_CHAN' then k=k+1
    if tnamrm[i] eq 'N_CHAN' then k=k+1
    if tnamrm[i] eq 'MATRIX' then k=k+1
    if tnamrm[i] eq 'FIRSTCHAN' then k=k+1
  endfor
  if k lt 11 then begin
    message,'RMSTR not a valid RMF structure; ignoring input',/informational
    nrs=0 & mrs=0
  endif
endif
if nrs eq 0 or mrs eq 0 then begin
  if nrf ne 0 then begin
    if szr[n_elements(szr)-2L] eq 7 then rmstr=rd_ogip_rmf(rmfile[0])
  endif
endif
nrs=n_elements(rmstr) & mrs=n_tags(rmstr)
if mrs ne 0 then begin
  ;	check whether RMF has ARF folded in
  ee=0.5*(rmstr.ELO+rmstr.EHI)
  alteffar=total(rmstr.MATRIX,1)
  altareff=interpol(alteffar,12.3985/ee,arwvl)
  if abs(mean(altareff)-1)/stddev(altareff) gt 0.05 and vv ge 5 then $
    message,'RMF probably contains ARF',/informational
  ;areff=areff*altareff	;this is commented out because the RMF is folded in later using conv_rmf
  ;	figure out input spectrum grid
  elo=rmstr.ELO & ehi=rmstr.EHI & os=sort(elo) & ee=[elo[os],max(ehi)]
  wgrid=12.3985/ee
endif
;
if n_elements(abund) lt 30 then abund=getabund('grevesse et al.')
;
nls=n_elements(lstr) & mls=n_tags(lstr)
if mls ne 0 then begin
  tnaml=tag_names(lstr) & k=0
  for i=0L,n_elements(tnaml)-1L do begin
    if tnaml[i] eq 'LINE_INT' then k=k+1
    if tnaml[i] eq 'LOGT' then k=k+1
    if tnaml[i] eq 'WVL' then k=k+1
    if tnaml[i] eq 'Z' then k=k+1
    if tnaml[i] eq 'ION' then k=k+1
    if tnaml[i] eq 'JON' then k=k+1
  endfor
  if k lt 6 then begin
    message,'LSTR not a valid line emissivity structure; reading in instead',$
	/informational
    nls=0 & mls=0
  endif
endif
if nls eq 0 or mls eq 0 then begin
  if not keyword_set(ldbdir) then ldbdir='$CHIANTI'
  lconf=rd_line(atom,wrange=[wmin,wmax],dbdir=ldbdir[0],n_e=1e9,fstr=lstr,$
	verbose=vv, _extra=e)
  cldb=strlowcase(ldbdir)
  if strpos(cldb,'aped',0) lt 0 and strpos(cldb,'atomdb',0) lt 0 and strpos(cldb,'apec',0) lt 0 then $
	lconf=fold_ioneq(lstr.LINE_INT,lstr.Z,lstr.JON,verbose=vv, _extra=e) else $
  	apedance,lconf,lstr.Z
  lstr.LINE_INT = lconf
endif
nls=n_elements(lstr) & mls=n_tags(lstr)
;
ncs=n_elements(cstr) & mcs=n_tags(cstr)
if mcs ne 0 then begin
  tnamc=tag_names(cstr) & k=0
  for i=0L,n_elements(tnamc)-1L do begin
    if tnamc[i] eq 'CONT_INT' then k=k+1
    if tnamc[i] eq 'LOGT' then k=k+1
    if tnamc[i] eq 'WVL' then k=k+1
    if tnamc[i] eq 'MIDWVL' then k=k+1
  endfor
  if k lt 4 then begin
    message,'CSTR not valid continuum emissivity structure; read in instead',$
	/informational
    ncs=0 & mcs=0
  endif
endif
if ncs eq 0 or mcs eq 0 then begin
  if not keyword_set(cdbdir) then cdbdir='$CONT'
  if not keyword_set(ceroot) then ceroot='cie'
  cconf=rd_cont(ceroot,wrange=[wmin,wmax],dbdir=cdbdir[0],n_e=1e9,$
  	abund=abund,fcstr=cstr,verbose=vv, _extra=e)
endif
ncs=n_elements(cstr) & mcs=n_tags(cstr)
;
NHcol=1d18 & if n_elements(NH) gt 0 then NHcol=NH[0]>0
if NHcol gt 10 and NHcol lt 30 then NHcol=10.D^(NHcol)

;	check to make sure LSTR and CSTR are on same T grid
ntl=n_elements(lstr.LOGT) & ntc=n_elements(cstr.LOGT)
logT=lstr.LOGT
ok='ok' & if ntl ne ntc then ok='LSTR and CSTR are incompatible' else $
 if mean(abs(logT-cstr.LOGT)) gt 0.01 then ok='LSTR and CSTR on different grids'
if ok ne 'ok' then begin
  message,ok,/informational
  message,'Ignoring continuum entirely',/informational
  mcs=0
endif

;	sundry initializations
os=sort(wgrid) & wgrid=wgrid[os] & wmid=0.5*(wgrid[1:*]+wgrid)
h=6.626176e-27	;[erg s]
c=2.9979e10	;[cm/s]
EM=1d23					;[cm^-5]
lemis=lstr.LINE_INT & cemis=cstr.CONT_INT
zab=abund[lstr.Z-1]
lwvl=abs(lstr.WVL) & cwvl=cstr.midWVL & cww=cstr.WVL & dcw=abs(cww[1:*]-cww)
nrgw=h*c*1e8/wmid & nrgl=h*c*1e8/lwvl & nrgc=h*c*1e8/cwvl	;[erg/ph]
areal=(interpol(areff,arwvl,lwvl)>0)<(max(areff))
areac=(interpol(areff,arwvl,cwvl)>0)<(max(areff))
ltrans=0.*areal & ctrans=0.*areac
if NHcol gt 0 then begin
  ltau=ismtau(lwvl,NH=NHcol,verbose=vv,abund=abund, _extra=e)
  ctau=ismtau(cwvl,NH=NHcol,verbose=vv,abund=abund, _extra=e)
  ltrans=exp(-ltau) & ctrans=exp(-ctau)
endif

;	define the output and the intermediates
fac=0.*logT & lspec=fac & cspec=fac & ctrate=fac & ergfx=ctrate

;	for each T, compute flux and compute counts
for i=0L,ntl-1L do begin		;{main loop

  ll=reform(lemis[i,*])*zab	;[1e23 ergs cm^3/s/line]
  lfx=ll*(EM/1d23)		;[ergs/s/cm^2/line]
  lspec[i]=total(lfx)		;[ergs/s/cm^2]
  if mcs gt 0 then begin	;(if CSTR is OK
    cc=reform(cemis[i,*])		;[1e23 ergs cm^3/s/AA]
    cfx=cc*dcw*(EM/1d23)		;[ergs/s/cm^2/bin]
    cspec[i]=total(cfx)			;[ergs/s/cm^2]
  endif				;MCS>0)

  lsp=lfx*ltrans*areal/nrgl
	;[ergs/s/cm^2/line]*[cm^2]*[ph/erg] ==[ph/s/line]
  lsp0=lfx*ltrans	;[ergs/s/cm^2/line]
  if mcs gt 0 then csp=cfx*ctrans*areac/nrgc
	;[ergs/s/cm^2/bin]*[cm^2]*[ph/erg] ==[ph/s/bin]
  if mcs gt 0 then csp0=cfx*ctrans	;[ergs/s/cm^2/bin]

  ;	estimate without RMF
  ctrate[i]=total(lsp)		;[ph/s]
  if mcs gt 0 then ctrate[i]=ctrate[i]+total(csp)	;[ph/s]
  ergfx[i]=total(lsp0)	;[ergs/s/cm^2]
  if mcs gt 0 then ergfx[i]=ergfx[i]+total(csp0)	;[ergs/s/cm^2]

  if vv ge 100 then print,$
	i,logT[i],lspec[i],cspec[i],ctrate[i],(lspec[i]+cspec[i])/ctrate[i]

  sp=hastogram(lwvl,wgrid,wts=lsp)	;[ph/s/bin]
  if mcs gt 0 then sp=sp+rebinw(csp,cww,wgrid,/perbin)	;[ph/s/bin]

  if ctrate[i] gt 0 then begin		;(bother only if it matters
    if mrs ne 0 then begin		;(if RMF is defined
      emid=12.3985/wmid & os=sort(emid) & emid=emid[os] & spe=sp[os]
      conv_rmf,emid,spe,chan,spec,rmstr

      if n_elements(chnrng) eq 2 then begin
        szc=size(chnrng) & nszc=n_elements(szc)
        if szc[nszc-2] lt 4 or max(chnrng) gt 15 then begin
	  ;	channel range in integers, assumed to be PHA or PI
	  ichan=lindgen(n_elements(spec))+rmstr.FIRSTCHAN
	  oo=where(ichan ge min(chnrng) and chan le max(chnrng),moo)
        endif else begin
	  ;	channel range in float, or max<15: assumed to be in keV
          oo=where(chan ge min(chnrng) and chan le max(chnrng),moo)
	  print,min(chnrng),max(chnrng),oo[0],oo[moo-1]
        endelse
        if moo eq 0 then begin
	  message,'CHNRNG missed out entire spectrum; returning',/informational
	  if vv gt 100 then stop,'HALTing.  type .CON to continue'
	  return,fac
        endif
        ctrate[i]=total(spec[oo])		;[ct/s]
      endif

      ;DEBUG ooo=where(emid ge min(chnrng) and emid le max(chnrng),mooo)
      ;DEBUG print,i,total(spe[ooo]),total(spec[oo]),total(spe[ooo])/total(spec[oo])
      ;DEBUG if total(spe[ooo])/total(spec[oo]) gt 2 then stop

    endif else begin			;yes RMF)(no RMF
      ;	channel range assumed to be in keV
      if n_elements(chnrng) eq 2 then $
	minchan=min(12.3985/chnrng,max=maxchan) else $
	minchan=min(wgrid,max=maxchan)
      oo=where(wgrid ge minchan and wgrid le maxchan,moo)
      ;DEBUG print,minchan,maxchan,moo,oo[0],oo[moo-1]
      if moo eq 0 then begin
	message,'CHNRNG missed out entire spectrum; returning',/informational
	if vv gt 100 then stop,'HALTing.  type .CON to continue'
	return,fac
      endif
      ctrate[i]=total(sp[oo])	;[ph/s]
    endelse				;no RMF)
  endif					;CTRATE>0)

  if vv ge 10 then print,$
	i,logT[i],lspec[i],cspec[i],ctrate[i],(lspec[i]+cspec[i])/ctrate[i]

endfor					;i=0,NTL-1}

;	so now we have lspec+cspec [ergs/s/cm^2] v/s ctrate [ct/s]
spec=lspec
if mcs ne 0 then spec=lspec+cspec
o0=where(ctrate gt 0,mo0)
if mo0 eq 0 then message,'All counts are zero',/informational
fac[o0]=spec[o0]/ctrate[o0]		;[ergs/cm^2/ct]

;	yeah, but how much of that is really reliable?
epsi=1d-6
if n_elements(eps) ne 0 then epsi=double(eps[0])
oy=where(spec lt epsi*max(spec) or ctrate lt epsi*max(ctrate),moy)
if moy ne 0 then begin
  if vv ge 5 then message,$
	strtrim(moy,2)+' points in output are being set to NaN',/informational
  fac[oy]=!values.F_NAN
endif

;	debugging step
if vv gt 200 then stop,'HALTing.  type .CON to continue'

return,fac
end

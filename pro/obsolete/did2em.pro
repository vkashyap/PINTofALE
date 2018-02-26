pro did2em,idstr,edens,fluxes=fluxes,DEM=DEM,logT=logT,ldir=ldir,dwvl=dwvl,$
	NH=NH,EMs=EMs,xtitle=xtitle,ytitle=ytitle,title=title, _extra=e
;+
;procedure	did2em
;	generate emission measure estimates at a variety of electron
;	density values and plot them
;
;syntax
;	did2em,idstr,edens,fluxes=fluxes,DEM=DEM,logT=logT,ldir=ldir,$
;	dwvl=dwvl,NH=NH,EMs=EMs,xtitle=xtitle,ytitle=ytitle,title=title,$
;	chifil=chifil,chidir=chidir,eqfile=eqfile,abund=abund,/noph,$
;	effar=effar,wvlar=wvlar,/ikev,tol=tol,fH2=fH2,He1=He1,HeII=HeII,/Fano
;
;obsolete
;	superceded by DID2EMIS2EM
;
;parameters
;	idstr	[INPUT; required] structure containing ID information
;		(see LINEID for description)
;	edens	[I/O] array of electron densities [cm^-3]
;		* default is [1e8,1e14]
;		* if not specified, insufficiently specified, or otherwise
;		  meaningless, uses default and overwrites input
;
;keywords
;	fluxes	[I/O] observed fluxes
;		* if size does not match the number of components or
;		  the number of IDs, will use IDSTR.().FLUX and overwrite
;	DEM	[INPUT] differential emission measure to use in computing
;		the fluxes
;		* default is to use 1e12 [cm^-5] at each specified LOGT
;	logT	[INPUT] temperatures at which DEM is defined
;		* if not specified, set to 6.5, unless DEM has more steps,
;		  in which case interpolated into the 4..8 range
;	ldir	[INPUT] string array of database directories to use to
;		search for the emissivities
;		* if size does not match either number of IDSTR components
;		  or IDs, only the first element is used
;	dwvl	[INPUT] search the line database for the wavelength of
;		interest with this assumed roundoff error
;		* default=5e-3
;	NH	[INPUT] H column density [cm^-2]
;	EMs	[OUTPUT] array of emission measures EMS(EDEN,WVL) that
;		are plotted.
;		* at any given point, EMS = \int DEM(LOGT)*dLOGT,
;		  unless DEM has one or two elements only, in which
;		  case EMS = SUM(DEM)
;		* what this implies is that the input DEM is scaled
;		  at each density to match the observed flux at each
;		  wavelength and this scaled DEM is represented in the
;		  plot by the summed measure.  This ensures that all
;		  wavelengths considered are compared over the same
;		  temperature range, and the ambiguity regarding the
;		  temperature of peak emissivity, etc. has been removed.
;	xtitle	[INPUT] passed to PLOT
;	ytitle	[INPUT] passed to PLOT
;	title	[INPUT] passed to PLOT
;	_extra	[INPUT] used to specify defined keywords to subroutines:
;		FOLD_IONEQ: CHIFIL
;		RD_IONEQ: CHIDIR, EQFILE
;		LINEFLX: ABUND, NOPH, EFFAR, WVLAR, IKEV
;		ARRAYEQ: TOL
;		ISMTAU: fH2, He1, HeII, FANO
;
;subroutines
;	UPDATID
;	LSD
;	    RD_LINE
;	    FOLD_IONEQ
;	        RD_IONEQ
;	            READ_IONEQ (CHIANTI subroutine)
;	    LINEFLX
;	    ARRAYEQ
;	    KABOOM
;	    INICON
;	ISMTAU
;
;side-effects
;	generates a plot which erases any plot that exists beforehand
;
;history
;	vinay kashyap (Dec98)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	allowed fluxes from multiple IDs to be added up (VK; FebMM)
;	modified behavior of handling multiple-matches to allow proper
;	  subset selection (VK; JulMM)
;	changed call to STR_2_ARR from STR2ARR (VK; AprMMV)
;-

message,'OBSOLETE; use DID2EMIS2EM instead',/informational

;	usage
nid=n_tags(idstr)
if nid eq 0 then begin
  print,'Usage: did2em,idstr,edens,fluxes=fluxes,DEM=DEM,logT=logT,ldir=ldir,$'
  print,'       dwvl=dwvl,NH=NH,EMs=EMs,xtitle=xtitle,ytitle=ytitle,title=title;
  print,'       chifil (FOLD_IONEQ); chidir,eqfile (RD_IONEQ);'
  print,'       abund,noph,effar,wvlar,ikev (LINEFLX); tol (ARRAYEQ);
  print,'       fH2,He1,HeII,Fano (ISMTAU)'
  print,'  generate emission measure estimates at a variety of electron'
  print,'  density values and plots them'
  return
endif

;	verify input
ok='ok' & tagid=tag_names(idstr)
if tagid(0) ne 'WVL' then ok='ID Structure not in right format' else $
 if tagid(0) eq 'WVL_COMMENT' then ok='ID Structure contains no data' else $
  if tagid(0) eq 'Z' then ok='ID structure appears to have been stripped'
if ok ne 'ok' then begin
  message,ok,/info & return
endif
;
;	how many components?  how many IDs?  which wavelengths?
obswvl=idstr.WVL & nwvl=n_elements(obswvl)
wvl=abs(idstr.(1).WVL) & for i=1L,nwvl-1L do wvl=[wvl,abs(idstr.(i+1).WVL)]
mwvl=n_elements(wvl)

;	check other inputs
nden=n_elements(edens)
if nden lt 2 then begin
  edens=[1e8,1e14] & nden=2L
endif

;	check and set keywords
	
strid=updatid(idstr,fluxes)			;fold in new fluxes
	
mD=n_elements(DEM) & mT=n_elements(logT)
if mT eq 0 then begin			;(set LOGT
  if mD le 1 then begin
    if mD eq 0 then DEM=1d12
    logT=6.5 & mT=1 & mD=1
  endif else begin
    mT=mD & logT=findgen(mT)*(4./float(mD-1.))+4.
  endelse
endif					;LOGT)
if mD eq 0 then begin			;(set DEM
  DEM=dblarr(mT)+1d12 & mD=mT
endif					;DEM)
if mD gt 0 then begin			;(convert log(DEM) to DEM
  if max(DEM) lt 100. then DEM=10.D^(DEM)
endif					;log(DEM)->DEM
	
ldbdir=strarr(mwvl)		;set line database directory for each ID
mdir=n_elements(ldir) & szd=size(ldir) & nszd=n_elements(szd)
if mdir gt 0 and szd(nszd-2) ne 7 then begin	;(LDIR is not a string
  mdir=1L
  if fix(ldir(0)) eq 1 then ldir='$SPEX' else $
   if fix(ldir(0)) eq 2 then ldir='$CHIANTI' else $
    ldir='$SCAR'
endif						;decoding LDIR)
if mdir gt 0 then begin			;(set LDIR
  ldbdir(*)=ldir(0)
  if mdir eq mwvl then ldbdir=ldir
  if mdir eq nwvl then begin
    k=0L
    for i=0L,nwvl-1L do begin
      mw=n_elements(strid.(i+1).WVL)
      ldbdir(k:k+mw-1L)=ldir([i])
      k=k+mw
    endfor
  endif
endif					;LDIR)
	
if not keyword_set(dWVL) then dWVL=5e-3
ndW=n_elements(dWVL)
	
if keyword_set(NH) then begin		;(set absorption
  tau=ismtau([wvl],NH=NH, _extra=e)
  oo=where(tau lt 69.,moo)
  etau=fltarr(mwvl) & if moo gt 0 then etau(oo)=exp(-tau(oo))
endif else etau=fltarr(mwvl)+1.		;ETAU)

;	initializations
atom=1 & rom=1 & inicon,atom=atom,roman=rom
atom=['X',atom] & rom=['',rom]				;atomic symbols
EMs=dblarr(nden,mwvl)					;output
labels=strarr(mwvl)					;plot labels

;	figure out the baseline, default EM
if mD gt 2 then begin
  mdT=alog(10.)*median(logT(1:*)-logT)
  tmp=DEM & tmp([0L,mD-1L])=0.5*tmp([0L,mD-1L])
  defEM=total(tmp*mdT)
endif else defEM=total(DEM)

;	and now fill up the output array
k=0L
for i=0L,nwvl-1L do begin			;{for each component
  tmpID=strid.(i+1) & idWVL=tmpID.WVL & idZ=tmpID.Z & idion=tmpID.ION
  idFLX=tmpID.FLUX & mw=n_elements(idWVL)
  ;
  for j=0L,mw-1L do begin			;{for each ID
    deltawvl=dWVL(0) & if ndW gt k then deltawvl=dWVL(k)
    wrange=idWVL(j)+deltawvl*[-1,1]
    elem=atom(idZ(j))+' '+rom(idION(j))
    labels(k)=atom(idZ(j))+rom(idION(j))+' !4k!3'+$
	strtrim(string(abs(idWVL(j)),'(g10.5)'),2)
    ;
    tmp='' & flx=0*edens		;default values
    if elem ne 'X' then tmp=$
      lsd(wrange,flx,elem=elem,edens=edens,DEM=DEM,tlog=logT,$
      dbdir=ldbdir(k),ceiling=1.,floor=0.,ratmax=ratmax,flxmax=flxmax,$
      _extra=e)
    ;
    ntmp=n_elements(tmp)
    if ntmp gt 1 then begin
      message,'Too many emissivities found!',/info
      for ii=0,ntmp-1 do print,strtrim(ii,2),':',tmp(ii),flxmax(ii),ratmax(ii)
      fmx=max(flxmax,jj) & c1=''
      read,prompt='Type indices of correct match (-1 for all) ['+$
	strtrim(jj,2)+'] > ',c1
      if strtrim(c1,2) ne '' then jj=str_2_arr(c1)
      if jj[0] eq -1 then jj=lindgen(ntmp) & njj=n_elements(jj)
      if njj eq 1 then flx=reform(flx(*,jj[0])) else begin
	fflx=0*edens
	for in_e=0,nden-1 do fflx(in_e)=total(flx(in_e,jj))
	flx=fflx
      endelse
    endif
    ;
    oo=where(flx gt 0,moo)
    if moo gt 0 then EMs(oo,k)=defEM*idFLX(j)/flx(oo)
    ;
    k=k+1L
  endfor					;J=0,MW-1}
endfor					;I=0,NWVL-1}

;	plot
;	first figure out plot ranges
oo=where(EMs gt 0,moo)
if moo eq 0 then begin
  message,'No lines found',/info & return
endif
emmin=min(EMs(oo),max=emmax)
if emmax/emmin gt 10 then begin
  ylog=1	;plot Y-axis in log scale
  emmin=emmin/5. & emmax=emmax*5.		;stretch the window
  drop=10.^((alog10(emmax)-alog10(emmin))/mwvl)	;label drops
endif else begin
  ylog=0	;plot Y-axis in normal scale
  emmin=emmin/2. > 0
  emmax=emmax*2.				;stretch the window
  drop=(emmax-emmin)/mwvl			;label drops
endelse
;
denmin=min(edens,max=denmax)
if denmax/denmin gt 100 then begin
  xlog=1	;plot X-axis in log scale
  denmin=denmin/10. & denmax=denmax*10.		;stretch the window
endif else begin
  xlog=0	;plot X-axis in normal scale
  denmin=denmin-10. > 0.
  denmax=denmax+10.				;stretch the window
endelse
;
;	set up the plot
plot,[1],/nodata,xrange=[denmin,denmax],yrange=[emmin,emmax],$
	xlog=xlog,ylog=ylog,xstyle=1,ystyle=1,$
	xtitle=xtitle,ytitle=ytitle,title=title
;
;	and plot
curvcol=fix(200*(!d.n_colors)/256.+1)
lablcol=fix(150*(!d.n_colors)/256.+1)
lablpos=fltarr(mwvl)
for i=0L,mwvl-1L do begin
  y=reform(EMs(*,i)) & oo=where(y gt 0,moo)
  if moo gt 0 then begin
    oplot,edens(oo),y(oo),color=curvcol
    lpos=max(y,ipos) & ok=where(lpos eq lablpos,mok)
    if mok eq 0 then lablpos(i)=lpos else lablpos(i)=lablpos(ok(mok-1L))-drop
    if ipos eq 0 then align=1 else align=0
    xyouts,edens(ipos),lpos,labels(i),align=align,col=lablcol
    if mok ne 0 then oplot,[max(edens),denmax],[max(y),lpos],col=lablcol
  endif
endfor

return
end

pro did2emis2em,idstr,edens,fluxes=fluxes,deblend=deblend,DEM=DEM,logT=logT,$
	ldir=ldir,NH=NH,EMs=EMs,FXs=FXs,xtitle=xtitle,ytitle=ytitle,$
	title=title,verbose=verbose,n_e=n_e, _extra=e
;+
;procedure	did2emis2em
;	for a specified line, read in emissivities at various densities,
;	compute fluxes, and invert to get emission measure estimates at
;	each density.
;
;syntax
;	did2emis2em,idstr,edens,fluxes=fluxes,/deblend,DEM=DEM,logT=logT,$
;	ldir=ldir,NH=NH,EMs=EMs,FXs=FXs,xtitle=xtit,ytitle=ytit,title=tit,$
;	verbose=v,dWVL=dWVL,eps=eps,/incieq,mapping=mapping,chifil=chifil,$
;	chidir=chidir,eqfile=eqfile,abund=abund,/noph,effar=effar,wvlar=wvlar,$
;	/kev,tol=tol,fH2=fH2,He1=He1,HeII=HeII,/Fano,/wam,/bam,/mam,/Zbeda,$
;	/Ibeda,Wbeda=Wbeda,/Lbeku,wform=wform,wstyle=wstyle,ziform=ziform,$
;	/nuthin,/noday,/notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$
;	stathk=stathk,charsize=charsize,align=align
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
;		* if size does not match the number of components
;		  or the number of IDs, will instead use IDSTR.(#).FLUX
;		  and overwrite the input
;		* see UPDATID() for details
;	deblend	[INPUT] if _not_ set (which is the default), squishes multiple
;		IDs of a single feature into one single flux and calculates
;		only a "composite" ID.  on the other hand, if set, leaves
;		the multiple IDs alone.
;	DEM	[INPUT] differential emission measure to use in computing
;		the fluxes
;		* default is to use 1e12 [cm^-5] at each specified LOGT
;	logT	[INPUT] temperatures at which DEM is defined
;		* if not specified, set to 6.5, unless DEM has more elements,
;		  in which case is interpolated into the 4..8 range
;	ldir	[INPUT] string array of database directories to use to
;		search for the emissivities
;		* default is '$CHIANTI'
;		* if size does not match either number of IDSTR components
;		  or IDs, only the first element is used
;	NH	[INPUT] H column density [cm^-2]
;	EMs	[OUTPUT] array of emission measures EMS[EDEN,WVL] that
;		are plotted.
;		* at any given point, EMS = \int DEM[LOGT]*dLOGT,
;		  unless DEM has one or two elements only, in which
;		  case EMS = SUM(DEM)
;		* what this implies is that the input DEM is scaled
;		  at each density to match the observed flux at each
;		  wavelength and this scaled DEM is represented in the
;		  plot by the summed measure.  This ensures that all
;		  wavelengths considered are compared over the same
;		  temperature range, and the ambiguity regarding the
;		  temperature of peak emissivity, etc. has been removed.
;	FXs	[OUTPUT] the predicted fluxes FXS[EDEN,WVL] which are used
;		to scale the EMs above.
;	xtitle	[INPUT] passed to PLOT
;	ytitle	[INPUT] passed to PLOT
;	title	[INPUT] passed to PLOT
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT] pass defined keyword values to subroutines:
;		ID2EMIS2ID: DWVL
;		RD_LIST: EPS, INCIEQ, MAPPING
;		FOLD_IONEQ: CHIFIL, EQFILE
;		RD_IONEQ: CHIDIR
;		SQUISHEM: ABUND
;		LINEFLX: ABUND, NOPH, EFFAR, WVLAR, KEV
;		ARRAYEQ: TOL
;		ISMTAU: fH2, He1, HeII, FANO, WAM, BAM, MAM
;		IDLABEL: ZBEDA, IBEDA, WBEDA, LBEKU, WFORM, WSTYLE, ZIFORM
;		STAMPLE: NUTHIN,NODAY,NOTIME,NOUSER,NOPACK,STACOL,STASIZ,STATHK
;		XYOUTS: CHARSIZE,ALIGN
;	n_e	[IGNORE] here only to trap boxing gloves
;
;restrictions
;	requires IDL5.3+ because of use of STRMATCH, etc.
;	requires subroutines
;	  BAMABS
;	  CAT_LN
;	  FOLD_IONEQ
;	  GETABUND
;	  ID2EMIS2ID
;	  IDLABEL
;	  INICON
;	  ISMTAU
;	  LAT2ARAB
;	  LINEFLX
;	  RDABUND
;	  RD_IONEQ
;	  RD_LINE
;	  RD_LIST
;	  READ_IONEQ
;	  SQUISHEM
;	  SYMB2ZION
;	  SYZE
;	  UPDATID
;	  WHEE
;	  ZION2SYMB
;
;side-effects
;	generates a plot which erases any plot that exists beforehand
;
;history
;	vinay kashyap (JanMMI; based on DID2EM)
;	improved color-scale setting for 24-bit consoles (VK; FebMMI)
;-

;	usage
nid=n_tags(idstr)
if nid eq 0 then begin
  print,'Usage: did2emis2em,idstr,edens,fluxes=fluxes,/deblend,DEM=DEM,$'
  print,'       logT=logT,ldir=ldir,NH=NH,EMs=EMs,FXs=FXs,xtitle=xt,ytitle=yt,$'
  print,'       title=tt,verbose=v, dWVL=dWVL,/incieq,mapping=mapping,$'
  print,'       chifil=chifil,chidir=chidir,eqfile=eqfile,abund=abund,/noph,$'
  print,'       effar=effar,wvlar=wvlar,/kev,tol=tol,fH2=fH2,He1=He1,$'
  print,'       HeII=HeII,/Fano,/wam,/bam,/mam,/Zbeda,/Ibeda,Wbeda=Wbeda,$'
  print,'       /Lbeku,wform=wform,wstyle=wstyle,ziform=ziform,/nuthin,$'
  print,'       /noday,/notime,/nouser,/nopack,stacol=stacol,stasiz=stasiz,$'
  print,'       stathk=stathk,charsize=charsize,align=align'
  print,'  generates emission measure estimates at a variety of electron'
  print,'  density values and plots them'
  return
endif

if float(strmid(!version.RELEASE,0,3)) lt 5.3 then begin
  message,'Requires IDL v5.3 or higher.  Returning',/info
  return
endif

;	verify input
ok='ok' & tagid=tag_names(idstr)
if tagid[0] ne 'WVL' then ok='ID Structure not in right format' else $
 if tagid[0] eq 'WVL_COMMENT' then ok='ID Structure contains no data' else $
  if tagid[0] eq 'Z' then ok='ID structure appears to have been stripped'
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
  message,'Setting density array to default values',/info
  edens=[1e8,1e14] & nden=2L
endif

;	check and set keywords

;	verbosity
v=0 & if keyword_set(verbose) then v=long(verbose[0]) > 1

;	deblend or no?
if not keyword_set(deblend) then begin		;(squish multiple IDs
  if v gt 3 then message,'squishing together multiple IDs into one flux',/info
  nflux=nwvl
  obsflx=fltarr(nwvl) & for i=0L,nwvl-1L do obsflx[i]=total(idstr.(i+1L).FLUX)
endif else begin				;)(leave them deblended
  if v gt 3 then message,'multiple IDs are deblended',/info
  nflux=mwvl
  obsflx=idstr.(1).FLUX & for i=1L,nwvl-1L do obsflx=[obsflx,idstr.(i+1).FLUX]
endelse						;DEBLEND)

;	fluxes
nf=n_elements(fluxes)
if nf eq 0 then fluxes=obsflx
if nf gt 0 and nf ne nflux then begin
  message,'Ignoring input fluxes; using IDSTR.(#).FLUX',/info
  fluxes=obsflx
endif

;fold in new fluxes
strid=updatid(idstr,fluxes,verbose=v)

;DEM and logT
mD=n_elements(DEM) & mT=n_elements(logT)
if mT eq 0 then begin			;(set LOGT
  if mD le 1 then begin
    if mD eq 0 then begin
      DEM=1d12
      if v ge 5 then message,'DEM undefined: setting to default',/info
    endif
    logT=6.5 & mT=1 & mD=1
    if v ge 5 then message,'LOGT undefined; defaulting to '+$
	strtrim(logT,2),/info
  endif else begin
    mT=mD & logT=findgen(mT)*(4./float(mD-1.))+4.
    if v ge 5 then message,'LOGT undefined; interpolating to ['+$
	strtrim(min(logT),2)+','+strtrim(max(logT),2)+']',/info
  endelse
endif					;LOGT)
if mD eq 0 then begin			;(set DEM
  DEM=dblarr(mT)+1d12 & mD=mT
  if v ge 5 then message,'DEM undefined: setting to default',/info
endif					;DEM)
if mD gt 0 then begin			;(convert log(DEM) to DEM
  if max(DEM) lt 100. then begin
    DEM=10.D^(DEM)
    message,'converting input DEM from log10',/info
    unlog=1	;check this keyword later to convert back if needed
  endif
endif					;log(DEM)->DEM

;database directories
ldbdir=strarr(mwvl)+'$CHIANTI'	;set line database directory for each ID
mdir=n_elements(ldir) & szd=size(ldir) & nszd=n_elements(szd)
if mdir gt 0 and szd[nszd-2] ne 7 then begin	;(LDIR is not a string
  mdir=1L
  if fix(ldir[0]) eq 1 then ldir='$SPEX' else $
   if fix(ldir[0]) eq 2 then ldir='$CHIANTI' else $
    ldir='$SCAR'
  if v gt 1 then message,'Setting emissivity database to '+ldir[0],/info
endif						;decoding LDIR)
if mdir gt 0 then begin			;(set LDIR
  ldbdir[*]=ldir[0]
  if mdir eq mwvl then ldbdir=ldir
  if mdir eq nwvl then begin
    k=0L
    for i=0L,nwvl-1L do begin
      mw=n_elements(strid.(i+1).WVL)
      ldbdir[k:k+mw-1L]=ldir[[i]]
      k=k+mw
    endfor
  endif
endif					;LDIR)

;	NH
etau=fltarr(nflux)+1.
if keyword_set(NH) then begin		;(set absorption
  if v ge 5 then message,'Computing ISM absorption',/info
  if keyword_set(deblend) then tau=ismtau(abs([wvl]),NH=NH, _extra=e) else $
    tau=ismtau(abs([obswvl]),NH=NH, _extra=e)
  oo=where(tau lt 69.,moo)
  if moo gt 0 then etau[oo]=exp(-tau[oo])
endif 					;ETAU)

;	initializations
atom=1 & rom=1		;need this for IDLv<5
inicon,atom=atom,roman=rom
atom=['X',atom] & rom=['',rom]				;atomic symbols
EMs=dblarr(nden,nflux)					;output
FXs=fltarr(nden,nflux)					;output
if keyword_set(deblend) then $
  labels=idlabel(idstr,idlabel_idx,/Wbeda,Lbeku=0, _extra=e) else $
  labels=idlabel(squishem(idstr,_extra=e),idlabel_idx,/Wbeda,Lbeku=0,_extra=e)


;	figure out the baseline, default EM
if mD gt 2 then begin
  mdT=alog(10.)*median(logT[1:*]-logT)
  tmp=DEM & tmp[[0L,mD-1L]]=0.5*tmp[[0L,mD-1L]]
  defEM=total(tmp*mdT)
endif else defEM=total(DEM)

;	the basic strategy is to use ID2EMIS2ID to read in
;	new emissivities at each density for each of the lines
;	in the input ID structure, and to compute fluxes in
;	each case.
for i=0L,nden-1L do begin			;{for each density
  n_e=edens[i]
  newid=id2emis2id(strid,ldbdir,verbose=v,n_e=n_e, _extra=e)
  ;
  predflx=fltarr(nflux) & k=0L
  for j=0L,nwvl-1L do begin			;{for all components
    tmpID=newid.(j+1L) & mw=n_elements(tmpID.wvl)
    tmpfx=lineflx(tmpID.EMIS,tmpID.LOGT,tmpID.WVL,tmpID.Z, _extra=e)
    if keyword_set(deblend) then begin		;(deblend
      predflx[k:k+mw-1L]=tmpfx & k=k+mw	
    endif else predflx[j]=total(tmpfx)		;squished)
  endfor					;J=0,NFLUX-1}
  predflx=predflx*etau	;ISM correction
  ;
  oo=where(predflx gt 0,moo)
  if moo gt 0 then EMs[i,oo]=defEM*fluxes[oo]/predflx[oo]
  FXs[i,*]=predflx
endfor						;I=0,NDEN-1}

;	plot
;	first figure out plot ranges
oo=where(EMs gt 0,moo)
if moo eq 0 then begin
  message,'No lines found',/info & stop & return
endif
emmin=min(EMs[oo],max=emmax)
if emmax/emmin gt 10 then begin
  ylog=1	;plot Y-axis in log scale
  emmin=emmin/5. & emmax=emmax*5.		;stretch the window
  drop=10.^((alog10(emmax)-alog10(emmin))/nflux)	;label drops
endif else begin
  ylog=0	;plot Y-axis in normal scale
  emmin=emmin/2. > 0
  emmax=emmax*2.				;stretch the window
  drop=(emmax-emmin)/nflux			;label drops
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
stample, _extra=e
;
;	and plot
dncolors=256. > !D.N_COLORS	;24-bit color screen temporary fix
curvcol=fix(200*(!d.n_colors)/dncolors+1)
lablcol=fix(150*(!d.n_colors)/dncolors+1)
lablpos=fltarr(nflux)
for i=0L,nflux-1L do begin
  y=reform(EMs[*,i]) & oo=where(y gt 0,moo)
  if moo gt 0 then begin
    oplot,edens[oo],y[oo],color=curvcol
    lpos=max(y,ipos) & ok=where(lpos eq lablpos,mok)
    if mok eq 0 then lablpos[i]=lpos else lablpos[i]=lablpos[ok[mok-1L]]-drop
    if ipos eq 0 then align=1 else align=0
    xyouts,edens[ipos],lpos,labels[i],align=align,col=lablcol, _extra=e
    if mok ne 0 then oplot,[max(edens),denmax],[max(y),lpos],col=lablcol
  endif
endfor

if keyword_set(unlog) then DEM=alog10(DEM)

return
end

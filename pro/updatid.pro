function updatid,idstr,fluxes,waves,emis,flxerr=flxerr,lambda=lambda,spec=spec,$
	bkg=bkg,digit=digit,verbose=verbose, _extra=e
;+
;function	updatid
;	update the values in an ID structure and return the updated structure
;
;syntax
;	newid=updatid(idstr,fluxes,waves,emis,flxerr=flxerr,lambda=lambda,$
;	spec=spec,bkg=bkg,/digit,verbose=v, width=width,sigma=sigma,$
;	nsigma=nsigma,sigmay=sigmay,xsize=xsize,ysize=ysize,wid=wid,$
;	dynrng=dynrng,markx=markx,marky=marky,markp=markp,markc=markc,$
;	marks=marks,markl=markl,marko=marko,/legg,/quiet)
;
;parameters
;	idstr	[INPUT; required] ID structure to be updated
;	fluxes	[I/O] if given, updates the IDSTR.(#).FLUX field
;		* size must match number of components or number of IDs
;		* if size does not match, but is given as a scalar or a
;		  single-element vector, then
;		  if +ve, treated as a global multiplicative factor
;		  if -ve, the strongest line is normalized to abs(FLUXES)
;		* if not single element, then gets overwritten on output
;		  using IDSTR.(#).FLUX, so will be different if
;		  DIGIT, LAMBDA, and SPEC are set
;	waves	[INPUT] if given, updates the values of the observed
;		wavelengths in IDSTR.WVL
;		* size must match the number of components
;		* if size does not match, but is given as a scalar or a
;		  single-element vector, then taken to be wavelength shift
;		* to update the wavelengths without messing with the fluxes,
;		  call UPDATID() as
;		  newid=updatid(idstr,1.,newwvl, ...)
;	emis	[INPUT] if given, updates the emissivities in IDSTR.(#).EMIS
;		* to update only the emissivities, call UPDATID() as
;		  newid=updatid(idstr,1.,idstr.WVL,newemis, ...)
;
;keywords
;	flxerr	[I/O] if given, updates the errors in IDSTR by these
;		values.  if not given, or is illegal, no changes are made.
;		* overwritten on output using IDSTR.(#).FLUXERR
;		* generally taken to be 1-sigma errors
;		* size must match that of FLUXES
;	lambda	[INPUT] wavelengths of base spectrum to be used to find
;		the fluxes in case FLUXES is not given
;	spec	[INPUT] spectrum to use to guess at the fluxes
;		* size of SPEC and LAMBDA must match, else both are ignored
;		* basically, at each IDSTR.WVL, find the total counts in a
;		  small region nearby, add up the counts, and that's it.
;		* SPEC and LAMBDA are used if
;		  1: FLUXES is not given, OR
;		  2: DIGIT is set
;	bkg	[INPUT] if size matches that of SPEC, subtract from SPEC
;		prior to finding the fluxes
;		* if single element, assumed to be global background value
;	digit	[INPUT] if set, then for each IDLSTR.WVL lets the user
;		interactively pick the range in which to compute the fluxes
;		* overrides FLUXES
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT] pass defined keywords to subroutines
;		GETLOCMAX: WIDTH,SIGMA,NSIGMA
;		PICKRANGE: SIGMA,XSIZE,YSIZE,WID,DYNRNG,MARKX,MARKY,MARKP,
;		  MARKC,MARKS,MARKL,MARKO,LEGG,QUIET
;
;subroutines
;	GETLOCMAX
;	PICKRANGE
;
;history
;	vinay kashyap (Dec98)
;	added keyword FLXERR (VK; MarMM)
;	added keyword VERBOSE; streamlined FLUXES/FLXERR output;
;	  converted to IDL5 notation (VK; JanMMI)
;	added parameter EMIS (VK; JulMMVIII)
;-

;	usage
nid=n_tags(idstr)
if nid eq 0 then begin
  print,'Usage: newid=updatid(idstr,fluxes,waves,flxerr=flxerr,lambda=lambda,$'
  print,'       spec=spec,bkg=bkg,/digit,verbose=v, width=width,sigma=sigma,$'
  print,'       nsigma=nsigma,sigmay=sigmay,xsize=xsize,ysize=ysize,wid=wid,$'
  print,'       dynrng=dynrng,markx=markx,marky=marky,markp=markp,markc=markc,$'
  print,'       marks=marks,markl=markl,marko=marko,/legg,/quiet)'
  print,'  update IDs obtained using LINEID with new fluxes and wavelengths'
  print,''
  return,lineid()
endif

;	verify input
ok='ok' & tagid=tag_names(idstr)
if tagid[0] ne 'WVL' then ok='ID Structure not in right format' else $
 if tagid[0] eq 'WVL_COMMENT' then ok='ID Structure contains no data' else $
  if tagid[0] eq 'Z' then ok='ID structure appears to have been stripped'
if ok ne 'ok' then begin
  message,ok,/info & return,idstr
endif

;	how many components?
obswvl=idstr.WVL & nwvl=n_elements(obswvl)
if keyword_set(verbose) then print,'there are ',NWVL,' components'

;	how many IDs?
mwvl=0L & for i=0L,nwvl-1L do mwvl=mwvl+n_elements(idstr.(i+1).WVL)
if keyword_set(verbose) then print,'there are a total of ',MWVL,' IDs'

;	check other inputs and see what to do
nf=n_elements(fluxes) & nw=n_elements(waves) & nfe=n_elements(flxerr)
nl=n_elements(lambda) & ns=n_elements(spec) & nb=n_elements(bkg)
nem=n_elements(emis) & szem=size(emis)
;
todo='nada'		;first, do no harm
todof=todo & todow=todo & todos=todo & todofe=todo & todoem=todo
;
if nf eq nwvl or nf eq mwvl then todof='fluxes'
if nf eq 1 and todof eq 'nada' then todof='fluxmult'
if nw eq nwvl then todow='waves'
if nw eq 1 and todow eq 'nada' then todow='wshift'
if nfe eq nf and todof ne 'nada' then todofe='fluxerrs'
if nem gt 0 then begin
  if szem[0] eq 1 then begin
    if mwvl eq 1 then todoem='emissivity' else $
     message,'EMIS must be 2D array; ignoring',/informational
  endif
  if szem[0] eq 2 then begin
    if mwvl eq szem[2] then todoem='emissivity' else $
     message,'EMIS is incompatible with IDSTR; ignoring',/informational
  endif
  if (size(idstr.(1).EMIS))[1] ne szem[1] then begin
    ;	don't touch it if the array sizes are different
    message,'IDSTR.(#).EMIS incompatible with EMIS; ignoring',/informational
    todoem='nada'
  endif
endif
;
if nl gt 2 and nl eq ns then todos='locmax'
if todof ne 'nada' then todos='nada'
;
dig=0 & if keyword_set(digit) then dig=1
if dig eq 1 then begin
  if todos eq 'nada' then begin
    message,'spectrum not defined; ignoring request for interactive use',/info
    dig=0
  endif else begin
    if todof ne 'nada' then message,'Ignoring input fluxes!',/info
    todos='digit'
    todof='nada'
  endelse
endif

;	define the output
strid=idstr

;	take care of the wavelengths -- they are easy
if todow ne 'nada' then begin
  if todow eq 'waves' then obswvl=waves
  if todow eq 'wshift' then obswvl=obswvl+waves[0]
  strid.WVL = obswvl
endif

;	handle the spectrum
if todos ne 'nada' then begin
  fluxes=fltarr(nwvl)
  sp=spec
  if nb eq nl then sp=spec-bkg else $
	if nb eq 1 then sp=spec-bkg[0]		;subtract background

  if todos eq 'locmax' then begin	;(best guess near peak
    ilm=getlocmax(lambda,sp,_extra=e)	;find all the local maxima
    if ilm[0] eq -1 then begin
      message,'no local maxima found!',/info
      ilm=lonarr(nwvl)
      for i=0L,nwvl-1L do begin
        tmp=min(abs(lambda-obswvl[i]),imn) & ilm[i]=imn
      endfor
    endif
    spp=fltarr(nwvl)		;find the nearest peaks
    if not keyword_set(wrange) then wrange=0.3*(max(obswvl)-min(obswvl))
    for i=0L,nwvl-1L do begin
      tmp=min(abs(lambda[ilm]-obswvl[i]),imn) & spp=sp[ilm[imn]]
      tflx=spp & k=imn & spk=sp[k]
      while spk ge 0.1*spp and k ge 0 do begin
	k=k-1L & spk=sp[k]
	tflx=tflx+sp[k]
      endwhile
      k=imn & spk=sp[k]
      while spk ge 0.1*spp and k lt ns do begin
	k=k+1L & spk=sp[k]
	tflx=tflx+sp[k]
      endwhile
      fluxes[i]=tflx
    endfor
  endif					;TODOS='LOCMAX')

  if todos eq 'digit' then begin	;(get the ranges manually
    if not keyword_set(wrange) then wrange=0.3*(max(obswvl)-min(obswvl))
    for i=0L,nwvl-1L do begin
      oxr=[obswvl[i]-wrange,obswvl[i]+wrange]
      xtitle='@'+strtrim(obswvl[i],2)
      tmp=cat_id(idstr,pick=[i])
      ow=pickrange(lambda,sp,oxr=oxr,oyr=oyr,xtitle=xtitle, _extra=e)
      if ow[0] eq -1 then begin
	message,'no range selected.  setting flux to max value in range',/info
	hw=where(lambda ge oxr[0] and lambda le oxr[1],mhw)
	if mhw gt 0 then fluxes[i]=max(sp[hw])
      endif else fluxes[i]=total(sp[ow])
    endfor
  endif					;TODOS='DIGIT'

  todof='fluxes' & nf=n_elements(fluxes)	;overwrite
endif

;	take care of the fluxes if that's all is needed
if todof ne 'nada' then begin
  if todof eq 'fluxmult' then begin
    fmult=fluxes[0]
    if keyword_set(verbose) then message,'multiplying the fluxes by '+$
	strtrim(fmult,2),/info
  endif else fmult=1.
  k=0L
  for i=0L,nwvl-1L do begin
    flx=idstr.(i+1).FLUX & flx=fmult*flx
    flxe=idstr.(i+1).FLUXERR & flxe=fmult*flxe
    mw=n_elements(flx)
    if todof eq 'fluxes' then begin
      if nf eq nwvl then begin
	;	split the given flux in the same ratio as original flux
	tflx=total(flx)
	if tflx ne 0 then rflx=flx/tflx else rflx=(fltarr(mw)+1.)/float(mw)
	flx=rflx*fluxes[i]
	;	split the square of given error in the same ratio as
	;	square of original errors, or if initially not given,
	;	then in the ratio of the fluxes
	if todofe eq 'fluxerrs' then begin
	  tflxe=total(flxe^2)
	  if tflxe ne 0 then rflxe=flxe^2/tflxe else rflxe=rflx
	  flxe=sqrt(rflxe*flxerr[i]^2)
	endif
      endif else begin
	;	overwrite each flux value
	flx = fluxes[k:k+mw-1L]
	if todofe eq 'fluxerrs' then flxe=flxerr[k:k+mw-1L]
      endelse
    endif
    strid.(i+1).FLUX = flx
    if todofe eq 'fluxerrs' then strid.(i+1).FLUXERR = flxe
    k=k+mw
  endfor
endif

;	update the emissivities if needed
if todoem ne 'nada' then begin
  if keyword_set(verbose) then message,'replacing emissivities',$
	/informational
  k=0L
  for i=0L,nwvl-1L do begin
    iwvl=n_elements(idstr.(i+1).WVL)
    tmpem=idstr.(i+1).EMIS
    for j=0L,iwvl-1L do begin
      tmpem[*,j]=emis[*,k] & k=k+1L
    endfor
    strid.(i+1).EMIS = tmpem
  endfor
endif

;	extract fluxes, and fluxerrors and overwrite the I/O keywords
k=0L
if nf eq 0 then begin & fluxes=fltarr(nwvl) & nf=nwvl & endif
if nfe eq 0 then begin & flxerr=fltarr(nwvl) & nfe=nwvl & endif
for i=0L,nwvl-1L do begin
  tmpid=strid.(i+1L) & mw=n_elements(tmpid.WVL)
  if nf eq nwvl then fluxes[i]=total(tmpid.FLUX) else $
   if nf eq mwvl then fluxes[k:k+mw-1L]=tmpid.FLUX else $
    if nf ne 1 then $
     message,'FLUXES: incompatible array sizes; returning zeros',/info
  if nfe eq nwvl then flxerr[i]=sqrt(total((tmpid.FLUXERR)^2)) else $
   if nfe eq mwvl then flxerr[k:k+mw-1L]=tmpid.FLUXERR else $
    message,'FLUXERR: incompatible array sizes; returning zeros',/info
  k=k+mw
endfor

return,strid
end

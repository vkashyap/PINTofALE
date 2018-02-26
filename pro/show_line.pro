pro show_line,wvl,flux,Z=Z,ion=ion,lambda=lambda,spec=spec,order=order,$
	wmark=wmark,fmark=fmark,lmark=lmark,cmark=cmark,markc=markc,$
	xo=xo,yo=yo,oxr=oxr,oyr=oyr, _extra=e
;+
;procedure	show_line
;	plot lines and label them
;
;syntax
;	show_line,wvl,flux,Z=Z,ion=ion,lambda=lambda,spec=spec,order=order,$
;	wmark=wmark,fmark=fmark,lmark=lmark,cmark=cmark,xo=xo,yo=yo,oxr=oxr,$
;	oyr=oyr,sep=sep,/squish,dynrng=dynrng,markp=markp,markc=markc,$
;	marks=marks,marko=marko,/quiet,xtitle=xtitle,ytitle=ytitle,$
;	title=title
;
;parameters
;	wvl	[INPUT; required] wavelengths of lines
;	flux	[INPUT; optional] fluxes of lines
;
;keywords
;	Z	[INPUT] atomic number of element producing WVL
;	ion	[INPUT] ionic state of element producing WVL
;		* Z and ION do not have to match the size of WVL:
;		  as many are used as available or necessary
;	spec	[INPUT] spectrum over which to plot the lines
;	lambda	[INPUT] as in SPEC(LAMBDA)
;		* SPEC and LAMBDA should match in size, or else they
;		  are both ignored
;	order	[INPUT] string of comma-separated integers denoting
;		grating orders to include in plot
;		* default is to use just the first order
;		* anything meaningless (including "0") defaults to first order
;	wmark	[INPUT] additional locations to mark (the output of SMUDGE,
;		for example)
;	fmark	[INPUT] fluxes at which to mark WMARKs
;	lmark	[INPUT] labels for WMARKs
;		* FMARK and LMARK are ignored unless they match WMARK
;	cmark	[INPUT] colors for LMARKs
;	markc	[INPUT] colors for WVLs
;		* default is 222
;	xo	[OUTPUT] output from PICKRANGE
;	yo	[OUTPUT] output from PICKRANGE
;	oxr	[OUTPUT] output from PICKRANGE
;	oyr	[OUTPUT] output from PICKRANGE
;	_extra	[INPUT] pass defined keywords to subroutines:-
;		STR_2_ARR: SEP,SQUISH
;		PICKRANGE: DYNRNG,MARKP,MARKS,MARKO,QUIET,XTITLE,YTITLE,TITLE
;
;subroutines
;	STR_2_ARR
;	INICON
;	PICKRANGE
;
;history
;	vinay kashyap (Nov98)
;	added keywords WMARK, FMARK, LMARK (VK; Dec98)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	added keywords MARKC and CMARK (VK; JanMMI)
;	bug fix: if FLUX not set, 1st 2 were being marked at [0,1] (VK; Feb01)
;	changed call to STR_2_ARR from STR2ARR (VK; AprMMV)
;-

;	usage
nw=n_elements(wvl) & nf=n_elements(flux)
if nw eq 0 then begin
  print,'Usage: show_line,wvl,flux,Z=Z,ion=ion,lambda=lambda,spec=spec,$'
  print,'       order=order,wmark=wmark,fmark=fmark,lmark=lmark,cmark=cmark,$'
  print,'       markc=markc,xo=xo,yo=yo,oxr=oxr,oyr=oyr,sep=sep,/squish,$'
  print,'       dynrng=dynrng,markp=markp,marks=marks,marko=marko,/quiet,$'
  print,'       xtitle=xtitle,ytitle=ytitle,title=title'
  print,'  plot lines and label them'
  return
endif

;	cast inputs in acceptable form
if nw eq 1 then ww=[abs(wvl)] else ww=abs(wvl(*))
if nf ne 0 then begin
  ff=0.*ww+flux(0)
  if nf gt 1 and nf le nw then ff(0L:nf-1L)=flux(*)
  if nf gt nw then ff(*)=flux(0L:nw-1L)
endif else ff=[0.,1.]

;	check keywords
nZ=n_elements(Z) & nI=n_elements(ion)
zz=intarr(nw) & jon=zz
if nZ gt 0 then begin
  if nZ gt 1 and nZ le nw then zz(0L:nZ-1L)=Z(*)
  if nZ gt nw then zz(*)=Z(0L:nw-1L)
endif
if nI gt 0 then begin
  if nI gt 1 and nI le nw then jon(0L:nI-1L)=ion(*)
  if nI gt nw then jon(*)=ion(0L:nw-1L)
endif
;
nx=n_elements(lambda) & ny=n_elements(spec)
ok='ok'
if nx lt 2 then ok='LAMBDA missing' else $
 if ny lt 2 then ok='SPEC(LAMBDA) missing' else $
  if nx ne ny then ok='SPEC and LAMBDA incompatible'
if ok eq 'ok' then begin
  x=[lambda(*)] & y=[spec(*)]
endif else begin
  minx=min(ww,max=maxx) & minf=min(ff,max=maxf)
  delx=maxx-minx
  x=[minx-0.1*delx,maxx+0.1*delx] & y=[0.5*minf,2*maxf]
endelse
;
if keyword_set(order) then gord=str_2_arr(order, _extra=e) else gord=0
if gord(0) eq 0 then gord=[1] & ngo=n_elements(gord)

;	initialize
atom=1 & rom=1 & inicon,atom=atom,roman=rom & atom=['',atom] & rom=['',rom]
labl=strarr(nw)
for i=0L,nw-1L do labl(i)='  '+atom(zz(i))+' '+rom(jon(i))+' !4k!3'+$
	strtrim(string(ww(i),'(g10.6)'),2)
markol=intarr(nw) & mc=n_elements(mc)
if mc eq 0 then markol(*)=222 else $
 if mc eq 1 then markol(*)=markc(0) else $
  if mc eq nw then markol=markc else $
   if mc eq ngo then begin
     markol=intarr(ngo*nw)
     for i=0L,ngo-1L do markol(i*nw:(i+1L)*nw-1L)=markc(i)
   endif

;	account for multiple grating orders
for i=0L,ngo-1L do begin
  if gord(i) le 0 then begin
    wg=abs(ww)
    fg=ff
    lg=labl
  endif else begin
    wg=abs(ww)*gord(i)
    fg=ff/gord(i)^2
    ;if nf ne 0 then fg=ff/gord(i)^2 else fg=ff
    if gord(i) gt 1 then lg=labl+' O='+strtrim(gord(i),2) else lg=labl
  endelse
  if i eq 0 then begin
    gw=wg
    if nf ne 0 then gf=fg
    gl=lg
    gc=markol
  endif else begin
    gw=[gw,wg]
    if nf ne 0 then gf=[gf,fg]
    gl=[gl,lg]
    gc=[gc,markol]
  endelse
endfor

;	tack on extra labels, if any
nmw=n_elements(wmark)
nmf=n_elements(fmark) & nml=n_elements(lmark) & nmc=n_elements(cmark)
if nmw gt 0 then begin
  ffmark=fltarr(nmw) & llmark=strarr(nmw) & ccmark=intarr(nmw)
  if nmf gt 0 then for i=0L,nmw-1L do ffmark(i)=fmark([i])
  if nml gt 0 then for i=0L,nmw-1L do llmark(i)=lmark([i])
  if nmc gt 0 then for i=0L,nmw-1L do ccmark(i)=cmark([i])
  gw=[wmark,gw]
  if nf ne 0 then gf=[ffmark,gf] else gf=[ffmark]
  gl=[llmark,gl]
  gc=[ccmark,gc]
endif

;	plot via PICKRANGE
oo=pickrange(x,y,markx=gw,marky=gf,markl=gl,markc=gc,/legg,$
	xo=xo,yo=yo,oxr=oxr,oyr=oyr, _extra=e)

return
end

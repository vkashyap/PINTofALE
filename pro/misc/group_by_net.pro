pro group_by_net,srcx,bkgx,outsp,outx,outdx,backsc=backsc,netthr=netthr,$
	verbose=verbose, _extra=e
;+
;procedure	group_by_net
;	groups photons until the background-subtracted value exceeds a threshold
;
;syntax
;	group_by_net,srcx,bkgx,outsp,outx,outdx,backsc=backsc,netthr=netthr,verbose=verbose
;
;parameters (all are required)
;	srcx	[INPUT] array of source photon attributes (wavelength or time) that must be grouped
;	bkgx	[INPUT] array of background photon attributes (wavelength or time) that must be grouped
;	outsp	[OUTPUT] net rate spectrum or lightcurve that results from the grouping
;	outx	[OUTPUT] bin boundaries for which OUTSP is computed
;	outdx	[OUTPUT] bin widths for each grouping; OUTSP*OUTDX should look flat
;
;keywords
;	backsc	[INPUT] backscale factor, the ratio of the geometric*effective areas of source to background
;		* if not given, assumed to be 1
;		* if set to 0, BKGX is ignored
;	netthr	[INPUT] threshold value for the net accumulated counts to stop the grouping
;		* default is 10
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2015dec)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(srcx) & nb=n_elements(bkgx)
if np lt 5 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SRCX is undefined' else $
  if nb eq 0 then ok='BKGX is undefined' else $
   if ns eq 1 then ok='nothing to group, is there?'
if ok ne 'ok' then begin
  print,'Usage: group_by_net,srcx,bkgx,outsp,outx,outdx,backsc=backsc,netthr=netthr,verbose=verbose'
  print,'  group event lists to make spectra or light curves with set net counts'
  if np ne 0 then message,ok,/informational
  return
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
bkgscal=1.0 & if n_elements(backsc) gt 0 then bkgscal=float(backsc[0])	;yeah, could be -ve or 0
bfac=0.0 & if bkgscal ne 0 then bfac=1./bkgscal
thrnet=10. & if keyword_set(netthr) then thrnet=float(netthr[0])	;yeah, could be -ve

;	sort the inputs
os=reverse(sort(srcx)) & sx=srcx[os]
ob=reverse(sort(bkgx)) & bx=bkgx[ob]

;	define the outputs
xbound=mid2bound(sx,/halfit) & xdel=xbound[1:*]-xbound
outx=xbound & outdx=xdel & outsp=fltarr(ns)-1.	;the -ves will be stripped away later

;	step through the events
i0=0L & netct=1L & kg=0L & ibank=0
;outsp[kg]=1./outdx[0] & outx[kg]=sx[0] & outx[kg+1L]=sx[1] & outdx[kg]=sx[1]-sx[0]
;x0=sx[i0]
for i=0L,ns-1L do begin
  if ibank eq 0 then begin
    x0=xbound[i0] & outx[kg]=x0
    x1=xbound[i+1L] & outx[kg+1L]=x1
    delx=abs(x1-x0) & if delx eq 0 then delx=xdel[i] & outdx[kg]=delx
  endif
  if bfac ne 0 then ob=where(bx ge x0 and bx le x1,mob) else mob=0L
  netct=float(i+1L-i0)-mob*bfac
  if delx gt 0 then outsp[kg]=netct/delx
  if netct ge thrnet then ibank=1
  if ibank eq 1 then begin
    if vv gt 100 then print,i0,i,kg,x0,x1,netct,mob,outsp[kg]
    kg=kg+1L & i0=i+1L
    ibank=0
    if vv gt 5 then kilroy
  endif
endfor

ok=where(outsp ge 0,mok)
if mok gt 0 then begin
  outsp=outsp[ok] & outdx=outdx[ok] & outx=outx[[ok,mok]]
endif else begin
  outsp=-1L & outdx=-1L & outx=-1L
endelse

;	debug
if vv gt 10 then begin
  plot,(outx[1:*]+outx)/2.,outsp,thick=2
  xmin=min(outx,max=xmax) & dx=(xmax-xmin)/21L
  hs=histogram(sx,min=xmin,max=xmax,bin=dx)
  hb=histogram(bx,min=xmin,max=xmax,bin=dx)
  xx=findgen(n_elements(hs))*dx+xmin
  oplot,xx,(hs-hb*bfac)/dx,line=2
  oplot,xx,hb*bfac/dx,line=1
endif
if vv gt 1000 then stop,'HALTing; type .CON to continue'

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

srcx=randomn(seed,150)
bkgx=randomu(seed,150)*10-5.
backsc=10.
if not keyword_set(netthr) then netthr=10.
if not keyword_set(verbose) then verbose=10

group_by_net,srcx,bkgx,outsp,outx,outdx,backsc=backsc,netthr=netthr,verbose=verbose

end

function haarline,y,scales,xerr=xerr,thrct=thrct,sclout=sclout,$
	sclmin=sclmin,sclmax=sclmax,verbose=verbose,hy=hy,wy=wy,$
	_extra=e
;+
;function	haarline
;	return the best guesses for the locations of lines in a
;	1D spectrum, in units of the position index, using Haar
;	wavelet transforms (hence the name)
;
;	gets the filtered transform at multiple scales, and starting
;	from the smallest scale, determines line locations by averaging
;	over contiguous regions weighted by the wavelet transform.
;	at succeedingly larger scales, the value of the computed line
;	location is updated if the max of the transform is greater
;	than before, and ignored if (a) if it is not, and (b) if the
;	support of the line at this scale contains more than one line
;	at smaller scales.
;
;syntax
;	xpos=haarline(y,scales,xerr=xerr,thrct=thrct,sclout=sclout,$
;	sclmin=sclmin,sclmax=sclmax,ysig=ysig,thrsig=thrsig,$
;	thrloc=thrloc,verbose=verbose,hy=hy,wy=wy)
;
;parameters
;	y	[INPUT; required] 
;	scales	[I/O] if given on input as an integer array, computes
;		the Haar transforms at these scales.
;		* if illegal in any way, is calculated using SCLMIN and SCLMAX
;
;keywords
;	xerr	[OUTPUT] the 1-sigma errors on the line locations, also
;		in units of array indices
;	thrct	[INPUT; default=10] minimum signal in the reconstructed
;		function before calling it a line
;	sclout	[OUTPUT] the scale at which the line is detected
;	sclmin	[INPUT; default=4, hardcoded minimum=1] lowest scale in
;		wavelet transform
;	sclmax	[INPUT; default=32, hardcoded maximum=N(Y)/6] highest scale
;		in wavelet transform
;		* scales go from SCLMIN to SCLMAX in powers of 2
;		* if SCLMIN and SCLMAX are not integers, the ceiling is used
;		* if SCLMIN and SCLMAX are not powers of 2, the nearest
;		  power of 2 smaller than them are used
;	verbose	[INPUT] controls chatter
;	wy	[OUTPUT] the wavelet correlation coefficients, output of HAARTRAN()
;	hy	[OUTPUT] the filtered coefficients, output of HAARTRAN()
;	_extra	[INPUT ONLY] pass defined keywords to subroutine HAARTRAN():
;		-- YSIG 	default is Poisson for counts, 0 otherwise
;		-- THRSIG 	default is 1.0
;		-- THRLOC 	set by default, explicitly set to 0 if unwanted
;
;subroutines
;	HAARTRAN
;	KILROY
;	IS_KEYWORD_SET
;
;history
;	vinay kashyap (Dec'02)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	some edge case bug fixes (VK; 15Nov2014)
;-

;	usage
ok='ok' & np=n_params()
ny=n_elements(y) & ns=n_elements(scales) & styp=size(scales,/tname)
if np eq 0 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='Y is undefined'
if ok ne 'ok' then begin
  print,'Usage: xpos=haarline(y,scales,xerr=xerr,thrct=thrct,sclout=sclout,$'
  print,'       sclmin=sclmin,sclmax=sclmax,ysig=ysig,thrsig=thrsig,$'
  print,'       thrloc=thrloc,verbose=verbose,wy=wy,hy=hy)'
  print,'  return line location indices'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	scales
sok='ok'
if ns eq 0 then sok='SCALES is undefined' else $
 if scales[0] eq 0 then sok='SCALES must be recomputed' else $
  if styp ne 'INT' and styp ne 'LONG' then sok='SCALES must be an integer array'

;	keywords
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
;
cthr=10. & if is_keyword_set(thrct) then cthr=thrct[0]
;
minscl=4 & if keyword_set(sclmin) then minscl=long(sclmin[0]+0.5)>1
maxscl=32 & if keyword_set(sclmax) then maxscl=long(sclmax[0]+0.5)<(long(ny/6))
minscl=2L^(long(alog(minscl)/alog(2)))
maxscl=2L^(long(alog(maxscl)/alog(2)))
if sok ne 'ok' then begin
  if vv gt 2 then begin
    message,sok+': using SCLMIN/SCLMAX',/informational
    if vv gt 5 then begin
      message,'SCLMIN = '+strtrim(minscl,2),/informational
      message,'SCLMAX = '+strtrim(maxscl,2),/informational
    endif
  endif
  nscl=(long(alog(maxscl)/alog(2))-long(alog(minscl)/alog(2)))+1L
  scl=2L^(lindgen(nscl))*minscl
endif else begin
  scl=[long(scales)] & ok=where(scales ge 1 and scales le long(ny/6),mok)
  if mok eq 0 then begin
    message,'all SCALES are beyond allowed range',/informational
    return,[-1L]
  endif
  scl=scales[ok]
endelse
ns=n_elements(scl)

;	get the filtered transform
hy=haartran(y,scl,scl2,/thrloc,verbose=vv,wy=wy, _extra=e)
scales=scl2	;overwrite for output

;	define intermediate arrays
py=lonarr(ns,ny)	;to hold the percolated indices
cy=fltarr(ns,ny)+0*hy[0]	;for the largest corval
				;the "0*hy[0]" to get correct array type
iy=fltarr(ns,ny) & iye=iy	;to hold the line location index
xx=lindgen(ny)+1L		;the x-index

;	step through the transforms from small to large scales
os=sort(scl)
for i=0L,ns-1L do begin			;{for each scale
  j=os[i]
  hyj=reform(hy[j+1L,*]) & pyj=lonarr(ny) & jy=pyj & cyj=fltarr(ny)+0*hy[0]
  oo=where(hyj gt 0,moo)
  if moo gt 0 then begin
    jy[oo]=1
    for k=1L,moo-1L do begin
      if oo[k]-oo[k-1L] gt 1 then jy[oo[k:*]]=jy[oo[k:*]]+1
    endfor
    py[j,*]=jy
  endif
  ;
  for k=1L,max(jy) do begin
    ok=where(py[j,*] eq k,mok) & if mok eq 0 then message,'BUG!'
    x=total(abs(hy[j+1L,ok])*xx[ok],/nan)/total(abs(hy[j+1L,ok]),/nan)
    xe=sqrt(total(abs(hy[j+1L,ok])*xx[ok]^2,/nan)/total(abs(hy[j+1L,ok]),/nan)-x^2)
    if finite(xe) eq 0 then xe=0.
    iy[j,ok]=x & iye[j,ok]=xe & cy[j,ok]=max(hy[j+1L,ok])
  endfor
endfor					;I=0,NS-1}

;	and now step through the scales to identify the lines
xpos=0 & xerr=0
for i=0L,ns-1L do begin			;{for each subsequent scale
  j=os[i]
  iyj=reform(iy[j,*]) & iyje=reform(iye[j,*]) & cyj=reform(cy[j,*])
  o0=where(iyj gt 0,mo0)
  if mo0 gt 0 then begin		;(there exists something look at

    if not keyword_set(xerr) then begin		;(first pass
      xpos=reform(iy[j,o0]) & xerr=reform(iye[j,o0])
      cx=reform(cy[j,o0])
      ipos=uniq(xpos,sort(xpos))
      xpos=xpos[ipos] & xerr=xerr[ipos]
      sclout=0*xpos+scl[j]
      cortot=cx[ipos]
    endif else begin				;first)(subsequent passes
      xpos2=reform(iy[j,o0]) & xerr2=reform(iye[j,o0])
      cx2=reform(cy[j,o0])
      ipos2=uniq(xpos2,sort(xpos2))
      xpos2=xpos2[ipos2] & xerr2=xerr2[ipos2]
      nx2=n_elements(xpos2)
      ;
      for k=0L,nx2-1L do begin		;{for each line at this scale
        ok=where(abs(iyj-xpos2[k]) lt 1e-3,mok)
	if mok eq 0 then message,'BUG!'
        xtmp=reform(iy[j-1L,ok]) & ctmp=reform(cy[j-1L,ok])
        oi=where(xtmp gt 0,moi)
        if moi eq 0 then begin		;(new line
	  xpos=[xpos,iyj[ok[0]]]
	  xerr=[xerr,iyje[ok[0]]]
	  cx=[cx,cyj[ok[0]]]
	  sclout=[sclout,scl[j]]
	  cortot=[cortot,cy[ok[0]]]
        endif else begin			;MOI=0)(updating line
          xxtmp=xtmp[oi] & xxtmp=xxtmp[uniq(xxtmp,sort(xxtmp))]
          nxxtmp=n_elements(xxtmp)
          if nxxtmp eq 1 then begin		;(update
	    jnk=min(abs(xxtmp[0]-xpos),imatch)
	    if cyj[ok[0]] gt cx[imatch] then begin
	      xpos[imatch]=iyj[ok[0]]
	      xerr[imatch]=iyje[ok[0]]
	      cx[imatch]=cyj[ok[0]]
	      sclout[imatch]=scl[j]
	    endif
	    cortot[imatch]=cortot[imatch]+cy[ok[0]]
          endif 				;NXXTMP=1)(ignore otherwise)
        endelse				;MOI>0)
      endfor				;K=0,NX2-1}
    endelse					;subsequent passes)

  endif					;MO0>0)

endfor					;I=0,NS-1}

;	impose signal thresholding
mok=n_elements(xpos) & ok=lindgen(mok) & hy0=reform(hy[0,*])
if mok eq 1 and xpos[0] eq 0 then begin
  message,'No lines found',/informational
  return,-1L
endif
for i=0L,mok-1L do begin
  i0=long(xpos[i]-xerr[i]+0.5)-1L & i1=long(xpos[i]+xerr[i]+0.5)-1L
  if max(hy0[i0:i1]) lt cthr then ok[i]=-1L
endfor
oy=where(ok ge 0,moy)
if moy gt 0 then begin
  xpos=xpos[oy] & xerr=xerr[os] & sclout=sclout[os]
endif else begin
  message,'No lines found',/informational
  return,-1L
endelse

;	the outputs
os=sort(xpos) & xpos=xpos[os] & xerr=xerr[os] & sclout=sclout[os]

return,xpos
end

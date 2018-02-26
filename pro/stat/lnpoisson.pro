function lnpoisson,d,m,chilike=chilike,verbose=verbose, _extra=e
;+
;function	lnpoisson
;	returns the natural log of the poisson likelihood of observing
;	the specified counts given the model source intensity,
;		p(D|M) = M^D * exp(-M) / D!
;
;syntax
;	lnp=lnpoisson(d,m,/chilike,verbose=verbose)
;
;parameters
;	d	[INPUT; required] data counts
;	m	[INPUT; required] model intensity
;		* output will be array of size [N(D),N(M)]
;		  unless CHILIKE is set, in which case output will be N(M)
;		  (beware that IDL automatically compresses trailing
;		  dimensions if any of them are equal to 1)
;
;keywords
;	chilike	[INPUT] if set, returns -2*SUM_{D}(ln(p(D|M)))
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to avoid crashing the program
;
;history
;	vinay kashyap (Aug01)
;	bug correction: m^d was not translated correctly (VK; Mar02)
;-

;	usage
ok='ok' & np=n_params() & nd=n_elements(D) & nm=n_elements(M)
if np lt 2 then ok='Insufficient parameters' else $
 if nd eq 0 then ok='D not defined' else $
  if nm eq 0 then ok='M not defined'
if ok ne 'ok' then begin
  print,'Usage: lnp=lnpoisson(d,m,/chilike,verbose=verbose)'
  print,'  returns log of Poisson likelihood, p(D|M)'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	keywords
chi=1. & if keyword_set(chilike) then chi=-2.
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1

;	compute
lnp=fltarr(nd,nm) & dd=[d[*]]
o0=where(dd lt 0,mo0) & ok=where(dd ge 0,mok)
for i=0L,nm-1L do begin		;{for each model intensity
  mm=m[i]
  ;
  if mm lt 0 then begin						;(M<0
    if vv gt 0 then message,'Model intensity -ve; doing nothing',/info
    lnp[*,i]=0.
  endif else begin						;M<0)(M.GE.0
    if mm eq 0 then lnp[*,i]=-(!VALUES.F_INFINITY) else begin	;(M>0
      if mo0 gt 0 then begin				;(D<0
        if vv gt 0 then message,$
		'observed counts cannot be -ve; doing nothing',/info
        lnp[o0,i]=0.
      endif else begin					;D<0)(D.GE.0
        ;
        tmp=0.*dd
        if mok gt 0 then begin				;(D.GE.0
          tmp[ok]=-mm+dd[ok]*alog(mm+0.)-lngamma(dd[ok]+1.)
          lnp[ok,i]=tmp[ok]
        endif						;D.GE.0)
      endelse						;D>0)
    endelse							;M>0)
  endelse							;M.GE.0)
endfor				;I=0,NM-1}

;	the "chi-square" allegory
if keyword_set(chilike) then begin
  lnp=chi*lnp
  if nd gt 1 then begin
    tmp=fltarr(nm) & for i=0L,nm-1L do tmp[i]=total(lnp[*,i])
    lnp=tmp
  endif
endif

return,lnp

end

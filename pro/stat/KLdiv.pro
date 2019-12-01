function KLdiv,p1,q2,zeros=zeros,eps=eps,Hp=Hp,Hq=Hq,useln=useln,$
	verbose=verbose, _extra=e
;+
;function	KLdiv
;	compute and return the Kullback-Leibler Divergence of the
;	probability distribution q from the probability distribution p
;
;	K-L divergence denotes the excess information needed to describe
;	the distribution p if q is used instead.  It is defined as
;		D_KL(p||q) = sum p log(q/p)
;	Here we will assume log is in base 2 unless specified otherwise.
;
;syntax
;	DKLpq=KLdiv(p1,q2,zeros=zeros,eps=eps,Hp=Hp,Hq=Hq,/useln,verbose=verbose)
;
;parameters
;	p1	[INPUT; required] probability distribution from which D_KL is computed
;	q2	[INPUT; required] probability distribution for which D_KL is computed
;		* NOTE: it is assumed here that
;		  (1) both P1 and Q2 are +ve and normalized, i.e., integrate/sum to 1
;		  (2) are both defined over the exact same range and binning,
;		      and thus are of the same vector length
;
;keywords
;	zeros	[INPUT; default=0] how to deal with zeros in either of P1 or Q2
;		* The default (=0, or not set) is to ignore all bins where
;		  either P1 or Q2 is zero.
;		* if set, then all 0's are replaced by EPS, and all non-zeros are
;		  decremented by EPS*(number of zeros)/(number of non-zeros)
;	eps	[INPUT] a small number, to be used if ZEROS is set
;		* default is to pick the smallest non-zero value of {P1,Q2}
;		  (if everything is 0, then uses 1e-6)
;		* if set to a -ve value, converts to abs before using
;		* if >1, converts to 1/EPS before using
;	Hp	[OUTPUT] the entropy of P1, computed as -(sum P1*log(P1))
;		* This is useful as a way to figure out how significant or
;		  relevant a computed DKL is.  The KL divergence is the
;		  extra entropy needed relative to Hp when logs are computed
;		  in base 2.  So one can check the fractional increase in
;		  entropy relative to the baseline to judge whether it is large
;		  or not.  See https://stats.stackexchange.com/a/111523 for a
;		  good description.
;		* If USELN is set, does NOT divide DKL and Hp by alog(2) before
;		  assessing entropy changes
;	Hq	[OUTPUT] the entropy of Q2, computed as -(sum Q2*log(Q2)),
;		provided here strictly for completeness
;	useln	[INPUT] if set, do the calculations using natural log
;		* default is to use base 2
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run KLdiv
;	(pops up windows to make plots if !D.NAME='X')
;	(uses PEASECOLR)
;
;history
;	Vinay Kashyap (2019feb)
;	changed name of physical file to all lowercase because for some reason
;	  IDL 8.7.2 does not know about mixed case (VK; 2019nov)
;-

;	usage
ok='ok' & np=n_params() & np1=n_elements(p1) & nq2=n_elements(q2)
if np lt 2 then ok='Insufficient parameters' else $
 if np1 eq 0 then ok='P1 is not defined' else $
  if nq2 eq 0 then ok='Q2 is not defined' else $
   if np1 ne nq2 then ok='P1 and Q2 must be of the same size'
if ok ne 'ok' then begin
  print,'Usage: DKLpq=KLdiv(p1,q2,zeros=zeros,eps=eps,Hp=Hp,Hq=Hq,/useln,verbose=verbose)'
  print,'       compute and return Kullback-Leibler divergence D_KL(p1||q2) of q2 from p1'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
pp=p1 & qq=q2
ok=where(pp lt 0,mok)
if mok gt 0 then begin
  message,'P1 cannot have -ve values; quitting',/informational
  return,-1L
endif
ok=where(qq lt 0,mok)
if mok gt 0 then begin
  message,'Q2 cannot have -ve values; quitting',/informational
  return,-1L
endif
o0p=where(pp eq 0,mo0p,complement=o1p,ncomplement=mo1p)
o0q=where(qq eq 0,mo0q,complement=o1q,ncomplement=mo1q)


;	define keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
alog2=alog(2) & if keyword_set(useln) then alog2=1.
;
if mo1p gt 0 then eps0p=min(pp[o1p]) else eps0p=1d-6
if mo1q gt 0 then eps0q=min(qq[o1q]) else eps0q=1d-6
eps0=eps0p<eps0q
;
if keyword_set(zeros) then begin	;(ZEROS is set, so replace zeros with EPS
  if keyword_set(eps) then begin
    eps1=abs(eps[0])
    if eps1 gt 1 then eps1=1.D/eps1
    if eps1 gt 0 then eps0=eps1
  endif
  if mo0p gt 0 then begin
    pp[o0p]=eps0
    pp=pp*total(p1)/total(pp)
    ;if mo1p gt 0 then pp[o1p]=p1[o1p]-eps0*float(mo1p)/float(mo0p)
  endif
  if mo0q gt 0 then begin
    qq[o0q]=eps0
    qp=qq*total(q2)/total(qq)
    ;if mo1q gt 0 then qq[o1q]=q2[o1q]-eps0*float(mo1q)/float(mo0q)
  endif
endif else begin			;ZEROS)(default is to use non-zeros intersection
  oy=where(pp eq 0 or qq eq 0,moy,complement=ok,ncomplement=mok)
  if mok eq 0 then begin	;(no bins present with both P1 and Q2 non-zero
    if vv gt 0 then message,'no intersection between P1 and Q2',/informational
    if mo1p gt 0 then Hp=total(pp[o1p]*alog(pp[o1p])/alog2) else Hp=!values.F_INFINITY
    if mo1q gt 0 then Hq=total(qq[o1q]*alog(qq[o1q])/alog2) else Hq=!values.F_INFINITY
    return,!values.D_INFINITY
  endif				;MOK=0)
  pp=pp[ok] & qq=qq[ok]
endelse					;ZEROS is not set)

;	compute entropies
Hp=-total(pp*alog(pp)/alog2)
Hq=-total(qq*alog(qq)/alog2)

;	compute D_KL(p||q)
DKLpq=total(pp*(alog(pp)-alog(qq))/alog2)

if vv gt 1000 then stop,'halting; type .CON to continue'

return,DKLpq
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	examples
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	calling sequence
DKLpq=KLdiv()

;	example dataset from https://stats.stackexchange.com/a/111523
p1=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.64]
q2=[0.002,0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.928]
DKLpq=KLdiv(p1,q2,Hp=Hp)

print,''
print,'P1=['+strjoin(strtrim(p1,2),',')+']'
print,'Q2=['+strjoin(strtrim(q2,2),',')+']'
print,'D_KL(p||q)='+strtrim(DKLpq,2)
print,'cf. H(p)='+strtrim(Hp,2)+' for an increase in entropy of '+string(100*DKLpq/Hp,'(i2)')+'%'
if !d.name eq 'X' then begin
  window,0,xsize=900,ysize=900,title='KLdiv() example'
  plot,q2,psym=10,/ylog,line=2,xstyle=5,$
  	title='KLdiv()',subtitle='P1 (solid), Q2 (dashed)',$
	thick=2,xthick=2,ythick=2,charthick=2,charsize=2
  oplot,p1,psym=10,thick=3
  xyouts,0.5,10.^(!y.crange[1]-0.1*(!y.crange[1]-!y.crange[0])),$
  	'D!dKL!n(P||Q)='+strtrim(DKLpq,2),charsize=2.5,charthick=2
  xyouts,0.5,10.^(!y.crange[1]-0.15*(!y.crange[1]-!y.crange[0])),$
  	'H(P)='+strtrim(Hp,2),charsize=2.5,charthick=2
endif

;	how does D_KL/Hp go for Gaussians?
npos=201L & pos=findgen(npos)*0.01-1.0   	;position ranges from -1 to +1
nsig=171L & sig=1.+findgen(nsig)*0.01-0.7	;width ranges from 0.3 to 2
DKLarr=dblarr(npos,nsig) & Hparr=dblarr(npos,nsig)
xx=findgen(601)*0.01-3. & gdef=mk_gauss(xx,0.,1.,/norm)
for i=0L,npos-1L do begin
  for j=0L,nsig-1L do begin
    gg=mk_gauss(xx,pos[i],sig[j],/norm)
    DKLarr[i,j]=KLdiv(gdef,gg,Hp=Hp,Hq=Hq,zeros=1) & Hparr[i,j]=Hp
  endfor
endfor
print,''
help,DKLarr,Hparr,pos,sig
if !d.name eq 'X' then begin
  print,'contour plot shows the ratio of DKL/Hp for N(x,sigma) relative to N(0,1)'
  print,'contours are percentage change in entropy relative N(0,1)'
  window,2,xsize=900,ysize=900,title='KLdiv() compare Gaussians example'
  peasecolr
  contour,100.*DKLarr/Hparr,pos,sig,/xs,/ys,xtitle='!4d!Xx',ytitle='!4dr!X',$
  	title='KLdiv(): % change in entropy for Gaussians',$
  	levels=100*[0.001,0.005,0.01,0.02,0.03,0.05,0.08,0.13,0.21,0.34,0.55,0.89,1.44,2.33],$
	/downhill,$
	c_labels=[1,1,1,1,1,1,1,1,1,1],c_charsize=1.5,$
	c_color=[1,2,3,4,5,6,7,8,9],$
	;/fill,$
	xtickinterval=0.2,ytickinterval=0.2,xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,$
  	thick=2,xthick=2,ythick=2,charthick=2,charsize=2
endif

end

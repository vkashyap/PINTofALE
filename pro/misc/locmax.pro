;+
;function	locmax
;		returns the local maxima of a 1D array
;
;syntax
;	fmax=locmax(f,sigma=sigma,thresh=thresh,width=width,prob=prob,$
;		/poisson,nubin=nubin,fmin=fmin,fmax=fmax)
;
;parameters 	f	[INPUT; required] input array
;
;keywords	sigma	[INPUT] errors on F (ignored if POISSON is set)
;			* ignored if POISSON is set AND F(*).GE.0
;			* default is sqrt(abs(F)+0.75)+1
;		thresh	[INPUT; default: 0.5] threshold probability below
;			which to ignore a calculation
;		width	[INPUT; default: 1] compare value in bin to values
;			on either side averaged over WIDTH bins
;			* NOTE: WIDTH is the *half*width
;		poisson	[INPUT] if set, assumes a poisson error distribution
;			* automatically UNSET if min(F) < 0
;		prob	[OUTPUT] returns the computed probabilities
;
;	the following are not to be changed unless the user understands
;	what's going on!
;
;		nubin	[INPUT; default: 100] size of uniform grid to
;			superpose on F
;		fmin	[INPUT] minimum value for integration.
;			default is 0.2*min(F) or if min(F)<0 then 5*min(F)
;		fmax	[INPUT] maximum value for integration.
;			default is 5*max(F) or if max(F)<0 then 0.2*max(F)
;
;description
;	at each point, compute the probability that given value is a
;	local maximum.
;	  p = \int dF p(f|f[i],s[i])*p(f[i-1]<f|s[i-1])*p(f[i+1]<f|s[i+1])
;	find all the peaks in the distribution of probabilities that
;	lie above a given threshold.  the positions of these peaks are
;	identified with positions of local maxima.
;
;requires
;	LOCMAX_CDF (included in file)
;	IGAMMA (not in old IDLs)
;	ERRORF (not in old IDLs?)
;
;history
;	vinay kashyap (Oct96)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function locmax_cdf,x,f0,sigma=sigma,poisson=poisson
;----------------------------------------------------------
;function	locmax_cdf
;		returns the probability that f0 < x for a set of x
;		called exclusively from LOCMAX
;parameters	x,f0
;keywords	sigma	error on f0
;		poisson	if set, ignores SIGMA and generates poisson cdf
;history
;	vinay kashyap (Oct96)
;----------------------------------------------------------

nn=n_elements(x)

f=0*x
if keyword_set(poisson) then begin
  for j=0,nn-1 do f(j)=1.-igamma(x(j)>1e-8,f0>1e-8)	;incomplete Gamma func.
endif else begin
  if sigma(0) eq 0 then begin
    oo=where(x ge f0) & f(oo)=1.			;delta function
  endif else begin
    xx=(x-f0)/abs(sigma(0))/sqrt(2)
    f=0.5*(1.+errorf(xx))				;gaussian error func.
  endelse
endelse

return,f
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function locmax,f,sigma=sigma,thresh=thresh,width=width,prob=prob,$
	poisson=poisson,nubin=nubin,fmin=fmin,fmax=fmax

np=n_params(0)
if np eq 0 then begin
  print,'Usage: fx=locmax(f,sigma=sigma,thresh=thresh,width=width,prob=prob,$'
  print,'            /poisson,nubin=nubin,fmin=fmin,fmax=fmax)'
  print,'  returns array of local maxima in distribution'
  return,-1L
endif

;if sigma is not set, initialize assuming Poisson errors
nf=n_elements(f) & ns=n_elements(sigma)
if ns ne nf then sigma=sqrt(abs(f)+0.75)+1.

;check keywords
if not keyword_set(width) then width=1
;
if not keyword_set(fmin) then begin
  fmin=min(f) & if fmin lt 0 then fmin=5*fmin else fmin=0.2*fmin
endif
;
if not keyword_set(fmax) then begin
  fmax=max(f) & if fmax lt 0 then fmax=0.2*fmax else fmax=5*fmax
endif
;
if not keyword_set(nubin) then nubin=100
;
;if any points are -ve, don't bother trying to use poisson statistics!
if fmin lt 0 then poisson=0

;handle the trivial
if nf lt 2 then return,f

;check threshold
if not keyword_set(thresh) then thresh=0.5

;grid for probability computation.  compute at each point in f, plus a
;uniform distribution superposed on it.
;
df=(fmax-fmin)/(nubin-1.) & f1=f(uniq(f,sort(f))) & f2=findgen(nubin)*df+fmin
;
ff=[f1,f2] & ff=ff(sort(ff)) & ff=ff(uniq(ff,sort(ff))) & nn=n_elements(ff)
delf=ff(1:*)-ff & fav=0.5*(ff+ff(1:*))

;output arrays
lmax=0.*f & prob=lmax

;append fillers at ends
y=[fltarr(width)+fmin,f,fltarr(width)+fmin]
s=[fltarr(width),sigma,fltarr(width)]

;intermediate arrays
p_pre=fltarr(nn) & p_mid=p_pre & p_pos=p_pre & p_pre(*)=1.

;initialize p_mid
p_mid=locmax_cdf(ff,f(0),sigma=sigma(0),poisson=poisson)

for i=0L,nf-1L do begin			;step through the distribution

  ;initialize
  k=i+width & fmid=y(k) & smid=s(k)
  ;smooth over WIDTH bins
  fpre=total(y(k-width:k-1))/width & fpos=total(y(k+1:k+width))/width
  spre=sqrt(total(s(k-width:k-1)^2))/width
  spos=sqrt(total(s(k+1:k+width)^2))/width

  ;shift previously computed arrays
  if i gt 0 then begin & p_pre=p_mid & p_mid=p_pos & endif
  ;compute p_pos
  p_pos=locmax_cdf(ff,fpos,sigma=spos,poisson=poisson)

  ;compute probabilities
  pmid=abs(p_mid(1:*)-p_mid)/delf
  ppre=(0.5*(p_pre+p_pre(1:*)) > 0.) < 1.
  ppos=(0.5*(p_pos+p_pos(1:*)) > 0.) < 1.
  ;
  tmp=pmid*ppre*ppos & prob(i)=total(delf*tmp)/total(pmid*delf)
  c1=strtrim(i,2)+'  '+strtrim(fpre,2)+'+-'+strtrim(spre,2)+'  '+$
	strtrim(fmid,2)+'+-'+strtrim(smid,2)+'  '+$
	strtrim(fpos,2)+'+-'+strtrim(spos,2)+'  '+strtrim(prob(i),2)
  if prob(i) gt thresh then begin
    print,c1
    ; plot,fav,ppre,col=100,xr=[min([fpre,fpos])-spre-spos,fmid+2*smid]
    ; oplot,fav,ppos,col=150 & oplot,fav,pmid & oplot,fav,tmp,line=1
  endif

endfor

;find the peaks in the PROB distribution
oo=where(prob gt shift(prob,1) and prob gt shift(prob,-1) and prob ge thresh)
if oo(0) ne -1 then lmax(oo)=f(oo)

return,lmax
end

;+
;function	post_srcgamma
;	return the posterior probability distribution of the source
;	strength, marginalized over background.
;
;	start from gamma distribution priors and integrate the resulting
;	2D joint posterior distribution, wrt the background intensities,
;	  p(lS,lB|Y,I) ~ exp(-lB*(bB+1)+lS*(bS+1)) * 
;			 (lB+lS)^(Y) * lB^(aB-1) * lS^(aS-1),
;	where Y are the observed counts, lS,lB are the model source and
;	background intensities, and aB,aS,bB,bS are the hyperparameters
;	of the gamma-distribution priors p(lS|I) and p(lB|I) (see eqn 5
;	of van Dyk et al., 2001, ApJ 548, 224).  aB is set to the number
;	of observed background counts + 1, bB is the ratio of the
;	background area (or exposure time) to the source area (or
;	exposure time), and aA and bA are set based on prior knowledge
;	of the source intensity.
;
;syntax
;	psrc=post_srcgamma(yct,ybkg,strpr,asrc=asrc,abkg=abkg,agsrc=agsrc,$
;	bgsrc=bgsrc,clev=clev,nsgrid=nsgrid,nbgrid=nbgrid,srcmax=srcmax,$
;	srcmin=srcmin,/qromb,/tabulat,verbose=verbose, jmax=jmax,/double,$
;	K=K,eps=eps)
;
;parameters
;	yct	[INPUT; required] observed source counts
;		* must be scalar -- if vector, only 1st element is used
;	ybkg	[INPUT; required] observed counts in the background
;		region or exposure
;		* must be scalar -- if vector, adds up all the elements
;	strpr	[OUTPUT] an anonymous structure containing the fields
;		  {
;		    GRID : the grid over which the posterior probability
;			   is computed (depends on keywords AGSRC,BGSRC,
;			   NSGRID,SRCMAX,SRCMIN)
;		    MODE : the mode of the distribution
;		    STEP : the step size of the grid, denoting error on MODE
;		    MEAN : the mean of the distribution
;		    VAR  : the variance on the mean
;		    U68  : the 68% upper bound on the MODE
;		    L68  : the 68% lower bound on the MODE
;		    U90  : the 90% upper bound on the MODE
;		    L90  : the 90% lower bound on the MODE
;		    U95  : the 95% upper bound on the MODE
;		    L95  : the 95% lower bound on the MODE
;		    CRED68 : credible region enclosing cdf 0.16:0.84
;		    CRED90 : credible region enclosing cdf 0.05:0.95
;		    CRED95 : credible region enclosing cdf 0.025:0.975
;		    CRED99 : credible region enclosing cdf 0.005:0.995
;		    CLEV : the value of the keyword CLEV
;		    ULEV : the CLEV level upper bound on the MODE
;		    LLEV : the CLEV level lower bound on the MODE
;		    CREDLEV : credible region enclosing cdf
;			      (0.5-CLE/2.):(0.5+CLEV/2.)
;		    }
;
;keywords
;	asrc	[INPUT; default=1.0] area (or exposure time) for the
;		(background-contaminated) source observation
;	abkg	[INPUT; default=1.0] area (or exposure time) for the
;		background observation
;	agsrc	[INPUT; default=1.0] alpha parameter for the source intensity
;		prior, a gamma distribution (see PROB_GAMMADIST)
;	bgsrc	[INPUT; default=0.] beta parameter for the source prior
;		* note that the default values of AGSRC,BGSRC make for a
;		  very non-informative, but improper, prior
;	clev	[INPUT] a nominal level at which to determine bounds on
;		the mode and get credible regions
;		* default is 0.99
;		* if < 0, abs value is used
;		* if > 1 and < 100, then assumed to be given as a percentage
;		* if > 100, then 1-1/CLEV is used as the true value
;	nsgrid	[INPUT; default=100] number of points in the output grid
;		* used only if LSRC is invalid
;	nbgrid	[INPUT; default=NSGRID] number of points in marginalization
;		integration
;		* used only if QROMB is not set, and even then, really
;		  matters only if TABULAT is not set
;	srcmax	[INPUT; default=(ALFA_B/BETA_B)+10*(ALFA_B/BETA_B^2)] maximum
;		in source intensity value to use in calculating the grid
;		* set to YCT+10*sqrt(YCT) if BETA_B is 0
;		* there is a hardcoded lower limit of 10
;	srcmin	[INPUT; default=((ALFA_B/BETA_B)-10*(ALFA_B/BETA_B^2))>0]
;		minimum in source intensity value to use in calculating the
;		grid
;		* set to (YCT-10*sqrt(YCT))>0 if BETA_B is 0
;		* using these YCT based limits on SRCMAX and SRCMIN is an
;		  acknowledgement of the fact that, yes, virginia, even if
;		  BGSRC is 0, there _is_ a prior
;	tabulat	[INPUT] if set, uses the Newton-Coates 5-point integration
;		method (see INT_TABULATED())
;	qromb	[INPUT] if set, uses Romberg integration routine
;		(see QROMB())
;		* the default is to use simple quadrature integration
;		* QROMB takes precedence over TABULAT
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutine
;		-- QROMB: JMAX, DOUBLE, EPS, K
;
;common blocks
;	vandyk	to communicate function parameters to MK_BKGGAMMA
;		ALFA_S,BETA_S,ALFA_B,BETA_B,Y_SB,LAMDA_S,KSTEP,
;		NORM1,NORM2,NORM3,NORM4
;
;example
;	p=post_srcgamma(1,48,s,asrc=10.,abkg=120.,verbose=21) & help,s,/str
;	p=post_srcgamma(1,48,s,asrc=10.,abkg=120.,agsrc=4,bgsrc=1,verbose=21)
;
;history
;	vinay kashyap (Aug01)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mk_bkggamma,lamda_b
; +
; function mk_bkggamma
; 	compute the projection of the 2D joint posterior probability
;	distribution of the source and background intensities,
;	p(lamda_s,lamda_s|Y,I) at a given value of lamda_s
;	(cf. van Dyk et al., 2001, ApJ 548, 224)
;
; restrictions
; 	because it is expected to be called only from POST_SRCGAMMA,
;	no checks (except the most rudimentary) are made on the variables.
;	they better be all OK.
;
; common blocks
; 	vandyk	ALFA_S,BETA_S,ALFA_B,BETA_B,Y_SB,LAMDA_S,KSTEP
; -
;
common vandyk,alfa_S,beta_S,alfa_B,beta_B,y_SB,lamda_S,kstep,$
	norm1,norm2,norm3,norm4
;
if kstep lt 0 then begin
  ;	need this to avoid infinities!
  norm1=0.0D & norm2=0.0D & norm3=0.0D & norm4=0.0D
endif
;
x=lamda_b[0] & kstep=kstep+1L
;
p=0. & tmp1=0. & tmp2=0. & tmp3=0. & tmp4=0. & tmp5=0.D
p0=0	;=1 if p is identically 0.
;
tmp1=x*(beta_B+1.)+lamda_S*(beta_S+1.)
;
if x+lamda_S gt 0 then tmp2=y_SB*alog(x+lamda_S) else if y_SB gt 0 then p0=1
;
if x eq 0 and alfa_B eq 1 then tmp3=0. else begin
  if x ne 0 then tmp3=(alfa_B-1.)*alog(x) else p0=1
endelse
;
if lamda_S eq 0 and alfa_S eq 1 then tmp4=0. else begin
  if lamda_S ne 0 then tmp4=(alfa_S-1.)*alog(lamda_S) else p0=1
endelse
;
;	;need the following for large ALFA_B, ALFA_S, etc
tmp1=tmp1-norm1
tmp2=tmp2-norm2
tmp3=tmp3-norm3
tmp4=tmp4-norm4
	;need this for large Y, else p blows up
;if y_SB gt 1 then tmp5=lngamma(y_SB-1.)+lngamma(alfa_B/beta_B)
;if (-norm1+norm2+norm3+norm4)-tmp5 lt 0. then tmp5=0.D
;
if p0 eq 0 and (-tmp1+tmp2+tmp3+tmp4) lt -69 then p0=1
if p0 ne 1 then p = exp(-tmp1 + tmp2 + tmp3 + tmp4 - tmp5)
;
if kstep eq 0 then begin
  norm1=tmp1 & norm2=tmp2 & norm3=tmp3 & norm4=tmp4
endif
;
;c='X='+strtrim(x,2)+' lamda_S='+strtrim(lamda_S,2)+' p='+strtrim(p,2)+$
;	' kstep='+strtrim(kstep,2) & print,c
;
return,p
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function post_srcgamma,yct,ybkg,strpr,asrc=asrc,abkg=abkg,agsrc=agsrc,$
	bgsrc=bgsrc,clev=clev,nsgrid=nsgrid,nbgrid=nbgrid,srcmax=srcmax,$
	srcmin=srcmin,qromb=qromb,tabulat=tabulat,verbose=verbose, _extra=e

message,'OBSOLETE: use PPD_SRC() instead',/informational

;	usage
ok='ok' & np=n_params() & ny=n_elements(yct) & nb=n_elements(ybkg)
if np lt 2 then ok='Insufficient parameters' else $
 if ny eq 0 then ok='YCT is undefined' else $
  if nb eq 0 then ok='YBKG is undefined'
if ok ne 'ok' then begin
  print,'Usage: psrc=post_srcgamma(yct,ybkg,strpr,asrc=asrc,abkg=abkg,$'
  print,'       agsrc=agsrc,bgsrc=bgsrc,nsgrid=nsgrid,nbgrid=nbgrid,$'
  print,'       srcmax=srcmax,srcmin=srcmin,/qromb,/tabulat,verbose=verbose,$'
  print,'       jmax=jmax,/double,K=K,eps=eps)'
  print,'  compute posterior probability distribution function for source'
  print,'  intensity, marginalized over background intensities'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	now define the common block
common vandyk,alfa_S,beta_S,alfa_B,beta_B,y_SB,lamda_S,kstep,$
	norm1,norm2,norm3,norm4

;	some special input cases
y_SB=yct[0]+0.0
if ny gt 1 then message,$
	'observed counts must be a scalar; using only first element',/info
if y_SB lt 0 then begin
  message,'cannot observe negative counts; returning immediately',/info
  return,-1L
endif
;
if nb gt 1 then begin
  ybg=total(ybkg)
  bkgar=fltarr(nb)+1. & nbg=n_elements(abkg)
  if nbg gt 0 then begin
    bkgar[*]=abkg[nbg-1L]
    if nbg lt nb then bkgar[0L:nbg-1L]=abkg[*] else bkgar[*]=abkg[0L:nb-1L]
  endif
  bkgar=total(bkgar)
endif else begin
  ybg=ybkg[0]+0.0
  bkgar=1. & if keyword_set(abkg) then bkgar=0.0+abkg[0]
endelse

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;
srcar=1. & if keyword_set(asrc) then srcar=0.0+asrc[0]
;
alfa_S=1. & if n_elements(agsrc) gt 0 then alfa_S=0.0+agsrc[0]
;
beta_S=0. & if n_elements(bgsrc) gt 0 then beta_S=0.0+bgsrc[0]
;
alfa_B=ybg[0]+1.0
;
beta_B=1. & if srcar ne 0 then beta_B=(bkgar/srcar)
;
msgrid=100L & if keyword_set(nsgrid) then msgrid=long(abs(nsgrid[0])) > 2
;
mbgrid=msgrid>10 & if keyword_set(nbgrid) then mbgrid=long(abs(nbgrid[0]))>10
;
;
smax=y_SB+10.*sqrt(y_SB) > 10.
if beta_S ne 0 then smax=(alfa_S/beta_S)+10.*sqrt(alfa_S/beta_S^2)
if keyword_set(srcmax) then smax=0.0+srcmax[0]
;
smin=y_SB-10.*sqrt(y_SB) > 0.
if beta_S ne 0 then smin=(alfa_S/beta_S)-10.*sqrt(alfa_S/beta_S^2) > 0
if keyword_set(srcmin) then smin=0.0+srcmin[0]
;
if smax lt smin then begin
  message,'SRCMAX and SRCMIN being interchanged',/info
  junk=smin & smin=smax & smax=junk
endif
;
bmax=(alfa_B/beta_B)+10.*sqrt(alfa_B/beta_B^2) > 10.
;
bmin=(alfa_B/beta_B)-10.*sqrt(alfa_B/beta_B^2) > 0.
;
credlev=0.99 & if keyword_set(clev) then credlev=0.0+clev[0]
if credlev lt 0 then credlev=abs(credlev)
if credlev gt 1 and credlev lt 100 then credlev=credlev/100.
if credlev gt 100 then credlev=1.0D - 1.0D/credlev

;	now define the lamda_S and lamda_B grids
dS=(smax-smin)/(msgrid-1L) & Sgrid=findgen(msgrid)*dS+smin
dB=(bmax-bmin)/(mbgrid-1L) & Bgrid=findgen(mbgrid)*dB+bmin

;	and for each value of Sgrid, integrate over grid_B
psrc=dblarr(msgrid)
lamda_S=y_SB & kstep=-1L & tmp=mk_bkggamma(alfa_B/beta_B) ;set normalizations
for i=0L,msgrid-1L do begin
  lamda_S = Sgrid[i] & kstep=0L
  if keyword_set(qromb) then begin
    psrc[i]=qromb('mk_bkggamma',bmin,bmax, _extra=e)
  endif else begin
    pb=dblarr(mbgrid)
    for j=0L,mbgrid-1L do pb[j]=mk_bkggamma(double(Bgrid[j]))
    if keyword_set(tabulat) then psrc[i]=int_tabulated(Bgrid,pb) else $
    	psrc[i]=total(pb*dB)
  endelse
  if vv gt 30 then print,'p('+strtrim(lamda_S,2)+'|'+strtrim(y_SB,2)+')='+$
	strtrim(psrc[i],2)
endfor

;	normalize the output
norm=total(psrc) & if norm gt 0 then psrc=psrc/norm

;	make the summary structure
tmp=max(psrc,im) & prbmode=Sgrid[im]
prbav=total(psrc*Sgrid)
prbvar=total(psrc*Sgrid^2)-prbav^2
;
cdf=fltarr(msgrid+1L) & for i=1L,msgrid do cdf[i]=cdf[i-1]+psrc[i-1L]
cdf=cdf/cdf[msgrid] & xcdf=[Sgrid,Sgrid[msgrid-1L]+dS]
ucdf=cdf-cdf[im] & ucdf=ucdf/max(ucdf) & uu=[smax,smax,smax,smax]
if im lt msgrid then uu=interpol(xcdf,ucdf,[0.68,0.90,0.95,credlev])
lcdf=cdf & lcdf[im:*]=cdf[im] & lcdf=lcdf/max(lcdf) & ll=[smin,smin,smin,smin]
if im gt 0 then ll=interpol(xcdf,lcdf,[0.32,0.10,0.05,1.0-credlev])
cred=interpol(xcdf,cdf,[0.005,0.025,0.05,0.16,0.84,0.95,0.975,0.995,$
	(0.5-credlev/2.),(0.5+credlev/2.)])
;
strpr=create_struct('GRID',Sgrid,'MODE',prbmode,'STEP',dS,$
	'MEAN',prbav,'VAR',prbvar,$
	'U68',uu[0],'L68',ll[0],'U90',uu[1],'L90',ll[1],$
	'U95',uu[2],'L95',ll[2],$
	'CRED68',[cred[3],cred[4]],$
	'CRED90',[cred[2],cred[5]],$
	'CRED95',[cred[1],cred[6]],$
	'CRED99',[cred[0],cred[7]],$
	'CLEV',credlev,'ULEV',uu[3],'LLEV',ll[3],'CREDLEV',[cred[8],cred[9]])

;	report
if vv gt 0 then begin
  tt='Y='+strtrim(y_SB,2)+'/'+strtrim(srcar,2)+$
    ' - B='+strtrim(ybg,2)+'/'+strtrim(bkgar,2)
  if vv gt 2 then message,tt+' ==> E(S)='+strtrim(prbav,2),/info
  if vv gt 3 then message,'ALPHA_S='+strtrim(alfa_S,2)+' ; '+$
	'BETA_S='+strtrim(beta_S,2)+' ; '+$
	'ALPHA_B='+strtrim(alfa_B,2)+' ; '+$
	'BETA_B='+strtrim(beta_B,2),/info,/noname
  if vv gt 20 then begin
    plot,sgrid,psrc,charsize=2,$
	xtitle='!4k!3!uS!n',ytitle='!3p(!4k!3!uS!n|Y,I)',$
	title=strtrim(strpr.MODE),subtitle=tt
    tmp=prob_gammadist(sgrid,agamma=alfa_S,bgamma=beta_S)
    oplot,sgrid,tmp/total(tmp),line=1
  endif
endif

return,psrc
end

function ppd_src,ns,nb,sgrid,smap=smap,scred=scred,clev=clev,ppdstr=ppdstr,$
	funit=funit,asrc=asrc,abkg=abkg,agamma=agamma,bgamma=bgamma,$
	priorbg=priorbg,nsgrid=nsgrid,srcmin=srcmin,srcmax=srcmax,$
	verbose=verbose, _extra=e
;+
;function	ppd_src
;	returns the posterior probability density for the source
;	intensity, marginalized over background intensity, assuming
;	gamma-function priors (see van Dyk et al., 2001, ApJ 548, 224)
;
;	D:ns,nb ; I:asrc,abkg,f ; M:s,b ; p(M|DI)=p(M|I)*p(D|MI)/p(D|I)
;
;	r = abkg/asrc & A = nb+alfa_B & B = f*(beta_B+r)
;
;	p(b|I) = ((f*beta_B)^alfa_B/Gamma(alfa_B)) * b^(alfa_B-1) * exp(-b*beta_B*f)
;	p(nb|b,I) = ((r*f*b)^nb/Gamma(nb+1)) * exp(-b*f*r)
;	p(b|nb,I) = (B^A/Gamma(A)) * b^(A-1) * exp(-b*B)
;
;	I_k = b^(ns+A-k-1) * s^(k+alfa_S-1) * exp(-s*f*(1+beta_S) -b*(f+B)) /
;		(Gamma(k+1)*Gamma(ns-k+1))
;	I_kb = (Gamma(ns+A-k)/Gamma(k+1)/Gamma(ns-k+1)) * (f+B)^(-(ns+A-k)) *
;		s^(k+alfa_S-1) * exp(-s*f*(1+beta_S))
;	I_kbs = (Gamma(ns+A-k)*Gamma(k+alfa_S)/Gamma(k+1)/Gamma(ns-k+1)) *
;		(f+B)^(-(ns+A-k)) * (f*(1+beta_S))^(k+alfa_S)
;
;	p(s|I) = ((f*beta_S)^alfa_S/Gamma(alfa_S)) * s^(alfa_S-1) * exp(-s*beta_S*f)
;	p(ns|M,I) = (f*(s+b))^ns/Gamma(ns+1) * exp(-(s+b)*f)
;	p(sb|D,I) = p(b|nb,I)*p(s|I)*p(ns|M,I) / p(ns|I)
;		  = (Sum_{k=0}^{ns} I_k) / (Sum_{k=0}^{ns} I_kbs)
;
;	p(ns|I) = f^ns * (B^A/Gamma(A)) * ((f*beta_S)^alfa_S/Gamma(alfa_S)) *
;		Sum_{k=0}^{ns} I_kbs
;	p(s|D,I) = (Sum_{k=0}^{ns} I_kb) / (Sum_{k=0}^{ns} I_kbs)
;
;syntax
;	ppd=ppd_src(ns,nb,sgrid,smap=smap,scred=scred,clev=clev,ppdstr=ppdstr,$
;	funit=funit,asrc=asrc,abkg=abkg,agamma=agamma,bgamma=bgamma,$
;	priorbg=priorbg,nsgrid=nsgrid,srcmin=srcmin,srcmax=srcmax,verbose=verbose)
;
;parameters
;	ns	[INPUT; required] observed counts in source region
;	nb	[INPUT; required] observed counts in background region
;		* NS and NB must be scalar -- if vector, adds up all the elements
;	sgrid	[OUTPUT] the output grid over which the PPD is calculated
;		* a regular grid that uses NSGRID, SRCMIN, and SRCMAX
;
;keywords
;	smap	[OUTPUT] the mode of the distribution, corresponding to the
;		maximum a posteriori value of the PPD
;	scred	[OUTPUT] a two-tailed credible range on the source intensity,
;		defined to include CLEV/2 of the area under the PPD above
;		and below SMAP
;		* NOTE: what this means is, if the distribution is skewed,
;		  then the total area under PPD as defined by SCRED is not
;		  necessarily equal to CLEV
;		* The one-tailed credible regions, that is, the range that
;		  encloses an area CLEV under the PPD, are returned in the
;		  structure PPDSTR.  these regions are determined such that
;		  the mode is always _included_ within the region, but is not
;		  necessarily at the center
;	clev	[INPUT] a nominal level at which to determine bounds and
;		get credible regions
;		* default is 0.68
;		* if < 0, abs(CLEV) is used
;		* if > 1 and < 100, then assumed to be given as a percentage
;		* if > 100, then 1-1/CLEV is used
;	ppdstr	[OUTPUT] an anonymous structure containing the fields
;		{
;		  GRID : same as SGRID
;		  MODE : mode, same as SMAP
;		  MEAN : mean of the distribution
;		  SIG : s.d. on the MEAN
;		  MEDIAN : the median of the distribution
;		  CLEV : the probability level defined by the keyword CLEV
;		  HPD1 : the highest posterior density "1-sigma" credible range,
;			 the interval that encloses 68% of the area and also
;			 the highest values of the probability density
;		  HPD2 : as HPD1, for "2-sigma", 95% of the area
;		  HPD3 : as HPD1, for "3-sigma", 99.7% of the area
;		  HPD0 : as HPD1, but corresponding to a fraction CLEV
;		  CRED1 : the "1-sigma" credible range, corresponding to 68%
;		  	  of the area under the PPD that _includes_ the MODE
;		  CRED2 : as CRED1, for "2-sigma", 95% of the area
;		  CRED3 : as CRED1, for "3-sigma", 99.7% of the area
;		  CRED0 : as CRED1, but corresponding to a fraction CLEV of the
;			  area under the PPD that includes the MODE
;		  RNG1 : the two-tailed "1-sigma" range, corresponding to a
;		  	 fractional area under the PPD of 0.34 above and below
;		  	 the MODE
;		  RNG2 : as for RNG1, but the "2-sigma" range, 0.475 of the area
;		  RNG3 : as for RNG1, but the "3-sigma" range, 0.4985 of the area
;		  RNG0 : as for RNG1, but CLEV/2 of the area, same as SCRED
;		  lnpD : the log of the predictive probability, ALOG(p(ns|I))
;		}
;	funit	[INPUT] a factor to _divide_ the counts to convert to different units
;		e.g., ph/cm^2/s
;		* assumed to be the same for both source and background counts
;		  (any differences in exposure times can be encoded in ASRC and ABKG)
;		* default is 1.0
;	asrc	[INPUT] area (or exposure time) for the background contaminated
;		source observation, NS
;		* default is 1.0
;	abkg	[INPUT] area (or exposure time) for the background observation, NB
;		* default is 1.0
;	agamma	[INPUT] alpha parameter for the source intensity prior, assumed to
;		be a gamma distribution (see PROB_GAMMADIST)
;		* default is 1.0
;		* hardcoded minimum is 0.001
;		* generally taken to be the number of *expected* source counts+1
;	bgamma	[INPUT] beta parameter for the source intensity prior
;		* default is 0.0 (also the hardcoded minimum)
;		* generally the ratio of that area in which the AGAMMA-1 counts are
;		  obtained in, to ASRC
;		* note that the default values of AGAMMA and BGAMMA make for a
;		  non-informative, but improper, prior
;	priorbg	[INPUT] the expected number of counts due to the background,
;		a number that is known from prior information and NOT derived
;		from NB itself.
;		* a Gamma function prior is assumed on the background strength, with
;		  alpha=priorbg+1; beta=abkg/asrc, but if priorbg=0, then beta=0
;	nsgrid	[INPUT] number of points in output grid
;		* default=101L
;	srcmin	[INPUT] minimum in source intensity value to consider
;	srcmax	[INPUT] maximum in source intensity value to consider
;		* beware that setting [SRCMIN,SRCMAX] to too small a range will cause
;		  it to act as a prior over and above the gamma prior.  this may not be
;		  desirable.
;		* the default range is the widest of:
;		  -- NS +- 10*sqrt(NS)
;		  -- (NS-(ASRC/ABKG)*NB) +- 10*sqrt(NS+(ASRC/ABKG)^2*NB) (if ABKG > 0)
;		  -- (agamma/bgamma) +- 10*sqrt(agamma/bgamma^2) (if BGAMMA > 0)
;		  as long as SRCMIN isn't -ve or SRCMAX is smaller than 20
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	PROB_GAMMADIST
;
;history
;	vinay kashyap (May02)
;	added keyword FUNIT (VK; Oct03)
;	added HPD1,HPD2,HPD3,HPD0 to the output (VK; Mar06)
;	changed hardcoded minimum of AGAMMA from 1 to 0.001 (VK; Feb16)
;	added GDL bypass for SPLINE_P (VK; Apr17)
;	forced HPD minimum to be always >0 (VK; Jun18)
;-

;	usage
ok='ok' & np=n_params() & nns=n_elements(ns) & nnb=n_elements(nb)
if np lt 2 then ok='Insufficient parameters' else $
 if nns eq 0 then ok='NS: undefined' else $
  if nnb eq 0 then ok='NB: undefined'
if ok ne 'ok' then begin
  print,'Usage: ppd=ppd_src(ns,nb,sgrid,smap=smap,scred=scred,clev=clev,ppdstr=ppdstr,$'
  print,'       funit=funit,asrc=asrc,abkg=abkg,agamma=agamma,bgamma=bgamma,$'
  print,'       priorbg=priorbg,nsgrid=nsgrid,srcmin=srcmin,srcmax=srcmax,verbose=verbose)'
  print,'  return posterior probability density of source intensity, marginalized'
  print,'  over background intensities'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	check input parameters
if nns gt 1 then Ds=total(ns) else Ds=ns[0]
if nnb gt 1 then Db=total(nb) else Db=nb[0]
if Ds lt 0 or Db lt 0 then begin
  message,'cannot handle -ve observed counts; quitting',/info & return,-1L
endif

;	check input keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;
crlev=0.68 & if keyword_set(clev) then crlev=0.0+clev[0]
if crlev lt 0 then crlev=abs(crlev)
if crlev ge 1 and crlev lt 100 then crlev=crlev/100.
if crlev ge 100 then crlev = 1.0D - 1.0D/crlev
;
facu=1.0 & if keyword_set(funit) then facu=1.0*funit[0]
;
srcar=1.0 & if keyword_set(asrc) then srcar=0.0+asrc[0]
if srcar le 0 then begin
  message,'Cannot handle -ve ASRC; quitting',/info & return,-1L
endif
bkgar=1.0 & if keyword_set(abkg) then bkgar=0.0+abkg[0]
if bkgar lt 0 then begin
  message,'Cannot handle -ve ABKG; quitting',/info & return,-1L
endif
;
alfa_S=1. & if keyword_set(agamma) then alfa_S=0.0+agamma[0]
if alfa_S lt 1 then begin
  if vv gt 0 then message,'AGAMMA cannot be .LE. 0; resetting to 0.001',/info
  alfa_S = 0.001
endif
;
beta_S=0. & if n_elements(bgamma) gt 0 then beta_S=0.0+bgamma[0]
if beta_S lt 0 then begin
  if vv gt 0 then message,'BGAMMA cannot be < 0; resetting',/info
  beta_S = 0.
endif
;
alfa_B=1. & beta_B=0.
if keyword_set(priorbg) then begin
  alfa_B=priorbg[0]+1.>1.
  beta_B=bkgar/srcar
endif
;
msgrid=101L & if keyword_set(nsgrid) then msgrid=long(abs(nsgrid[0])) > 2
;
s_ML=Ds/facu & s_MLe=(sqrt(s_ML+0.75)+1.)/facu
smin1=0 & smax1=20/facu & smin=smin1 & smax=smax1
smin2=s_ML-10*s_MLe > smin1
smax2=s_ML+10*s_MLe > smax1
smin=smin2 & smax=smax2
if bkgar gt 0 then begin
  s_ML=Ds/facu-(srcar/bkgar)*Db/facu
  s_MLe=(sqrt((sqrt(Ds+0.75)+1.)^2+(srcar/bkgar)^2*(sqrt(Db+0.75)+1.)^2))/facu
  smin3=s_ML-10*s_MLe > smin1
  smax3=s_ML+10*s_MLe > smax1
  if smin3 lt smin2 then smin=smin3
  if smax3 gt smax2 then smax=smax3
endif
if beta_S gt 0 then begin
  smin4=(alfa_S/beta_S) - 10*sqrt(alfa_S/beta_S^2) > smin1
  smax4=(alfa_S/beta_S) + 10*sqrt(alfa_S/beta_S^2) > smax1
  if smin4 lt smin2 then smin=smin4
  if smax4 gt smax2 then smax=smax4
endif
if keyword_set(srcmin) then begin
  if srcmin[0] gt smin then begin
    if vv gt 0 then message,'SRCMIN > default; presumably you know what you are doing',/info
  endif
  smin=srcmin[0]
endif
if keyword_set(srcmax) then begin
  if srcmax[0] lt smax then begin
    if vv gt 0 then message,'SRCMAX < default; presumably you know what you are doing',/info
  endif
  smax=srcmax[0]
endif
if smin gt smax then begin
  message,"SRCMIN > SRCMAX?  can't be havin' wi' that",/info
  jnk=smax & smax=smin & smin=jnk
endif
deltS=(smax-smin)/(msgrid-1L) & sgrid=findgen(msgrid)*deltS+smin

;	output
lnppd=fltarr(msgrid) & ppd=dblarr(msgrid)

;	now compute the coefficients
kk=lindgen(Ds+1L)
lnI_kbs = lngamma(Ds+Db+(alfa_B-1)-kk+1L) +$
	lngamma(kk+alfa_S) -$
	lngamma(kk+1L) -$
	lngamma(Ds-kk+1L) -$
	(Ds+Db+(alfa_B-1)-kk+1L)*alog(1.0+(bkgar/srcar)+beta_B) -$
	(Ds+Db+(alfa_B-1)-kk+1L)*alog(facu) -$
	(kk+alfa_S)*alog(1.+beta_S) -$
	(kk+alfa_S)*alog(facu)
normkbs=max(lnI_kbs) & dkbs=lnI_kbs-normkbs
sumkbs=total(exp(dkbs)) & sumkbs=alog(sumkbs)+normkbs
lnI_kb = lngamma(Ds+Db+(alfa_B-1)-kk+1L) -$
	lngamma(kk+1L) -$
	lngamma(Ds-kk+1L) -$
	(Ds+Db+(alfa_B-1)-kk+1L)*alog(1.0+(bkgar/srcar)+beta_B) -$
	(Ds+Db+(alfa_B-1)-kk+1L)*alog(facu)

;	r = abkg/asrc & A = nb+alfa_B & B = f*(beta_B+r)
;	I_kbs = (Gamma(ns+A-k)*Gamma(k+alfa_S)/Gamma(k+1)/Gamma(ns-k+1)) *
;		(f+B)^(-(ns+A-k)) * (f*(1+beta_S))^(k+alfa_S)
;	p(D|I) = f^ns * (B^A/Gamma(A)) * ((f*beta_S)^alfa_S/Gamma(alfa_S)) *
;		Sum_{k=0}^{ns} I_kbs
if beta_S gt 0 then p_D_I=alfa_S*alog(beta_S) else p_D_I=0.
p_D_I = p_D_I + ns*alog(facu) +$
	(Db+alfa_B)*alog(beta_B+(bkgar/srcar)) +$
	(Db+alfa_B)*alog(facu) -$
	lngamma(Db+alfa_B) +$
	alfa_S*alog(facu) -$
	lngamma(alfa_S) +$
	sumkbs
;p_D_I = exp(p_D_I)

;	for each source intensity grid point, compute the ppd
for i=0L,msgrid-1L do begin
  s=sgrid[i]
  if s eq 0 then begin
    tmp=lnI_kb[0] & normkb=tmp & sumkb=1.
    lnppd[i]=lnI_kb[0]-sumkbs ;lnI_kbs[0]
  endif else begin
    tmp=lnI_kb+(kk+alfa_S-1L)*alog(s)
    normkb=max(tmp) & dkb=tmp-normkb & sumkb=total(exp(dkb))
    lnppd[i]=alog(sumkb)+normkb-sumkbs-(1.+beta_S)*s*facu
  endelse
  if vv ge 150 then begin
    print,'ln{p('+strtrim(s,2)+'|'+strtrim(Ds,2)+','+strtrim(Db,2)+',I)}='+strtrim(lnppd[i],2)
    if vv gt 200 then stop,s
  endif
endfor
lnormppd=max(lnppd) & ppd=exp(lnppd-lnormppd)
normppd=total(ppd*deltS) & ratnorm=exp(-lnormppd)/normppd	;just in case
ppd=ppd*exp(lnormppd)
;
if vv gt 90 then print,'normalization (should be v/s is) =',exp(-lnormppd),normppd
if ratnorm gt 1.5 or ratnorm lt 1/1.5 and vv gt 0 then begin
  message,'	grid appears to be non optimal;',/info
  message,'	increase NSGRID, or fiddle with SRCMIN and SRCMAX',/info
endif

;	find the mode
defsysv,'!GDL',exists=igdl
if not keyword_set(igdl) then begin
  if facu eq 1 then spline_p,sgrid,ppd,ss,pp,interval=deltS/100. else begin
    if facu lt 100. then spline_p,sgrid,ppd,ss,pp else begin
      ss=sgrid & pp=ppd
    endelse
  endelse
endif else begin
  ss=sgrid & pp=ppd
endelse
pp=pp>0
jnk=max(pp,imx) & smap=ss[imx]

;	find the highest posterior density intervals
os=reverse(sort(pp)) & cpp=total(double(pp[os]),/cumul) & cpp=cpp/max(cpp) & xx=ss[os]
om=where(xx le smap,mom) & hpdm=fltarr(4)+smap
op=where(xx ge smap,mop) & hpdp=fltarr(4)+smap
if mom gt 1 then hpdm=interpol(xx[om],cpp[om],[crlev,0.68,0.95,0.997]) > 0
if mop gt 1 then hpdp=interpol(xx[op],cpp[op],[crlev,0.68,0.95,0.997])

;	find the mean
;smean=total(ss*pp*deltS/100.) & ssig=sqrt(total(ss^2*pp*deltS/100.)-smean^2)
smean=total(sgrid*ppd*deltS) & ssig=sqrt(total(sgrid^2*ppd*deltS)-smean^2)
if vv gt 90 and finite(ssig) le 0 then stop,smean,ssig

;	make cumulative distribution
sp=dindgen(msgrid+1L)*deltS+smin-deltS/2.>0 & cp=dblarr(msgrid+1L)
for i=1L,msgrid do cp[i]=cp[i-1L]+ppd[i-1L] & cp=cp/max(cp)
smedian=interpol(sp,cp,0.5)

;	make cumulative 2-tailed distributions
cmode=interpol(cp,sp,smap)
ucp=cp-cmode & ou=where(ucp ge 0,mou)
if mou gt 1 then begin
  usp=[smap,sp[ou]] & ucp=[0.,ucp[ou]] & ucp=ucp/max(ucp)
endif else begin
  usp=[smap,sp[msgrid]] & ucp=[0.,1.]
endelse
lcp=-(cp-cmode) & ol=where(lcp ge 0,mol)
if mol gt 1 then begin
  lsp=[sp[ol],smap] & lcp=[lcp[ol],0.] & lcp=lcp/max(lcp)
endif else begin
  lsp=[sp[0],smap] & lcp=[1.,0.]
endelse
levs=[crlev,0.68,0.95,0.997]
uu=interpol(usp,ucp,levs) & ll=interpol(lsp,lcp,levs)
scred=[ll[0],uu[0]]

;	find credible regions that have the right area AND
;	include the mode in all cases
ulev=0.5D + levs/2.D - (0.5D - cmode)
llev=0.5D - levs/2.D - (0.5D - cmode)
oo=where(llev lt 0,moo)
if moo gt 0 then begin
  ulev[oo]=ulev[oo]-llev[oo] & llev[oo]=0.
endif
oo=where(ulev gt 1,moo)
if moo gt 0 then begin
  llev[oo]=llev[oo]-(ulev[oo]-1.) & ulev[oo]=1.
endif
ureg=interpol(sp,cp,ulev) < smax
lreg=interpol(sp,cp,llev) > smin

;	make the summary structure
ppdstr=create_struct('GRID',sgrid,'MODE',smap,'MEAN',smean,'SIG',ssig,$
	'MEDIAN',smedian,'CLEV',crlev,$
	'HPD1',[hpdm[1],hpdp[1]],'HPD2',[hpdm[2],hpdp[2]],$
	'HPD3',[hpdm[3],hpdp[3]],'HPD0',[hpdm[0],hpdp[0]],$
	'CRED1',[lreg[1],ureg[1]],'CRED2',[lreg[2],ureg[2]],$
	'CRED3',[lreg[3],ureg[3]],'CRED0',[lreg[0],ureg[0]],$
	'RNG1',[ll[1],uu[1]],'RNG2',[ll[2],uu[2]],$
	'RNG3',[ll[3],uu[3]],'RNG0',[ll[0],uu[0]],$
	'lnpD',p_D_I)

;	report
tt0='S = '+strtrim(smap,2)+$
	' +'+strtrim(scred[1]-smap,2)+',-'+strtrim(smap-scred[0],2)+$
	' @'+strtrim(crlev,2)
tt1='MaxLik='+strtrim(s_ML,2)+'+-'+strtrim(s_MLe,2)+$
	' ; MEAN='+strtrim(smean,2)+'+-'+strtrim(ssig,2)+$
	' ; MEDIAN='+strtrim(smedian,2)
if vv gt 0 then begin
  message,tt0,/info
  if vv gt 1 then message,tt1,/info,/noname
  if vv gt 3 then message,'AGAMMA='+strtrim(alfa_S,2)+' ; '+$
	'BGAMMA='+strtrim(beta_S,2)+' ; '+$
	'ASRC='+strtrim(srcar,2)+' ; '+$
	'ABKG='+strtrim(bkgar,2)+' ; '+$
	'FUNIT='+strtrim(facu,2),/informational,/noname
  if vv gt 20 then begin
    plot,sgrid,ppd,charsize=2,$
	xtitle='!4k!3!uS!n',ytitle='!3p(!4k!3!uS!n|n!dS!n,n!dB!n,I)',$
	title=strtrim(smap,2),subtitle=tt0
    tmp=prob_gammadist(sgrid,agamma=alfa_S,bgamma=beta_S)
    oplot,sgrid,tmp*total(ppd)/total(tmp),line=1,thick=2
  endif
  if vv ge 100 then stop,'HALTing.  type .CON to continue'
endif else begin
  if np eq 2 and not arg_present(smap) and not arg_present(ppdstr) then begin
    message,tt0,/info
    message,tt1,/info,/noname
  endif
endelse

return,ppd
end

;+
;low_counts_significance.com
;
;description:
;	script to estimate the significance of some small number of counts
;	in a source region, given a measurement of the background in a
;	non-overlapping background region.  See
;		http://groundtruth.info/AstroStat/slog/2008/significance-of-5-counts/
;	for the genesis of this script.
;
;usage
;	set the number of counts observed in the source region, OBS_COUNTS
;	set the number of counts observed in the background region, BG_COUNTS
;	set the area of the background region relative to the source region, BG_AREA
;	@low_counts_significance.com
;
;subroutines
;	LNPOISSON()
;
;vinay kashyap (Apr2008)
;-

;	check inputs
if n_elements(obs_counts) eq 0 then obs_counts=5	;default is 5
if not keyword_set(bg_counts) then bg_counts=6	;can't accept BG_COUNTS=0
if not keyword_set(bg_area) then bg_area=50.0	;can't accept BG_AREA=0

;	verify inputs
ok='ok'
if n_elements(obs_counts) gt 1 then ok='OBS_COUNTS must be scalar' else $
 if size(obs_counts,/type) gt 4 then ok='OBS_COUNTS must be a number' else $
  if n_elements(bg_counts) gt 1 then ok='BG_COUNTS must be scalar' else $
   if size(bg_counts,/type) gt 4 then ok='BG_COUNTS must be a number' else $
    if n_elements(bg_area) gt 1 then ok='BG_AREA must be scalar' else $
     if size(bg_area,/type) gt 4 then ok='BG_AREA must be a number'
if ok ne 'ok' then message,ok
help,obs_counts,bg_counts,bg_area

; compute the probability of seeing IB counts in the bg area
nb=long(bg_counts+10*sqrt(bg_counts)+0.5) & ib=findgen(nb+1)
; 	one way is to use IGAMMA() to compute the cumulative probability
; 	and take the difference to get the differential, like so --
; cpb=igamma(ib,BG_COUNTS) & dpb=abs(cpb[1:*]-cpb)
; 	or use the Poisson likelihood directly, like so --
dpb=exp(lnpoisson(ib,BG_COUNTS))
;	these IB counts translate to a background intensity in the source area,
xb=ib/float(BG_AREA) & dxb=xb[1:*]-xb
;	for each such background, compute the likelihood that
;	as many as OBS_COUNTS counts are obtained from the background,
;	which is the same as 1-Prob(getting more than OBS_COUNTS)
pp=fltarr(nb) & for i=0,nb-1 do pp[i]=total(exp(lnpoisson(dindgen(10)+obs_counts,xb[i])))
 ; so the p-value is,
pval = pp
 ; weighted for the probability of actually seeing that particular XB,
pval = pval * dpb/dxb
 ; integrating,
print,'probability that the observed counts can be obtained as a fluctuation from the background,'
print,'p=',total(pval*dxb)

 ; and also try the Monte Carlo way
if not keyword_set(nsim) then nsim=1000L
help,nsim
rb=randomu(seed,nsim,poisson=bg_counts) & rxb=rb/float(bg_area) & rpp=dblarr(nsim)
for i=0L,nsim-1L do rpp[i]=total(exp(lnpoisson(dindgen(10)+obs_counts,rxb[i])))
print,'<p>=',mean(rpp)

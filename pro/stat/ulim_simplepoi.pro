function ulim_simplepoi,bkg,alpha,beta,alfa=alfa,Nthr=Nthr,$
	verbose=verbose, _extra=e
;+
;function	ulim_simplepoi
;	compute the upper limit to the intensity of an undetected source,
;	given the background intensity, for a specified false detection
;	probability alpha and detection probability beta, in the Poisson
;	regime.
;
;syntax
;	ul=ulim_simplepoi(bkg,alpha,beta,alfa=alfa,Nthr=Nthr,verbose=v)
;
;parameters
;	bkg	[INPUT; required] background intensity over the source region
;	alpha	[INPUT] maximum false detection probability
;		* if not given, assumed to be 0.0027
;		* if .GE.1, assumed to be given in units of equivalent
;		  Gaussian sigma
;		* if .LT.0, assumed to be given as 1./abs(ALPHA)
;		* the actual false detection probability for the selected
;		  counts threshold may be less than the specified value.
;	beta	[INPUT] probability of detecting the source of given intensity
;		* note that some authors use the notation beta to denote the
;		  Type II error, the probability of a false negative, which is
;		  1-beta in our notation.
;		* if not given, assumed to 0.5
;		* if .LT.0 or .GE.1, assumed to be given as 1.D-1.D/abs(BETA)
;
;keywords
;	alfa	[OUTPUT] ALPHA, as noted above, is the maximum allowed
;		false detection probability.  Because of the discreteness
;		of counts, the actual false det. probability that is used
;		to figure out BETA will be different (though always less
;		than ALPHA).  This computed probability is returned in this
;		array.
;	Nthr	[OUTPUT] the counts detection threshold chosen such that the
;		false detection probability is < ALPHA
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines (all standard IDL)
;	IGAMMA
;	INTERPOL
;	ERRORF
;
;description
;	Kashyap, van Dyk, Connors, Freeman, Siemiginowska, Xu, Zezas, 2010,
;	"On Computing Upper Limits to Source Intensities", ApJ, 719, 900
;
;	This encodes the simple Poisson case (footnotes 7 and 8) where the
;	background intensity is known.  First, ALPHA, the false detection
;	probability, is used to determine a counts threshold for detection,
;	and that source intensity where the source would be detected with
;	probability BETA is returned as the upper limit.
;
;example usage
;	bkg=2 & alpha=-1024.^2 & beta=0.5
;	.run ulim_simplepoi
;
;history
;	Vinay Kashyap (Feb.MMX)
;-

;	usage
ok='ok' & np=n_params()
nbg=n_elements(bkg) & na=n_elements(alpha) & nb=n_elements(beta)
if np eq 0 then ok='Insufficient parameters' else $
 if nbg eq 0 then ok='BKG is not defined' else $
  if nbg gt 1 then begin
    if na gt 1 and na ne nbg then ok='BKG and ALPHA arrays incompatible' else $
    if nb gt 1 and nb ne nbg then ok='BKG and BETA arrays incompatible'
  endif else $
   if na gt 1 then begin
     if nb gt 1 and nb ne na then ok='ALPHA and BETA arrays incompatible'
   endif
if ok ne 'ok' then begin
  print,'Usage: ul=ulim_simplepoi(bkg,alpha,beta,alfa=alfa,Nthr=Nthr,verbose=v)'
  print,'  compute upper limit to intensity of an undetected source,'
  print,'  given background intensity, for specified false detection'
  print,'  probability alpha and detection probability beta, in the'
  print,'  Poisson regime.'
  print,'  (see Kashyap, van Dyk, Connors, Freeman, Siemiginowska, Xu, Zezas, 2010,'
  print,'  "On Computing Upper Limits to Source Intensities", submitted to ApJ)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

;	figure out inputs
nn=nbg > na > nb
if na eq 0 then aa0=0.0027 else aa0=alpha[0]
if nb eq 0 then bb0=0.5 else bb0=beta[0]
if nn eq 1 then begin
  bg=double(bkg[0]) & aa=double(aa0) & bb=double(bb0)
endif else begin
  bg=dblarr(nn) & aa=bg & bb=bg
  if nbg gt 1 then bg=double(bkg) else bg[*]=bkg[0]
  if na gt 1 then aa=double(alpha) else aa[*]=aa0
  if nb gt 1 then bb=double(beta) else bb[*]=bb0
endelse
for i=0L,nn-1L do begin
  if vv gt 10+i*10 then begin
    if aa[i] lt 0 then message,'['+strtrim(i,2)+']'+strtrim(aa[i],2)+$
    	': assuming false detection probability = 1./abs(ALPHA)',/informational
    if aa[i] gt 1 then message,'['+strtrim(i,2)+']'+strtrim(aa[i],2)+$
    	': assuming ALPHA given as Gaussian equivalent sigma',/informational
    if bb[i] lt 0 or bb[i] ge 1 then message,'['+strtrim(i,2)+']'+strtrim(bb[i],2)+$
    	': computing Type II error as 1./abs(BETA)',/informational
	;note: Type II error = 1-BETA in our notation
  endif
  if aa[i] eq 0 then aa[i]=aa0
  if aa[i] lt 0 then aa[i]=1.D/abs(aa[i])
  if aa[i] ge 1 then aa[i]=1.D - errorf(aa[i]/sqrt(2.))
  ;
  if bb[i] eq 0 then bb[i]=bb0
  if bb[i] lt 0 or bb[i] ge 1 then bb[i]=1.D - 1.D/abs(bb[i])
endfor

;	figure out outputs
if nn eq 1 then begin
  alfa=0.D & Nthr=0L & ul=0.D
endif else begin
  alfa=dblarr(nn) & Nthr=lonarr(nn) & ul=dblarr(nn)
endelse

;	loop through each case
if not keyword_set(smin) then smin=0.
if not keyword_set(smax) then smax=100.
if not keyword_set(ds) then ds=0.1
nsbin=long(abs(smax-smin)/ds)+1 & ss=findgen(nsbin)*ds+smin
for i=0L,nn-1L do begin
  bgi=bg[i] & aai=aa[i] & bbi=bb[i]
  aak=1.0 & k=0L
  while aak gt aai do begin	;{compute alpha for different counts thresholds
    k=k+1L
    aak=igamma(k,bgi,/double)
  endwhile			;AAK>AAI}
  alfa[i]=aak & Nthr[i]=k-1L
  bbk=igamma(Nthr[i]+1,ss+bgi,/double) & ssi=smax
  				;compute BETA for various source intensities
  while max(bbk) lt bbi do begin	;{in case max(SS) isn't enough
    ssi=smax+10.*ds & ss=[ss,ssi] & smax=ssi
    bbk=[bbk,igamma(Nthr[i]+1,ssi+bgi,/double)]
    if vv gt 2 then print,'.',form='($,a)'
  endwhile				;MAX(BBK)<BBI}
  ul[i]=interpol(ss,bbk,bbi)	;compute upper limit for given alpha,beta
  ;
  if vv gt 0 then begin
    if vv gt 1 then print,'i,bg[i],alpha[i],beta[i]'
    if vv gt 1 then print,i,bgi,aai,bbi
    print,'aa[i],Nthr[i],ul[i]'
    print,alfa[i],Nthr[i],ul[i]
  endif
  if vv gt 10000 then stop,'HALTing; type .CON to continue'
endfor

if vv gt 1000 then stop,'HALTing; type .CON to continue'

if nn eq 1 then return,ul[0] else return,ul
end
;..............................................................................

;	example calling sequence

print,''
ul=ulim_simplepoi()
print,''

if not keyword_set(bkg) then bkg=2
if not keyword_set(alpha) then alpha=[6d-7,5,-1024.^2]
if not keyword_set(beta) then beta=[0.5d,0.9d,0.99d]
if n_elements(verbose) eq 0 then verbose=51

ul=ulim_simplepoi(bkg,alpha,beta,alfa=alfa,Nthr=Nthr,verbose=verbose)

print,''
print,'upper limit to source intensity = '+strtrim(ul,2)
print,' given that background='+strtrim(bkg,2)
print,' for alpha='+strtrim(alpha,2)+' and beta='+strtrim(beta,2)
print,''

end

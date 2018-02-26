function mcmc_abund,simdem,logt,wvl,flx,fsigma,emis,z,asig,zout,abweight=abweight,$
elogt=elogt,nosol=nosol,olog=olog,crctn_factor=crctn_factor,abund=abund,_extra = e 

;+ 
;function    mcmc_abund 
;        
;        An abundance-independent DEM-reconstruction can be achieved 
;        by running a Monte Carlo Markov-Chain DEM reconstruction
;        mcmc_dem() on a set of He-like to H-like flux ratios. 
;
;        A DEM realized in this fashion together with the observed 
;        line fluxes give implied abundances. This function returns 
;        fractional abundnaces relative to solar 'grevesse & sauval'
;        (not logarithmic values) given a set of MCMC realized DEMS 
;        (mcmc_dem ouput), fluxes, and line details. 
; 
;parameters 
;        simdem [INPUT;required]  an array of size
;		(NT,NSIM+1), with last column containing the best-fit
;	 logt   [INPUT] log_10(Temperature [K]) at which DEMS are given
;        wvl	[INPUT; required] wavelengths [Ang] 
;	 flx	[INPUT; required] flux [ph/s] at each WVL ([ph/s/cm^2]
;        fsigma [INPUT; required] errors for flx
;	 emis	[INPUT; required] line emissivities [1e-23 ergs cm^3/s]
;		* EMIS==EMIS(LOGT,WVL)
;        z      [INPUT; required] atomic numbers of the elements
;               generating the	lines at WVL.
;        asig   [OUTPUT] error in the weighted mean abundance
;        zout   [OUTPUT] atomic numbers corresponding to output abundances
;keywords
;        abund  [INPUT] set  of abundances. see GETABUND(). Output 
;               will be relative to this set. Default is Grevesse &
;               Sauval (1998)
;        abweight [INPUT] if a different weighting is required to compute 
;               abundances then specify abweight as an array of weights of 
;               the same size as flx (e.g. one can use fsigma) 
;        elogt  [INPUT] log_10(Temperature[K]) at which EMIS are
;                       defined if different from logt 
;        nosol  [INPUT] if set, then output abundances will not be
;                       relative to solar 
;        olog   [INPUT] if set, outputs abundances and asig will be
;               logarithmic values.         
;        crctn_factor [INPUT] if set, correctin factor is applied to
;                     theoretical fluxes before calculation of abundnaces
;        _extra [INPUT ONLY] use this to pass defined keywords to subroutines
;           
;subroutines
;        MC_EROR (/home/llin/mc_eror.pro) not PoA standard yet 
;        LINEFLX 
;history 
;        Liwei Lin 5/03 
;           ADDED keyword log for logarithmic output (LL 7/03) 
;           BUG FIX elogt is actually taken into account when set (LL/WB 8/03)
;        LL/CA 9/03 BUGFIX Costanza points out error calculation bug
;           which gives overestimation in errors 
;       BUG FIX elogt handler (call to rebinx) should have logt
;           not dlogt as input (LL/WB 9/03) 
;        LL 01/04 
;           ADDED keyword crctn_factor for use with mixie
;        LL 01/04 
;           H abundance cannot be updated 
;        LL/SC 2/04 
;           BUGFIX HISTOGRAM will give a one element array if given 
;           just one atomic number to bin up.
;        LL ADDED _extra keyword to lineflx 
;        LL 5/05 Change the asig determination to only include the errors 
;           as determined from the DEMs.THIS SHOULD HAVE BEEN FIXED 9/03 AS 
;           stated above. (do not include the error in the 
;           fluxes a second time). 
;           ADDED keyword abweight 
;        LL 6/05 ADDED keyword abund
;           BUGFIX previous error calculation would fail for highly
;           assymetric cases
;        LL/DG 6/05 BUGFIX crash when MC_EROR returns two-element 
;           arrays 
;- 

nl=n_elements(flx)
ntd=n_elements(logt) & nsim=n_elements(SIMDEM[0,*])
em_rbnd=fltarr(ntd,nl)
mod_flxs=fltarr(nl,nsim)
mod_flx_erru = fltarr(nl)
mod_flx_errl = fltarr(nl)
if n_elements(abund) ne 30 then abund = getabund('grevesse & sauval')

; if correction factors set, see if array matches
if keyword_set(crctn_factor) then begin 
  ncf = n_elements(crctn_factor) 
  if ncf ne nl then begin 
    message, 'correction factor array does not match flux array', /info
  return, -1L
  endif 
endif 

; if another weighting scheme is specified, then 
if keyword_set(abweight) then begin 
  nabw = n_elements(abweight)  
  if nabw ne nl then begin 
    message, 'weighting scheme array does not match flux array', /info 
  return, -1L 
  endif
endif

;if emis is on different T grid than DEM then do this
if keyword_set(elogt) then begin 
 for i=0, nl-1 do em_rbnd(*,i) = rebinx(emis(*,i),elogt,logt)
 emis = em_rbnd
endif 

;calculate fluxes for all lines for each dem realization
 for q=0, nsim-1 do begin
    if keyword_set(nosol) then mod_flx = lineflx(emis,logt,wvl,z,DEM=SIMDEM[*,q],/noabund,_extra=e) else $
    mod_flx = lineflx(emis,logt,wvl,z, DEM=SIMDEM[*,q], abund=abund,_extra=e)
    if keyword_set(crctn_factor) then mod_flxs[*,q]=mod_flx/crctn_factor else mod_flxs[*,q]=mod_flx
 endfor 

;get confidence limits by examining distribution of fluxes about best set
 for w=0, nl-1 do begin 
    mc_eror, mod_flx(w), mod_flxs[w,*], bndu,bndl 
    mod_flx_erru[w] =  bndu[0]-mod_flx[w]
    mod_flx_errl[w] =  mod_flx[w]-bndl[0]
 endfor 

;the last iteration should correspond to the best mod_flx
abnd   = flx/mod_flx ; X/H relative to solar if nosol not set 
;abnerru = sqrt([abnd^2] * [(mod_flx_erru/mod_flx)^2])
;abnerrl = sqrt([abnd^2] * [(mod_flx_errl/mod_flx)^2])
mod_flx_err = (mod_flx_erru+mod_flx_errl)/2 
abnsig = mod_flx_err*flx/(mod_flx^2) 

;take the weighted mean for each element
zcount = histogram([0,z])  ; count how many of each element
zcount = zcount[1:*]       ; get rid of dummy

zwhere=where((zcount gt 0)); identify which elements 
zcount = zcount(zwhere)    ; zcount curtailed

wtd_mean_abnd = fltarr(n_elements(zwhere)) 
wtd_mean_abnde= fltarr(n_elements(zwhere))

;determine weights 
if keyword_set(abweight) then wght=abweight else wght = abnsig

for q = 0, n_elements(zwhere)-1 do begin 
nmatches = zcount(q) ; number of lines for this element
zelement = zwhere(q) ; remeber this is IDL indexed i.e. O will be 7  
ii = where(z eq zelement+1)

tmp_abnd  = total(abnd(ii)/(wght(ii)^2))/total(1/(wght(ii)^2))
tmp_abnde = sqrt(1/total(1/(wght(ii)^2)))

wtd_mean_abnd(q)  = tmp_abnd
wtd_mean_abnde(q) = tmp_abnde
endfor

asig = wtd_mean_abnde
zout = zwhere+1
Hw = where(zout eq 1) 
if Hw[0] ne -1L then wtd_mean_abnd(Hw) = 1d

if keyword_set(olog) then begin 
asig = abs(asig/(wtd_mean_abnd*alog(10)))
wtd_mean_abnd=alog10(wtd_mean_abnd)
endif 
return, wtd_mean_abnd
end 

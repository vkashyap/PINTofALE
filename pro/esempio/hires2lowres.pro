;+
; procedure     hires2lowres
;         given filenames containing HETG HEG and MEG PHA data,
;         calculated ARFs, and ACIS ARFs and RMFs, predict what would
;         be observed with ACIS-S or ACIS-I (without gratings) 
;
; syntax
;         make_lowres,hphafile,mphafile,hgarffile,mgarffile,$
;         hhetgwvl,hhetgcts,mhetgwvl,mhetgcts,hetgwvl,hhetgflx,$
;         hhetgflxerr,mhetgflx,mhetgflxerr,tothetgflx,tothetgflxerr,$
;         aciswvl,aciscts,aciserr,hphahdr,mphahdr,hgarfhdr,mgarfhdr,$
;         acisarffile,acisrmffile,debug=debug
;
; parameters
;         hphafile     [INPUT; required] PHA type I spectral file
;                      containing the HEG HETG spectrum for a given
;                      order
;         mphafile     [INPUT; required] PHA type I spectral file
;                      containing the MEG HETG spectrum for a given
;                      order
;         hgarffile    [INPUT; required] Gratings ARF corresponding to
;                      the HEG PHA file
;         mgarffile    [INPUT; required] Gratings ARF corresponding to
;                      the MEG PHA file
;         hhetgwvl     [OUTPUT; required] The midpoints of the HEG
;                      wavelength bins
;                      * units are Angstroms
;         hhetgcts     [OUTPUT; required] The background-subtracted
;                      counts spectrum extracted from the HEG PHA file
;                      * units are [cts/s/AA]
;         mhetgwvl     [OUTPUT; required] The midpoints of the MEG
;                      wavelength bins
;                      * units are Angstroms
;         mhetgcts     [OUTPUT; required] The background-subtracted
;                      counts spectrum extracted from the MEG PHA file
;                      * units are [cts/s/AA]
;         hetgwvl      [OUTPUT; required] The midpoints of the common
;                      HEG, MEG, and HEG+MEG flux spectrum's
;                      wavelength grid 
;                      * grid extends ~0-30 AA at HEG resolution
;         hhetgflx     [OUTPUT; required] Adaptively-smoothed,
;                      background-subtracted HEG flux spectrum
;                      * units are [ph/s/cm^2/AA]
;         hhetgflxerr  [OUTPUT; required] The error on the HEG flux
;                      spectrum
;         mhetgflx     [OUTPUT; required] Adaptively-smoothed,
;                      background-subtracted MEG flux spectrum
;                      * units are [ph/s/cm^2/AA]
;         mhetgflxerr  [OUTPUT; required] The error on the MEG flux
;                      spectrum
;         tothetgflx   [OUTPUT; required] Adaptively-smoothed,
;                      background-subtracted coadded HEG+MEG flux
;                      spectrum
;                      * units are [ph/s/cm^2/AA]
;         tothetgflxerr [OUTPUT; required] The error on the coadded
;                      HEG+MEG flux spectrum
;         aciswvl      [OUTPUT; required] The midpoints of the
;                      predicted ACIS spectrum's wavelength grid
;         aciscts      [OUTPUT; required] The predicted ACIS spectrum
;                      * units are [cts/s/AA]
;         aciserr      [OUTPUT; required] The error on the predicted
;                      ACIS spectrum
;         hphahdr      [OUTPUT; required] The header of the HEG PHA
;                      file
;         mphahdr      [OUTPUT; required] The header of the MEG PHA
;                      file
;         hgarfhdr     [OUTPUT; required] The header of the HEG gARF
;                      file
;         mgarfhdr     [OUTPUT; required] The header of the MEG gARF
;                      file
;         acisarffile  [INPUT; required] The ACIS ARF
;         acisrmffile  [INPUT; required] The ACIS RMF
;
; keywords
;         debug        [INPUT] controls debugging chatter
;
; restrictions
;         requires conv_rmf, smoothie, and rebinw from PINTofALE,
;         available at http://hea-www.harvard.edu/PINTofALE/
;
; history
;         Owen Westbrook, Peter Mendygral, and John Slavin (Mar07)
;-
pro hires2lowres,hphafile,mphafile,hgarffile,mgarffile,$
    hhetgwvl,hhetgcts,mhetgwvl,mhetgcts,hetgwvl,hhetgflx,$
    hhetgflxerr,mhetgflx,mhetgflxerr,tothetgflx,tothetgflxerr,$
    aciswvl,aciscts,aciserr,hphahdr,mphahdr,hgarfhdr,mgarfhdr,$
    acisarffile,acisrmffile,debug=debug

; show the usage if no parameters passed
if n_params(0) eq 0 then begin
   print, ' hires2lowres, hphafile, mphafile, hgarffile, mgarffile, '
   print, ' hhetgwvl, hhetgcts, mhetgwvl, mhetgcts, hetgwvl, hhetgflx, '
   print, ' hhetgflxerr, mhetgflx, mhetgflxerr, tothetgflx, tothetgflxerr, '
   print, ' aciswvl, aciscts, aciserr, hphahdr, mphahdr, hgarfhdr, mgarfhdr, '
   print, ' acisarffile, acisrmffile, debug=debug'
endif

if (not keyword_set(debug)) then debug=0 else debug=1

; read HETG spectral data from HEG pha2 file
hhetgsp    = mrdfits(hphafile,1,hphahdr)
hphabin_lo = hhetgsp.bin_lo			; wavelength (A)
hphabin_hi = hhetgsp.bin_hi			; wavelength (A)
hhetgcts   = hhetgsp.counts			; counts
hhetgbgu   = hhetgsp.background_up		; counts
hhetgbgd   = hhetgsp.background_down		; counts

; read HETG spectral data from MEG pha2 file
mhetgsp    = mrdfits(mphafile,1,mphahdr)
mphabin_lo = mhetgsp.bin_lo			; wavelength (A)
mphabin_hi = mhetgsp.bin_hi			; wavelength (A)
mhetgcts   = mhetgsp.counts			; counts
mhetgbgu   = mhetgsp.background_up		; counts
mhetgbgd   = mhetgsp.background_down		; counts

; HEG background scaling - units are counts
hhetgbg = (hhetgbgu+hhetgbgd)/(sxpar(hphahdr,'BACKSCUP')+sxpar(hphahdr,'BACKSCDN'))

; MEG background scaling - units are counts
mhetgbg = (mhetgbgu+mhetgbgd)/(sxpar(mphahdr,'BACKSCUP')+sxpar(mphahdr,'BACKSCDN'))

; exposure time for observation (s)
exptim = sxpar(hphahdr,'EXPOSURE')		; time (s)
en2wvl = 12.3985				; AA*keV conversion factor

; sort HEG data by wavelength
hwsrt      = sort(hphabin_lo)
hphabin_lo = hphabin_lo[hwsrt]			; wavelength (A)
hphabin_hi = hphabin_hi[hwsrt]			; wavelength (A)
hhetgcts   = hhetgcts[hwsrt]			; counts
hhetgbg    = hhetgbg[hwsrt]			; counts
hhetgwvl   = 0.5*(hphabin_lo+hphabin_hi)	; midpts of wvl grid (A)
hwgrid     = [hphabin_lo,max(hphabin_hi)]	; endpts of wvl grid (A)
hgwvl_sp   = 0.0025                             ; wavelength spacing of HEG grid
; sort MEG data by wavlength
mwsrt      = sort(mphabin_lo)
mphabin_lo = mphabin_lo[mwsrt]			; wavelength (A)
mphabin_hi = mphabin_hi[mwsrt]			; wavelength (A)
mhetgcts   = mhetgcts[mwsrt]			; counts
mhetgbg    = mhetgbg[mwsrt]			; counts
mhetgwvl   = 0.5*(mphabin_lo+mphabin_hi)	; wavelength (A)
mgwvl_sp   = 0.005                              ; wavelength spacing of MEG grid
; create a common wavelength grid from 0-30 AA spaced at HEG resolution
wgrid = hwgrid
maxwvl = 30. < max(mhetgwvl)
while (max(wgrid) le (maxwvl-hgwvl_sp)) do wgrid = [wgrid,max(wgrid)+hgwvl_sp]
wbin_lo = wgrid[0:n_elements(wgrid)-2]
wbin_hi = wgrid[1:n_elements(wgrid)-1]
hetgwvl = 0.5*(wbin_lo+wbin_hi)

; read HETG HEG gARF effective area data
hgarf       = mrdfits(hgarffile,1,hgarfhdr)
hgarfbin_lo = hgarf.BIN_LO			; energy (keV)
hgarfbin_hi = hgarf.BIN_HI			; energy (keV)
hgeffar_h   = hgarf.SPECRESP			; area (cm^2)

; read HETG HEG gARF effective area data
mgarf       = mrdfits(mgarffile,1,mgarfhdr)
mgarfbin_lo = mgarf.BIN_LO			; energy (keV)
mgarfbin_hi = mgarf.BIN_HI			; energy (keV)
mgeffar_h   = mgarf.SPECRESP			; area (cm^2)

; sort HEG gARF data by wavlength
hgwsrt      = sort(hgarfbin_lo)			
hgarfbin_lo = hgarfbin_lo[hgwsrt]		; energy (keV)
hgarfbin_hi = hgarfbin_hi[hgwsrt]		; energy (keV)
hgeffar_h   = hgeffar_h[hgwsrt]			; area (cm^2)
hgwvl_h     = 0.5*(hgarfbin_lo + hgarfbin_hi)	; wavelength (A)

; sort HEG gARF data by wavlength
mgwsrt      = sort(mgarfbin_lo)			
mgarfbin_lo = mgarfbin_lo[mgwsrt]		; energy (keV)
mgarfbin_hi = mgarfbin_hi[mgwsrt]		; energy (keV)
mgeffar_h   = mgeffar_h[mgwsrt]			; area (cm^2)
mgwvl_h     = 0.5*(mgarfbin_lo + mgarfbin_hi)	; wavelength (A)

; interpolate HETG HEG gARF effective area onto HEG spectral grid 
; (or could use hdwvl = hphabin_hi - hphabin_lo)
; units of hhetgeffar are cm^2
hhetgeffar = (interpol(hgeffar_h,hgwvl_h,hhetgwvl) > 0.) < (max(hgeffar_h))
hdwvl      = abs(hhetgwvl - shift(hhetgwvl,1))	; wavelength (A)
hdwvl[0]   = hdwvl[1]				; wavelength (A)

; interpolate HETG MEG gARF effective area onto MEG spectral grid 
; (or could use mdwvl = mphabin_hi - mphabin_lo)
; units of mhetgeffar are cm^2
mhetgeffar = (interpol(mgeffar_h,mgwvl_h,mhetgwvl) > 0.) < (max(mgeffar_h))
mdwvl      = abs(mhetgwvl - shift(mhetgwvl,1))	; wavelength (A)
mdwvl[0]   = mdwvl[1]				; wavelength (A)

; check total counts to make sure smoothie won't reject spectra
sqtoth=sqrt(total(hhetgcts))
sqtotm=sqrt(total(mhetgcts))
sqtothbg=sqrt(total(hhetgbg))
sqtotmbg=sqrt(total(mhetgbg))

; determine the signal to noise threshholds
hstn = 4.
if ((sqtoth lt hstn) or (sqtotm lt hstn)) then begin
   hstn = floor((sqtoth < sqtotm))
endif
bgstn = 5.
if ((sqtothbg lt bgstn) or (sqtotmbg lt bgstn)) then begin
   bgstn = floor((sqtothbg < sqtotmbg))
endif

; smooth the hi-res. source cnts/bin arrays with smoothie
smoothie,hhetgcts,shhetgcts,shc_err,igroup,yerr=sqrt(hhetgcts),snrthr=hstn
smoothie,mhetgcts,smhetgcts,smc_err,igroup,yerr=sqrt(mhetgcts),snrthr=hstn

; smooth the hi-res. background cnts/bin arrays with smoothie
smoothie,hhetgbg,shhetgbg,shbg_err,igroup,yerr=sqrt(hhetgbg),snrthr=bgstn
smoothie,mhetgbg,smhetgbg,smbg_err,igroup,yerr=sqrt(mhetgbg),snrthr=bgstn

; subtract background from smoothed data
shhetgcts = ((shhetgcts - shhetgbg) > 0.)	; counts/bin
smhetgcts = ((smhetgcts - smhetgbg) > 0.)	; counts/bin

; propagate errors for smoothed data
hc_err = sqrt((shc_err^2)+(shbg_err^2))
mc_err = sqrt((smc_err^2)+(smbg_err^2))

; HEG: hhetgcts units are counts bin^-1
hhetgcts = ((hhetgcts-hhetgbg) > 0.)		; counts/bin
hhetgflx = 0.*hhetgcts				; dimensionless

; determine HEG gARF range:
; set effective area cutoff
; limit gARF range to regions where there is actually flux
hnz = where(hhetgcts gt 0.)
hfph = hnz[1]
hlph = hnz[n_elements(hnz)-2]
hnonzero = where(hhetgeffar gt 1.)
hphrng = where((hnonzero ge hfph) and (hnonzero le hlph))
hnonzero = hnonzero[hphrng]

; MEG: hhetgcts units are counts bin^-1
mhetgcts = ((mhetgcts-mhetgbg) > 0.)		; counts/bin
mhetgflx = 0.*mhetgcts				; dimensionless

; determine MEG gARF range:
; set effective area cutoff
; limit gARF range to regions where there is flux
mnz = where(mhetgcts gt 0.)
mfph = mnz[1]
mlph = mnz[n_elements(mnz)-2]
mnonzero = where(mhetgeffar gt 1.)
mphrng = where((mnonzero ge mfph) and (mnonzero le mlph))
mnonzero = mnonzero[mphrng]

; divide by HETG effective area for HEG and MEG
; shhetgflx and smhetgflx units are counts bin^-1 cm^-2
hhetgflx = 0.*shhetgcts
mhetgflx = 0.*smhetgcts
hhetgflx[hnonzero] = shhetgcts[hnonzero]/hhetgeffar[hnonzero]
mhetgflx[mnonzero] = smhetgcts[mnonzero]/mhetgeffar[mnonzero]

; calculate the approximate Poisson error on each wvl bin
hf_err = 100.*hc_err	    			; HEG flux error
mf_err = 100.*mc_err				; MEG flux error
hf_err[hnonzero] = hc_err[hnonzero]/hhetgeffar[hnonzero]
mf_err[mnonzero] = mc_err[mnonzero]/mhetgeffar[mnonzero]
hhetgflxerr = hf_err
mhetgflxerr = mf_err

; calculate where the MEG data overlaps the combined HEG+MEG range
mwvl_crg = where(mhetgwvl ge hetgwvl[n_elements(hetgwvl)-1])
mhetgwvl_c = mhetgwvl[0:mwvl_crg[0]]
smhetgcts_c = smhetgcts[0:mwvl_crg[0]]
mhetgflx_c = mhetgflx[0:mwvl_crg[0]]
mf_err = mf_err[0:mwvl_crg[0]]

; rebin MEG counts and errors onto combined wavelength grid
if (debug eq 1) then print,'DEBUG: total MEG flux before rebinning =',total(mhetgflx_c)
mhetgflx_c = rebinw(mhetgflx_c,mhetgwvl_c,wgrid,/perbin)
mf_err = rebinw(mf_err,mhetgwvl_c,wgrid,/perbin)
if (debug eq 1) then print,'DEBUG: total MEG flux after rebinning =',total(mhetgflx_c)

; get ACIS-S (or ACIS-I or Zero-Order) effective area & wavelength
acisarf = mrdfits(acisarffile,1,aahdr)
enrg_a  = (acisarf.energ_lo + acisarf.energ_hi)/2.	; energy (keV)
effar_a = acisarf.specresp				; area (cm^2)
wvl_a   = en2wvl/enrg_a					; wavelength (A)

; sort effective area & wavelength by wavlength
arfsort = sort(wvl_a)
wvl_a   = wvl_a[arfsort]				; wavelength (A)
effar_a = effar_a[arfsort]				; area (cm^2)

; get ACIS-S/ACIS-I/Zero-Order effective area at HETG HEG resolution
effar_hires = (interpol(effar_a,wvl_a,hetgwvl) > 0.) < (max(effar_a))

; combine the HEG and MEG flux arrays, weighting by their errors
; combine the HEG and MEG errors (if z=a*x+b*y, sigz=sqrt(a2*sigx2+b2*sigy2))
hnum = n_elements(hhetgflx)
tothetgflx = fltarr(n_elements(mhetgflx_c))
tothetgflxerr = 0.*tothetgflx
tothetgflx[0:hnum-1] = (hhetgflx/hf_err+mhetgflx_c[0:hnum-1]/mf_err[0:hnum-1])/(1/hf_err+1/mf_err[0:hnum-1]) ; cts/(bin*cm^2)
tothetgflx[hnum:*] = mhetgflx_c[hnum:*]
a = mf_err[0:hnum-1]/(mf_err[0:hnum-1]+hf_err)	; weighting factor on HEG
b = hf_err/(mf_err[0:hnum-1]+hf_err)		; weighting factor on MEG
tothetgflxerr[0:hnum-1] = sqrt((a^2)*hf_err^2+(b^2)*(mf_err[0:hnum-1])^2); combined HEG+MEG error
tothetgflxerr[hnum:*] = mf_err[hnum:*]

; multiply the combined fluxes and errors by the aimpoint or zero-order effective areas
totaciscts = tothetgflx*effar_hires			;cts/bin
totaciserr = tothetgflxerr*effar_hires

; get ACIS-S/ACIS-I/Zero-Order RMF
rmf = rd_ogip_rmf(acisrmffile)

; convert high-resolution wvl to energy
enrg_h = en2wvl/hetgwvl				; energy (keV)

; convolve combined HEG+MEG hi-res. data and error with response matrix
; input spectrum is in counts/bin
conv_rmf,enrg_h,totaciscts,enrg_lr,phbin_lr,rmf
conv_rmf,enrg_h,totaciserr^2,enrg_lre,phbin_lre,rmf
if (debug eq 1) then print,'DEBUG: total counts before convolving with rmf =',total(totaciscts)
if (debug eq 1) then print,'DEBUG: total counts after convolving with rmf =',total(phbin_lr)

; phbin_lr returned here has units of (cts/bin)
; enrg_lr is bin boundaries -> must convert to midbin values
m_enrg_lr = 0.*phbin_lr
for k=0,n_elements(phbin_lr)-1 do begin
    m_enrg_lr[k] = 0.5*(enrg_lr[k]+enrg_lr[k+1])
endfor
aciserr = sqrt(phbin_lre)
aciswvl = en2wvl/m_enrg_lr				; wavelength (A)
dwlr = abs(aciswvl - shift(aciswvl,1))
dwlr[0] = dwlr[1]

aciscts = phbin_lr/(dwlr*exptim)			; cts/(s*A)
if (debug eq 1) then print,'DEBUG: total(aciscts*dwlr)*exptim =',total(aciscts*dwlr)*exptim
srt = sort(aciswvl)
dwlr = dwlr[srt]

aciswvl = aciswvl[srt]					; wavelength (A)
aciscts = aciscts[srt]					; cts/(s*A)
aciserr = aciserr[srt]					; error (cts/bin)
hhetgcts = hhetgcts/(hdwvl*exptim)			; cts/(s*A)
mhetgcts = mhetgcts/(mdwvl*exptim)			; cts/(s*A)

return
end

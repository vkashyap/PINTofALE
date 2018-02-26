;+
; procedure     make_lowres
;         wrapper script to run hires2lowres and put the simulated
;         ACIS spectral data into a ASCII tables and fits files
;
; syntax
;         make_lowres,hphafile,mphafile,hgarffile,mgarffile,$
;         hhetgwvl,hhetgcts,mhetgwvl,mhetgcts,hetgwvl,hhetgflx,$
;         hhetgflxerr,mhetgflx,mhetgflxerr,tothetgflx,tothetgflxerr,$
;         aciswvl,aciscts,aciserr,hphahdr,mphahdr,hgarfhdr,mgarfhdr,$
;         acisarffile,acisrmffile,debug=debug,txtfile=txtfile,$
;         fitsfile=fitsfile,flxfits=flxfits,flxtxt=flxtxt
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
;         txtfile      [INPUT; required] Output text file for the HETG
;                      and predicted ACIS spectra
;         fitsfile     [INPUT; required] Output FITS file for the HETG
;                      and predicted ACIS spectra
;         flxfits      [INPUT; required] Output FITS file for the HEG,
;                      MEG, and combined HEG+MEG flux spectra
;         flxtxt       [INPUT; required] Output text file for the HEG,
;                      MEG, and combined HEG+MEG flux spectra
;
; restrictions
;         requires conv_rmf, smoothie, and rebinw from PINTofALE,
;         available at http://hea-www.harvard.edu/PINTofALE/
;
; history
;         Owen Westbrook, Peter Mendygral, and John Slavin (Mar07)
;-
pro make_lowres,hphafile,mphafile,hgarffile,mgarffile,$
    hhetgwvl,hhetgcts,mhetgwvl,mhetgcts,hetgwvl,hhetgflx,$
    hhetgflxerr,mhetgflx,mhetgflxerr,tothetgflx,tothetgflxerr,$
    aciswvl,aciscts,aciserr,hphahdr,mphahdr,hgarfhdr,mgarfhdr,$
    acisarffile,acisrmffile,debug=debug,txtfile=txtfile,$
    fitsfile=fitsfile,flxfits=flxfits,flxtxt=flxtxt

; show the usage if no parameters passed
if n_params(0) eq 0 then begin
   print, ' make_lowres, hphafile, mphafile, hgarffile, mgarffile, '
   print, ' hhetgwvl, hhetgcts, mhetgwvl, mhetgcts, hetgwvl, hhetgflx, '
   print, ' hhetgflxerr, mhetgflx, mhetgflxerr, tothetgflx, tothetgflxerr, '
   print, ' aciswvl, aciscts, aciserr, hphahdr, mphahdr, hgarfhdr, mgarfhdr, '
   print, ' acisarffile, acisrmffile, debug=debug, txtfile=txtfile, '
   print, ' fitsfile=fitsfile, flxfits=flxfits, flxtxt=flxtxt '
endif

; run the hires2lowres routine 
hires2lowres,hphafile,mphafile,hgarffile,mgarffile,$
    hhetgwvl,hhetgcts,mhetgwvl,mhetgcts,hetgwvl,hhetgflx,$
    hhetgflxerr,mhetgflx,mhetgflxerr,tothetgflx,tothetgflxerr,$
    aciswvl,aciscts,aciserr,hphahdr,mphahdr,hgarfhdr,mgarfhdr,$
    acisarffile,acisrmffile,debug=debug

; open the output file for writing the hi and low-res counts data
openw,unit,txtfile,/get_lun

; print the column names
printf,unit,'   Hi-Res. HEG data          Hi-Res. MEG data     Sim. Low-Res. data'
printf,unit,'  wavelength    flux       wavelength    flux     wavelength    flux    Error'
printf,unit,'     (A)      (cts/s/A)       (A)      (cts/s/A)     (A)      (cts/s/A) cts/bin'

; loop through the elements of the arrays returned by hires2lowres and output the data
nhhi = n_elements(hhetgwvl)
nmhi = n_elements(mhetgwvl)
nhi = nhhi < nmhi
nlo = n_elements(aciswvl)

for i = 0, nhi - 1 do begin
    if (i lt nlo) then printf,unit,format='(F9.4,E15.4,F9.4,E15.4,F9.4,E15.4,F9.4)',hhetgwvl[i],hhetgcts[i], $
        mhetgwvl[i],mhetgcts[i],aciswvl[i],aciscts[i],aciserr[i] $
    else printf,unit,format='(F9.4,E15.4,F9.4,E15.4)',hhetgwvl[i],hhetgcts[i],mhetgwvl[i],mhetgcts[i]
endfor
free_lun,unit

data=create_struct('HEG_WVL',hhetgwvl,'HEG_CTS',hhetgcts,'MEG_WVL',mhetgwvl,$
   'MEG_CTS',mhetgcts,'LOW_WVL',aciswvl,'LOW_CTS',aciscts,'LOW_ERR',aciserr)
fits_header = hphahdr
sxaddpar, fits_header, 'TFIELDS', 7
sxaddpar, fits_header, 'CREATOR', 'make_lowres.pro', $
  'tool that created this output'
mwrfits,data,fitsfile,fits_header,/CREATE

; write the hi-res. flux data & errors
openw,unit2,flxtxt,/get_lun

; print the column names
printf,unit2,'           HETG HEG+MEG data			HETG HEG data				  HETG MEG data
printf,unit2,'    wvl         flux	        error	   wvl	       flux           error        wvl         flux           error'
printf,unit2,'    (A)     (ph/s/A/cm2)                   (A)     (ph/s/A/cm2)                    (A)    (ph/s/A/cm2)'

; loop through the flux arrays returned by hires2lowres & output the data
nhflx = n_elements(hhetgflx)
nmflx = n_elements(mhetgflx)
ntflx = n_elements(tothetgflx)

for j = 0, ntflx - 1 do begin
   if (j lt nhflx) then printf,unit2,format='(F9.4,E15.4,E15.4,F9.4,E15.4,E15.4,F9.4,E15.4,E15.4)',hetgwvl[j],tothetgflx[j],tothetgflxerr[j],hhetgwvl[j],hhetgflx[j],hhetgflxerr[j],mhetgwvl[j],mhetgflx[j],mhetgflxerr[j] $ 
   else printf,unit2,format='(F9.4,E15.4,E15.4)',hetgwvl[j],tothetgflx[j],tothetgflxerr[j]
endfor

data2 = create_struct('HEG_WVL',hhetgwvl,'HEG_FLX',hhetgflx,'HEG_ERR',hhetgflxerr,'MEG_WVL',mhetgwvl,'MEG_FLX',mhetgflx,'MEG_ERR',mhetgflxerr,'TOT_WVL',hetgwvl,'TOT_FLX',tothetgflx,'TOT_ERR',tothetgflxerr)
header2 = hphahdr
sxaddpar, header2, 'TFIELDS', 9
sxaddpar, header2, 'CREATOR', 'make_lowres.pro', $
   'tool that created this output'
mwrfits,data2,flxfits,header2,/CREATE

end

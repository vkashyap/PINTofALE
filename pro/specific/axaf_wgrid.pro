function axaf_wgrid,wmin,leg=leg,meg=meg,heg=heg,$
	detlen=detlen,PixSize=PixSize,OverSampl=OverSampl,$
	TG_period=TG_period,X_Rowland=X_Rowland,MaxOrd=MaxOrd,$
	WavCtr=WavCtr,WavBinSize=WavBinSize
;+
;function	axaf_wgrid
;	return wavelength grid of bin-beginning values and the final bin-ending
;	value in which the sampling follows the optimal binning for AXAF
;	transmission grating spectra to accurately reflect higher resolution
;	in higher orders.
;
;syntax
;	ww=axaf_wgrid(wmin,/LEG,/MEG,/HEG,detlen=detlen,PixSize=PixSize,$
;	OverSampl=OverSampl,TG_period=TG_period,X_Rowland=X_Rowland,$
;	MaxOrd=MaxOrd,WavCtr=WavCtr,WavBinSize=WavBinSize)
;
;parameters
;	wmin	[INPUT; required] minimum wavelength for output array [Ang]
;
;keywords
;	leg		[INPUT] if set, sets all keyword defaults to the
;			Low-Energy Transmission Grating specific numbers.
;			* default
;	meg		[INPUT] if set, sets all keyword defaults to the
;			Medium-Energy Transmission Grating specific numbers.
;			* overrides LEG
;	heg		[INPUT] if set, sets all keyword defaults to the
;			High-Energy Transmission Grating specific numbers.
;			* overrides MEG and LEG
;	detlen		[INPUT] Approximate length of detector array,
;			from zero-order to outer limit [mm]
;	PixSize		[INPUT] Pixel size of the array (not necessarily a
;			resolution element) [mm]
;	OverSampl	[INPUT] Amount of oversampling of pixels of size
;			PIXSIZE (e.g., if dithering gives half-pixel
;			resolution, then specify OverSampl of 2); value
;			is factor to divide PixSize by.
;	TG_period	[INPUT] Period of the diffraction grating [Ang]
;	X_Rowland	[INPUT] Rowland diameter of the grating [mm]
;			(canonically 8635mm, but larger at XRCF)
;	MaxOrd		[INPUT] Maximum order to consider
;	WavCtr		[OUTPUT] Array of bin central wavelengths [Ang]
;	WavBinSize	[OUTPUT] Array of bin full widths, corresponding
;			to the WAVCTR bins [Ang]
;
;description
;	In order to accurately model the higher diffraction orders of a
;	predicted spectrum, the input model spectrum must have sufficient
;	resolution to be sampled by the detector.  A grid which matches the
;	first order spectrum will be undersampled in higher orders, since
;	higher orders are increasingly dispersed by a factor of m, the
;	diffraction order.  If the grid were set to match the highest order
;	to be considered, then the input spectrum would have a very large
;	number of grid points.  This is unnecessary, since much of the
;	spectrum is diffracted distances greater than the length of the
;	detector.  This procedure constructs a grid which has increasingly
;	finer resolution at progressively shorter wavelengths.  The grid
;	spacing is discontinuous at boundaries where the next higher order
;	will fall off the detector.
;
;	Gridding is as follows:
;	* From DetLen to DetLen/2, grid at the first order resolution
;	  (essentially the detector pixel size, possibly oversampled a bit).
;	* From DetLen/2 to DetLen/3, grid at half the first order.
;	* From DetLen/3 to DetLen/4, grid at one-third the first order.
;	* ... and so forth.
;
;	The general form is:
;	  n_m = the number of grid points for order m,
;	  n = the number of grid points required to span from wavelength=0
;	    (zero-order centroid) to DetLen at the first order resolution. 
;	  n_m = [ 1/m - 1/(m+1) ] * m  * n  =   n / (m+1) 
;	Thus the total number of grid-points required is: 
;	  N = sum_{m=1..MaxOrd} { [1/m - 1/(m+1)] * m  * n} 
;	    = n * sum_{m=1..MaxOrd} {1 / (m+1)} 
;	    = n * { Psi(2+MaxOrd) - 1 + gamma}, 
;	where Psi(x) is the digamma function, and gamma is Euler's constant.
;
;	This is a diverging but slowly increasing function.  For MaxOrd of
;	10, 20, and 30, N is 2.0*n, 2.6*n, and 3.0*n, respectively.  (If the
;	spectrum were gridded at the highest resolution required, N would
;	scale directly with MaxOrd.)
;
;	Given some reasonable MaxOrd (such as 30 for LETGS), the minimum
;	wavelength grid point calculated as outlined above will still be
;	above the minimum needed in model spectra.  Thus, we have added a
;	minimum wavelength as a parameter, and the grid will be filled from
;	the minimum grid point to the specified minimum at the gridding of
;	the MaxOrd segment.
;
;history
;	v1.0 as WAVE_GRID.PRO by Dave Huenemoerder (2 december 1995)
;	%=====================================================================
;	% Time-stamp: <97/04/24 10:00:48 dph>
;	% MIT Directory: ~dph/h1/Model-Spec/pro
;	% CfA Directory: /dev/null
;	% Mac Directory: :Bruno:science:Model-Spec:pro
;	% File: line_spectrum.pro
;	% Author: D. Huenemoerder
;	% original version: 951202
;	%=====================================================================
;	modified to make it SCAR-friendly, changed WMIN to mean minimum of
;	  bin-beginning values rather than minimum of bin-center values
;	  (VK;Jan98)
;	corrected spelling mistake w. OVERSAMPL (VK;Jun99)
;	corrected spike at WAVE_MIN at end of MAXORD (VK; Jun99)
;	changed default value of X_ROWLAND (VK; FebMM)
;-

;	usage
nn=n_elements(wmin)
if nn ne 1 then begin
  print,'Usage: ww=axaf_wgrid(wmin,/LEG,/MEG,/HEG,detlen=detlen,PixSize=PixSize,$'
  print,'       OverSampl=OverSampl,TG_period=TG_period,X_Rowland=X_Rowland,$'
  print,'       MaxOrd=MaxOrd,WavCtr=WavCtr,WavBinSize=WavBinSize)'
  print,'  return wavelength grid for AXAF transmission gratings'
  return,-1L
endif

;	set keyword defaults
;if keyword_set(leg) then begin		;(LETG is default
  dlen=160.0	;[mm] approx length of one side of HRC-S
  pix=0.0064	;[mm] set by electronic readout
  osamp=1.0	;HRC readout oversamples the FWHM of the charge cloud (0.025 mm)
  tg_per=9921.0	;[Ang] LETG baseline period
  x_row=8635	;[mm] Rowland diameter, flight config
  x_row=8632.48	;[mm] DPH 16FebMM
  mxord=30	;should be plenty -- strong lines in 27th order can contribute
  		;a few percent of counts at some positions
;endif					;LETG)
if keyword_set(meg) then begin		;(METG
  dlen=80.0	;[mm]
  pix=24e-3	;[mm]
  osamp=2.0	;maybe change to 4?
  tg_per=4001.41	;[Ang]
  x_row=8633.69	;[mm] acc. to Obs Guide
  x_row=8632.48	;[mm] DPH 16FebMM
  mxord=7
endif					;METG)
if keyword_set(heg) then begin		;(HETG
  dlen=80.0	;[mm]
  pix=24e-3	;[mm]
  osamp=4.0	;maybe change to 2?
  tg_per=2000.81	;[Ang]
  x_row=8633.69	;[mm] acc. to Obs Guide
  x_row=8632.48	;[mm] DPH 16FebMM
  mxord=7
endif					;HETG)

;	override defaults
if keyword_set(detlen) then dlen=detlen
if keyword_set(pixsize) then pix=pixsize
if keyword_set(oversampl) then osamp=oversampl
if keyword_set(tg_period) then tg_per=tg_period
if keyword_set(x_rowland) then x_row=x_rowland
if keyword_set(maxord) then mxord=maxord

;; determine the number of grid points for full first order:
;
   dlen = float(dlen)             ; make sure it's real
   n_first_order = long(dLen / pix * osamp)

;; determine maximum wavelength - wavelength at dLen:
;
   wave_max = tg_per * dLen / x_row

; Construct the first-order part of the grid, from dLen/2 (maximum
; wavelength/2) to dLen (maximum wavelength)
   
   this_m = 1
   
   wave_min =  wave_max / (this_m + 1)
   delta_wave = wave_max / (this_m * (this_m + 1))

   n_this_m = n_first_order / (this_m + 1)
   grid_center = dindgen(n_this_m)/n_this_m * delta_wave + wave_min
   grid_bins = replicate(delta_wave/n_this_m, n_this_m)

;; form and concatenate grids up to highest order specified
;
   FOR this_m = 2, MxOrd DO BEGIN
      
      wave_min =  wave_max / (this_m + 1)
      delta_wave = wave_max / (this_m * (this_m + 1))
      n_this_m = n_first_order / (this_m + 1)

      grid_center = [dindgen(n_this_m)/n_this_m * delta_wave +$
	wave_min, grid_center]
      grid_bins = [replicate(delta_wave/n_this_m, n_this_m), grid_bins]
      
   ENDFOR

;; extend down to the minimum wavelength requested
;
   IF wmin LT wave_min THEN BEGIN
      
      delta_wave = wave_min-wmin
      n_left = delta_wave / grid_bins(0)
      
      ;{VK:	this was the original version, and caused spikes
      ;		at the end of WAVE_MIN
      ;grid_center = [dindgen(n_left)/n_left * delta_wave + wmin, $
      ;               grid_center]
      ;}{VK:	and this is the altered version, which "goes the
      ;		other way" and shouldn't have any spikes.
      dbin=delta_wave/long(n_left+1)
      grid_center=[reverse(dindgen(long(n_left+1))*(-dbin)+wave_min-dbin),$
	grid_center]
      ;VK: NOTE: same problem occurs at all other WAVE_MINs -- correct later}
      grid_bins = [replicate(grid_bins(0), n_left), grid_bins]
      
   ENDIF
   
   WavCtr = grid_center
   WavBinSize = grid_bins

nw=n_elements(WavCtr)

;WW=[WavCtr-0.5*WavBinSize,WavCtr(nw-1)+0.5*WavBinSize(nw-1)]
WW=[WavCtr(0),WavCtr+WavBinSize]

return,ww
end

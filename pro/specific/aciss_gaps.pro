pro aciss_gaps,hgap,mgap,lgap,onchip=onchip,offset=offset,order=order,$
	ikeV=ikeV,mm=mm,pix=pix, _extra=e
;+
;procedure	aciss_gaps
;	calculate the wavelength (or energy) locations of the CXO chip gaps
;
;	the outputs are all in the following sequence:
;	[left edge of S0, right edge of S0,
;	 left edge of S1, right edge of S1,
;	 left edge of S2, right edge of S2,
;	 left edge of S3, right edge of S3,
;	 left edge of S4, right edge of S4,
;	 left edge of S5, right edge of S5]
;
;syntax
;	aciss_gaps,hgap,mgap,lgap,onchip=onchip,offset=offset,$
;		order=order,/ikeV,/mm,/pix
;
;parameters
;	hgap	[OUTPUT] HEG gaps (ACIS-S/HETG) [Ang]
;	mgap	[OUTPUT] MEG gaps (ACIS-S/HETG) [Ang]
;	lgap	[OUTPUT] LEG gaps (ACIS-S/LETG) [Ang]
;
;keywords
;	onchip	[OUTPUT] the number of the chip on which the 0th order
;		position lies.
;		* 0=>S0, 1=>S1, .. 5=>S5
;		* 0.5=>in gap between S0 & S1, etc.
;	offset	[INPUT (default=0)] offset of source position from
;		nominal aim point [arcmin]
;		* nominal offset is 6.04-0.024 mm from first active column
;		* +ve is towards center of S3
;	order	[INPUT (default=1] grating order
;	ikeV	[INPUT] if set, outputs are in [keV] (default is [Ang])
;	mm	[INPUT] if set, OFFSET is assumed to be in [mm]
;	pix	[INPUT] if set, OFFSET is assumed to be in [pix]
;		* PIX overrides MM
;	_ref_extra	[JUNK] here only to prevent crashing the program
;
;history
;	9/9/97 dd Written for Fred Seward for PG
;	11/27/99 dd Updated values for NRA2
;	MIM.XII -- vk copied from
;	  http://space.mit.edu/HETG/technotes/chip_gaps/hetgs_gaps.pro
;	  and modified the interface to be more useful to acisgarf.pro
;	changed value of Rowland diameter (VK; FebMM)
;	changed keyword KEV to IKEV (VK; JanMMI)
;	changed focal length and Rowland diameter (VK; Sep01)
;	changed leg_ang to -0.07 (VK; Jan04)
;-

;	usage
if n_params() eq 0 then begin
  print,'Usage: aciss_gaps,hgap,mgap,lgap,onchip=onchip,$'
  print,'       offset=offset,order=order,/ikeV,/mm,/pix'
  print,'  compute grating chip gaps for ACIS-S'
  return
endif

;	keywords
noff=n_elements(offset)
if noff gt 1 then message,'can only handle one offset value at a time',/info
if noff ne 0 then offs=offset(0) else offs=0.0
;
keV=0 & if keyword_set(ikeV) then ikeV=1	;output in [keV]
imm=0 & if keyword_set(mm) then imm=1		;OFFSET in [mm]
ipix=0 & if keyword_set(pix) then ipix=1	;OFFSET in [pix]
if ipix eq 1 then imm=0			;PIX overrides MM
;
gord=1 & if keyword_set(order) then gord=abs(fix(order(0))) > 1

;	{the following is used as is from original code

; ACIS
; everybody knows this:
; * Change chip size to be 1022 pixels since first and
;   last columns don't report events and are part of the gap!
;   Increase the gap size by these two pixels also...
chip_size = 1022.0 * 0.024  ; mm

; From ICD diagram:
; gap_size = 0.43 ; mm
;
; Measured TDET Y S0-S5 corrections:
;   -4.97, -3.58, -2.42, 0.0, -0.08, 0.33
; and corresponding gap size changes:
;    1.39, 1.16, 2.42, 0.0, -0.08, 0.41
; * Make gap one extra pixel wide for TDET corr of
;   -3.0,  -2.0,  -1.0,  0.0,  1.0,  2.0
; so gap errors are:
;   -1.97  -1.58  -1.42  0.  0.92  1.67
; --> all within 2 pixels
;
; total of 3 pixels extra gap size:
gap_size = 0.43 + 3.0*0.024; mm

; From diagram in ICD
nom_offset = 6.04         ; mm from S3 edge (first column)
nom_offset = 6.04 - 0.024 ; mm from S3's first active column (#2)

; HETG
; values from HETG CIP file HETGperiod.rdb
heg_p = 2000.81      ; A
heg_ang = -5.235      ; degrees
meg_p = 4001.41      ; A
meg_ang = 4.725       ; degrees
; include LEG values as well
leg_p = 9912.16      ; A
leg_ang = -0.04   ; degrees; very approximate
leg_ang = -0.07   ; degrees; still very approximate (1/04)

; Flight Rowland spacing
rc_spacing = 8632.65 ; mm	nominal design value
rc_spacing = 8632.48 ; mm	DPH 16FebMM
rc_spacing = 8635.50 ; mm	Jeremy told me on May 11, 2000
rc_spacing = 8637.00 ; mm	acc. to POG 3.0, p226

; HRMA
; new value...
focal_len = 10061.62 ; mm
focal_len = 10070.00 ; mm	as of Sep 2001

; hc
hc = 12.3985
; ......................................

;	VK's modifications begin from here on down}

;	compute offset in [arcmin], [mm], and [pix]
if imm eq 0 and ipix eq 0 then begin	;input in [arcmin]
  offs_arcmin=offs
  offs_mm=focal_len*!DTOR*(offs/60.0)
  offs_pix=long(offs_mm/0.024)
endif
if imm eq 1 then begin			;input in [mm]
  offs_mm=offs
  offs_pix=long(offs_mm/0.024)
  offs_arcmin=60.*(offs_mm/focal_len)/!DTOR
endif
if ipix eq 1 then begin			;input in [pix]
  offs_pix=offs
  offs_mm=offs_pix*0.024
  offs_arcmin=60.*(offs_mm/focal_len)/!DTOR
endif

; Arrays for the chip edges w.r.t. zero-order in mm
low_edge = ([0.,1.,2.,3.,4.,5.]-3.0) * (chip_size + gap_size) - $
	(nom_offset + offs_mm)
high_edge = chip_size + low_edge

; on which chip does the source lie?
onchip=-0.5
if low_edge(0) ge 0 then onchip=-0.5	;off scale
if low_edge(0) lt 0 and high_edge(0) ge 0 then onchip=0.	;S0
if high_edge(0) lt 0 and low_edge(1) ge 0 then onchip=0.5
if low_edge(1) lt 0 and high_edge(1) ge 0 then onchip=1.	;S1
if high_edge(1) lt 0 and low_edge(2) ge 0 then onchip=1.5
if low_edge(2) lt 0 and high_edge(2) ge 0 then onchip=2.	;S2
if high_edge(2) lt 0 and low_edge(3) ge 0 then onchip=2.5
if low_edge(3) lt 0 and high_edge(3) ge 0 then onchip=3.	;S3 usually
if high_edge(3) lt 0 and low_edge(4) ge 0 then onchip=3.5
if low_edge(4) lt 0 and high_edge(4) ge 0 then onchip=4.	;S4
if high_edge(4) lt 0 and low_edge(5) ge 0 then onchip=4.5
if low_edge(5) lt 0 and high_edge(5) ge 0 then onchip=5.	;S5
if high_edge(5) lt 0 then onchip=5.5	;off scale

print,'Calculating edges for offset of:'
print,	string(offs_arcmin,'(f8.2)')+' [arcmin] == '+$
	string(offs_mm,'(f8.2)')+' [mm] == '+$
	string(long(offs_pix),'(i7)')+' [pix]'+$
	'	==> Chip S'+strtrim(string(onchip,'(f4.1)'),2)

; Make a list of gaps
gap_low = [-1.E6, high_edge]
gap_high = [low_edge, 1.E6]

; Convert these to Energy for meg and heg
;   energy = grating order * hc / ( period * dispersion_distance/rc_spacing )
;    and dispersion_distance = gap_distance/COS(dispersion_angle)

eheg_low = gord*hc/( heg_p * (ABS(gap_low)/COS(!DTOR*heg_ang))/rc_spacing  )
eheg_high = gord*hc/( heg_p * (ABS(gap_high)/COS(!DTOR*heg_ang))/rc_spacing  )

emeg_low = gord*hc/( meg_p * (ABS(gap_low)/COS(!DTOR*meg_ang))/rc_spacing  )
emeg_high = gord*hc/( meg_p * (ABS(gap_high)/COS(!DTOR*meg_ang))/rc_spacing  )

; and LETG...
eleg_low = gord*hc/( leg_p * (ABS(gap_low)/COS(!DTOR*leg_ang))/rc_spacing  )
eleg_high = gord*hc/( leg_p * (ABS(gap_high)/COS(!DTOR*leg_ang))/rc_spacing  )

;	convert to [Ang]
if keV eq 0 then begin
  eheg_low=hc/eheg_low & eheg_high=hc/eheg_high
  emeg_low=hc/emeg_low & emeg_high=hc/emeg_high
  eleg_low=hc/eleg_low & eleg_high=hc/eleg_high
endif

;	form the outputs
hgap=fltarr(12) & mgap=hgap & lgap=hgap
for i=0,5 do begin
  hgap(2*i)=eheg_high(i) & hgap(2*i+1)=eheg_low(i+1)
  mgap(2*i)=emeg_high(i) & mgap(2*i+1)=emeg_low(i+1)
  lgap(2*i)=eleg_high(i) & lgap(2*i+1)=eleg_low(i+1)
endfor

return
end

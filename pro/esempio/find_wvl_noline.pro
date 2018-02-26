;+
;find_wvl_noline.pro
;	find wavelength bins which are line-free in a given spectrum
;
;	reads in the line emissivity database, makes a frequency histogram
;	of the lines in each wavelength bin, and identifies bins where the
;	number of lines is 0.
;
;requires
;	WGRID == a wavelength grid
;		(if not given, a standard MEG grid of 1:41.96:0.005 is assumed)
;	LDBDIR == location of line emissivity database
;		(if not given, assumed to be $APED)
;	LSFWDT == the expected spread in the line.  This is basically the
;		window half-width of a boxcar smoothing function.
;		If integer, assumed to be in pixels, and if float, taken to be
;		in units of bin width and converted to pixels accordingly.
;		(if not given, assumed to be 2)
;
;vinay k (2011may26)
;-

;	initialize
peasecolr & loadct,3 & peasecolr
;
if n_elements(wgrid) lt 2 then begin
  mwvl=8192L & if n_elements(nwvl) ge 1 then mwvl=nwvl[0]
  wvlmin=1.0 & if n_elements(wmin) ge 1 then wvlmin=abs(wmin[0])
  wvlmax=41.96 & if n_elements(wmax) ge 1 then wvlmax=abs(wmax[0])
  dwvl=(wvlmax-wvlmin)/mwvl & wgrid=findgen(mwvl+1L)*dwvl+wvlmin
endif
wvlmin=min(wgrid,max=wvlmax) & dwvl=median(wgrid[1:*]-wgrid)
;
if n_elements(ldbdir) eq 0 then begin
  ldbdir='$APED'
  if not keyword_set(n_e) then n_e=1e9
  if n_tags(lstr) eq 0 then lemis=rd_line(atom,n_e=n_e,wrange=[wvlmin,wvlmax],dbdir=ldbdir,fstr=lstr,verbose=verbose)
  ;if strpos(strupcase(ldbdir),'APED') lt 0 then begin
  ;  lemis=fold_ioneq(lstr.LINE_INT,lstr.Z,lstr.JON,eqfile=ioneqf,chifil=chifil,chidir=chidir,verbose=verbose)
  ;  lstr.LINE_INT=lemis
  ;endif
endif
;
hwidth=2
if n_elements(lsfwdt) ne 0 then begin
  slsf=size(lsfwdt,/type)
  if slsf lt 4 or slsf eq 12 then hwidth=lsfwdt[0]
  if slsf eq 4 or slsf eq 5 then hwidth=long(lsfwdt[0]/dwvl+0.5)>1
endif
;
if not keyword_set(linethresh) then linethresh=0.

hw=histogram(abs(lstr.WVL),min=wvlmin,bin=dwvl) & nhw=n_elements(hw) & xhw=findgen(nhw)*dwvl+wvlmin
shw=smooth(float(hw),2*hwidth+1L)
og=lindgen(nhw) & for k=1L,nhw-1L do if hw[k] eq hw[k-1] and shw[k] le linethresh then og[k]=og[k-1L]
dog=og[1:*]-og & odog=where(dog gt 0) & og[odog+1L]=-1
ocont=where(og ge 0,mocont)

plot,xhw,hw,psym=10,xtitle='wavelength',ytitle='# lines',title='lines in each wavelength bin',/xs
if mocont gt 0 then begin
  oplot,xhw[ocont],hw[ocont],psym=1,col=2
  xyouts,!x.crange[0],!y.crange[1]-0.05*(!y.crange[1]-!y.crange[0]),'line free bins are marked with +',col=2,align=-0.1
endif
print,''
message,'WARNING: bins identified as line-free are conditional on the completeness of the atomic line database',/informational

end

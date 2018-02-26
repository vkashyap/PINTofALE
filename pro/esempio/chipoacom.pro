;+ 
;script chipoacom
;
;       Compares emissivities in PoA database to emissivities
;       calculated in real time with CHIANTI.
; 
;	The CHIANTI calculation is done with CH_SYNTHETIC.PRO,
;	a major component of CHIANTI's flagship spectral-synthesis
;	program CH_SS.PRO.  CH_SYNTHETIC is used to calculate G(T),
;	after which PoA's LINEFLX() is used to calculate intensities.
;	These intensities are compared to intensities generated
;	using PoA's RD_LINE(), FOLD_IONEQ(), and LINEFLX().
;	Comparisons of the binned spectra are plotted.
;	VOORSMOOTH() is used to smooth binned spectra.
;
;syntax
;        .run chipoacom
; 
;inputs 
;	!IONEQF	ion-balance file
;	!LDBDIR	atomic database to use default =!LDBDIR 
;	!EDENS	electron density
;	V_EM	emission measure components corresponding 
;	        to each T [cm^5/logK]
;	        default is [6.1d11, 6.1d11, 7.1e11]         
;	V_LOGT	temperatures at which emission measures are defined 
;	        default is [6.1, 6.8, 7.2] 
;	V_WMIN	minimum wavelength default = 5.0
;	V_WMAX	maximum wavelength default = 15.0
;	V_NWBIN	number of bins in output spectra default = 8000
;	V_Np	number of plots with which the resulting spectra 
;	        are displayed. The wavelength range in each plot
;	        is determined by signal strength. 
;	V_SMOOT	smoothing scale used to smooth binned spectra 
;outputs 
;	V_OUTPS	name of ps output file 
;        
;how to use 
;        1. initialize PINTofALE using INITALE
;        2. parameters set with initale are: 
;           -- ATOMIC: !LDBDIR, !IONEQF, !CHIDIR
;           -- STELLAR: !EDENS
;        3. set the wavelength range of interest 
;              V_WMIN = 5.0 
;              V_WMAX = 180.0
;        4. set the emission measures and corresponding temperatures
;        5. set the output filename if output ps is desired
;              V_OUTPS ='CHI4PoAcomparison.ps'
;        6. run this script
;
;history 
;	Liwei Lin (Jun03)  
;	bug fix noclip=0 and polyfill could result in 'inverse' polyfill
;	  (LL; Jul03)
;	bug fix using floating point precission when defining wgrid
;	  is inadequate and results in incongruous specrtal binning
;	  results via hastogram(). use double precission instead (LL; Jul03)
;	couple of bug fixes (VK; Apr04)
;       add check for CHIANTI and IDL version compatiblity (LL; Dec05) 
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;       bug fix extraneous /2 after call to lineflx (LL; Aug08)  
;- 

if not keyword_set(V_WMIN)  then V_WMIN = 5.0 
if not keyword_set(V_WMAX)  then V_WMAX = 15.0
if not keyword_set(V_LOGT)  then V_LOGT = [6.1, 6.8, 7.2]
if not is_keyword_set(V_EM) then V_EM   = [6.1d11, 6.1d11, 7.1e11] 
if not keyword_set(V_SMOOT) then V_SMOOT = 0.02
if not keyword_set(V_Np)    then V_Np = 5
if not keyword_set(V_NWBIN) then V_NWBIN = 8000

defsysv,'!PoA',exists=ivar
if ivar eq 0 then message,'Must initialize PINTofALE first via INITALE'

;check versions of IDL and CHIANTI to ensure compatibility 
;(CHIANTI 5.0,5.1 are not compatible with IDL versions older than 5.6) 
; because of the use of the keyword DIMENSION in built in function MAX() 
idl_ver = !version.release 
;lifted from use_chianti() v8:
ff = findfile(concat_dir(!xuvtop,'VERSION'))

if ff(0) ne '' then begin 
   openr, 1, ff(0) 
   readf, 1, version & close, 1 
   if version ge 5.1 and idl_ver lt 5.6 then begin
      print, 'WARNING: IDL 5.6 or greater might be required for this CHIANTI version.' 
   stop
   endif 
endif  

;read line-cooling emissivites & fold in ion equilibrium
 lconf = rd_line(atom, n_e = !EDENS, wrange = [V_WMIN,V_WMAX]$
          ,wvl=LWVL,logT=LLOGT,Z=Z, ion=ION,jon=JON,fstr=lstr) 
 lconf=  fold_ioneq(lconf,Z,JON,$
          logT=LLOGT,eqfile=!IONEQF,verbose=!VERBOSE)

;construct dem and compute intensities
 lwvl = abs(lwvl);theoretical lines are annotated w/ negative wavelengths

 !DEM=mk_dem('delta', logT = !LOGT, pardem=V_LOGT, indem=V_EM)
 int   = lineflx(lconf, !LOGT, lwvl, z,dem=!DEM)
 
;compute contribution functions  using CHIANTI's CH_SYNTHETIC
 IONEQ_NAME = !CHIDIR+!IONEQF ;get full path to default IONEQ file

 ch_synthetic,V_WMIN,V_WMAX, output = synthout,$
             pressure = 1e15, density = !EDENS, goft = 1,$
             IONEQ_NAME=IONEQ_NAME, /all

;construct and compute CHIANTI intensities
 CHI_DEM=mk_dem('delta',logT= synthout.ioneq_logt, pardem=V_LOGT, indem=V_EM)

 chi_int = lineflx(4*!pi*synthout.lines.goft/1e-23, synthout.ioneq_logt, $ 
 synthout.lines.wvl, synthout.lines.iz, dem=CHI_DEM)

;commented out by VK (Apr04)
;;filter out weak lines which which may be included in PoA and not
;;CHIANTI or vice versa which can cause large oscillations
;
; cutof = (1e-5)*max(int)
; aa = where(int lt cutof) 
; if total(aa) gt 0 then int(aa) = 1e-24
; ee = where(chi_int lt cutof) 
; if total(ee) gt 0 then chi_int(ee) = 1e-24
; ;aa=where(int lt max(int)/1e5,maa)
; ;if maa gt 0 then int(aa)=max(int)/1e10
; ;ee=where(chi_int lt max(int)/1e5,mee)
; ;if mee gt 0 then chi_int(ee)=max(int)/1e10

;bin the two intensity sets onto a common grid
 dwvl=double((V_WMAX-V_WMIN)/V_NWBIN)
 wgrid=(dindgen(V_NWBIN+1L)*dwvl+V_WMIN)
 cwvl = float(synthout.lines.wvl)
 chispc0 = hastogram(cwvl,wgrid,wts=chi_int)
 poaspc0 = hastogram(lwvl, wgrid, wts = int)
 
 chispc = voorsmooth(chispc0,V_SMOOT,wgrid)
 poaspc = voorsmooth(poaspc0,V_SMOOT,wgrid)
 chispc = chispc >1e-24
 poaspc = poaspc >1e-24

 !p.multi = [0,1,2]
if keyword_set(V_OUTPS) then begin 
 set_plot, 'PS'
 device, filename = string(V_OUTPS), /color
endif else set_plot, 'x'

 w=wgrid & y = chispc & Np = V_Np
 oW=sort(W) & F=Y[oW]-min(Y) & X=W[oW] & Xb=mid2bound(X)
 NX=n_elements(X) & CF=fltarr(NX+1L)

 for i=1L,Nx do CF[i]=CF[i-1L]+F[i-1L] & CF=CF/max(CF)
 o0=max(where(CF eq 0)) & o1=min(where(CF eq 1)) & CF=CF[o0:o1] & Xb=Xb[o0:o1]
 xx=interpol(Xb,CF,findgen(Np+1L)/Np)
 wrange=fltarr(2,Np)
 wrange[0,*]=xx[0:Np-1L]
 wrange[1,*]=xx[1:*]+(xx[1:*]-xx[0:Np-1L])*0.1 < max(Xb)

for ii=0,Np-1 do begin 
 oo  = where((wgrid ge wrange[0,ii]) and (wgrid le wrange[1,ii]))
 if total(oo) gt 0 then begin 
  npoaspc = poaspc(oo)/max([chispc(oo),poaspc(oo)])
  nchispc = chispc(oo)/max([chispc(oo),poaspc(oo)])
  nwgrid  = wgrid(oo)

  plot, nwgrid, npoaspc, title = 'PoA(blue) & CHIANTI(red)',xtickname = replicate(' ',10), $
    /nodata, ytitle='Normalized Flux',xtitle=' ', xrange=wrange[*,ii],ymargin = [0,4], xstyle = 1, $
    yrange=[0,1], ystyle=1
  oplot, nwgrid, nchispc
  ;NOTE: oo should not be necessary given noclip = 0 but a nasty
  ;glitch with polyfill and noclip results in 'inverse' polyfilling'
  bottom = fltarr(n_elements(nwgrid))
  polyfill,[nwgrid,reverse(nwgrid)],[npoaspc,bottom],col=1, noclip=0
    oplot, nwgrid, nchispc, thick = 3, color = 2

  plot,wgrid, (poaspc-chispc)/(poaspc+chispc),xrange=wrange[*,ii], $
    xtitle = 'Angstroms', ytitle = '(PoA-Chianti)/(PoA+Chianti)',ymargin = [4,0], yrange=[-1,1], xstyle = 1, /nodata
    oplot, wgrid, (poaspc-chispc)/(poaspc+chispc)
    oplot, wrange[*,ii],[0,0], color = 3, thick =1, noclip=0
  wait,1 
  plot,wgrid, (poaspc-chispc)/poaspc,xrange=wrange[*,ii], $
    xtitle = 'Angstroms', ytitle = '(PoA-Chianti)/PoA',ymargin = [4,0], yrange=[-1,1], xstyle = 1, /nodata
    oplot, wgrid, (poaspc-chispc)/poaspc
    oplot, wrange[*,ii],[0,0], color = 3, thick =1, noclip=0
  plot,wgrid, (poaspc-chispc)/max(poaspc),xrange=wrange[*,ii], $
    xtitle = 'Angstroms', ytitle = '(PoA-Chianti)/max(PoA)',ymargin = [4,0], yrange=[-1,1], xstyle = 1, /nodata
    oplot, wgrid, (poaspc-chispc)/max(poaspc)
    oplot, wrange[*,ii],[0,0], color = 3, thick =1, noclip=0
  wait,1 
 endif
endfor

if keyword_set(V_OUTPS) then begin 
 device, /close_file
 set_plot, 'x'
endif
 !p.multi = 0

end

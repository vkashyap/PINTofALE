;+
;script	eg_haarline
;	a script to demonstrate the use of the wavelet based
;	line detection program HAARLINE
;
;usage
;	pha2fil='/data/snafu/kashyap/Capella/o1248_pha2.fits'
;	pharow=1
;	.run eg_haarline
;
;vinay kashyap (Dec'02)
;-

;	read in spectrum
if not keyword_set(pha2fil) then $
	pha2fil='/data/snafu/kashyap/Capella/o1248_pha2.fits'
if n_elements(pharow) eq 0 then pharow=1
if n_tags(pha) eq 0 then begin
  pha=mrdfits(pha2fil,1,h)
  bgscl=sxpar(h,'BACKSCUP')+sxpar(h,'BACKSCDN')
endif
ct=pha[pharow].COUNTS
wlo=pha[pharow].BIN_LO & whi=pha[pharow].BIN_HI & wvl=0.5*(wlo+whi)
bgu=pha[pharow].BACKGROUND_UP & bgd=pha[pharow].BACKGROUND_DOWN
y=ct-(bgu+bgd)/bgscl
ysig=sqrt((sqrt(ct+0.75)+1.)^2+(sqrt(bgu+bgd+0.75)+1)^2/bgscl^2)

;	call HAARLINE
xpos=haarline(y,scales,xerr=xerr,thrct=10,sclout=sclout,$
	sclmin=2,sclmax=64,ysig=ysig,thrsig=1,thrloc=1,verbose=10,$
	hy=hy,wy=wy)

;	convert indices to wavelength grid
wpos=interpolate(wvl,xpos-1) & werrp=interpolate(wvl,xpos+xerr-1) & werrm=interpolate(wvl,xpos-xerr-1)
cpos=interpol(ct,wvl,wpos)

;	display in triumph
window,xsize=1200,ysize=1000 & peasecolr & pmulti=!p.multi & !p.multi=[0,1,2]
plot,wvl,ct,xr=[12,18],/xs,psym=10,thick=2,xtitle='Wavelenght [Ang]',ytitle='[ct]',subtitle=pha2fil
for i=0L,n_elements(wpos)-1L do oplot,wpos[i]*[1,1],cpos[i]*[0,1],thick=2,color=2
for i=0L,n_elements(wpos)-1L do oplot,[werrm[i],werrp[i]],cpos[i]*[1,1],thick=2,color=3
plot,wvl,ct,xr=[6,12],/xs,psym=10,thick=2,xtitle='Wavelenght [Ang]',ytitle='[ct]',subtitle=pha2fil
for i=0L,n_elements(wpos)-1L do oplot,wpos[i]*[1,1],cpos[i]*[0,1],thick=2,color=2
for i=0L,n_elements(wpos)-1L do oplot,[werrm[i],werrp[i]],cpos[i]*[1,1],thick=2,color=3
!p.multi=pmulti

end

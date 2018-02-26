function vradimage,vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$
	vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$
	dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$
	imgwdt=imgwdt,imgpix=imgpix,imgout=imgout,imgroot=imgroot,$
	imglog=imglog,r2p2b=r2p2b,verbose=verbose, _extra=e
;+
;function	vradimage
;	wrapper for RADPROJECT(), takes a given radial vector field (like,
;	say velocity), computes the projections onto the YZ plane at \infty,
;	bins up the projections along X to return a velocity profile, and
;	makes a projected image of the velocity field.  returns a structure
;	of the form
;	{velocity grid, velocities histogram, mid-bins of grid,
;	 weighted velocity profile, velocity images, particle weights,
;	 RGB image, RGB image of blue-shifts, RGB image of red-shifts}
;
;syntax
;	vstr=vradimage(vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$
;	vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$
;	dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$
;	imgwdt=imgwdt,imgpix=imgpix,imgout=imgout,imgroot=imgroot,/imglog,/r2p2b,$'
;	verbose=verbose, inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$
;	opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0)
;
;parameters
;	vfunct	[INPUT; required] a string denoting what kind of radial
;		vector (velocity) profile to use.  current options are:
;		- BETA : v = VINF*(1-BCONST*(RPEG/R))^VBETA
;		         BCONST = 1-(VPEG/VINF)^(1/VBETA)
;		- CONST : v = VINF
;		- USER : v = execute(UVFUNCT)
;	rfunct	[INPUT; required] a string denoting what kind of radial
;		number density weighting to use.  current options are:
;		- EXP : n = NPEG * exp(-(R/RSCL))
;		- INVSQ : n = NPEG * (RSCL/R)^2
;		- CONST : n = NPEG
;		- USER : n = execute(URFUNCT)
;
;keywords
;	latrng	[INPUT; default=[-90,90]] range of latitudes to consider [deg]
;	lonrng	[INPUT; default=[0,360]] range of longitudes to consider [deg]
;		* (LAT,LON)=(0,0) projects to (Y,Z)=(0,0)
;		* if only one element is given, then the extremum value
;		  nearest to the given value is reset
;	dlat	[INPUT; default=1.] step size to go through latitudes [deg]
;		* longitude steps are corrected according to cos(lat)
;		* ideally, to avoid gridding artifacts in the projected images,
;		  what we should have is that dlat should be small enough to
;		  seamlessly cover the image pixels set via IMGPIX and IMGWDT,
;		  and to avoid double counting rays, this should also match
;		  the opening angle ARCONE of the rays.  Because these combine
;		  the input to different routines, they are not automatically
;		  computed here.  But here is a useful scaling guide:
;		  first compute the maximum recommended subtended angle of
;		  image pixels at MAXRAD,
;		  	DLAT=((IMGWDT/IMGPIX)/MAXRAD)*(180./!pi)
;		  For small angles, the pencil beam thickness ARCONE=DLAT^2
;	vinf	[INPUT; default=1000] value of radial vector field at \infty
;		* e.g., if velocity, 1000 km/s
;	vpeg	[INPUT; default=0] value of radial vector field at RPEG
;	vbeta	[INPUT; default=1] index for radial vector beta profile
;	rpeg	[INPUT; default=1] radius value at which VPEG is specified
;	npeg	[INPUT; default=1e9] value of density weighting at default
;		location
;		* e.g., if number density, 1e9 cm^-3
;	rscl	[INPUT; default=1] radius parameter that defines the radial
;		density function
;	vmin	[INPUT] minimum value of projected vector
;	vmax	[INPUT] maximum value of projected vector
;		* default is to take them from output of RADPROJECT()
;	dvel	[INPUT] projected vector value binning width
;		* default is to set it to produce 101 bins in [VMIN,VMAX]
;	uvfunct	[INPUT] functional form f(RR) to use if VFUNCT='USER'
;	urfunct	[INPUT] functional form f(RR) to use if RFUNCT='USER'
;		* examples of these user defined functions are --
;		  	uvfunct='VINF/(1+(RR-RPEG)^2/RPEG^2)'
;		  	urfunct='NPEG^2*exp(-2*RR/RSCL)'
;		  	urfunct='NPEG*RR'
;		  etc. Goes without saying that this is a powerful tool
;		  and can hurt badly if misused.  Please make sure you
;		  understand how this works before using it.
;	imgwdt	[INPUT; default=10] the range of projected Y,Z over which
;		to construct the velocity image
;		* covers the range +-IMGWDT in both Y and Z
;	imgpix	[INPUT; default=128] number of pixels along each axis of the image
;	imgout	[OUTPUT] the weighted velocity image, of size (3,IMGPIX,IMGPIX)
;		with the first slice being the total, the second the blue-shifted,
;		and the third the red-shifted
;	imgroot	[INPUT; default='vradimage'] root name of output jpg file
;		* three files are written out, one for the full range of
;		  the velocities, and one each for just the red-shifted and
;		  blue-shifted components
;	imglog	[INPUT] if set, returns RGB images in log scale of intensities
;		* be warned that this can highlight small numerical irregularities
;	r2p2b	[INPUT] decides how to let hue go from red to green
;		via green (clockwise, the default, R2P2B=0), or
;		via pink (anti-clockwise, /R2P2B)
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to RADPROJECT()
;		INRAD,MAXRAD,RRNG,RDEL,OPAQ,ARCONE,THETA0,PHI0
;
;example usage
;	.run vradimage
;
;history
;	vinay kashyap (2015dec; supersedes vradprofile)
;-

;	usage
ok='ok' & np=n_params()
nv=n_elements(vfunct) & sv=size(vfunct,/type)
nr=n_elements(rfunct) & sr=size(rfunct,/type)
if np lt 2 then ok='Insufficient parameters' else $
 if nv eq 0 then ok='VFUNCT is undefined' else $
  if nr eq 0 then ok='RFUNCT is undefined' else $
   if sv ne 7 then ok='VFUNCT must be a string constant' else $
    if sr ne 7 then ok='RFUNCT must be a string constant' else $
     if nv gt 1 then ok='VFUNCT must be a scalar' else $
      if nr gt 1 then ok='RFUNCT must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: vstr=vradimage(vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$'
  print,'       vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$'
  print,'       dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$'
  print,'       imgwdt=imgwdt,imgpix=imgpix,imgout=imgout,imgroot=imgroot,/imglog,/r2p2b,$'
  print,'       verbose=verbose, inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$'
  print,'       opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0)'
  print,'  wrapper for radproject, returns binned values of'
  print,'  radial vector field projected onto plane at \infty'
  print,'  and greyscale and color images of the projected image'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	figure out VFUNCT and RFUNCT
;note: uvf and urf are used only for user defined functions, but are
;constructed here as illustrative examples.  The built-in functions
;are explictly called at runtime, without going through execute, for
;speed reasons.

vf='unk'
if strpos(strupcase(vfunct[0]),'BETA') ge 0 then vf='beta'
if strpos(strupcase(vfunct[0]),'CONST') ge 0 then vf='const'
if strpos(strupcase(vfunct[0]),'USER') ge 0 then vf='user'
if vf eq 'unk' then begin & message,vfunct[0]+': not understood; returning',/informational & return,-1L & endif
if vf eq 'beta' then begin
  pegv=0. & if keyword_set(vpeg) then pegv=float(vpeg[0])
  infv=1e3 & if keyword_set(vinf) then infv=float(vinf[0])
  pegR=1. & if keyword_set(rpeg) then pegR=float(rpeg[0])
  betav=2. & if keyword_set(vbeta) then betav=vbeta[0]
  bconst = 1.D - (pegv/infv)^(1./betav) & uvf='infv*(1.-bconst*(pegr/rr))^(betav)'
endif
if vf eq 'const' then begin
  infv=1e3 & if keyword_set(vinf) then infv=float(vinf[0])
  uvf='infv+0.*rr'
endif
if vf eq 'user' then begin
  if not keyword_set(uvfunct) then begin & message,"set UVFUNCT='f(RR)'",/informational & return,-1L & endif
  uvf=uvfunct[0]
endif

rf='unk'
if strpos(strupcase(rfunct[0]),'CONST') ge 0 then rf='const'
if strpos(strupcase(rfunct[0]),'INVSQ') ge 0 then rf='invsq'
if strpos(strupcase(rfunct[0]),'EXP') ge 0 then rf='exp'
if strpos(strupcase(rfunct[0]),'USER') ge 0 then rf='user'
if rf eq 'unk' then begin & message,rfunct[0]+': not understood; returning',/informational & return,-1L & endif
if rf eq 'const' then begin
  pegN=1e9 & if keyword_set(npeg) then pegN=float(npeg[0])
  urf='pegN+0.*rr'
endif
if rf eq 'invsq' then begin
  pegN=1e9 & if keyword_set(npeg) then pegN=float(npeg[0])
  sclR=1. & if keyword_set(rscl) then sclr=float(rscl[0])
  urf='pegN*(sclr/rr)^2'
endif
if rf eq 'exp' then begin
  pegN=1e9 & if keyword_set(npeg) then pegN=float(npeg[0])
  sclR=1. & if keyword_set(rscl) then sclr=float(rscl[0])
  urf='pegN*exp(-rr/sclr)'
endif
if rf eq 'user' then begin
  if not keyword_set(urfunct) then begin & message,"set URFUNCT='f(RR)'",/informational & return,-1L & endif
  urf=urfunct[0]
endif

;	figure out image for output
iwdt=10. & if keyword_set(imgwdt) then iwdt=abs(imgwdt[0]) > 0.1
npix=128L & if keyword_set(imgpix) then npix=long(abs(imgpix[0])) > 1L
xminy=-iwdt & xmaxy=iwdt & yminz=-iwdt & ymaxz=iwdt & xbiny=2.*iwdt/(npix-1.) & ybinz=xbiny
ox=findgen(npix)*xbiny+xminy & oy=findgen(npix)*ybinz+yminz
imgout=fltarr(3,npix,npix) & imgrgb=bytarr(3,npix,npix) & imghue=imgout & imgval=imgout
imgdenom=lonarr(3,npix,npix) & imgwt=fltarr(3,npix,npix)
iroot='vradimage' & if keyword_set(imgroot) then iroot=strtrim(imgroot[0],2)

;	now let's see what other inputs are needed before we begin
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
bmin=-90. & bmax=90.
if keyword_set(latrng) then begin
  if n_elements(latrng) eq 1 then begin
    dmin=latrng[0]-bmin & dmax=bmax-latrng[0]
    if dmin lt dmax then bmin=latrng[0] else bmax=latrng[0]
  endif
  if n_elements(latrng) eq 2 then begin & bmin=min(latrng) & bmax=max(latrng) & endif
  if vv gt 0 then begin
    if bmin lt -90 then message,'resetting minimum latitude to -90 deg',/informational
    if bmax gt 90 then message,'resetting maximum latitude to +90 deg',/informational
  endif
  bmin = bmin > (-90.)
  bmax = bmax < (90.)
endif
;
lmin=0. & lmax=360.
if keyword_set(lonrng) then begin
  if n_elements(lonrng) eq 1 then begin
    dmin=lonrng[0]-lmin & dmax=lmax-lonrng[0]
    if dmin lt dmax then lmin=lonrng[0] else lmax=lonrng[0]
  endif
  if n_elements(lonrng) eq 2 then begin & lmin=lonrng[0] & lmax=lonrng[1] & endif
  if vv gt 0 then begin
    if lmin lt -360 then message,'resetting minimum longitude',/informational
    if lmax gt 90 then message,'resetting maximum latitude',/informational
  endif
  if lmin lt -360 then lmin=(lmin mod 360)
  if lmax gt 360 then lmax=(lmax mod 360)
endif
db=1.0 & if keyword_set(dlat) then db=abs(dlat[0])

nb=long((bmax-bmin)/db+0.5)>1 & bb=findgen(nb)*db+bmin
zvarr=0 & zwarr=0 & zyarr=0 & zzarr=0
kilroy,dot='['+strtrim(nb,2)+'] '
for ib=0L,nb-1L do begin			;{for all latitudes
  xb=bb[ib] & dbcos=abs(cos(xb*!pi/180.))
  if dbcos ne 0 then begin
    dl=db/dbcos < (lmax-lmin)
  endif else dl=(lmax-lmin)
  nl=long((lmax-lmin)/dl+0.5)>1 & ll=findgen(nl)*dl+lmin
  if vv gt 0 then kilroy,dot=strtrim(ib,2)+'('+strtrim(nl,2)+')'

  zvarrb=0 & zwarrb=0 & zyarrb=0 & zzarrb=0
  for il=0L,nl-1L do begin			;{for all longitudes
    xl=ll[il]
    rstr=radproject(xb,xl,verbose=vv, _extra=e)
    rr=rstr.R & pp=rstr.P & vol=rstr.VOL & yy=rstr.Y & zz=rstr.Z
    o0=where(pp ne 0,mo0)
    if mo0 gt 0 then begin
      rr=rr[o0] & pp=pp[o0] & vol=vol[o0]

      case vf of
        'user': jnk=execute('zv='+uvf)
        'beta': zv=infv*(1.-bconst*(pegr/rr))^(betav)
        'const': zv=infv+0.*rr
      endcase
      case rf of
        'user': jnk=execute('zw='+urf)
        'const': zw=pegN+0.*rr
        'invsq': zw=pegN*(sclr/rr)^2
        'exp': zw=pegN*exp(-rr/sclr)
      endcase

      zvp=-zv*pp 	;projected velocities
      ;	-ve sign to ensure approaching velocities (blue shifts) have -ve velocities
      zvw=zw*vol	;weighting
      if not keyword_set(zvarrb) then begin
        zvarrb=zvp & zwarrb=zvw & zyarrb=yy & zzarrb=zz
      endif else begin
        zvarrb=[zvarrb,zvp] & zwarrb=[zwarrb,zvw] & zyarrb=[zyarrb,yy] & zzarrb=[zzarrb,zz]
      endelse

    endif
  endfor					;IL=0,NL-1}

  if not keyword_set(zvarr) then begin
    zvarr=zvarrb & zwarr=zwarrb ;& zyarr=zyarrb & zzarr=zzarrb
  endif else begin
    zvarr=[zvarr,zvarrb] & zwarr=[zwarr,zwarrb] ;& zyarr=[zyarr,zyarrb] & zzarr=[zzarr,zzarrb]
  endelse

  ;	the following is an obfuscated way of computing the weighted velocities in each image pixel
  ;it is obfuscated because it is way faster than the naive code, which would look like
  ;  for i=0L,npix-1L do begin & for j=0L,npix-1L do begin & ok=where(ix eq i and iy eq j,mok) ... & endfor & endfor
  ;or even the slightly optimized version
  ;  oz=where(tmpimg gt 0,moz) & ii=array_indices(tmpimg,oz) & for k=0,moz-1 do begin & i=ii[0,k] & j=ii[1,k] & ok=where(ix eq i and iy eq j,mok) & if mok gt 0 then begin & yzv=zvarrb[ok] & ... & okb=where(yzv le 0,mokb) & if mokb gt 0 then begin & ... & endif & okr=where(yzv gt 0,mokr) & ... & endif & endfor

  ix=long((zyarrb-xminy)/xbiny) & iy=long((zzarrb-yminz)/ybinz)
  ;
  o0=where(zyarrb ge xminy and zyarrb lt xmaxy and zzarrb ge yminz and zzarrb lt ymaxz)
  yzy0=zyarrb[o0] & yzz0=zzarrb[o0] & yzv0=zvarrb[o0] & yzw0=zwarrb[o0]
  htmp=reform(histogram(ix[o0]+npix*iy[o0],min=0,max=npix*npix-1L,bin=1,reverse_indices=ri),npix,npix) & oz=where(htmp gt 0,moz)
  ;	NOTE: this is exactly identical to
  ;	alt_htmp=hist_2d(zyarrb[o0],zzarrb[o0],min1=xminy,min2=yminz,max1=xmaxy,max2=ymaxz,bin1=xbiny,bin2=ybinz)
  ;	but the 1D version gets the reverse indices, which is what gives the speed bump
  ii=array_indices(htmp,oz)
  for k=0L,moz-1L do begin
    jk=oz[k] & ok=ri[ri[jk]:ri[jk+1L]-1L] & mok=ri[jk+1L]-ri[jk] & yzv=yzv0[ok] & yzw=yzw0[ok] & i=ii[0,k] & j=ii[1,k]
    imgout[0,i,j]=imgout[0,i,j]+total(yzv*yzw) & imgwt[0,i,j]=imgwt[0,i,j]+total(yzw) & imgdenom[0,i,j]=imgdenom[0,i,j]+mok
  endfor
  if vv gt 10 then window,0 & loadct,0 & tvscl,rebin(reform(imgout[0,*,*]),npix*3,npix*3)
  ;
  o0b=where(zyarrb ge xminy and zyarrb lt xmaxy and zzarrb ge yminz and zzarrb lt ymaxz and zvarrb le 0,mo0b)
  if mo0b gt 0 then begin
    yzy0=zyarrb[o0b] & yzz0=zzarrb[o0b] & yzv0=zvarrb[o0b] & yzw0=zwarrb[o0b]
    htmpb=reform(histogram(ix[o0b]+npix*iy[o0b],min=0,max=npix*npix-1L,bin=1,reverse_indices=ri),npix,npix) & ozb=where(htmpb gt 0,mozb)
    ii=array_indices(htmpb,ozb)
    for k=0L,mozb-1L do begin
      jk=ozb[k] & ok=ri[ri[jk]:ri[jk+1L]-1L] & mok=ri[jk+1L]-ri[jk] & yzv=yzv0[ok] & yzw=yzw0[ok] & i=ii[0,k] & j=ii[1,k]
      imgout[1,i,j]=imgout[1,i,j]+total(yzv*yzw) & imgwt[1,i,j]=imgwt[1,i,j]+total(yzw) & imgdenom[1,i,j]=imgdenom[1,i,j]+mok
    endfor
    if vv gt 10 then window,1 & loadct,1 & tvscl,rebin(reform(imgout[1,*,*]),npix*3,npix*3)
  endif
  ;
  o0r=where(zyarrb ge xminy and zyarrb lt xmaxy and zzarrb ge yminz and zzarrb lt ymaxz and zvarrb gt 0,mo0r)
  if mo0r gt 0 then begin
    yzy0=zyarrb[o0r] & yzz0=zzarrb[o0r] & yzv0=zvarrb[o0r] & yzw0=zwarrb[o0r]
    htmpr=reform(histogram(ix[o0r]+npix*iy[o0r],min=0,max=npix*npix-1L,bin=1,reverse_indices=ri),npix,npix) & ozr=where(htmpr gt 0,mozr)
    tmpimg=hist_2d(zyarrb[o0],zzarrb[o0],min1=xminy,min2=yminz,max1=xmaxy-xbiny,max2=ymaxz-ybinz,bin1=xbiny,bin2=ybinz)
    ii=array_indices(htmpr,ozr)
    for k=0L,mozr-1L do begin
      jk=ozr[k] & ok=ri[ri[jk]:ri[jk+1L]-1L] & mok=ri[jk+1L]-ri[jk] & yzv=yzv0[ok] & yzw=yzw0[ok] & i=ii[0,k] & j=ii[1,k]
      imgout[2,i,j]=imgout[2,i,j]+total(yzv*yzw) & imgwt[2,i,j]=imgwt[2,i,j]+total(yzw) & imgdenom[2,i,j]=imgdenom[2,i,j]+mok
    endfor
    if vv gt 10 then window,2 & loadct,3 & tvscl,rebin(reform(imgout[2,*,*]),npix*3,npix*3)
  endif

  if vv gt 1000 then stop,il,' of ',nl,ib,' of ',nb

endfor						;IB=0,NB-1}

;	and now make the velocity histogram 
minv=min(zvarr,max=maxv,/nan)
if keyword_set(vmin) then minv=vmin
if keyword_set(vmax) then maxv=vmax
dv=(maxv-minv)/100. & if keyword_set(dvel) then dv=abs(dvel[0])
hv=histogram(zvarr,min=minv,max=maxv,binsize=dv,reverse_indices=ri)
nhv=n_elements(hv) & vgrid=findgen(nhv+1L)*dv+minv
vmid=findgen(nhv)*dv+minv+0.5*dv
;	and the weighted velocity profile
vpp=dblarr(nhv)
for i=0L,nhv-1L do begin
  if ri[i] ne ri[i+1] then begin
    ok=ri[ri[i]:ri[i+1]-1]
    vpp[i]=total(zwarr[ok])
  endif
endfor

;iout=imgout	;& imgout=imgout/(imgdenom>1)	;no need to do this because ARCONE is supposed to handle the number of particles correctly

;	compute the Hue
imghue=fltarr(npix,npix) & imghuer=imghue & imghueb=imghue
tmp11=reform(imgout[1,*,*]) & tmp21=reform(imgwt[1,*,*]) & ok1=where(tmp21 gt 0,mok1)
if mok1 gt 0 then begin
  if keyword_set(r2p2b) then imghueb[ok1]=((tmp11[ok1]/tmp21[ok1])/vinf)*60.+300. else imghueb[ok1]=abs((tmp11[ok1]/tmp21[ok1])/vinf)*120.+120.
endif
tmp12=reform(imgout[2,*,*]) & tmp22=reform(imgwt[2,*,*]) & ok2=where(tmp22 gt 0,mok2)
if mok2 gt 0 then begin
  if keyword_set(r2p2b) then imghuer[ok2]=((tmp12[ok2]/tmp22[ok2])/vinf)*60.+300. else imghuer[ok2]=((tmp12[ok2]/tmp22[ok2])/vinf)*120.
endif
tmp10=reform(imgout[0,*,*]) & tmp20=reform(imgwt[0,*,*]) & ok0=where(tmp20 gt 0,mok0)
if mok0 gt 0 then begin
  if keyword_set(r2p2b) then imghue[ok0]=((tmp10[ok0]/tmp20[ok0])/vinf)*60.+300. else imghue[ok0]=-((tmp10[ok0]/tmp20[ok0])/vinf)*120.+120.
endif

;	saturation is useless, set to 1
imgsat=0.*imghue+1.

;	set the brightness based on how much material there is along each pixel
imgval=reform(imgwt[0,*,*]) & ok=where(imgval gt 0,mok) & if keyword_set(imglog) then imgval[ok]=alog10(imgval[ok]) & imgval[ok]=imgval[ok]/max(imgval[ok])
imgvalb=reform(imgwt[1,*,*]) & ok=where(imgvalb gt 0,mok) & if keyword_set(imglog) then imgvalb[ok]=alog10(imgvalb[ok]) & imgvalb[ok]=imgvalb[ok]/max(imgvalb[ok])
imgvalr=reform(imgwt[2,*,*]) & ok=where(imgvalr gt 0,mok) & if keyword_set(imglog) then imgvalr[ok]=alog10(imgvalr[ok]) & imgvalr[ok]=imgvalr[ok]/max(imgvalr[ok])

;	convert to RGB
imgrgb=reform(hsv2rgb(imghue,imgsat,imgval),npix,npix,3)
bimgrgb=reform(hsv2rgb(imghueb,imgsat,imgvalb),npix,npix,3)
rimgrgb=reform(hsv2rgb(imghuer,imgsat,imgvalr),npix,npix,3)

;	write to files
write_jpeg,iroot+'.jpg',imgrgb,true=3
write_jpeg,iroot+'_blue.jpg',bimgrgb,true=3
write_jpeg,iroot+'_red.jpg',rimgrgb,true=3
spawn,'ls -l '+iroot+'*.jpg'

;	output structure
vstr=create_struct('VGRID',vgrid,'VHIST',hv,'VMID',vmid,'VPROFILE',vpp,'VIMAGE',imgout,'VWEIGHT',imgwt,$
	'RGBIMAGE',imgrgb,'BLUERGB',bimgrgb,'REDRGB',rimgrgb)

if vv gt 900 then stop,'HALTing; type .CON to continue'

return,vstr
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example calling sequence

peasecolr & loadct,3 & peasecolr

if not keyword_set(vfunct) then vfunct='BETA'
	if not keyword_set(vinf) then vinf=1e3
	if not keyword_set(vpeg) then vpeg=0.
	if not keyword_set(rpeg) then rpeg=1.
	if not keyword_set(vbeta) then vbeta=2.
if not keyword_set(rfunct) then rfunct='INVSQ'
	if not keyword_set(npeg) then npeg=1e9
	if not keyword_set(rscl) then rscl=1.
if not keyword_set(verbose) then verbose=1

if n_elements(latrng) eq 0 then latrng=[-3.,3.]
if n_elements(lonrng) eq 0 then lonrng=[0.,360.]

if not keyword_set(imgwdt) then imgwdt=5.
if not keyword_set(imgpix) then imgpix=128
if not keyword_set(imglog) then imglog=0
if not keyword_set(r2p2b) then r2p2b=0

if not keyword_set(inrad) then inrad=2.
if not keyword_set(maxrad) then maxrad=20.
if not keyword_set(rrng) then rrng=[1.0,20.0]
if not keyword_set(rdel) then rdel=0.1
if n_elements(opaq) eq 0 then opaq=1
theta0=0.
phi0=0.

dlat=((float(imgwdt)/imgpix)/maxrad)*(180./!pi)*2.
arcone=dlat^2
print,dlat,arcone

vstr=vradimage(vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$
	vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$
	dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$
	imgwdt=imgwdt,imgpix=imgpix,imgout=imgout,imgroot=imgroot,$
	imglog=imglog,r2p2b=r2p2b,$
	verbose=verbose, inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$
	opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0)

window,0,retain=2,title='weighted velocity profile'
plot,vstr.VMID,vstr.VPROFILE,xtitle='velocity [km/s]',ytitle='velocity profile (n*vol)'

window,1,retain=2,xsize=imgpix*4,ysize=imgpix*4,title='blue-shifts'
loadct,1 & tvlct,rr,gg,bb,/get & rr=reverse(rr) & gg=reverse(gg) & bb=reverse(bb) & tvlct,rr,gg,bb
tvscl,rebin((reform(vstr.VIMAGE[1,*,*])),imgpix*4,imgpix*4)

window,2,retain=2,xsize=imgpix*4,ysize=imgpix*4,title='red-shifts'
loadct,3
tvscl,rebin((reform(vstr.VIMAGE[2,*,*])),imgpix*4,imgpix*4)

end

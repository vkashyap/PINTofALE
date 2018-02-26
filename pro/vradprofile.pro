function vradprofile,vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$
	vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$
	dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$
	verbose=verbose, _extra=e
;+
;function	vradprofile
;	wrapper for RADPROJECT(), takes a given radial vector field (like, say
;	velocity), computes the projections onto the YZ plane at \infty, and
;	bins up the projections along X to return a velocity profile.
;	returns a structure of the form
;	{velocity grid, velocities histogram, mid-bins, weighted velocity profile}
;
;syntax
;	vstr=vradprofile(vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$
;	vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$
;	dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$
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
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to RADPROJECT()
;		INRAD,MAXRAD,RRNG,RDEL,OPAQ,ARCONE,THETA0,PHI0
;		WARNING: THETA0,PHI0 are only valid for small angles
;		in RADPROJECT(), and don't work at all when called
;		from this program.
;
;example usage
;	.run vradprofile
;	@vradprofile.par
;
;history
;	vinay kashyap (2012dec10; supersedes sphrofile)
;	bug fixes (theta0,phi0=0; latrng,lonrng reversal; latrng,lonrng
;	  not ranges; dvel) (VK; 2013mar)
;	bug fix with user defined functions (JJD/VK; 2014feb/apr)
;	now approaching velocities get negative values (VK; 2015dec)
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
  print,'Usage: vstr=vradprofile(vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$'
  print,'       vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$'
  print,'       dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$'
  print,'       verbose=verbose, inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$'
  print,'       opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0)'
  print,'  wrapper for radproject, returns binned values of'
  print,'  radial vector field projected onto plane at \infty'
  print,''
  print,'  This only makes integrated velocity profiles.  For images, see vradimage()'
  print,''
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
zvarr=0 & zwarr=0
for ib=0L,nb-1L do begin			;{for all latitudes
  xb=bb[ib] & dbcos=abs(cos(xb*!pi/180.))
  if dbcos ne 0 then begin
    dl=db/dbcos < (lmax-lmin)
  endif else dl=(lmax-lmin)
  nl=long((lmax-lmin)/dl+0.5)>1 & ll=findgen(nl)*dl+lmin
  if vv gt 0 then kilroy,dot=strtrim(ib,2)+'('+strtrim(nl,2)+')'

  zvarrb=0 & zwarrb=0
  for il=0L,nl-1L do begin			;{for all longitudes
    xl=ll[il]
    rstr=radproject(xb,xl,verbose=vv, _extra=e)
    rr=rstr.R & pp=rstr.P & vol=rstr.VOL
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
      ;	-ve sign to mesh with convention of blue-shifts being -ve
      zvw=zw*vol	;weighting
      if not keyword_set(zvarrb) then begin
        zvarrb=zvp & zwarrb=zvw
      endif else begin
        zvarrb=[zvarrb,zvp] & zwarrb=[zwarrb,zvw]
      endelse

    endif
  endfor					;IL=0,NL-1}
  if not keyword_set(zvarr) then begin
    zvarr=zvarrb & zwarr=zwarrb
  endif else begin
    zvarr=[zvarr,zvarrb] & zwarr=[zwarr,zwarrb]
  endelse

endfor						;IB=0,NB-1}

;	and now make the histogram 
minv=min(zvarr,max=maxv,/nan)
if keyword_set(vmin) then minv=vmin
if keyword_set(vmax) then maxv=vmax
dv=(maxv-minv)/100. & if keyword_set(dvel) then dv=abs(dvel[0])
hv=histogram(zvarr,min=minv,max=maxv,binsize=dv,reverse_indices=ri)
nhv=n_elements(hv) & vgrid=findgen(nhv+1L)*dv+minv
vmid=findgen(nhv)*dv+minv+0.5*dv

vpp=dblarr(nhv)
for i=0L,nhv-1L do begin
  if ri[i] ne ri[i+1] then begin
    ok=ri[ri[i]:ri[i+1]-1]
    vpp[i]=total(zwarr[ok])
  endif
endfor

if vv gt 900 then stop,'HALTing; type .CON to continue'

vstr=create_struct('VGRID',vgrid,'VHIST',hv,'VMID',vmid,'VPROFILE',vpp)

return,vstr
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example calling sequence

if not keyword_set(vfunct) then vfunct='BETA'
	if not keyword_set(vinf) then vinf=1e3
	if not keyword_set(vpeg) then vpeg=0.
	if not keyword_set(rpeg) then rpeg=1.
	if not keyword_set(vbeta) then vbeta=2.
if not keyword_set(rfunct) then rfunct='INVSQ'
	if not keyword_set(npeg) then npeg=1e9
	if not keyword_set(rscl) then rscl=1.
if not keyword_set(verbose) then verbose=1

if n_elements(latrng) eq 0 then latrng=[-90.,90.]
if n_elements(lonrng) eq 0 then lonrng=[0.,360.]

if not keyword_set(inrad) then inrad=1.
if not keyword_set(maxrad) then maxrad=20.
if not keyword_set(rrng) then rrng=[19.5,20.0]
if not keyword_set(rdel) then rdel=0.1
if n_elements(opaq) eq 0 then opaq=1
if not keyword_set(arcone) then arcone=1.
theta0=0.
phi0=0.

vstr=vradprofile(vfunct,rfunct,latrng=latrng,lonrng=lonrng,dlat=dlat,$
	vinf=vinf,vpeg=vpeg,vbeta=vbeta,rpeg=rpeg,npeg=npeg,rscl=rscl,$
	dvel=dvel,vmin=vmin,vmax=vmax,uvfunct=uvfunct,urfunct=urfunct,$
	verbose=verbose, inrad=inrad,maxrad=maxrad,rrng=rrng,rdel=rdel,$
	opaq=opaq,arcone=arcone,theta0=theta0,phi0=phi0)

plot,vstr.VMID,vstr.VPROFILE,xtitle='velocity [km/s]',ytitle='velocity profile (n*vol)'

end

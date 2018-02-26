pro fidgit,x,y,lstr,cstr,ysigma=ysigma,dem=dem,logt=logt,abund=abund,$
	abfrac=abfrac,ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,lsf=lsf,$
	wdem=wdem,wspc=wspc,verbose=verbose, _extra=e
;+
;procedure	fidgit
;	digitally (i.e., using Darwin-given fingers) fit a DEM
;	to a given spectrum
;
;syntax
;	fidgit,x,y,lstr,cstr,ysigma=ysigma,dem=dem,logt=logt,abund=abund,$
;	abfrac=abfrac,ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,lsf=lsf,$
;	wdem=wdem,wspc=wspc,verbose=verbose,$
;	pres=pres,logP=logP,n_e=n_e,desig=desig,econf=econf,/allah,$
;	chidir=chidir,chifil=chifil,chidir=chidir,eqfile=eqfile,$
;	effar=effar,wvlar=wvlar
;
;parameters
;	x	[INPUT; required] wavelengths or channels of observed
;		spectrum
;		* if channels, must specify an RMF via the keyword LSF
;		  that matches the binning of the spectrum
;		  BEWARE: no checks are performed to verify that the
;		  supplied RMF is consistent with the input data
;	y	[INPUT; required] spectrum
;	lstr	[I/O; required] line emissivity structure of the sort read in
;		by RD_LINE()
;		* if not given on input, will call RD_LINE() and generate it
;	cstr	[I/O; required] continuum emissivity structure of the sort
;		read in by RD_CONT()
;		* if not given on input, will call RD_LINE() and generate it
;
;keywords
;	ysigma	[INPUT] errors on Y; default is 1+SQRT(0.75+ABS(Y))
;	dem	[I/O] DEM to start out with, will contain the final
;		result on output
;	logt	[I/O] log(T [K]) at which DEM are defined.  ignored on
;		input if it doesn't match the size of DEM
;	abund	[INPUT] abundances
;		* if not set, assumed to be from Anders & Grevesse 1989
;	abfrac	[I/O] multiplicative factor by which to modify ABUND
;		* if scalar, assumed to be metallicity
;		* if vector, must match size of ABUND, else gets reset to 1
;		* default is 1
;	ldbdir	[INPUT] path to line emissivity database to read in if LSTR
;		is not given
;		* default is '$CHIANTI'
;		* uses !LDBDIR if set
;	cdbdir	[INPUT] path to continuum emissivity database to read in if
;		CSTR is not given
;		* default is '$CONT'
;		* uses !CDBDIR if set
;	ceroot	[INPUT] prefix for continuum emissivity files
;		* default is 'cie'
;		* uses !CEROOT if set
;	lsf	[INPUT] line spread function
;		* if scalar, the width of the line in bins used as a
;		  halfwidth for boxcar smoothing
;		* if vector, assumed to be the kernel with which to smooth
;		* if RMF structure (of the type read in from RD_OGIP_RMF())
;		  then convolves the theoretical spectrum with this RMF
;	wdem	[INPUT] window number to display DEM in
;	wspc	[INPUT] window number to display spectrum in
;	verbose	[INPUT] controls chatter
;		* uses !VERBOSE if set
;	_extra	[INPUT] allows specifying defined keywords to subroutines
;		called by this program
;		* RD_LINE: PRES, LOGP, N_E, DESIG, ECONF, ALLAH
;		* RD_CONT: PRES, LOGP, N_E
;		* FOLD_IONEQ: CHIFIL, VERBOSE
;		* RD_IONEQ: CHIDIR, EQFILE
;		* LINEFLX: EFFAR, WVLAR
;
;restrictions
;	* requires subroutines
;	  -- DEMACS
;	  -- WHEE
;	  -- RD_LINE [FOLD_IONEQ, RD_IONEQ, READ_IONEQ]
;	  -- RD_CONT
;	  -- LINEFLX
;	  -- CONV_RMF
;	  -- INICON
;
;history
;	vinay kashyap (Feb97)
;	modified ion balance calcs (VK; 99May)
;	completely rewritten (VK; Aug04)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y)
nl=n_elements(lstr) & ml=n_tags(lstr)
nc=n_elements(cstr) & mc=n_tags(cstr)
if np lt 4 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X: undefined' else $
  if ny eq 0 then ok='Y: undefined' else $
   if nx ne ny and nx ne ny+1L then ok='X and Y: incompatible' else $
    if nx lt 2 then ok='Insufficient bins'
if ok ne 'ok' then begin
  print,'Usage: fidgit,x,y,lstr,cstr,ysigma=ysigma,dem=dem,logt=logt,abund=abund,$'
  print,'       abfrac=abfrac,ldbdir=ldbdir,cdbdir=cdbdir,ceroot=ceroot,lsf=lsf,$'
  print,'       wdem=wdem,wspc=wspc,verbose=verbose,$'
  print,'       pres=pres,logP=logP,n_e=n_e,desig=desig,econf=econf,/allah,$'
  print,'       chidir=chidir,chifil=chifil,chidir=chidir,eqfile=eqfile,$'
  print,'       effar=effar,wvlar=wvlar'
  print,'  digital fitting of DEMs to line spectra'
  if np ne 0 then message,ok,/informational
  return
endif

;	check keywords
nsig=n_elements(ysigma)
;
ysig=1.+sqrt(abs(y)+0.75) & if nsig eq nx then ysig=0.+ysigma
;
ivar=0 & defsysv,'!LDBDIR',exists=ivar
if ivar ne 0 then setsysval,'LDBDIR',ldb,/getval
if keyword_set(ldbdir) then ldb=strtrim(ldbdir[0],2)
if not keyword_set(ldb) then ldb='$CHIANTI'
;
ivar=0 & defsysv,'!CDBDIR',exists=ivar
if ivar ne 0 then setsysval,'CDBDIR',cdb,/getval
if keyword_set(cdbdir) then cdb=strtrim(cdbdir[0],2)
if not keyword_set(cdb) then cdb='$CONT'
;
ivar=0 & defsysv,'!CEROOT',exists=ivar
if ivar ne 0 then setsysval,'CEROOT',cer,/getval
if keyword_set(ceroot) then cer=strtrim(ceroot[0],2)
if not keyword_set(cer) then cer='cie'
;
ivar=0 & defsysv,'!ABUND',exists=ivar
if ivar ne 0 then setsysval,'ABUND',abnd,/getval
if keyword_set(abund) then abnd=long(abund[0])>1
if not keyword_set(abnd) then abnd=getabund('anders & grevesse')
nab=n_elements(abnd) & mab=n_elements(abfrac)
if mab eq 1 then begin
  tmp=fltarr(nab)+1. & tmp[2:*]=abfrac[0] & abfrac=tmp
endif else begin
  if mab ne nab then begin
    message,'ABFRAC incompatible with ABUND; assuming metallicity=1',/informational
    abfrac=fltarr(nab)+1.
  endif
endelse
;
ivar=0 & defsysv,'!VERBOSE',exists=ivar
if ivar ne 0 then setsysval,'VERBOSE',vv,/getval
if keyword_set(verbose) then vv=long(verbose[0])>1
if not keyword_set(vv) then vv=0L
;
if not keyword_set(wdem) then wdem=1
if not keyword_set(wspc) then wspc=0

;	is there an RMF?
nrmf=n_elements(lsf) & mrmf=n_tags(lsf)
if nrmf eq 0 then begin
  if vv gt 0 then message,$
	'assuming that X are wavelengths in [Angstroms]',/informational
  wrange=minmax(x)
endif else begin
  if mrmf eq 0 then begin
    if vv gt 0 then message,$
	'assuming that X are wavelengths in [Angstroms]',/informational
    wrange=minmax(x)
    if nrmf eq 1 then linwdt=abs(long(lsf[0]))>2L else linwdt=0
    if nrmf gt 1 and nrmf lt nx then linkern=float(lsf[*]) else linkern=0
    if nrmf ge nx then message,$
	'LSF array too long; ignoring',/informational
  endif else begin
    if vv gt 0 then message,$
	'assuming that X channels',/informational
    wrange=12.3985/[max(lsf.EMN),min(lsf.EMX)]
  endelse
endelse

;	read in the emissivities
if ml eq 0 then begin
  emis=rd_line(atom,wrange=wrange,dbdir=ldb,verbose=vv,fstr=lstr, _extra=e)
  emis=fold_ioneq(lstr.LINE_INT,lstr.Z,lstr.JON,logt=lstr.LOGT,verbose=vv,$
	_extra=e)
  lstr.LINE_INT=emis
endif
if mc eq 0 then begin
  emis=rd_cont(cer,wrange=wrange,dbdir=cdb,abund=abfrac*abnd,fcstr=cstr,$
	_extra=e)
endif

;	initialize
go_on=1						;terminate while loop if 0
window,wdem					;window for DEM editing
loadct,3 & peasecolr & inicon,atom=atom
if n_elements(logt) le 1 then logT=lstr.logt
if n_elements(dem) ne n_elements(logt) then dem=logt
log_T=logt
xrange=[min(log_T),max(log_T)] & yrange=[min(dem),max(dem)]
dem=demacs(lstr.LOGT,dem0=dem,logt0=log_t,group=grp,igroup=igrp,$
	xrange=xrange,yrange=yrange)
d0=dem & d1=d0 & d2=d0 & d3=d0 & d4=d0 		;DEM save buffers
d5=d0 & d6=d0 & d7=d0 & d8=d0 & d9=d0
window,wspc					;window for showing spectra
lfx=lineflx(lstr.LINE_INT,lstr.LOGT,lstr.WVL,lstr.Z,DEM=DEM,$
	abund=abfrac*abnd, _extra=e)
cfx=lineflx(cstr.CONT_INT,cstr.LOGT,cstr.midWVL,DEM=DEM, _extra=e)
cww=mid2bound(cstr.midWVL) & cdw=cww[1:*]-cww
if mrmf eq 0 then wgrid=mid2bound(x) else begin
  emn=lsf.EMN & emx=lsf.EMX
  wmn=12.3985/emx & wmx=12.3985/emn
  os=sort(wmn) & wgrid=[wmn[os],max(wmx)]
endelse
lspec=hastogram(lstr.WVL,wgrid,wts=lfx)
cspec=rebinw(cfx*cdw,cww,wgrid,/perbin)
spc=lspec+cspec
if mrmf ne 0 then begin
  nrg=12.3985/wgrid & os=sort(nrg) & nrg=nrg[os] & tmp=spc[os]
  conv_rmf,nrg,tmp,xx,spc,lsf
  xx=12.3985/xx & os=sort(xx) & xx=xx[os] & spc=spc[os]
  if total(abs(xx-x)) gt 1 then begin
    message,'possible that LSF is incompatible with X',/informational
    if vv ge 5 then stop,'HALTing.. type .CON to continue'
  endif
endif
chi2=total((y-spc)^2/ysig^2)/nx
wset,wspc
sxrange=minmax(x) & syrange=minmax(y)
plot,x,y,title='Reduced Chi^2:'+strtrim(chi2,2),xrange=sxrange,yrange=syrange
for i=0L,nx-1L do oplot,x[i]*[1,1],y(i)+ysig[i]*[1,-1],col=1
oplot,x,y,psym=10,col=2

;	begin loop of setting DEM and comparing observed spectrum to
;	computed spectrum
c1=''
while go_on eq 1 do begin		;{terminate if go_on=0

  ;	get DEM
  wset,wdem
  if strlowcase(c1) ne 'z' and $
   strlowcase(c1) ne 'n' and $
   strlowcase(c1) ne 'a' and $
   c1 ne 'X' and c1 ne 'Y' then $
   dem=demacs(lstr.logT,dem0=dem,group=grp,igroup=igrp,$
	xrange=xrange,yrange=yrange)

  ;	get spectrum
  wset,wspc
  lfx=lineflx(lstr.LINE_INT,lstr.LOGT,lstr.WVL,lstr.Z,DEM=DEM,$
	abund=abfrac*abnd, _extra=e)
  cfx=lineflx(cstr.CONT_INT,cstr.LOGT,cstr.midWVL,DEM=DEM, _extra=e)
  lspec=hastogram(abs(lstr.WVL),wgrid,wts=lfx)
  cspec=rebinw(cfx*cdw,cww,wgrid,/perbin)
  spc=lspec+cspec
  if keyword_set(linwdt) then spc=smooth(spc,linwdt>2,/edge)
  if keyword_set(linkern) then spc=convol(spc,linkern,/edge_wrap)
  xx=x
  if mrmf ne 0 then begin
    nrg=12.3985/wgrid & os=sort(nrg) & nrg=nrg[os] & tmp=spc[os]
    conv_rmf,nrg,tmp,xx,spc,lsf
    xx=12.3985/xx & os=sort(xx) & xx=xx[os] & spc=spc[os]
  endif
  chi2=total((y-spc)^2/ysig)/nx
  plot,x,y,title='Reduced Chi^2:'+strtrim(chi2,2),$
	xrange=sxrange,yrange=syrange
  oplot,xx,spc,psym=10,col=2

  ;	what to do next?
  c1='q to quit, z to halt, n to "renorm", a to change abundances' & print,c1
  c1='x,y to set Xrange, Yrange for DEM plot' & print,c1
  c1='X,Y to set Xrange, Yrange for spectrum plot' & print,c1
  c1='s,0-9 to save in buffer, r,<shift>0-9 to restore from buffer'
  c1='any other key to continue' & print,c1
  c1=get_kbrd(1)
  if c1 ne 'X' and c1 ne 'Y' then c1=strlowcase(c1)
  case c1 of				;(options
    'q': go_on=0			;	quit
    'z': stop,'HALTing.. type .CON to continue'
    'n': begin
      message,'multiplying DEM by '+strtrim(total(y)/total(spc),2),/info
      dem=dem*total(y)/total(spc)	;	renormalize
      yrange=yrange*total(y)/total(spc)
    end
    'a': begin				;(change abundances
      print,'changing abundance' & iz=-1 & az=1.
      read,prompt='type atomic number and factor by which to change: ',iz,az
      print,'changing the abundance of '+atom[iz-1]+' from '+$
	strtrim(abnd[iz-1],2)+' to '+strtrim(abnd[iz-1]*az,2)
	abfrac[iz-1]=az
    end					;abundances)
    'x': begin
      read,prompt='type in X-range values for DEM plot ['+$
	strtrim(xrange[0],2)+','+strtrim(xrange[1],2)+']> ',$
	xrange
    end
    'y': begin
      read,prompt='type in Y-range values for DEM plot ['+$
	strtrim(yrange[0],2)+','+strtrim(yrange[1],2)+']> ',$
	yrange
    end
    'X': begin
      read,prompt='type in X-range values for spectral plot ['+$
	strtrim(sxrange[0],2)+','+strtrim(sxrange[1],2)+']> ',$
	sxrange
    end
    'Y': begin
      read,prompt='type in Y-range values for spectral plot ['+$
	strtrim(syrange[0],2)+','+strtrim(syrange[1],2)+']> ',$
	syrange
    end
    's': d0=dem				;	save
    '0': d0=dem
    '1': d1=dem
    '2': d2=dem
    '3': d3=dem
    '4': d4=dem
    '5': d5=dem
    '6': d6=dem
    '7': d7=dem
    '8': d8=dem
    '9': d9=dem
    '(': dem=d9				;	restore
    '*': dem=d8
    '&': dem=d7
    '^': dem=d6
    '%': dem=d5
    '$': dem=d4
    '#': dem=d3
    '@': dem=d2
    '!': dem=d1
    ')': dem=d0
    'r': dem=d0
    else: 				;	nothing
  endcase				;case c1 of)

endwhile				;go_on=1}

return
end

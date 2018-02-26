function rd_coco,cocofil,wrange=wrange,nbin=nbin,abund=abund,$
	noerg=noerg,noang=noang,cemis=cemis,wvl=wvl,wmid=wmid,$
	logT=logT,verbose=verbose
;+
;function	rd_coco
;	reads in the ATOMDB continuum emissivities and returns it in
;	an IDL structure
;
;syntax
;	cocostr=rd_coco(cocofil,wrange=wrange,nbin=nbin,abund=abund,$
;	/noerg,/noang,cemis=cemis,wvl=wvl,wmid=wmid,logT=logT,$
;	verbose=verbose)
;
;parameters
;	cocofil	[INPUT; required] full path to the location of the
;		continuum emissivities database file
;		* if set to '$COCO', looks for $APECROOT/apec_coco.fits
;		* if set to '$oldAPED', looks for /fubar/atomdb/atomdb_v2.0.2/apec_coco.fits
;		* if set to '$APED', looks for /fubar/atomdb/atomdb_v3.0.2/apec_coco.fits
;		* if the environment variable APECROOT is not set, uses
;		  getpoadef('TOPDIR')+'/atomdb/'
;
;keywords
;	wrange	[INPUT] 2-element array to define the range over which
;		the continuum is extracted
;		* default is [0.01,1300]
;	nbin	[INPUT] number of bins in the output array
;		* if -ve, the bins will be in log scale, with longer
;		  widths at longer wavelengths
;		* default is 2000
;	abund	[INPUT] if supplied, replaces the hardcoded abundances
;		with new ones
;		* if not set, keeps it as Anders & Grevesse
;	noerg	[INPUT] if set, refrains from converting the units
;		from [ph ...] to [erg ...]
;	noang	[INPUT] if set, refrains from converting the units
;		from [... /keV] to [... /Ang]
;	cemis	[OUTPUT] emissivity array of size (TEMPERATURE,WAVELENGTH)
;		collated from all the various extensions in COCOFIL,
;		interpolated onto a standard wavelength grid, and with
;		PoA style units [1e23 ergs cm^3/s/Ang] (unless NOERG and
;		NOANG are set)
;	wvl	[OUTPUT] the standardized output wavelength grid
;	wmid	[OUTPUT] the mid-bin values of the wavelength grid
;	logT	[OUTPUT] the temperature grid over which CEMIS is defined
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Apr2009)
;	changed default $APED to point to v2.0.0 (VK; May2011)
;	changed default $APED to point to v3.0.2, now looks for the
;	  environment variable APECROOT (VK; Dec2015)
;-

;	usage
ok='ok' & np=n_params() & nc=n_elements(cocofil) & szc=size(cocofil,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if nc eq 0 then ok='COCOFIL is undefined' else $
  if szc ne 7 then ok='COCOFIL cannot be understood' else $
   if nc gt 1 then ok='COCOFIL must be a scalar string'
if ok ne 'ok' then begin
  print,'Usage: cocostr=rd_coco(cocofil,wrange=wrange,nbin=nbin,abund=abund,$'
  print,'  /noerg,/noang,cemis=cemis,wvl=wvl,wmid=wmid,logT=logT,verbose=verbose)'
  print,'  reads in the ATOMDB continuum emissivities'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
mbin=2000L & if keyword_set(nbin) then mbin=nbin[0]
;
wmin=0.01 & wmax=1300. & nw=n_elements(wrange)
if nw eq 1 then begin
  wmin = wmin < wrange[0]
  wmax = wmax > wrange[0]
endif
if nw eq 2 then begin
  wmin=min(wrange)
  wmax=max(wrange)
endif
if wmin lt 0 then begin
  message,'wavelength range cannot go negative; resetting to 0.01',/informational
  wmin=0.01
endif
if mbin lt 0 and wmin eq 0 then begin
  message,'wavelength minimum cannot be zero; resetting to 0.01',/informational
  wmin=0.01
endif
;
defabund=getabund('anders & grevesse')
newabund=defabund & if n_elements(abund) ge 30 then newabund=abund

;	wavelength grid
y1=wmax & if mbin lt 0 then y1=alog10(wmax)
y0=wmin & if mbin lt 0 then y0=alog10(wmin)
delw=(y1-y0)/(abs(mbin))
wvl=y1-findgen(abs(mbin)+1)*delw & if mbin lt 0 then wvl=10.^(wvl)
wmid=0.5*(wvl[1:*]+wvl) & dW=abs(wvl[1:*]-wvl)

;	input file
filnam=cocofil[0]
apecroot=getenv('APECROOT')
if not keyword_set(apecroot) then apecroot=filepath('atomdb',root=getpoadef('TOPDIR'))
if strmatch(strupcase(filnam),'$COCO') then filnam=filepath('apec_coco.fits',root=apecroot)
;if strmatch(strupcase(filnam),'$oldAPED') then filnam='/fubar/SCAR/atomdb_1.3.1/apec_coco.fits'
if strmatch(strupcase(filnam),'$oldAPED') then filnam='/fubar/atomdb/atomdb_v2.0.2/apec_coco.fits'
if strmatch(strupcase(filnam),'$APED') then filnam='/fubar/atomdb/atomdb_v3.0.2/apec_coco.fits'

;	read in first extension and figure what we've got
t1=mrdfits(filnam,1,ht1)
logT=alog10(t1.kT/0.08617)+6. & edens=t1.EDENSITY & zz=t1.NELEMENT & ncont=t1.NCONT & npseudo=t1.NPSEUDO
tagnam=tag_names(t1) & tmp=create_struct('header',ht1) & for i=0,n_elements(tagnam)-1 do tmp=create_struct(tmp,tagnam[i],t1.(i))
nT=n_elements(logT) & cocostr=create_struct('APEC_COCO',tmp)

;	output emissivity
cemis=dblarr(nT,abs(mbin))

;	now read in all the different temperature extensions
for iT=0L,nT-1L do begin
  if vv gt 0 then kilroy
  t=mrdfits(filnam,iT+2,ht)
  hdunam=sxpar(ht,'HDUNAME')
  while strpos(hdunam,'.') ge 0 do begin & strput,hdunam,'_',strpos(hdunam,'.') & endwhile
  tagnam=tag_names(t) & tmp=create_struct('header',ht) & for i=0,n_elements(tagnam)-1 do tmp=create_struct(tmp,tagnam[i],t.(i))
  cocostr=create_struct(cocostr,hdunam,tmp)

  tZ=t.Z & n_cont=t.N_CONT & E_cont=t.E_CONT
  n_pseudo=t.N_PSEUDO & E_pseudo=t.E_pseudo
  coc=t.CONTINUUM & cop=t.PSEUDO

  for iZ=0L,zz[iT]-1L do begin
    zab=newabund[tZ[iZ]-1]/defabund[tZ[iZ]-1]
    ;
    nZc=n_cont[iZ] & nrgc=1.6021892e-12*E_cont[0L:nZc-1L,iZ]*1e3*1d23 & w_cont=12.398521/E_cont[0L:nZc-1L,iZ]
    coc[*,iZ]=coc[*,iZ]*zab	;[ph cm^3/s/keV]
    if not keyword_set(noerg) then coc[0L:nZc-1L,iZ]=coc[0L:nZc-1L,iZ]*nrgc	;[1e23 ergs cm^3/s/keV]
    tmp=(interpol(coc[0L:nZc-1L,iZ],w_cont,wmid)>0)
    if not keyword_set(noang) then tmp=tmp*12.398521/wmid^2	;[1e23 ergs cm^3/s/Ang]
    cemis[iT,*]=cemis[iT,*]+tmp
    ;
    nZp=n_pseudo[iZ] & nrgp=1.6021892e-12*E_pseudo[0L:nZp-1L,iZ]*1e3*1d23 & w_pseudo=12.398521/E_pseudo[0L:nZp-1L,iZ]
    cop[*,iZ]=cop[*,iZ]*zab	;[ph cm^3/s/keV]
    if not keyword_set(noerg) then cop[0L:nZp-1L,iZ]=cop[0L:nZp-1L,iZ]*nrgp	;[1e23 ergs cm^3/s/keV]
    tmp=(interpol(cop[0L:nZp-1L,iZ],w_pseudo,wmid)>0)
    if not keyword_set(noang) then tmp=tmp*12.398521/wmid^2	;[1e23 ergs cm^3/s/Ang]
    cemis[iT,*]=cemis[iT,*]+tmp
  endfor
endfor

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,cocostr
end

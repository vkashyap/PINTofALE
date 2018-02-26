function ionabs,w,abund=abund,ionfrac=ionfrac,vfkydir=vfkydir,ikev=ikev,$
	nconsec=nconsec,DEM=DEM,logt=logt,verbose=verbose,noHeH=noHeH,$
	icrstr=icrstr, _extra=e
;+
;function	ionabs
;       returns the photoelectric absorption cross-sections for a given
;	photon energy for specified chemical composition, and ion fractions.
;	Uses ground state photoionisation cross-sections computed using the
;	fortran code of Verner, D.A., Ferland, G.J., Korista, K.T., & Yakovlev,
;	D.G., 1996, ApJ, 465, 487 and is supplemented with high resolution
;	cross-sections in the vicinity of 1s-2p resonances for
;	O (from Garcia et al., 2005, ApJS, 158, 68),
;	Ne (from Juett et al., 2006, ApJ, 648, 1066), and
;	C (from Hasogle et al, 2010, ApJ, 724, 1296).
;
;syntax
;	sigabs=ionabs(w,abund=abund,ionfrac=ionfrac,vfkydir=vfkydir,/ikeV,$
;	nconsec=nconsec,DEM=DEM,logt=logt,verbose=verbose,noHeH=noHeH,$
;	icrstr=icrstr,chidir=chidir,eqfile=eqfile)
;
;parameters
;	w	[INPUT; required] photon energy values at which to compute
;		the absorption cross-section
;		* default units are [Ang], unless IKEV is set, when they are
;		  assumed to be [keV]
;
;keywords
;	abund	[INPUT] abundances relative to H=1
;		* default is to use Grevesse & Sauval
;	ionfrac	[I/O] the ion fraction for all elements
;		* must be a 2D array of size (NZ,NION), where
;		  NION=NZ+1 and usually NZ=30
;		* note that in normal cases, total(IONFRAC,2)
;		  is an array of 1s.
;		* if not given, assumed to be all neutral, IONFRAC[*,0]=1
;		  - unless DEM is given, see below
;		* if non-zero scalar, assumed to be all neutral
;		  - DEM is ignored
;		* if given as an array of size NZ, each element
;		  is assumed to have a total of that fraction,
;		  and all the ions have the same fraction
;		  - this is clearly an artificial case, to be used mainly
;		    for debugging and the like
;		  - DEM is studiously ignored
;		* if given as a 1D array of size not matching NZ,
;		  is ignored with a warning and the default is used
;		  - unless DEM is given, see below
;		* if given as a 2D array of size not matching (NZ,NZ+1),
;		  program quits with an error
;		  - unless DEM is given, see below
;		* if DEM (see below) is given, is read in as a 3D array
;		  of size (NT,NZ+1,NZ1) from the database (using the PINTofALE
;		  function RD_IONEQ(), which calls the CHIANTI procedure
;		  READ_IONEQ), weights the fractions according to the DEM
;		  and is converted into a 2D array of size (NZ,NZ+1)
;	vfkydir	[INPUT] where to find the save files that contain the
;		cross-sections calculated for each ionization state
;		* default is '$ARDB'
;		* NOTE: the cross-sections are not read in if the keyword
;		  ICRSTR is properly defined
;	ikeV	[INPUT] if set, assumes that W are in units of keV
;	nconsec	[INPUT] number of grid points that define the cross-sections
;		that must fall within one bin of the user defined grid before
;		the switch from averaging to interpolating occurs
;		* default=1
;	DEM	[INPUT] the Differential Emission Measure that will be
;		used to weight the different temperatures to compute the
;		ionization fraction
;		* looked at only if IONFRAC is not supplied as input
;		* size must match LOGT, otherwise ignored
;	logT	[INPUT] the temperature grid over which DEM is defined
;		* size must match DEM, otherwise DEM is ignored
;	noHeH	[INPUT] if set to any number other than 1 or 2, excludes
;		H and He from the cross-sections
;		* if set to 1, only excludes H
;		* if set to 2, only excludes He
;	icrstr	[I/O] the cross-sections calculated for each ionization
;		state, stored in a structure of structures of the form
;		ICRSTR.(Z) = {EVSPLAC, CROSSI}
;		* if given on input, will use these numbers
;		* if the structure does not contain enough tags, or if
;		  any one of the substructures does not contain enough
;		  tags, will be read in using the data in VFKYDIR
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined variables to subroutines
;		-- RD_IONEQ : CHIDIR, EQFILE
;
;subroutines
;	GETABUND()
;	RD_IONEQ()
;	READ_IONEQ
;	INICON
;	GETPOADEF()
;
;history
;       written as IONTAU by Jeremy Drake (Apr07)
;       further shoehorned into PoA format by Vinay Kashyap (Apr07)
;	bug fixes and corrections and changed behavior of IONFRAC
;	  (VK; May09)
;	updated references (JJD; Aug11)
;	BUG: DEM was not being normalized when weighting IONFRAC (JJD; Nov11)
;	added call to GETPOADEF (VK; Aug15)
;-

;       usage
ok='ok' & np=n_params() & nw=n_elements(w)
if np eq 0 then ok='Insufficient parameters' else $
 if nw eq 0 then ok='energy/wavelength grid is undefined' else $
  if nw eq 1 then ok='energy/wavelength grid must be an array'
if ok ne 'ok' then begin
  print,'Usage: sigabs=ionabs(w,abund=abund,ionfrac=ionfrac,vfkydir=vfkydir,$'
  print,'       /ikeV,nconsec=nconsec,logt=logt,verbose=verbose,noHeH=noHeH,$'
  print,'       icrstr=icrstr,chidir=chidir,eqfile=eqfile)'
  print,'  compute absorption cross-section given abundance and ion fractions'
  print,'  based on Verner et al. (1996) calculations'
  if np ne 0 then message,ok,/informational
  return,0.
endif

;message,'*** MEN AT WORK *** DO NOT USE ***'

;	need this
inicon,atom=atom,aname=aname

;       check keywords

;VERBOSE
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;ABUND
abun=getabund('Grevesse & Sauval')
if n_elements(abund) gt 2 then abun=abund & nab=n_elements(abun)
atom=strlowcase(atom[0L:nab-1L]) & Z=indgen(nab)+1 & nZ=n_elements(Z)

;VFKYDIR
;vdir='/data/letg10/gorczyca/iontau'
;if size(vfkydir,/type) eq 7 then vdir=vfkydir[0]

;NCONSEC
mconsec=1 & if keyword_set(nconsec) then mconsec=fix(nconsec)

;IONFRAC, DEM(LOGT)
;	if IONFRAC not given and DEM not given: neutral
;	if IONFRAC not given and DEM is given: read in
;	if IONFRAC is non-zero scalar: neutral, and always ignore DEM
;	if IONFRAC is array of size NZ: flat across ions, and ignore DEM
;	if IONFRAC is array of size !NZ: neutral, but read in if DEM given
;	if IONFRAC is 2D array of size [NZ,NION]: use as is and ignore DEM
;	if IONFRAC is 2D but ![NZ,NION]: quit with error, but read if DEM given
;	anything else: read if DEM is given, quit with error otherwise
ok='neutral' & nion=n_elements(ionfrac) & szion=size(ionfrac)
aye='nay' & nDEM=n_elements(DEM) & nT=n_elements(logT)
if nDEM ne 0 and nT ne 0 and nDEM eq nT then aye='aye'
if not keyword_set(ionfrac) and aye eq 'aye' then ok='read'
if keyword_set(ionfrac) and nion eq 1 then ok='neutral'
if szion[0] eq 1 then begin
  if nion eq nZ then ok='flat' else begin
    if aye eq 'aye' then ok='read' else ok='neutral'
  endelse
endif
if szion[0] eq 2 then begin
  if szion[1] eq nZ and szion[2] eq nZ+1 then ok='ok' else begin
    if aye eq 'aye' then ok='read' else ok='quit'
  endelse
endif
if szion[0] gt 2 then begin
  if aye eq 'aye' then ok='read' else ok='quit'
endif
;
case ok of
  'ok': ifraction=ionfrac
  'neutral': begin	;(default
    ifraction=fltarr(nZ,nZ+1) & ifraction[*,0]=1.
  end			;NEUTRAL)
  'flat': begin		;(same for all ions
    tmp=fltarr(nZ,nZ+1)
    for iab=0,nZ-1 do tmp[i,0:i]=ionfrac[i]/(i+1.)
    ifraction=tmp
  end			;FLAT)
  'read': begin		;(read in from database and fold in DEM
    ionfrac=rd_ioneq(indgen(nZ)+1,logT, _extra=e)
    tmp=ionfrac
    for i=0L,nT-1L do tmp[i,*,*]=DEM[i]*ionfrac[i,*,*]/total(DEM)	;had to divide by total(DEM) to ensure proper fractional weighting
    ifraction=total(tmp,1);/nT		;this division by nT seems to be a residue of previous attempts to normalize for the DEM, and so, given above correction by total(DEM), has now been commented out
    ;	the output of RD_IONEQ is in [NT,NION,NZ],
    ;	so the remaining axes must be transposed
    ifraction=transpose(ifraction)
  end			;READ)
  else: begin		;(whatchewtakinaboot, willis?
    message,'IONFRAC is not understandable; quitting',/informational
    return,-1L
  end			;QUIT)
endcase
ionfrac=ifraction

;ok='ok' & nion=n_elements(ionfrac) & szion=size(ionfrac)
;aye='flatDEM' & nDEM=n_elements(DEM) & nT=n_elements(logT)
;if nDEM eq 0 then aye='flatDEM' else $
; if nT eq 0 then aye='flatDEM' else $
;  if nDEM ne nT then begin
;    if nDEM eq 81 or nDEM eq 41 then aye='specialDEM'
;  endif else aye='DEM'
;case szion[0] of
;  0: begin
;    if nDEM gt 0 then ok='read' else ok='flat'
;    mZ=nab
;  end
;  1: begin
;    if szion[1] eq nab then ok='flatZ' else $
;     if szion[1] eq nab+1 then ok='flatION' else ok='flat'
;  end
;  2: begin
;    mZ=szion[1] & mION=szion[2]
;    if mZ ne nab or mION ne nab+1 then ok='read'
;  end
;  3: begin
;    mT=szion[1] & mZ=szion[2] & mION=szion[3]
;    if aye eq 'DEM' or aye eq 'specialDEM' then begin
;      if mZ eq nab and mION eq nab+1 then ok='DEM' else $
;       ok='read'
;    endif
;  end
;  else: ok='read'
;endcase
;;
;ifraction=fltarr(nab,nab+1) & ifraction[*,0]=1.
;if ok eq 'ok' then ifraction=ionfrac else $
; if ok eq 'read' then ifraction=rd_ioneq(indgen(nZ)+1,logT, _extra=e)
; if ok eq 'flatZ' then for i=0L,nab-1L do ifraction[i,*]=ionfrac[i] else $
;  if ok eq 'flatION' then for i=0L,nab do ifraction[*,i]=ionfrac[i] else $
;   if ok eq 'DEM' or ok eq 'DEM' then begin
;     xDEM=DEM
;     if aye eq 'specialDEM' then begin
;       if nDEM ne mT then begin
;         xlogT=findgen(41)*0.1+4.
;	 if nDEM eq 81 then xlogT=findgen(81)*0.05+4.
;         xDEM=interpol(DEM,xlogT,findgen(mT)*4./float(mt-1L)+4.)
;       endif
;     endif else xlogT=logT
;     tmp=ionfrac
;     for i=0L,mT-1L do tmp[i,*,*]=xDEM[i]*ionfrac[i,*,*]
;     ifraction=total(tmp,1)/mT
;     ionfrac=ifraction
;   endif
;ionfrac=ifraction

;IKEV
if not keyword_set(ikev) then begin
  ev=reverse(12398.5/w)	;input is in [Ang] -- convert to [eV]
endif else begin
  ev=w*1e3		;input is in [keV] -- convert to [eV]
endelse

;ICRSTR
nICR=n_tags(icrstr)

; cross-sections are stored in save files in the form of total
; cross-section for each ion for each element with Z=1-30 
; energy grid is highly non-uniform so as to capture edge energies
; and structure accurately.  
; structure of each array is (nion,nenergy)

;if not keyword_set(vfkydir) then vfkydir='/data/letg10/gorczyca/iontau/'
;zTOPDIR='/data/fubar/SCAR/'
;ivar=0 & defsysv,'!TOPDIR',exists=ivar  ;if !TOPDIR exists
;if ivar ne 0 then setsysval,'TOPDIR',zTOPDIR,/getval else begin
;  tmp=routine_info('IONABS',/source,/functions) & scardir=tmp.PATH
;  jvar=strpos(scardir,'pro/ionabs.pro')		;we know where IONABS.PRO is
;  zTOPDIR=strmid(scardir,0,jvar-1)
;  if vv ge 10 then message,'PoA directory is: '+zTOPDIR,/info
;endelse
;zARDB=filepath('',root_dir=filepath('ardb',root_dir=zTOPDIR))
zTOPDIR=getpoadef()
zARDB=getpoadef('ARDB')
;if not keyword_set(vfkydir) then vfkydir=zARDB
if size(vfkydir,/type) ne 7 then vdir=zARDB else vdir=vfkydir[0]

elem=strlowcase(atom[0:nab-1])
elemnam=aname[0:nab-1]

; set up ev array bin boundaries and the binsize for rebin when user
; grid is more coarse than save file data (needed for conservation of
; cross-section over fine resonances)

dev=ev[1:*]-ev & dev=[dev[0],dev]
bev=[ev[0]-dev[0]/2.,ev+dev/2.]

; set up interpolation and rebin array to contain the cumulative cross-section

totcrossrb=fltarr(n_elements(ev))*0.

if vv gt 10000L then stop,'HALTing; type .CON to continue'

; main loop to read in data for each element, then loop over each 
; element to sum up the cross-sections according to weighting
; by ion fractions.  Weight total element cross-section by abundance,
; where abundance of H=1

for i=0,n_elements(elem)-1 do begin

  splacfil=filepath(elem[i]+'splac.cross',root_dir=vdir)
  if nICR eq 0 then begin
    restore,splacfil
    tmp=create_struct(elemnam[i],create_struct('EVSPLAC',EVSPLAC,'CROSSI',CROSSI))
    if i eq 0 then icrstr=tmp else icrstr=create_struct(icrstr,tmp)
    if vv gt 0 then print,splacfil
    if vv gt 100000L then stop,'HALTing; type .CON to continue'
  endif else begin
    evsplac=icrstr.(i).EVSPLAC & crossi=icrstr.(i).CROSSI
  endelse

; compute total cross-section:
; need total(ionfrac*cross-secton)*abun for each element
; then rebin onto supplied grid
; In order to sample properly when original data grid is finer than
; user requested, we need rebin
; To deal with regions in which user grid is finer than data grid,
; need to use interpolation.
; Here, both methods are used and then combined according to which of
; these regimes the grids are in.

  totcross=abun[i]*reform(ifraction[i,0:i]#crossi)
  totcrossi=interpol(totcross,evsplac,ev)

  devsplac=evsplac[1:*]-evsplac & devsplac=[devsplac[0],devsplac]
  devsplaci=interpol(devsplac,evsplac,ev)
  oo=where(devsplaci lt dev,count)
  if (count gt mconsec) then begin
    doo=oo[1:*]-oo & doo=smooth([doo[0],doo],mconsec)
    oo=where(doo lt 1.5,count)
    if (count gt mconsec) then totcrossi[oo]=(rebinw(totcross,evsplac,bev))[oo]
  endif
  totcrossrb=totcrossrb+totcrossi

endfor

; convert to optical depth using column density.  Recall that units
; from Verner routines are Mb (10^-18 cm^2)

;totcrossrb=nh*totcrossrb*1.e-18
totcrossrb=totcrossrb*1.e-18
if not keyword_set(ikev) then totcrossrb=reverse(totcrossrb)

return,totcrossrb

end

function acisgarf0,wgrid,onchip,arm=arm,caldb=caldb,hrmaea=hrmaea,$
	acissqe=acissqe,grfH=grfH,grfM=grfM,grfL=grfL,grfcol=grfcol,$
	_extra=e
;+
;function	acisgarf0
;	returns 0th order effective area as a function of wavelength for
;	specified CXC grating on the ACIS-S
;
;	note that with the ACIS-S/HETG combination, the effective area for
;	only one of the components (HEG/MEG) is returned at a time.  Run
;	twice and add to get total:
;		arf=acisgarf0(wgrid,arm='meg')+acisgarf0(wgrid,arm='heg')
;
;syntax
;	arf=acisgarf0(wgrid,onchip,arm=arm,caldb=caldb,hrmaea=hrmaea,$
;	acissqe=acissqe,grfH=grfH,grfM=grfM,grfL=grfL,grfcol=grfcol)
;
;parameters
;	wgrid	[I/O; required] wavelengths at which effective area is computed
;		* if not given on input, will be computed for the same
;		  grid as for the 1st order spectrum.
;	onchip	[INPUT] the chip on which the source is located
;		* if not given, or is illegal, assumed to be the S3 chip.
;
;keywords
;	arm	[INPUT] string specifying which grating arm to be used.
;		* 'MEG', 'HEG', or 'LEG', in decreasing precedence
;		* if illegal or undecipherable, MEG is assumed
;	caldb	[INPUT] root $CALDB directory
;		* default: /data/draft_caldb/draft_caldb/
;		* CALDBhrma=CALDB/local/axaf/cal/hrma/cip/
;		* CALDBacis=CALDB/local/axaf/cal/acis/cip/
;		* CALDBhetg=CALDB/local/axaf/cal/hetg/cip/
;		* CALDBletg=CALDB/local/axaf/cal/letg/cip/
;	hrmaea	[INPUT] name of file containing HRMA effective areas
;		* if not given,
;			CALDBhrma/hrmaD1996-12-20effareaN0005.rdb
;		* must be RDB files with columns ENERGY, EA_[1,3,4,6], EA_HRMA
;	acissqe	[INPUT] array of files containing ACIS QEs
;		* if not given,
;			CALDBacis/aciss{ONCHIP}D1997-04-17qeobfN0002.dat
;		* must be ascii files, one for each chip, and must contain
;		  columns "Energy(keV)", "QE", "Filt_trans", "QE*Filt_trans"
;		  with the first line being the column headers.
;	grfH	[INPUT] name of file containing HEG grating efficiencies
;		* if not given,
;			CALDBhetg/hetghegD1996-11-01greffpr001N0003.rdb
;	grfM	[INPUT] name of file containing MEG grating efficiencies
;		* if not given,
;			CALDBhetg/hetgmegD1996-11-01greffpr001N0003.rdb
;	grfL	[INPUT] name of file containing LEG grating efficiencies
;		* if not given,
;			CALDBletg/letgD1996-11-01greffpr001N0003.rdb
;		* GRF[H,M,L] must be RDB files with columns ENERGY, GRFCOL
;	grfcol	[INPUT] the column name in the grating efficiencies files
;		* default behavior: if GRF[H,M,L] are set, GRFCOL is set
;		  to "O_0" else to "OZ"
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	pay careful attention to the formats of HRMAEA,ACISQE,GRF*
;	handles only the 0th order
;
;subroutines
;	RDB [STR_2_ARR, CREATE_STRUCT, KILROY]
;
;see also
;	ACISGARF,ACISS_GAPS,AXAF_WGRID,ARIA,EUVE_DS,EUVE_SW,EUVE_MW,EUVE_LW,
;	RD_PIMMS_FILE,RDARF,SMUDGE
;
;history
;	vinay kashyap (MarMM)
;	changed default for HRMAEA (VK; SepMM)
;-

message,'OBSOLETE!  Use mkacisgarf(...,order=0) instead',/informational

;	usage
if n_params() eq 0 then begin
  print,'effar=acisgarf0(wgrid,onchip,arm=arm,caldb=caldb,hrmaea=hrmaea,$'
  print,'      acissqe=acissqe,grfH=grfH,grfM=grfM,grfL=grfL,grfcol=grfcol)'
  print,'  return grating effective areas'
  return,0.
endif

;	check input
nw=n_elements(wgrid)
chipid=3 & if n_elements(onchip) gt 0 then chipid=fix(onchip(0))
if chipid lt 0 or chipid gt 5 then chipid=3	;default is ACIS-S3

;	keywords
garm='MEG' & sza=size(arm) & nsza=n_elements(sza)
if sza(nsza-2) eq 7 then begin
  if strpos(strlowcase(strtrim(arm(0),2)),'leg',0) ge 0 then garm='LEG'
  if strpos(strlowcase(strtrim(arm(0),2)),'heg',0) ge 0 then garm='HEG'
  if strpos(strlowcase(strtrim(arm(0),2)),'meg',0) ge 0 then garm='MEG'
endif
dbcal='/data/draft_caldb/draft_caldb/' & szc=size(caldb) & nszc=n_elements(szc)
if szc(nszc-2) eq 7 then dbcal=strtrim(caldb(0),2)

;{	initialize

;	$CALDB
CALDBaxafdata=dbcal+'/data/axaf/'
CALDBaxafpcf=dbcal+'/data/axaf/pcf/'
CALDBaxafcal=dbcal+'/local/axaf/cal/'
CALDBhrmacal=CALDBaxafcal+'hrma/'
CALDBaspectcal=CALDBaxafcal+'aspect/'
CALDBaciscal=CALDBaxafcal+'acis/'
CALDBhrccal=CALDBaxafcal+'hrc/'
CALDBhetgcal=CALDBaxafcal+'hetg/'
CALDBletgcal=CALDBaxafcal+'letg/'
CALDBxrcfcal=CALDBaxafcal+'xrcf/'

legmin=3.0 & legmax=180.0	;LEG 0th order wavelength range
legbin=0.0075			;LEG wavelength grid spacing
megmin=1.0 & megmax=42.0	;MEG 0th order wavelength range
megbin=0.005			;MEG wavelength grid spacing
hegmin=0.5 & hegmax=21.0	;HEG 0th order wavelength range
hegbin=0.0025			;HEG wavelength grid spacing
if nw eq 0 then begin
  legmin=legmin & legmax=legmax & legbin=legbin
  megmin=megmin & megmax=megmax & megbin=megbin
  hegmin=hegmin & hegmax=hegmax & hegbin=hegbin
  ;	reference wavelength grid
  nleg=(legmax-legmin)/legbin & wleg=findgen(nleg)*legbin+legmin
  nmeg=(megmax-megmin)/megbin & wmeg=findgen(nmeg)*megbin+megmin
  nheg=(hegmax-hegmin)/hegbin & wheg=findgen(nheg)*hegbin+hegmin
endif else begin
  nleg=nw & wleg=wgrid
  nmeg=nw & wmeg=wgrid
  nheg=nw & wheg=wgrid
endelse

;	mirror effective area
sz=size(hrmaea) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then hrmafil=strtrim(hrmaea(0),2) else $
	hrmafil=CALDBhrmacal+'cip/hrmaD1996-12-20effareaN0005.rdb'

;	QE data
sz=size(acissqe) & nsz=n_elements(sz) & nacisqe=n_elements(acissqe)
if sz(nsz-2) eq 7 then acisqefil=acissqe else acisqefil=$
	findfile(CALDBaciscal+'cip/aciss'+strtrim(chipid,2)+$
	'D1997-04-17qeobfN0002.dat',count=nacisqe)

;	grating efficiencies
sz=size(grfH) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then hegeff=strtrim(grfH(0),2) else $
	hegeff=CALDBhetgcal+'cip/hetghegD1996-11-01greffpr001N0003.rdb'
sz=size(grfM) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then hegeff=strtrim(grfM(0),2) else $
	megeff=CALDBhetgcal+'cip/hetgmegD1996-11-01greffpr001N0003.rdb'
sz=size(grfL) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then hegeff=strtrim(grfL(0),2) else $
	legeff=CALDBletgcal+'cip/letgD1996-11-01greffpr001N0003.rdb'

case garm of
  'LEG': begin
    ngrid=nleg & wgrid=wleg & greffil=legeff
    hrmacols=['ENERGY','EA_HRMA'] & addhrmacols=0
    if keyword_set(grfL) then gcol='o_0' else gcol='oz'
  end
  'MEG': begin
    ngrid=nmeg & wgrid=wmeg & greffil=megeff
    hrmacols=['ENERGY','EA_1','EA_3'] & addhrmacols=1
    if keyword_set(grfM) then gcol='o_0' else gcol='oz'
  end
  'HEG': begin
    ngrid=nheg & wgrid=wheg & greffil=hegeff
    hrmacols=['ENERGY','EA_4','EA_6'] & addhrmacols=1
    if keyword_set(grfH) then gcol='o_0' else gcol='oz'
  end
  else: begin
    message,'watchewtakinabout, willis?',/info
    return,0.
  end
endcase
if keyword_set(grfcol) then gcol=strtrim(grfcol(0),2)

;	initialzed}

;	file existence checks
ok='ok'
tmp=findfile(greffil,count=ngref)
tmp=findfile(hrmafil,count=nhrma)
if nacisqe eq 0 then ok='ACIS-S QE files not found'
if ngref eq 0 then ok=greffil+': Grating efficiency file missing'
if nhrma eq 0 then ok=hrmafil+': Mirror effective areas missing'
if ok ne 'ok' then begin
  message,ok,/info & return,0*wgrid
endif

;	read in CALDB files
;	HRMA EA
print,'Extracting columns: '+strjoin(hrmacols,', ')
hrma=rdb(hrmafil,cols=hrmacols)
;	Grating Efficiencies
print,'Extracting columns: ENERGY,'+gcol
gref=rdb(greffil,cols=['ENERGY',gcol])
;	QEs
nqelin=wc(acisqefil(0))
sNnrg=fltarr(nqelin-1L) & sNqef=fltarr(nqelin-1L)
print,'' & print,'Reading file '+acisqefil(0)
var=dblarr(4,nqelin-1L)
openr,uqe,acisqefil(0),/get_lun
line='' & readf,uqe,line & readf,uqe,var
;
sNnrg(0:nqelin-2L)=reform(var(0,*))
sNqef(0:nqelin-2L)=reform(var(3,*))
;
close,uqe & free_lun,uqe

;	interpolate HRMA to correct grid
hrmaw=12.3985/hrma.energy & hrmav=hrma.(2)
if addhrmacols then hrmav=hrmav+hrma.(3)
oo=sort(hrmaw) & hrmaw=hrmaw(oo) & hrmav=hrmav(oo)
vhrma=interpol(hrmav,hrmaw,wgrid)
;	interpolate Grating Efficiencies to correct grid
grefw=12.3985/gref.energy & grefv=gref.(2) & oo=sort(grefw)
grefw=grefw(oo) & grefv=grefv(oo)
vgref=interpol(grefv,grefw,wgrid)
;	interpolate QE to correct grid
qew=12.3985/sNnrg & qev=sNqef
oo=sort(qew) & qew=qew(oo) & qev=qev(oo)
vqe=interpol(qev,qew,wgrid)

;	make ARF
garf=vhrma*vgref*vqe

return,garf
end

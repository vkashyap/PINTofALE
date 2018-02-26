function acisgarf,offset,wgrid,order=order,arm=arm,caldb=caldb,$
	hrmaea=hrmaea,acissqe=acissqe,grfH=grfH,grfM=grfM,grfL=grfL,$
	grcol=grcol, _extra=e
;+
;function	acisgarf
;	returns effective area as a function of wavelength for specified
;	CXC grating and order on the ACIS-S
;
;syntax
;	arf=acisgarf(offset,wgrid,order=order,arm=arm,caldb=caldb,$
;	hrmaea=hrmaea,acissqe=acissqe,grfH=grfH,grfM=grfM,grfL=grfL,grcol=grcol)
;
;parameters
;	offset	[INPUT; required] offset of source position from the
;		nominal aim point on S3 [arcmin]
;		* +ve is towards the center of S3
;		* someday, after I figure out how to do it, I'll change
;		  the required inputs to (SRC_X,SRC_Y), and possibly ROLL.
;		  meanwhile, this parameter goes straight to ACISS_GAPS.
;	wgrid	[OUTPUT] wavelengths at which effective area is computed
;
;keywords
;	order	[INPUT] diffraction order.  if not set, assumes +ve 1st order.
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
;			CALDBacis/aciss?D1997-04-17qeobfN0002.dat
;		* must be ascii files, one for each chip, and must contain
;		  columns "Energy(keV)", "QE", "Filt_trans", "QE*Filt_trans"
;		  with the first line being the column headers.
;		* the files are matched to the chip simply by the order
;		  in which the files are sorted.  So there >must< be
;		  6 files, one for each ACIS-S chip, and in the right order,
;		  with ACISSQE[0] corresponding to S0, etc.
;	grfH	[INPUT] name of file containing HEG grating efficiencies
;		* if not given,
;			CALDBhetg/hetghegD1996-11-01greffpr001N0003.rdb
;	grfM	[INPUT] name of file containing MEG grating efficiencies
;		* if not given,
;			CALDBhetg/hetgmegD1996-11-01greffpr001N0003.rdb
;	grfL	[INPUT] name of file containing LEG grating efficiencies
;		* if not given,
;			CALDBletg/letgD1996-11-01greffpr001N0003.rdb
;		* GRF[H,M,L] must be RDB files with columns ENERGY, GRCOL+0..N
;		  where N.GE.ORDER
;	grcol	[INPUT] the prefix to the column names in the grating
;		efficiencies files
;		* default behavior: if GRF[H,M,L] are set,
;		- GRCOL is set to "O_" else to
;		- "OM" or "OP" depending on +ve or -ve ORDER.
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	pay careful attention to the formats of HRMAEA,ACISQE,GRF*
;	does not handle the 0th order
;
;subroutines
;	ACISS_GAPS
;	RDB [STR_2_ARR, CREATE_STRUCT, KILROY]
;
;see also
;	ACISGARF0,ACISS_GAPS,AXAF_WGRID,ARIA,EUVE_DS,EUVE_SW,EUVE_MW,EUVE_LW,
;	RD_PIMMS_FILE,RDARF,SMUDGE
;
;history
;	vinay kashyap (MIM.XI)
;	corrected for HEG and MEG's mirror shells (VK; MM.I)
;	bug correction -- added ORDER keyword to ACISS_GAPS call
;	  changed CALDB default (VK; FebMM)
;	added keywords HRMAEA,ACISSQE,GRFH,GRFM,GRFL,GRCOL (VK; MarMM)
;	changed HRMAEA default (VK; SepMM)
;-

message,'OBSOLETE!  Use mkacisgarf() instead',/informational

;	usage
noff=n_elements(offset) & ok='ok'
if noff eq 0 then ok='missing input' else $
 if noff gt 1 then ok='cannot handle array of OFFSETs'
if ok ne 'ok' then begin
  print,'Usage: effar=acisgarf(offset,wvlar,order=order,arm=arm,caldb=caldb,$'
  print,'  hrmaea=hrmaea,acissqe=acissqe,grfH=grfH,grfM=grfM,grfL=grfL,grcol=grcol)'
  print,'  return grating effective areas'
  if n_params() ne 0 then message,ok,/info
  return,0.
endif
srcpos=offset(0)

;	keywords
gord=1 & if n_elements(order) gt 0 then gord=fix(order(0))
if gord eq 0 then begin
  message,'cannot handle the 0th order.  look to ACISGARF0',/info
  return,0.
endif
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
if keyword_set(gord) then begin
  legmin=legmin/abs(gord) & legmax=legmax/abs(gord) & legbin=legbin/abs(gord)
  megmin=megmin/abs(gord) & megmax=megmax/abs(gord) & megbin=megbin/abs(gord)
  hegmin=hegmin/abs(gord) & hegmax=hegmax/abs(gord) & hegbin=hegbin/abs(gord)
endif

;	reference wavelength grid
nleg=(legmax-legmin)/legbin & wleg=findgen(nleg)*legbin+legmin
nmeg=(megmax-megmin)/megbin & wmeg=findgen(nmeg)*megbin+megmin
nheg=(hegmax-hegmin)/hegbin & wheg=findgen(nheg)*hegbin+hegmin

;	mirror effective area
sz=size(hrmaea) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then hrmafil=strtrim(hrmaea(0),2) else $
	hrmafil=CALDBhrmacal+'cip/hrmaD1996-12-20effareaN0005.rdb'

;	QE data
sz=size(acissqe) & nsz=n_elements(sz) & nacisqe=n_elements(acissqe)
if sz(nsz-2) eq 7 then acisqefil=acissqe else acisqefil=$
	findfile(CALDBaciscal+'cip/aciss?D1997-04-17qeobfN0002.dat',$
	count=nacisqe)

;	grating efficiencies
sz=size(grfH) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then hegeff=strtrim(grfH(0),2) else $
	hegeff=CALDBhetgcal+'cip/hetghegD1996-11-01greffpr001N0003.rdb'
sz=size(grfM) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then megeff=strtrim(grfM(0),2) else $
	megeff=CALDBhetgcal+'cip/hetgmegD1996-11-01greffpr001N0003.rdb'
sz=size(grfL) & nsz=n_elements(sz)
if sz(nsz-2) eq 7 then legeff=strtrim(grfL(0),2) else $
	legeff=CALDBletgcal+'cip/letgD1996-11-01greffpr001N0003.rdb'

;	chip gap locations
aciss_gaps,hgap,mgap,lgap,onchip=onchip,offset=srcpos,order=gord
;	input OFFSET in [arcmin], output ?GAP in [Ang]

case garm of
  'LEG': begin
    ngrid=nleg & wgrid=wleg & greffil=legeff & gaps=lgap
    hrmacols=['ENERGY','EA_HRMA'] & addhrmacols=0
    if keyword_set(grfL) then grfcol='o_' else begin
      if gord lt 0 then grfcol='om' else grfcol='op'
    endelse
  end
  'MEG': begin
    ngrid=nmeg & wgrid=wmeg & greffil=megeff & gaps=mgap
    hrmacols=['ENERGY','EA_1','EA_3'] & addhrmacols=1
    if keyword_set(grfM) then grfcol='o_' else begin
      if gord lt 0 then grfcol='om' else grfcol='op'
    endelse
  end
  'HEG': begin
    ngrid=nheg & wgrid=wheg & greffil=hegeff & gaps=hgap
    hrmacols=['ENERGY','EA_4','EA_6'] & addhrmacols=1
    if keyword_set(grfH) then grfcol='o_' else begin
      if gord lt 0 then grfcol='om' else grfcol='op'
    endelse
  end
  else: begin
    message,'watchewtakinabout, willis?',/info
    return,0.
  end
endcase
if keyword_set(grcol) then grfcol=strtrim(grcol(0),2)

;	initialzed}

;	file existence checks
ok='ok'
tmp=findfile(greffil,count=ngref)
tmp=findfile(hrmafil,count=nhrma)
if nacisqe lt 6 then ok='ACIS-S QE files not found'
if ngref eq 0 then ok=greffil+': Grating efficiency file missing'
if nhrma eq 0 then ok=hrmafil+': Mirror effective areas missing'
if ok ne 'ok' then begin
  message,ok,/info & return,0*wgrid
endif

;	read in CALDB files
;	HRMA EA
print,'Extracting columns: ENERGY,EA_HRMA
hrma=rdb(hrmafil,cols=hrmacols)
;	Grating Efficiencies
col2=grfcol+strtrim(abs(gord),2)
print,'Extracting columns: ENERGY,'+col2
gref=rdb(greffil,cols=['ENERGY',col2])
;	QEs
nqelin=0L & for i=0,5 do nqelin=nqelin > wc(acisqefil(i))
sNnrg=fltarr(6,nqelin-1L) & sNqef=fltarr(6,nqelin-1L)
print,''
for i=0,5 do begin
  print,'Reading file '+acisqefil(i)
  nlin=wc(acisqefil(i)) & var=dblarr(4,nlin-1L)
  openr,uqe,acisqefil(i),/get_lun
  line='' & readf,uqe,line & readf,uqe,var
  ;
  sNnrg(i,0:nlin-2L)=reform(var(0,*))
  sNqef(i,0:nlin-2L)=reform(var(3,*))
  ;
  close,uqe & free_lun,uqe
endfor

;	interpolate HRMA to correct grid
hrmaw=12.3985/hrma.energy & hrmav=hrma.(2)
if addhrmacols then hrmav=hrmav+hrma.(3)
oo=sort(hrmaw) & hrmaw=hrmaw(oo) & hrmav=hrmav(oo)
vhrma=interpol(hrmav,hrmaw,wgrid)
;	interpolate Grating Efficiencies to correct grid
grefw=12.3985/gref.energy & grefv=gref.(2) & oo=sort(grefw)
grefw=grefw(oo) & grefv=grefv(oo)
vgref=interpol(grefv,grefw,wgrid)
;;	interpolate QEs to correct grid
;qef=fltarr(6,ngrid)
;for i=0,5 do begin
;  qew=12.3985/reform(sNnrg(i,*)) & qev=reform(sNqef(i,*)) & oo=sort(qew)
;  qew=qew(oo) & qev=qev(oo)
;  vqe=interpol(qev,qew,wgrid) & qef(i,*)=vqe
;endfor

;	figure out which chips we have to deal with and multiply the
;	appropriate chip QE with HRMA effar and grating efficiency
garf=0.*wgrid
if gord lt 0 and onchip lt 0 then return,garf	;spectrum off of array
if gord gt 0 and onchip gt 5 then return,garf	;spectrum off of array
if gord lt 0 then begin			;(-ve orders
  for i=0,fix(onchip) do begin		;{for each chip..
    w1=gaps(2*i) & w0=0. & if onchip gt i then w0=gaps(2*i+1)
    ow=where(wgrid ge w0 and wgrid le w1,mow)
    if mow gt 0 then begin
      ; interpolate QE to correct grid
      qew=12.3985/reform(sNnrg(i,*)) & qev=reform(sNqef(i,*)) & oo=sort(qew)
      qew=qew(oo) & qev=qev(oo)
      vqe=interpol(qev,qew,wgrid)
      ; add to ARF
      garf(ow)=garf(ow) + vhrma(ow)*vgref(ow)*vqe(ow)
    endif
  endfor				;I=0,ONCHIP}
endif else begin			;)(+ve orders
  for i=fix(onchip),5 do begin		;{for each chip..
    w1=gaps(2*i+1) & w0=0. & if onchip lt i then w0=gaps(2*i)
    ow=where(wgrid ge w0 and wgrid le w1,mow)
    if mow gt 0 then begin
      ; interpolate QE to correct grid
      qew=12.3985/reform(sNnrg(i,*)) & qev=reform(sNqef(i,*)) & oo=sort(qew)
      qew=qew(oo) & qev=qev(oo)
      vqe=interpol(qev,qew,wgrid)
      ; add to ARF
      garf(ow)=garf(ow) + vhrma(ow)*vgref(ow)*vqe(ow)
    endif
  endfor				;I=ONCHIP,5}
endelse					;GORD)

return,garf
end

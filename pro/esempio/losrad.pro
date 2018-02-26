;+
;script	losrad
;	compute radiative loss function
;
;syntax
;	.run losrad
;
;warning
;	all internal program variables have names that begin as "v_" or "jnk_"
;	if you have any preexisiting variables with these names, they will be
;	overwritten.
;
;inputs
;	!LDBDIR	 line emissivity database
;	!CDBDIR	 continuum emissivity database directory
;	!CEROOT	 continuum emissivity database root
;	!IONEQF	 ion-balance file
;	!EDENS	 electron density
;	!ABUND	 abundances
;	!VERBOSE controls chatter
;	V_WMIN	 [1.0 Ang] minimum wavelength to consider
;	V_WMAX	 [3000.0 Ang] maximum wavelength to consider
;
;outputs
;	V_LLOGT	 temperatures at which radiative loss is computed
;	V_LOSS	 radiative loss function P(T)
;	V_LLOS	 loss function due to lines
;	V_CLOS	 loss function due to continuum processes
;	V_LLOSE	 loss function due to excitation transitions only
;	V_LLOSR	 loss function due to radiative transitions only
;	V_LLOSI	 loss function due to inner-shell ioniszation transitions only
;
;how to use
;	1. initialize the PINTofALE environment using INITALE
;	2. the parameters that will be set with INITALE are:
;	   --	ATOMIC: !LDBDIR, !IONEQF, !CDBDIR, !CEROOT, !CHIDIR
;	   --	STELLAR: !EDENS, !ABUND, !NH, !FH2, !HE1, !HEII
;	   verify that they are initialized to the right values
;	3. set the wavelength range of interest
;		V_WMIN = 1.0
;		V_WMAX = 3000.
;	4. run this script
;
;history
;	vinay kashyap (May99; modified from radloss.pro)
;	major cosmetic overhaul (VK; DecMM)
;	minor cosmetic changes (VK; Jul01)
;	changed default WMAX to 3000 AA (VK; Sep01)
;	plot now handles the 24-bit color case (VK; Apr03)
;	updated for IDL5.6 keyword_set([0]) behavior change for vectors
;	  (VK; 20Mar2006)
;	allowed for APED emissivities (VK; 6Aug2009)
;-

;	{	check that all inputs are defined

;2:
v='LDBDIR' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_LDBDIR,/getval
v='IONEQF' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_IONEQF,/getval
v='CDBDIR' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_CDBDIR,/getval
v='CEROOT' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_CEROOT,/getval
v='CHIDIR' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_CHIDIR,/getval
v='EDENS' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_EDENS,/getval
v='ABUND' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_ABUND,/getval
v='NH' & defsysv,'!'+v,exists=i
if i eq 0 then message,v+' undefined' else setsysval,v,v_NH,/getval
v='FH2' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'FH2',0.26,verbose=2
v='HE1' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'HE1',0.,verbose=2
v='HEII' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'HEII',0.,verbose=2
v='VERBOSE' & defsysv,'!'+v,exists=i
if i eq 0 then setsysval,'VERBOSE',5,verbose=2

;3:
if not keyword_set(v_WMIN) then begin
  v_WMIN=1.0
  message,'Setting minimum wavlenegth to '+strtrim(v_WMIN,2),/info
endif
if not keyword_set(v_WMAX) then begin
  v_WMAX=3000.0
  message,'Setting maximum wavlenegth to '+strtrim(v_WMAX,2),/info
endif

;	OK, all set!	}

;;	this routine cannot be used for APED databases, which already
;;	include BOTH ion balance and Anders & Grevesse abundances.
;if strpos(v_LDBDIR,'aped',0) ge 0 then begin
;  message,v_LDBDIR+': This appears to be an APED database.  Program cannot run.',/informational
;  message,'(ion balance and abundances already in emissivity tables)'
;endif

;	now define some useful variables
;	(these are used to figure out if databases must be reread)
if not keyword_set(jnk_LDBDIR) then jnk_LDBDIR=v_LDBDIR+' old'
if not keyword_set(jnk_CDBDIR) then jnk_CDBDIR=v_CDBDIR+' old'
if not keyword_set(jnk_CEROOT) then jnk_CEROOT=v_CEROOT+' old'
if not keyword_set(jnk_IONEQF) then jnk_IONEQF=v_IONEQF+' old'
if not is_keyword_set(jnk_EDENS) then jnk_EDENS=v_EDENS+1e10
if not keyword_set(jnk_ABUND) then jnk_ABUND=v_ABUND*0.9
if not keyword_set(jnk_WMIN) then jnk_WMIN=v_WMIN+100
if not keyword_set(jnk_WMAX) then jnk_WMAX=v_WMAX+100

;	obtain line contribution function
if (v_LDBDIR ne jnk_LDBDIR) or (v_IONEQF ne jnk_IONEQF) or $
   (v_EDENS ne jnk_EDENS) or (v_WMIN ne jnk_WMIN) or $
   (v_WMAX ne jnk_WMAX) then begin
  ;	read in line database
  v_lconf=rd_line(atom,n_e=v_EDENS[0],wrange=[v_WMIN,v_WMAX],$
	dbdir=v_LDBDIR[0],verbose=v_VERBOSE,$
	wvl=v_LWVL,logT=v_LLOGT,Z=v_Z,ion=v_ION,jon=v_JON,fstr=v_lstr) > 0.
  ;	if APED, remove Anders & Grevesse abundances from the emissivities
  if strpos(strlowcase(v_LDBDIR),'aped',0) ge 0 then apedance,v_lconf,v_Z
  ;	include ion balance if NOT APED
  if strpos(strlowcase(v_LDBDIR),'aped',0) lt 0 then begin
    if strpos(strlowcase(v_IONEQF),'none',0) lt 0 then begin
      v_lconf=fold_ioneq(v_lconf,v_Z,v_JON,chidir=v_CHIDIR,$
	logT=v_LLOGT,eqfile=v_IONEQF,verbose=v_VERBOSE) > 0.
      v_lstr.LINE_INT = v_lconf
    endif
  endif
endif

;	figure out how many lines are due to DR/radiative recombination,
;	excitation, and inner-shell ionization transitions
v_OOR=where(v_ION-v_JON lt 0,v_moor)	;recombination
v_OOE=where(v_ION-v_JON eq 0,v_mooe)	;excitation
v_OOI=where(v_ION-v_JON gt 0,v_mooi)	;iii

;	obtain continuum contribution function
if (v_CDBDIR ne jnk_CDBDIR) or (v_CEROOT ne jnk_CEROOT) or $
   (v_EDENS ne jnk_EDENS) or (v_WMIN ne jnk_WMIN) or $
   (v_WMAX ne jnk_WMAX) or (not arrayeq(v_ABUND,jnk_ABUND)) then begin
  ;	read in continuum database
  v_cconf=rd_cont(v_CEROOT[0],n_e=v_EDENS[0],wrange=[v_WMIN,v_WMAX],$
	dbdir=v_CDBDIR[0],abund=v_ABUND,verbose=v_VERBOSE,$
	wvl=v_CWW,logT=v_ClogT)
endif

;	P(T)
v_NT=n_elements(v_LLOGT)
v_LLOS=dblarr(v_nT) & v_CLOS=v_LLOS
v_LLOSR=v_LLOS & v_LLOSE=v_LLOS & v_LLOSI=v_LLOS
v_zab=v_ABUND(v_Z-1) & v_CDW=abs(v_CWW[1:*]-v_CWW)
for v_IT=0L,v_NT-1L do begin
  v_LLOS[v_IT]=total( v_lconf[v_IT,*] * v_zab )
  if v_MOOE gt 0 then v_LLOSE[v_IT]=total( v_lconf[v_IT,v_OOE] * v_zab[v_OOE] )
  if v_MOOR gt 0 then v_LLOSR[v_IT]=total( v_lconf[v_IT,v_OOR] * v_zab[v_OOR] )
  if v_MOOI gt 0 then v_LLOSI[v_IT]=total( v_lconf[v_IT,v_OOI] * v_zab[v_OOI] )
  v_CLOS[v_IT]=total( v_cconf[v_IT,*] * v_CDW )
endfor
v_LOSS = v_LLOS + v_CLOS

;	set the junk variables
jnk_LDBDIR=v_LDBDIR
jnk_CDBDIR=v_CDBDIR
jnk_CEROOT=v_CEROOT
jnk_IONEQF=v_IONEQF
jnk_EDENS=v_EDENS
jnk_WMIN=v_WMIN
jnk_WMAX=v_WMAX

;	make plots
print,"plot,v_LLOGT,v_LOSS,/xs,/ys,/nodata,/ylog,$"
print,"	yrange=max(v_LOSS)*[1e-4,1.5],$"
print,"	xtitle='[log!d10!n (Temperature [K])',$"
print,"	ytitle='P(T) [10!u-23!n ergs cm!u3!n s!u-1!n]',$"
print,"	title='Radiative Loss Function'"
print,"oplot,v_LLOGT,v_LOSS,color=v_COL1 & print,v_C1+': Radiative Loss'"
print,"oplot,v_LLOGT,v_LLOS,color=v_COL2 & print,v_C2+': lines only'"
print,"oplot,v_LLOGT,v_CLOS,color=v_COL3 & print,v_C3+': continuum only'"
print,"oplot,v_LLOGT,v_LLOSR,color=v_COL4 & print,v_C4+': DR/rad recombination'"
print,"oplot,v_LLOGT,v_LLOSI,color=v_COL5 & print,v_C5+': inner-shell ionization'"
;
print,"oplot,v_LLOGT,v_LLOSE,color=v_COL6 & print,v_C6+': excitation'
plot,v_LLOGT,v_LOSS,/xs,/ys,/nodata,/ylog,yrange=max(v_LOSS)*[1e-4,1.5],$
	xtitle='[log!d10!n (Temperature [K])]',$
	ytitle='P(T) [10!u-23!n ergs cm!u3!n s!u-1!n]',$
	title='Radiative Loss Function'
v_COL1=((!D.N_COLORS-2) > 1) mod 256L & v_C1='yellow' & setkolor,v_C1,v_COL1
v_COL2=((!D.N_COLORS-3) > 1) mod 256L & v_C2='red' & setkolor,v_C2,v_COL2
v_COL3=((!D.N_COLORS-4) > 1) mod 256L & v_C3='green' & setkolor,v_C3,v_COL3
v_COL4=((!D.N_COLORS-5) > 1) mod 256L & v_C4='blue' & setkolor,v_C4,v_COL4
v_COL5=((!D.N_COLORS-6) > 1) mod 256L & v_C5='gold' & setkolor,v_C5,v_COL5
v_COL6=((!D.N_COLORS-7) > 1) mod 256L & v_C6='magenta' & setkolor,v_C6,v_COL6
oplot,v_LLOGT,v_LOSS,color=v_COL1 & print,v_C1+': Radiative Loss'
oplot,v_LLOGT,v_LLOS,color=v_COL2 & print,v_C2+': lines only'
oplot,v_LLOGT,v_CLOS,color=v_COL3 & print,v_C3+': continuum only'
oplot,v_LLOGT,v_LLOSR,color=v_COL4 & print,v_C4+': DR/rad recombination'
oplot,v_LLOGT,v_LLOSI,color=v_COL5 & print,v_C5+': inner-shell ionization'
oplot,v_LLOGT,v_LLOSE,color=v_COL6 & print,v_C6+': excitation'

end

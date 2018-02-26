pro merge_line,dir1,dir2,outdir,atom=atom,logp=logp,n_e=n_e,prec=prec,$
	doubt=doubt,tdex=tdex,frat=frat,mish=mish,mash=mash,$
	omatch=omatch,batch=batch,_extra=e
;+
;procedure	merge_line
;	merges the spectral line information available in RD_LINE
;	readable format in two directories
;
;	there is no problem if none of the lines are common to the
;	two databases; the databases are simply concatenated.  if
;	there are common lines, then the following rules apply to
;	I. recognizing that two spectral lines are identical
;	*  first of all, the atomic number and ionic state of the 
;	   two lines in question must be identical.  if they are, and
;	   o if the wavelengths of the two lines are identical within a
;	     user specified precision then they are considered "strong"
;	     candidates for a match;
;	   o if the lines differ by greater than the stated precision, but
;	     "not by much" (controlled by another user specified keyword),
;	     the lines are "weak" candidates for a match.
;	*  the temperature response of the candidate lines, in particular
;	   the temperatures of maximum response and the intensity of a line
;	   summed over all temperatures are then compared.
;	   o if the temperature maxima are within a user specified fraction
;	     of each other AND the intensities summed over all temperature
;	     bins are within a factor 4,
;	     x "strong" candidates are accepted as identical
;	     x "weak" candidates are "provisionally" accepted as identical
;	   o if the temperature maxima are beyond the specified fraction AND
;	     the summed intensities are beyond a factor 4,
;	     x "strong" candidates are "provisionally" accepted as identical
;	     x "weak" candidates are considered distinct lines
;	*  for "provisionally" accepted matches, all the relevant data
;	   are displayed and the user is given a choice of accepting or
;	   rejecting the match.  the default action is to reject the
;	   match unless the total line intensities summed over temperature
;	   (where both profiles lie above a threshold=1e-10*max) are within
;	   a user specified fraction.
;	II. what to do with a matched set of spectral lines
;	*  the final database should not contain duplicates, so the
;	   information from one of the matched lines must be discarded.
;	   this is done arbitrarily by assuming that the first directory
;	   named in the parameter list contains the better information.
;	*  it may turn out that many lines from DIR2 match many lines
;	   from DIR1.  in that case, only the first of the best candidates
;	   are accepted.
;
;syntax
;	merge_line,dir1,dir2,outdir,atom=atom,logP=logP,n_e=n_e,$
;	prec=prec,doubt=doubt,tdex=tdex,frat=frat,mish=mish,mash=mash,$
;	omatch=omatch,/batch,wrange=wrange,/allah
;
;restrictions
;	* assumes that the temperature gridding of the two datasets are
;	  identical, and does *no corrections* (quits on inconsistency)
;	* requires subroutines
;	  -- KILROY
;	  -- RD_LINE [SYMB2ZION [LAT2ARAB]]
;	  -- CREATE_STRUCT
;
;side-effects
;	* if plotting device is X, then displays the temperature profile
;	  of all lines that are "weak" candidates.
;	* if OUTDIR is set, writes (possibly large) files to disk.
;
;parameters
;	dir1	[INPUT; required] directory which contains the database files
;	dir2	[INPUT; required] directory which contains the database files
;	outdir	[INPUT] output directory to hold merged files -- if not
;		specified, will not write anything to disk.
;
;keywords
;	atom	[INPUT; default: ALL] confine attention to specified atom
;		* passed w/o comment to RD_LINE
;		* from a practical viewpoint, always set this keyword!
;	logP	[INPUT; default: 15] log10(Pressure [cm^-3 K])
;	n_e	[INPUT] electron density [cm^-3].  if set,
;		* overrides LOGP
;		* appends a "D" to DIR1, DIR2, and OUTDIR
;	prec	[INPUT; default: 0.01] expected precision of the database
;		wavelength information; lines whose wavelengths lie within
;		PREC of each other are considered identical if they pass
;		all the other tests.
;	doubt	[INPUT; default: 1] room for doubt; creates a fuzzy
;		boundary for PREC -- any lines with wavelengths that
;		differ by > PREC but < (1+DOUBT)*PREC are also considered
;		identical (if they pass other tests)
;	tdex	[INPUT; default: 0.5] same as PREC, but for the temperatures
;		of maximum response.
;	frat	[INPUT; default: 3] fraction within which the summed
;		intensities of candidate lines must be for "tendency to
;		accept as match" (see description above)
;	batch	[INPUT] if set, all default actions are performed w/o
;		user input -- i.e., runs in batch mode.
;	mish	[I/O] index matching DIR1 to DIR2
;	mash	[I/O] position indices of extra lines in DIR2
;		* MISH and MASH may get OVERWRITTEN during the program
;		* if they are properly defined, then the matching
;		  algorithm is skipped.
;		* they must <<both>> be correctly defined or else they
;		  are ignored (and worse, overwritten!)
;		* use them to avoid redoing the matching (especially all
;		  those "weak" candidates -- see below)
;		* the quality of the match is uniformly set to 1 ("weak")
;		  unless OMATCH is correctly defined at input, in which
;		  case OMATCH.Q is left as is.
;		* don't even <<think>> of trying to define these arrays
;		  from scratch.  to use, first do
;		  IDL> merge_line,dir1,dir2,mish=mish,mash=mash,logp=lp1 ...
;		  at which time they get initialized, and may be used
;		  over and over,
;		  IDL> merge_line,dir1,dir2,mish=mish,mash=mash,logp=lp2 ...
;	omatch	[OUTPUT] structure that contains all the matched wavelengths
;		and corresponding fluxes
;	_extra	[INPUT] used to pass defined keywords to RD_LINE
;		(WRANGE, ALLAH, VERBOSE; no point in setting PRES!)
;
;history
;	vinay kashyap (Jan97)
;	added keyword N_E (VK; Jun97)
;	removed unwanted 'DD's going into RD_LINE (VK; Jul97)
;-

;	usage
np=n_elements(dir1)+n_elements(dir2)
if np ne 2 then begin
  print,'Usage: merge_line,dir1,dir2,outdir,atom=atom,logP=logP,n_e=n_e,$'
  print,'       prec=prec,doubt=doubt,tdex=tdex,frat=frat,mish=mish,mash=mash,$'
  print,'       omatch=omatch,/batch'
  print,'  merges spectral line information from DIR2 into DIR1 and writes'
  print,'  it into OUTDIR.  Also accepts defined keywords for RD_LINE:'
  print,'  WRANGE, ALLAH'
  return
endif

;	check inputs
if not keyword_set(logP) then logP=15
if not keyword_set(prec) then prec=0.01
if not keyword_set(doubt) then doubt=1
if not keyword_set(tdex) then tdex=0.5
if not keyword_set(frat) then frat=3. & if frat lt 1. then frat=1./abs(frat)
fuzz=prec*(1.+doubt)

;	read in info from DIR1 and DIR2
fdr1=1 & fdr2=1 & dr1=dir1 & dr2=dir2
tmp=rd_line(atom,dbdir=dir1,/desig,/econf,fstr=fdr1,logp=logp,n_e=n_e,_extra=e)
tmp=rd_line(atom,dbdir=dir2,/desig,/econf,fstr=fdr2,logp=logp,n_e=n_e,_extra=e)
;don't need next line *before* RD_LINE because RD_LINE takes care of "+'D'"s
if keyword_set(n_e) then begin & dr1=dr1+'D' & dr2=dr2+'D' & endif

;	unpack
;	MISH:j1 keeps track of whether wvl1 has a match
;	MASH:j2 stores the index of wvl1 at appropriate place
logt=fdr1.logt
z1=fdr1.z & ion1=fdr1.ion & wvl1=fdr1.wvl & fx1=fdr1.line_int
z2=fdr2.z & ion2=fdr2.ion & wvl2=fdr2.wvl & fx2=fdr2.line_int
nt=n_elements(logt) & n1=n_elements(z1) & n2=n_elements(z2)

;	initialize
k=0					;number distinct lines in DIR2
j1=lonarr(n1)-1				;index of matching line from DIR2
j2=lonarr(n2)-1				;index of unmatched new line from DIR2
s1=intarr(n1)				;match "strength" (0:no,1:weak,2:strong)
elem=['H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al','Si','P',$
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
					;elements from 1-30
rom=['I','II','III','IV','V','VI','VII','VIII','IX','X']
rom=[rom,'X'+rom,'XX'+rom,'XXX'+rom]	;roman numerals from 1-40
tx1=lonarr(n1) & tx2=lonarr(n2)		;index of maximum temperature responses
for i=0L,n1-1L do tx1(i)=(where(fx1(*,i) eq max(fx1(*,i))))(0)
for i=0L,n2-1L do tx2(i)=(where(fx2(*,i) eq max(fx2(*,i))))(0)
fxA=fltarr(n1) & fxB=fltarr(n2)		;intensities summed over temperatures
for i=0L,n1-1L do fxA(i)=total(fx1(*,i))
for i=0L,n2-1L do fxB(i)=total(fx2(*,i))
iskip=0					;skip matching algo (0) or not (1)?
fnum=logP(0) & If keyword_set(n_e) then fnum=alog10(n_e(0))	;suffix
;	output files
wfil='wvl' & tfil='tem' & ifil='ion' & sfil='src'
lfil='lvl' & efil='ecn'
if fnum ge 10 or fnum le 0 then ffil=string(fnum,'(f4.1)') else $
	ffil='0'+string(fnum,'(f3.1)')

;	catch trivial errors
	c1='No data found in '+dr1
if z1(0) eq 0 and z2(0) ne 0 then begin
  message,c1,/info
  iskip=1 & mish=[-1L] & mash=lindgen(n2)
  logt=fdr2.logt & nt=n_elements(logt)
endif
	c1='No data found in '+dr2
if z1(0) ne 0 and z2(0) eq 0 then begin
  message,c1,/info
  iskip=1 & mish=lonarr(n1)-1 & mash=[0L]
  logt=fdr1.logt & nt=n_elements(logt)
endif
	c1='No data found in either directory'
if z1(0) eq 0 and z2(0) eq 0 then begin
  message,c1,/info
  iskip=1 & mish=[-1L] & mash=[0L]
endif
if z1(0) ne 0 and z2(0) ne 0 then begin
	c1='Temperature grid mismatch'
  if n_elements(fdr1.logt) ne n_elements(fdr2.logt) then message,c1
	c1='Temperature mismatch'
  if total(abs(fdr1.logt-fdr2.logt)) gt 1e-7 then message,c1
endif

;	figure out if matching algorithm may be skipped
if keyword_set(mish) and keyword_set(mash) then begin
  iskip=1				;in good faith...
  if n_elements(mish) ne n1 then iskip=0	;MISH not right
  if n_elements(mash) ne n2 then iskip=0	;MASH not right
endif else iskip=0
;
if iskip eq 1 then begin
  j1=mish & j2=mash
  oo=where(j1 gt 0) & if oo(0) ne -1 then q=0*oo+1 else q=intarr(n1)
  szom=size(omatch) & nszom=n_elements(szom)
  if szom(nszom-2) eq 8 then begin
    if (tag_names(omatch))(6) eq 'Q' then q=omatch.q
  endif
  k2=where(j2 ge 0,k)
  goto, skipalgo	;ok, ready to skip
endif

;	step through each spectral line, checking
for i=0L,n2-1L do begin				;{check-chock, check-chock

  ;	initialize
  w1=abs(wvl1) & w2=abs(wvl2(i)) & dw=abs(w1-w2) & i1=-1L
  ;	the wavelength check
  ow=where(z2(i) eq z1 and ion2(i) eq ion1 and dw le prec)
  od=where(z2(i) eq z1 and ion2(i) eq ion1 and dw gt prec and dw le fuzz)
  ;	the temeperature check
  dt=abs(logt(tx1)-logt(tx2(i))) & rat=fxA/(fxB(i)>1e-8) & ss=lonarr(n1)
  if ow(0) ne -1 then begin
    ss(ow)=2				;"strong"
    ot=where(dt(ow) gt tdex or (rat(ow) lt 0.25 or rat(ow) gt 4.))
    if ot(0) ne -1 then ss(ow(ot))=1	;downgrade to "weak"
  endif
  if od(0) ne -1 then begin
    ss(od)=1				;"weak"
    ot=where(dt(od) gt tdex or (rat(od) lt 0.25 or rat(od) gt 4.))
    if ot(0) ne -1 then ss(od(ot))=0	;downgrade to "nope"
  endif
  ;	no polygamy please...
  op=where(ss eq 1 and s1 eq 2,mop)
  if mop gt 0 then ss(op)=0		;downgrade to "nope"

  o2=where(ss eq 2,mo2) & o1=where(ss eq 1,mo1)
  if mo2 gt 0 then mo1=0		;why bother with iffy when sure?

  ;	if sure of many, pick nearest
  if mo2 ge 1 then i1=(o2(sort(dw(o2))))(0)

  ;	if doubtful, ask
  if mo1 ge 1 then begin		;{doubt exists?
    i1=-1L		;index of match
    mrat=fltarr(mo1)	;array of emissivity ratios
    irat=intarr(mo1)	;accepted as match or not? 
    for io1=0,mo1-1 do begin		;{ask for each possible match
      ;	profile check over temperatures with non-zero emissivities ONLY
      tmp=[fx1(*,o1(io1)),fx2(*,i)] & yr=[0.,max(tmp)>1] & ymin=1e-10*yr(1)
      hh=where(fx1(*,o1(io1)) gt ymin and fx2(*,i) gt ymin,nhh)
      if nhh gt 0 then begin
        ratnum=total(fx1(hh,o1(io1))) & ratden=total(fx2(hh,i))
	rat=ratnum/ratden
      endif else rat=1e8
      if rat lt 1. then mrat(io1)=1./rat else mrat(io1)=rat	;min(mrat).GE.1

      ;	what info to display to ease selection?
      c1='1: '+elem(z1(o1(io1))-1)+rom(ion1(o1(io1))-1)+' ['+$
	strtrim(wvl1(o1(io1)),2)+'] Tmax='+string(logt(tx1(o1(io1))),'(f5.2)')
      c2='2: '+elem(z2(i)-1)+rom(ion2(i)-1)+' ['+strtrim(wvl2(i),2)+$
	'] Tmax='+string(logt(tx2(i)),'(f5.2)')
      c3=strtrim(i+1,2)+'('+strtrim(n2,2)+'): '+$
	(fdr1.desig)(0,o1(io1))+'('+(fdr2.desig)(0,i)+')'+' >><< '+$
	(fdr1.desig)(1,o1(io1))+'('+(fdr2.desig)(1,i)+')'
      c4=strtrim(io1+1,2)+'['+strtrim(mo1,2)+']'
      c5=(fdr1.config)(0,o1(io1))+'('+(fdr2.config)(0,i)+')'+' -X- '+$
        (fdr1.config)(1,o1(io1))+'('+(fdr2.config)(1,i)+')'+$
	' {'+strtrim(rat,2)+'}'

      ;	show and ask
      print,c3 & print,'('+strtrim(rat,2)+'); '+c1+' -x- '+c2
      ;;;
      ;;print,'remember to uncomment the if statement below!'
      ;;;
      if !d.name eq 'X' then begin
	plot,logt,fx1(*,o1(io1)),title=c1+' '+c2,subtitle=c3,yr=yr,/ys,$
	  xtitle='log!d10!n(T [K])',ytitle='10!u-23!n ergs cm!u3!n s!u-1!n'
	oplot,logt,fx2(*,i),psym=2
	if mo1 gt 1 then xyouts,logt(1),0.9*yr(0)+0.1*yr(1),c4,align=2
	xyouts,logt(1),0.1*yr(0)+0.9*yr(1),c5
      endif
      ;
      if rat ge 1./frat and rat le frat then begin	;accept match?
	accept='y' & c1=accept
	if not keyword_set(batch) then read,prompt='accept match? [y/n] ',c1
      endif else begin					;reject match?
	accept='n' & c1=accept
	if not keyword_set(batch) then read,prompt='accept match? [n/y] ',c1
      endelse

      if keyword_set(batch) then print,form="($,2a)",accept,string("15b)
      ;if keyword_set(batch) then wait,batch

      c2=strlowcase(strmid(strtrim(c1,2),0,1))
      if c2 eq 'n' or c2 eq 'y' then accept=c2
      if accept eq 'y' then irat(io1)=1

    endfor				;io1=0,mo1-1}
    o2=where(irat eq 1)
    if o2(0) ne -1 then begin		;pick the first of the best!
      minrat=min(mrat(o2),jrat) & i1=o1(o2(jrat))
    endif
  endif					;mo1 ge 1}

  if i1 ne -1 then begin				;{match found
    if s1(i1) gt 0 then begin			;{already matched
      ji=j1(i1)
      oldw=wvl2(ji) & dwo=abs(wvl1(i1)-oldw) & dwn=abs(wvl1(i1)-wvl2(i))
      ro=fxA(i1)/(fxB(ji)>1e-8) & rn=fxA(i1)/(fxB(i)>1e-8)
      if ro lt 1e-8 then ro=1e8 & if rn lt 1e-8 then rn=1e8
      if ro lt 1 then ro=1./ro & if rn lt 1 then rn=1./ro
      if s1(i1) eq 1 then begin			;(old match was "weak"
        if dwo gt dwn or ro gt rn then begin	;(worth exchanging?
          kilroy,dot='%'; made an exchange
          j1(i1)=i			;reset match to current
          j2(ji)=k & k=k+1		;append old "match" to line list
        endif else begin		;)(don't exchange, append
          kilroy,dot='&'; added to growing line list
          j2(i)=k & k=k+1
        endelse				;dwo.LE.dwn.AND.ro.LE.rn)
      endif else begin				;)(old match was "strong"
        kilroy,dot='+'; added to growing line list
        j2(i)=k & k=k+1
      endelse					;end s1(i1)=2)
    endif else begin				;}{new match
      kilroy,dot='*'; made a mark.
      j1(i1)=i & s1(i1)=ss(i1)
    endelse					;}
  endif	else begin					;}{no matches
    kilroy; was here.
    j2(i)=k & k=k+1
  endelse						;i1=-1}

endfor						;i=0,n2-1}

skipalgo: mish=j1 & mash=j2	;reset

;	so which are the common elements?
o1=where(j1 ge 0,m1) & o2=where(j2 lt 0,m2)
if m1 ne m2 then message,'bug in program!'
if m1 ne 0 then begin
  jj=j1(o1)
  w1=wvl1(o1) & w2=wvl2(jj) & g1=fx1(*,o1) & g2=fx2(*,jj)
  tmp1=fdr1.desig(*,o1)+' ('+fdr1.config(*,o1)+')'
  tmp2=fdr2.desig(*,jj)+' ('+fdr2.config(*,jj)+')'
  comm=tmp1+' / '+tmp2
  atm=z1(o1) & jon=ion1(o1) & ss=s1(o1)
  if iskip eq 1 then ss=q	;so that the quality of match is not strained
  omatch=create_struct('WVL1',w1,'INT1',g1,'WVL2',w2,'INT2',g2,$
  	'ELEM',elem(atm-1),'ION',rom(jon-1),'Q',ss,'DESCR',comm)
endif else omatch=create_struct('WVL1',0.,'INT1',0.,'WVL2',0.,'INT2',0.,$
	'ELEM',elem([z1-1]),'ION',rom([ion1-1]),'Q',s1,'DESCR','None found')

;	merge
k2=where(j2 ge 0,m2) & nw=n1+m2
if m2 ne k then message,'missing entries?'		;error check
;
wvl=fltarr(nw)
z=lonarr(nw) & ion=z & src=z
fx=dblarr(nt,nw)
tmp_trans=strarr(2,nw) & tmp_econf=tmp_trans
;
wvl(0:n1-1L)=wvl1 & z(0:n1-1L)=z1 & ion(0:n1-1L)=ion1 & src(0:n1-1L)=fdr1.src
;for i=0,n1-1 do begin
;  fx(*,i)=(fdr1.line_int)(*,i)
;  tmp_trans(*,i)=(fdr1.desig)(*,i)
;  tmp_econf(*,i)=(fdr1.config)(*,i)
;endfor
fx(*,0:n1-1L)=fdr1.line_int
tmp_trans(*,0:n1-1L)=fdr1.desig & tmp_econf(*,0:n1-1L)=fdr1.config
;
if m1 ne 0 then begin				;accumulate descriptions
  ;src(o1)=src(o1)+src(o2)
  for i=0,1 do begin
    tmp_trans(i,o1)=tmp_trans(i,o1)+' {'+(fdr2.desig)(i,j1(o1))+'}'
    tmp_econf(i,o1)=tmp_econf(i,o1)+' {'+(fdr2.config)(i,j1(o1))+'}'
  endfor
;  for i=0,m1-1 do begin
;    tmp_trans(*,o1(i))=tmp_trans(*,o1(i))+' {'+(fdr2.desig)(*,j1(o1)(i))+'}'
;    tmp_econf(*,o1(i))=tmp_econf(*,o1(i))+' {'+(fdr2.config)(*,j1(o1)(i))+'}'
;  endfor
endif
;
if k2(0) ne -1 then begin
  ii=lindgen(m2)+n1
  wvl(ii)=wvl2(k2) & z(ii)=z2(k2) & ion(ii)=ion2(k2) & src(ii)=fdr2.src(k2)
  fx(*,ii)=(fdr2.line_int)(*,k2)
  tmp_trans(*,ii)=(fdr2.desig)(*,k2) & tmp_econf(*,ii)=(fdr2.config)(*,k2)
endif
;if k2(0) ne -1 then begin
;  ii=lindgen(m2)+n1
;  wvl(ii)=wvl2(k2) & z(ii)=z2(k2) & ion(ii)=ion2(k2) & src(ii)=fdr2.src(k2)
;  for i=0,m2-1 do begin
;    fx(*,ii(i))=(fdr2.line_int)(*,k2(i))
;    tmp_trans(*,ii(i))=(fdr2.desig)(*,k2(i))
;    tmp_econf(*,ii(i))=(fdr2.config)(*,k2(i))
;  endfor
;endif

print,''
if n_elements(outdir) eq 0 then return		;no disk output

;if OUTDIR doesn't exist, create it
;	clearly, this bit works only on UNIX
odir=outdir & if keyword_set(n_e) then odir=odir+'D'
spawn,'csh -c "if \!( -e ' + odir +' ) mkdir '+string(odir)+'"'

zu=z(uniq(z,sort(z))) & nz=n_elements(zu)	;number of unique elements
for iz=0,nz-1 do begin				;{for each element
  ;	if any lines were found...
  oo=where(z eq zu(iz),moo)
  if zu(iz) eq 0 then moo=0L		;never mind
  if moo gt 0 then begin
    fpre=strcompress(elem(zu(iz)-1),/remove_all)	;file prefixes
    message,'	writing to '+odir+'/'+fpre+'_*',/info
    ;	open files
    openw,uw,odir+'/'+fpre+'_'+wfil,/get_lun		;wavelengths
    openw,ut,odir+'/'+fpre+'_'+tfil,/get_lun		;temperatures
    openw,ui,odir+'/'+fpre+'_'+ifil,/get_lun		;ionic state
    openw,uf,odir+'/'+fpre+'_'+ffil,/get_lun		;intensities
    openw,us,odir+'/'+fpre+'_'+sfil,/get_lun		;data source
    ;        odir+'/'+fpre+'_'+lfil			;level designations
    ;        odir+'/'+fpre+'_'+efil			;e configurations
    ;	write
    writeu,uw,moo,wvl(oo)
    writeu,ut,nt,logT
    writeu,ui,moo,ion(oo)
    writeu,uf,nt,moo,fx(*,oo)
    writeu,us,moo,src(oo)
    tmp=tmp_trans & tmp_trans=tmp_trans(*,oo)
    save,file=odir+'/'+fpre+'_'+lfil,tmp_trans & tmp_trans=tmp
    tmp=tmp_econf & tmp_econf=tmp_econf(*,oo)
    save,file=odir+'/'+fpre+'_'+efil,tmp_econf & tmp_econf=tmp
    ;	close files
    close,uw & free_lun,uw & close,ut & free_lun,ut
    close,ui & free_lun,ui & close,uf & free_lun,uf
    close,us & free_lun,us
  endif
endfor						;iz=0,nz-1}

return
end

function piffle,parfil,var,pardir=pardir,defdir=defdir,follow=follow,$
	v=v, _extra=e
;+
;function	piffle
;	return a structure containing the values of parameters listed
;	in specified CIAO compatible parameter file
;
;syntax
;	p=piffle(parfil,var,pardir=pardir,defdir=defdir,/follow,v=v)
;
;parameters
;	parfil	[INPUT; required] name of task or parameter file
;		* the ".par" extension is not necessary
;	var	[INPUT] if scalar string, returns only the value of
;		matching parameter.
;		* ignored if illegal
;
;keywords
;	pardir	[INPUT] director[y|ies] in which to look for PARFIL
;		* default: $PDIRS, $PFILES, $HOME/cxcds_param/
;	defdir	[INPUT] director[y|ies] in which to look as a last resort
;		* default: $PDIRS, $PFILES, /home/ascds/DS.daily/param/
;	follow	[INPUT] follow linked parameters (defined by -- not '(' --
;		a `)' in the parameter list)
;	v	[INPUT] if set, prints lots of stuff to the screen
;	_extra	[JUNK] here only to prevent crashing the program
;
;restriction
;	PARFIL must have parameter declarations in the following format:
;		name, type, mode, value, minimum, maximum, prompt
;	with "," being the field separator, and the following TYPEs
;		b: boolean, i: integer, r: real, s: string, f: filename
;	are recognized.  (see, e.g., /iraf/irafx/irafx/doc/clman.ms )
;
;subroutines
;	STR_2_ARR
;	CREATE_STRUCT
;
;history
;	vinay kashyap (MIM.XI.XXIX)
;	changed call to STR_2_ARR from STR2ARR (VK; MMV.IV)
;-

;	usage
ok='ok' & szp=size(parfil) & nszp=n_elements(szp)
if szp(nszp-2) ne 7 then ok='illegible parameter file' else $
 if szp(0) gt 0 then ok='cannot handle array of parameter files'
if ok ne 'ok' then begin
  print,'Usage: p=piffle(parfil,var,pardir=pardir,defdir=defdir,/follow)'
  print,'  return structure containing parameters and values'
  return,-1L
endif

;	inputs
paramf=strtrim(parfil(0),2)
c=str_sep(paramf,'.') & nc=n_elements(c) & xtn=''
if nc gt 1 then begin		;(is there an extension?
  cc='' & for i=0,nc-2 do cc=cc+c(i)+'.'	;reconstruct sans extension
  paramf=cc & xtn=c(nc-1L)
endif else xtn='.par'		;NC>1)
paramf=paramf+xtn
;
nvar=n_elements(var)
if nvar gt 1 then begin
  nvar=0	;ignore arrays
  message,'VAR is an array; cannot handle arrays, returning everything',/info
endif

;	keywords
pfiles=getenv('PFILES') & sc=str_sep(pfiles,';') & nsc=n_elements(sc)
cc=sc(0) & for i=1,nsc-1 do cc=cc+':'+sc(i) & pfiles=cc
c=str_sep(pfiles,':') & pfiles=c
pdirs=getenv('PDIRS') & sc=str_sep(pdirs,';') & nsc=n_elements(sc)
cc=sc(0) & for i=1,nsc-1 do cc=cc+':'+sc(i) & pdirs=cc
c=str_sep(pdirs,':') & pdirs=c
home=getenv('HOME')
pdir=[pdirs,pfiles,home+'/cxcds_param']
ddir=[pdirs,pfiles,'/home/ascds/DS.daily/param']
if keyword_set(pardir) then pdir=[strtrim(pardir,2),pdir]
if keyword_set(defdir) then ddir=[strtrim(defdir,2),ddir]
oo=where(pdir ne '',moo) & if moo gt 0 then pdir=pdir(oo)
oo=where(ddir ne '',moo) & if moo gt 0 then ddir=ddir(oo)
;
verbose=0 & if keyword_set(v) then verbose=fix(v)>1

;	find the parameter file
ifil=0 & kp=0L & kd=0L
while ifil eq 0 do begin		;{keep on finding..
 ; first look in PDIR
 if kp lt n_elements(pdir) then begin	;(first in PDIR..
   if verbose ge 3 then print,'looking for '+paramf+' in: '+pdir(kp)
   pfile=findfile(pdir(kp)+'/'+paramf,count=ifil) & kp=kp+1L
 endif else begin			;)(else in DDIR..
   if verbose ge 3 then print,'looking for '+paramf+' in: '+ddir(kd)
   pfile=findfile(ddir(kd)+'/'+paramf,count=ifil) & kd=kd+1L
   if kd eq n_elements(ddir) and not keyword_set(pfile) then begin
     message,paramf+': file not found',/info & return,-1L
   endif
 endelse				;KP)
endwhile				;IFIL v/s 0}
paramf=pfile(0)
message,'Reading from file: '+paramf,/info

;	read from parameter file
openr,upar,paramf,/get_lun		;{open parfil
while not eof(upar) do begin		;{read from parfil
  line='' & readf,upar,line
  c0=strmid(strtrim(line,2),0,1)
  if strtrim(line,2) eq '' then c0='#'	;empty line -- treat as comment
  if keyword_set(whybother) then c0='#'	;why bother if VAR's been found?
  if c0 eq '#' then goto,skipcomm	;(I know I can use an IF instead.

  ;	remove the `"'s
  iq=strpos(line,'"',0)
  while iq ge 0 do begin
    strput,line,' ',iq & iq=strpos(line,'"',0)
  endwhile

  ;	parse
  flds=str_sep(line,',') & nflds=n_elements(flds)
  if nflds ge 1 then pname=strtrim(flds(0),2)
  if nflds ge 2 then ptype=strtrim(flds(1),2)
  if nflds ge 3 then pmode=strtrim(flds(2),2)
  if nflds ge 4 then pval=strtrim(flds(3),2)
  if nflds ge 5 then pmin=strtrim(flds(4),2)
  if nflds ge 6 then pmax=strtrim(flds(5),2)
  if nflds ge 7 then pprompt=strtrim(flds(6),2)

  ;	if PVAL is listed in some other parameter file, go look for it
  if keyword_set(follow) then begin
    ;this here `(' is here for absolutely no reason whatsoever
    i0=strpos(pval,')',0)	;this tells variable is defined elsewhere
    if i0 ge 0 then begin	;(link to another file..
      i1=strpos(pval,'.',i0) & i2=strlen(pval)-1
      if i1 gt i0 then begin
	ffil=strmid(pval,i0+1,i1-i0-1)
	fval=strmid(pval,i1+1,i2-i1)
	fpar=piffle(ffil,fval,pardir=pardir,defdir=defdir,/follow)
	if n_tags(fpar) eq 1 then begin
	  ptype=fpar.(0).type
	  pmode=fpar.(0).mode
	  pval=fpar.(0).value
	  pmin=fpar.(0).min
	  pmax=fpar.(0).max
	  pprompt=fpar.(0).prompt
	endif else message,'could not follow '+pname+' to '+ffil+'.'+fval,/info
      endif else message,'could not follow to find '+pval,/info
    endif			;I0.GE.0)
  endif

  ;	decipher
  p0=strmid(ptype,0,1)	;treat fe/fn/fr/fw as simply f
  case p0 of
    'b': begin					;byte
      val='' & min='' & max=''
      if strlowcase(pval) eq 'yes' then val='yes'
      if strlowcase(pval) eq 'no' then val='no'
      if strlowcase(pmin) eq 'yes' then min='yes'
      if strlowcase(pmin) eq 'no' then min='no'
      if strlowcase(pmax) eq 'yes' then max='yes'
      if strlowcase(pmax) eq 'no' then max='no'
      pval=val & pmin=min & pmax=max
    end
    'i': begin					;integer
      val=str_2_arr(pval,/i4) & pval=val(0)
      min=str_2_arr(pmin,/i4) & pmin=min(0)
      max=str_2_arr(pmax,/i4) & pmax=max(0)
    end
    'r': begin					;real
      val=str_2_arr(pval) & pval=val(0)
      min=str_2_arr(pmin) & pmin=min(0)
      max=str_2_arr(pmax) & pmax=max(0)
    end
    'f': begin					;filename
      i0=strpos(pval,'$',0)
      if i0 ge 0 then begin
	i1=strlen(pval)-i0-1
	ii=strpos(pval,'/',i0) & if ii gt 0 then i1=ii-1
	env=getenv(strmid(pval,i0+1,i1-i0))
	rest=strmid(pval,i1,strlen(pval)-i1-1)
	if env ne '' then pval=env+rest
      endif
    end
    else: begin  				;assume string
      i0=strpos(pval,'$',0)
      if i0 ge 0 then begin
	i1=strlen(pval)-i0-1
	ii=strpos(pval,'/',i0) & if ii gt 0 then i1=ii-1
	env=getenv(strmid(pval,i0+1,i1-i0))
	rest=strmid(pval,i1,strlen(pval)-i1-1)
	if env ne '' then pval=env+rest
      endif
    end
  endcase

  if nvar eq 1 then begin
    if pname eq var(0) then whybother=1 else $	;why bother with rest of file?
	goto,skipcomm				;don't add to PARSTR
  endif

  ;	make structure
  tmp=create_struct('type',ptype, 'mode',pmode, 'value',pval,$
	'min',pmin, 'max',pmax, 'prompt',pprompt)
  if n_tags(parstr) eq 0 then parstr=create_struct(pname,tmp) else $
	parstr=create_struct(parstr,pname,tmp)

  skipcomm:				;GOTO comes here)
endwhile				;EOF(UPAR)}
close,upar & free_lun,upar		;close parfil}

;	report
if verbose ge 1 then begin
  pname=tag_names(parstr) & np=n_elements(pname)
  for i=0,np-1 do begin
    cc=pname(i)+' = '+strtrim(parstr.(i).value,2)
    if verbose ge 2 then cc=cc+' ('+parstr.(i).prompt+')'
    print,cc
  endfor
endif

return,parstr
end

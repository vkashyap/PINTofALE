pro setsysval,sysv,val,getval=getval,verbose=verbose,Read_Only=Read_Only,$
	_extra=e
;+
;procedure	setsysval
;	sets a system variable (the kind preceded by a "!") to the
;	specified value
;
;syntax
;	setsysval,sysv,val,/getval,verbose=v,/Read_Only
;
;parameters
;	sysv	[INPUT; required] a string scalar containing the
;		name of the system variable, written out without the
;		obligatory "!" (e.g., 'PATH')
;		* if !SYSV does not exist, creates it using defsysv
;		  inside an execute()
;	val	[I/O] if defined on input, will set the value of
;		!SYSV to this value, if !SYSV allows overwriting;
;		if not defined, or is illegal (e.g., wrong type of
;		input) then on output contains the value.
;		* if GETVAL is set, then regardless of all above considerations
;		  VAL contains the value of !SYSV on _output_
;
;keywords
;	getval   	[INPUT] if set, returns the existing value of
;			!SYSV in VAL
;	verbose 	[INPUT; default=0] controls chatter
;	Read_Only	[INPUT] if set, defines new !SYSV as read only
;	_extra		[JUNK] here only to prevent crashing the program
;
;subroutines
;	LEGALVAR
;
;history
;	vinay kashyap (OctMM)
;-

;	usage
ok='ok'
np=n_params() & ns=n_elements(sysv) & szs=size(sysv) & nszs=n_elements(szs)
if np eq 0 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SYSVAR undefined' else $
  if ns gt 1 then ok='cannot handle array of SYSVARs' else $
   if szs[nszs-2] ne 7 then ok='SYSVAR must be a string'
if ok ne 'ok' then begin
  print,'Usage: setsysval,sysvar,value,/getval,verbose=v,/Read_Only'
  print,'  sets the value of a system variable'
  if np ne 0 then message,ok,/info
  return
endif

;	inputs
v=0 & if keyword_set(verbose) then v=long(verbose[0]) > 1
;
ss=strupcase(strtrim(sysv[0],2)) & ifact=strpos(ss,'!',0)
if ifact eq 0 then begin
  if v gt 0 then message,'The prefixing "!" is unnecessary',/info
  ss=strmid(ss,1,strlen(ss)-1)
endif
ileg=legalvar(ss)
if ileg eq 0 then begin
  message,ss+': watchewtalkinaboot, willis?',/info & return
endif

;	is VAL given?
nv=n_elements(val)
if keyword_set(getval) then nv=0

;	check if !SYSV exists
i=0 & defsysv,'!'+ss,exists=i

;	if not..
if i eq 1 then begin	;(does exist
  if nv eq 0 then j=execute('val=!'+ss) else begin
    j=execute('!'+ss+'=val')
    if j eq 0 then begin
      if v gt 0 then message,'assignment to !'+ss+' failed',/info
      message,'returning existing value',/info
      j=execute('val=!'+ss)
    endif
  endelse
endif else begin	;cogito..)(does not exist
  if nv eq 0 then begin		;(VAL not given
    if v gt 0 then begin
      message,'!'+ss+': does not exist',/info
      message,'Nothing given to initialize it with; ignoring',/info
    endif
  endif else begin		;)(VAL is given
    if keyword_set(Read_Only) then defsysv,'!'+ss,val,Read_Only else $
	defsysv,'!'+ss,val
  endelse			;nV)
endelse			;vegetato..)

return
end

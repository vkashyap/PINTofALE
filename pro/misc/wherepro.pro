function wherepro,name
;+
;function	wherepro
;	returns the full path name to the source file that defines
;	the specified IDL function or procedure.  the function or
;	procedure must have been previously compiled, or else it
;	is not found.
;
;syntax
;	path=wherepro(name)
;
;parameters
;	name	[INPUT; required] scalar string containing the name
;		of the IDL function or procedure whose location must
;		be found
;
;keywords	NONE
;
;restrictions
;	because it depends on reading off the HELP output, the
;	function or procedure must have been already compiled,
;	or else it will not be found.  other procedures that do
;	the same job in different ways are:
;	  ROUTINE_INFO (built-in)
;	  FINDPRO (IDL-Astro)
;	  WHICH (PoA)
;	  WHEREIS (PoA; deprecated)
;
;history
;	vinay kashyap (FebMMVIII)
;-

;	usage
ok='ok' & np=n_params() & mp=n_elements(name) & szp=size(name,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if mp eq 0 then ok='NAME is undefined' else $
  if szp[0] ne 7 then ok='NAME is not a string'
if ok ne 'ok' then begin
  print,'Usage: path=wherepro(name)
  print,'  returns the full path name to the source file that defines
  print,'  the specified IDL function or procedure'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check the help
help,/source,output=temp

;	find it
cc=strtrim(strupcase(name[0]),2) & ncc=strlen(cc)
ii=where(strcmp(temp,cc,ncc) eq 1,mii)
if mii eq 0 then begin
  message,cc+': could not be found; compile and try again',/informational
  return,''
endif

;	extract it
ss=strarr(mii)
for i=0,mii-1 do ss[i]=(strsplit(temp[ii[i]],[' '],/extr))[1]

return,ss
end

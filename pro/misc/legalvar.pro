function legalvar,var,verbose=verbose, _extra=e
;+
;function	legalvar
;	returns 1 or 0 depending on whether VAR is a legal IDL variable
;	or not.  simply calls EXECUTE.  need to put this inside a subroutine
;	to protect possible existing variables of same name in main routine.
;
;syntax
;	l=legalvar(vars,verbose=verbose)
;
;parameters
;	var	[INPUT; required] string of candidate variable names
;		* if array, output will also be an array of 1s and 0s
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (MM.I)
;	bug fix if input is scalar (VK; MM.X)
;	added keyword verbose (VK; IMVIM.VIII)
;-

;	usage
ok='ok' & np=n_params()
nv=n_elements(var) & szv=size(var) & nszv=n_elements(szv)
if np eq 0 then ok='no' else $
 if nv eq 0 then ok='input undefined' else $
  if szv[nszv-2] ne 7 then ok='input should be a character string'
if ok ne 'ok' then begin
  print,'Usage: l=legalvar(vars,verbose=verbose)'
  print,'  checks whether a given string fragment is a legal IDL variable'
  if ok ne 'no' then message,ok,/info
  return,0
endif

;	define output
l=0 & if szv[0] ne 0 then l=make_array(size=szv,/int,value=0)

;	check for legality
for i=0L,nv-1L do begin
  l[i]=execute(var[i]+'=0')
  if l[i] eq 0 and keyword_set(verbose) then $
	jnk=execute('message,"<'+var[i]+'> : not a legal variable",/informational')
endfor

return,l
end

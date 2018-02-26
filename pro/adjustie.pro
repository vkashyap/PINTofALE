pro adjustie,a,ties=ties,vname=vname,$
	ifit=ifit,itie=itie,tieto=tieto,tieops=tieops,$
	_extra=e
;+
;procedure	adjustie
;	adjust the parameters to account for any ties among them, and
;	return the corrected set
;
;syntax
;	adjustie,a,ties=ties,vname=vname
;
;parameters
;	a	[I/O; required] parameters to be updated
;
;keywords
;	ties	[INPUT] string array describing how a parameter is tied to
;		other parameters.
;		* must be a mathematical expression of the form
;		  "a0 = a1/10.", "a2=a3+a4", "a5 = a5 < 0", etc.
;	vname	[INPUT; default='a'] in case TIES contains variable names
;		such as "p0=p2^2", etc, then set VNAME='p'
;
;		NOTE:  the following keywords are generated in situ the first
;		time the function is called.  on subsequent calls they are
;		used instead of EXECUTEing TIES in order to improve speed.
;
;		NOTE2: actually, they are not implemented yet.  someday.
;
;	ifit	[I/O] integer array describing which parameters are frozen
;		(0's) and which are not (1's).  updated upon exit to include
;		information given in TIES
;	itie	[I/O] integer array showing which parameters are tied.
;	tieto	[I/O] linked list showing which parameters each of the
;		ITIE are tied to which parameter
;	tieops	[I/O] as TIETO, but containing information on which operator
;		is used, e.g., "+", "-", "*", "/", "<", ">", etc.
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Oct98)
;	converted to IDL5 (VK; OctMM)
;-

;	usage
na=n_elements(a)
if na eq 0 then begin
  print,'Usage: adjustie,a,ties=ties,vname=vname'
  print,'  adjust parameters to include inter-parameter ties'
  return
endif

nt=n_elements(ties) & if nt eq 0 then return		;nothing to do!
if strtrim(ties[0],2) eq '0' then return		;nothing to do!

varn='a' & nv=n_elements(vname)
if nv eq 1 then begin
  if strtrim(vname[0],2) ne '' then varn=strtrim(vname[0],2)
endif

;	initialize
for i=0,na-1 do tmp=execute(varn+strtrim(i,2)+'=a[i]')

for i=0,nt-1 do begin			;{step through and decipher
  c=strlowcase(strtrim(ties[i],2)) & lc=strlen(c) & bc=byte(c)
  j=1 & if c ne '' then j=execute(c)
  if j ne 1 then message,'expression: '+c+'  not understood',/info
endfor					;I=0,NT-1}

;	output
for i=0,na-1 do tmp=execute('a[i]='+varn+strtrim(i,2))

return
end

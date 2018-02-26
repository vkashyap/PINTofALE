function whereis,file,lookin=lookin,dironly=dironly,nospawn=nospawn,$
	verbose=verbose, _extra=e
;+
;function	whereis
;	apes the csh command `which` in the IDL environment, and
;	finds the path to the specified file by looking through
;	all the directories specified in the !PATH variable.
;
;	NOTE: it is possible to use ROUTINE_INFO to get the path
;	to a compiled subroutine directly:
;	  help,routine_info(routine,/source,/function),/structure
;	but clearly you need to know the exact name of the routine
;	(can't use wildcards) as well as whether what you need is
;	a procedure or a function beforehand.  also, they have to
;	be compiled procedures, not, for example, command files or
;	scripts or datafiles lurking under the horizon.  and what
;	if there are duplicates in different directories?
;
;	btw, FILEPATH() is completely uselss for this purpose and
;	FINDFILE() is of limited use because you still need to
;	explicitly say which directory to look in.
;
;syntax
;	fpath=whereis(file,lookin=lookin,/dironly,/nospawn,verbose=verbose)
;
;parameters
;	file	[INPUT; required] name(s) of file to search for
;		* may include shell wildcards, regexps, etc.
;		  (but remember to prepend the escape character "\")
;		* may also be an array.
;	path	[OUTPUT] full path to the requested file(s)
;		* there is no reason to expect that the array sizes of
;		  FILE and PATH will match.
;
;keywords
;	lookin	[INPUT] path specification in which to look for FILE
;		* this will be passed without comment to EXPAND_PATH()
;	dironly	[INPUT] if set, returns only the name(s) of the
;		directory(ies) containing the specified file(s)
;	nospawn	[INPUT] if set, acts as wrapper for FINDFILE()
;		even on UNIX systems, i.e., does spawn find
;	verbose	[INPUT] if set, spits out informational messages
;	_extra	[JUNK] here only to prevent crashing the program
;
;notes
;	spawns "find" on UNIX systems, and calls FINDFILE() otherwise
;
;history
;	vinay kashyap (Dec2001)
;-

forward_function filepath	;for pre 5.2 versions

;	usage
ok='ok' & np=n_params() & nf=n_elements(file) & path=''
szf=size(file) & nszf=n_elements(szf)
if np eq 0 then ok='Missing filename' else $
 if nf eq 0 then ok='Input filename undefined' else $
  if szf(nszf-2) ne 7 then ok='Filename not character string' else $
   if !version.OS_FAMILY ne 'unix' and $
   float(strmid(!version.release,0,3)) lt 5 then ok='UNIX only: spawns find'
if ok ne 'ok' then begin
  print,'Usage: fpath=whereis(file,path,/dironly,/nospawn,verbose=verbose)'
  print,'  search through !PATH to find file'
  print,'  see also: which, findfile(), routine_info()'
  if np ne 0 then message,ok,/info
  return,''
endif

;	recast inputs
ff=[file(*)]
vv=0 & if keyword_set(verbose) then vv=long(verbose(0)) > 1
tospawn=0 & if !version.OS_FAMILY eq 'unix' then tospawn=1
if keyword_set(nospawn) then tospawn=0

;	expand !path
path=!path
if keyword_set(lookin) then begin
  path=lookin(0)
  for i=1L,n_elements(lookin)-1 do path=path+':'+lookin(i)
endif
dirs=expand_path(path,/array) & ndirs=n_elements(dirs)

;	define the output
path=''

for i=0L,nf-1L do begin		;{for each file
  if vv ge 1 then kilroy,dot=strtrim(ff(i),2)+': '
  if keyword_set(tospawn) then begin		;(UNIX; spawn find
    cmd='find' & for j=0L,ndirs-1L do cmd=cmd+' '+dirs(j)
    cmd=cmd+' -name "'+strtrim(ff(i),2)+'"'
    if keyword_set(dironly) then cmd=cmd+' -exec dirname {} \;' else $
      cmd=cmd+' -print'
    spawn,cmd,cc
    if vv ge 2 then print,cc
    if keyword_set(path) then path=[path,cc] else path=cc
  endif else begin				;UNIX)(use FINDFILE()
    for j=0L,ndirs-1L do begin
      c1=filepath(strtrim(ff(i),2),root_dir=dirs(j))
      cc=findfile(c1,count=nfil)
      if nfil gt 0 then begin
	if keyword_set(dironly) then begin
	  if keyword_set(path) then path=[path,dirs(j)] else path=dirs(j)
	endif else begin
	  if keyword_set(path) then path=[path,cc] else path=cc
	endelse
      endif
    endfor
  endelse					;no spawn)
endfor				;I=0,NF-1}

;	get rid of duplicates
path=path(uniq(path,sort(path)))

if np eq 1 then begin
  if keyword_set(path) then print,path else print,file+': Not found'
endif

return,path
end

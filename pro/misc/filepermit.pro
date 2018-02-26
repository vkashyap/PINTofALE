pro filepermit,dir,plus=plus,world=world,group=group,exclude=exclude,$
	opts=opts,dironly=dironly,nouser=nouser,verbose=verbose, _extra=e
;+
;procedure	filepermit
;	a wrapper for find.  removes (or adds) group and world rwx
;	permissions as appropriate (i.e., duplicates u permissions)
;	to all files in the hierarchy
;
;syntax
;	filepermit,dir,plus=plus,world=world,group=group,exclude=exclude,$
;	opts=opts,/dironly,/nouser,verbose=verbose
;
;parameters
;	dir	[INPUT; required] topmost level of directory structure
;		* if array of directories, runs the same commands on each
;
;keywords
;	plus	[INPUT] if set, _adds_ group rwx permissions to all files
;		as appropriate
;		* the default is to _remove_ group rwx permission
;	world	[INPUT] if set, _adds_ world rwx permissions to all file
;		as appropriate
;		* the default is to _remove_ setting of world permissions
;	group	[INPUT] if set to a scalar string, changes group name
;		to the specified GROUP
;	exclude	[INPUT] if set to a scalar string, excludes the given
;		file pattern from the searches
;		* wildcards must be escaped with a "\"
;	opts	[INPUT] if set to a scalar string, assumes that what is
;		given are extra options to FIND and tacks them on just
;		prior to the -type command
;		* unless you know exactly what you are doing, avoid
;		  -user, -perm, -exec, -type
;		* example:
;		  to avoid traversing links, include OPTS='! -type l'
;		  to echo each file found to STDOUT, use OPTS='-ls'
;		  to check only files modified within the past week,
;		  	use OPTS='-mtime -7'
;	dironly	[INPUT] if set, operates on directories only
;	nouser	[INPUT] if set, does not include "-user WHOAMI" in call
;	verbose	[INPUT] controls chatter
;
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	spawns the following instances of find:
;		find DIR \
;		-user WHOAMI \
;		[ OPTS ] \
;		\( ! -name EXCLUDE \)  \
;		[ -type d ] \
;		-perm -u+[rwx] \
;		-exec chmod g[-+][rwx],o[-+][rwx] {} \; \
;		-exec chgrp GROUP {} \;
;
;restrictions
;	works only on UNIX
;	cannot EXCLUDE directories -- must use
;	  /NOUSER,OPTS='\( -name DIR_TO_EXCLUDE -prune -o -print \)'
;	or something of that sort
;
;history
;	vinay kashyap (MMIVfeb; yes, a shell or perl script would
;	have been more appropriate, but I found it easier to handle
;	the flexibility better in IDL)
;-

;	usage
ok='ok' & np=n_params() & ndir=n_elements(dir) & szd=size(dir,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if szd ne 7 then ok='DIR must be a character string' else $
  if strlowcase(!version.os_family) ne 'unix' then ok='works only on UNIX'
if ok ne 'ok' then begin
  print,'Usage: filepermit,dir,plus=plus,world=world,group=group,$'
  print,'       exclude=exclude,opts=opts,/dironly,/nouser,$'
  print,'       verbose=verbose'
  print,'  wrapper for find.  removes/adds group and world rwx permissions'
  print,'  as appropriate (i.e., duplicates u permissions)'
  if np ne 0 then message,ok,/informational
  return
endif

;	figure out keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
groupname='' & szg=size(group,/type) & if szg eq 7 then groupname=group[0]
;
pattern='' & szx=size(exclude,/type)
if szx eq 7 then begin
  nx=n_elements(exclude) & for i=0L,nx-1L do pattern=pattern+' ! -name '+exclude[i]
endif
;
onlydir=0 & if keyword_set(dironly) then onlydir=1
;
opsfind='' & szf=size(opts,/type) & nopts=n_elements(opts)
if szf eq 7 then begin
  if nopts eq 1 then opsfind=' '+opts[0] else $
	opsfind=' '+strjoin(opts,' ',/single)
endif

;	who am I?
spawn,'whoami',whoami

;	construct the commands
cmd='find '+dir[0]			;find DIR
if not keyword_set(nouser) then $
	cmd=cmd+' -user '+whoami[0]	;-user WHOAMI
if keyword_set(pattern) then $
	cmd=cmd+' \( '+pattern+' \) '	;! -name EXCLUDE
cmd=cmd+opsfind				;OPTS
if keyword_set(onlydir) then $
	cmd=cmd+' -type d'		;-type d
cmdr=cmd+' -perm -u+r'			;-perm -u+r
cmdw=cmd+' -perm -u+w'			;-perm -u+w
cmdx=cmd+' -perm -u+x'			;-perm -u+x
if not keyword_set(plus) then begin
  cmdz=' -exec chmod g-w,o-w {} \;'	;-exec chmod g-w,o-w {} \;
endif else begin
  cmdzr=' -exec chmod g+r'
  if keyword_set(world) then cmdzr=cmdzr+',o+r'
  cmdzr=cmdzr+' {} \;'			;-exec chmod g+r[,o+r] {} \;
  cmdzw=' -exec chmod g+w'
  if keyword_set(world) then cmdzw=cmdzw+',o+w'
  cmdzw=cmdzw+' {} \;'			;-exec chmod g+w[,o+w] {} \;
  cmdzx=' -exec chmod g+x'
  if keyword_set(world) then cmdzx=cmdzx+',o+x'
  cmdzx=cmdzx+' {} \;'			;-exec chmod g+x[,o+x] {} \;
endelse
cmdg='' & if keyword_set(groupname) then $	;-exec chgrp GROUP {} \;
	cmdg=' -exec chgrp '+groupname+' {} \;'

;	spawn the commands
if vv gt 0 then begin
  print,'spawning ...'
  if not keyword_set(plus) then begin
    print,cmd+cmdz+cmdg
  endif else begin
    print,cmdr+cmdzr+cmdg
    print,cmdw+cmdzw
    print,cmdx+cmdzx
  endelse
endif
if vv gt 1 and vv lt 10 then wait,5
if vv ge 10 then begin
  c='type any key to continue, q or z to halt' & print,c
  c=get_kbrd(1)
  if strlowcase(c) eq 'q' or strlowcase(c) eq 'z' then $
	stop,'HALTing.  type .CON to continue'
endif
if not keyword_set(plus) then begin
  spawn,cmd+opsfind+cmdz+cmdg
endif else begin
  spawn,cmdr+opsfind+cmdzr+cmdg
  spawn,cmdw+opsfind+cmdzw
  spawn,cmdx+opsfind+cmdzx
endelse

return
end

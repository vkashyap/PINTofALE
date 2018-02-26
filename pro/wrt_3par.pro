pro wrt_3par,x,y,pardir=pardir,root=root,add=add,missing=missing,$
	fwhm=fwhm,norm=norm, _extra=e
;+
;function	wrt_3par
;	choose parameters describing 3-parameter curve interactively
;	and write/modify parameter files
;
;syntax
;	wrt_3par,x,y,root=root,pardir=pardir,add=add,missing=missing,$
;	type=type,/asis,/fwhm,/norm,beta=betap
;
;parameters
;	x	[INPUT; required] points at which the "reference
;		curve" is defined
;	y	[INPUT] the reference curve (usually data)
;		if not given, X is taken to be Y, and
;		the indices of X are taken to be X!
;
;keywords
;	pardir	[INPUT; default: .] directory in which to look for
;		and write parameter files
;	root	[INPUT; default: gauss] prefix for parameter files
;	add	[INPUT] full path name with prefix of a set of
;		previously existing parameter files that must be
;		appended to the interactively chosen points
;		* PARDIR will not be prefixed to ADD
;	missing	[INPUT] 3-element array [p,s,h] containing defaults
;		to override hardcoded defaults to use in case any
;		of the values and/or files are missing
;	_extra	[INPUT] pass defined keywords to subroutines:-
;		MK_3MODEL: TYPE, ASIS, ALLCOMP
;		MK_GAUSS: MISSING, FWHM, NORM
;		MK_LORENTZ: BETAP, MISSING, NORM
;
;usage summary
;	* call as a procedure
;	* requires a reference curve (y[x]) as input
;			
;parameter file format
;	SEE file RD_3PAR.PRO
;	briefly:  3 files in pardir --
;		root_pos.par
;		root_wdt.par
;		root_hgt.par
;	with each line describing an independent component.
;	multiple entries on single line indicates linked components.
;
;side effects
;	creates 3 files pardir/root_{pos,wdt,hgt}.par and overwrites them
;	  if they already exist
;	opens a plot window and takes cursor input, and is set to work
;	  for X-windows only!
;
;subroutines
;	RD_3PAR
;	MK_3MODEL
;	MK_GAUSS
;	MK_LORENTZ
;
;history
;	vinay kashyap (Oct 96)
;	added _EXTRA, removed keywords FWHM, NORM (VK; Oct98)
;	button press status now stored in !MOUSE, not !ERR (VK; Apr09)
;-

np=n_params(0)
if np eq 0 then begin
  print, 'Usage: wrt_3par,x,y,root=root,pardir=pardir,add=add,missing=missing,$'
  print, '       type=type,/asis,/fwhm,/norm,betap=betap'
  print, '  choose parameters for 3-parameter curve components interactively'
  print, '  and write to parameter files'
  return
endif

;check keywords
if not keyword_set(root) then root='gauss' else root=strtrim(root,2)
nx=n_elements(x) & x0=x(nx/2) & mxx=max(x,min=mnx)
if not keyword_set(missing) then missing=[x0,0.1*(mxx-mnx),1.]

;if Y is not given, take x to be Y
xx=x & if np eq 1 then begin & y=x & xx=lindgen(nx) & endif

if !d.name ne 'X' then begin
  print,'sorry, not implemented for anything but X windows'
  return
endif

;plot
plot,xx,y,psym=10,xr=[x(0),x(nx-1)],/xs,yr=[min(y),1.1*max(y)],/ys

;read in from previously set parameters (if any) and overplot
if keyword_set(add) then begin
  n=rd_3par(p,w,h,group,delp,root=add,missing=missing, _extra=e)
  if n gt 0 then begin
    m=mk_3model(xx,p,w,h,group,delp,missing=missing, _extra=e)
    oplot,xx,m
  endif
endif

idlvers=float(!VERSION.release) & if idlvers le 0 then idlvers=7
;now read in positions, heights, and widths interactively
pp=[0.] & ww=[0.] & hh=[0.] & go_on=1
print,'click mouse-button to denote location and peak, drag to denote width'
print,'	<left> continue; <middle> undo previous; <right> exit'
while go_on eq 1 do begin
  cursor,x0,y0,/down,/data
  if idlvers lt 5 then mbutton=!err else mbutton=!MOUSE.button
  cursor,x1,y1,/up,/data
  if mbutton eq 2 then begin			;middle button
    nn=n_elements(pp)
    if nn gt 1 then begin			;eliminate last point
      oplot,[pp(nn-1)],[hh(nn-1)],psym=2,symsize=2
      oplot,pp(nn-1)+ww(nn-1)*[1.,-1.],hh(nn-1)*[1,1],col=0
      pp=pp(0:nn-2) & ww=ww(0:nn-2) & hh=hh(0:nn-2)
    endif
  endif else begin
    if mbutton eq 4 then begin			;right button
      go_on=0					;quit
    endif else begin				;left button
      pp=[pp,x0] & hh=[hh,y0] & ww=[ww,abs(x1-x0)]
      oplot,[x0],[y0],psym=6 & oplot,x0+abs(x1-x0)*[1.,-1.],y0*[1.,1.]
      c1='Position = '+strtrim(x0,2)
      c1=c1+' Width = '+strtrim(abs(x1-x0),2)
      c1=c1+' Height = '+strtrim(y0,2) & print,c1
    endelse
  endelse
endwhile
nn=n_elements(pp)-1
if nn lt 1 then begin
  c1='oh ok.  here i am, slaving away at keeping track of *everything*,'
  c2='and computing away like crazy...do you know there is a terrible pain'
  c3='in all the diodes down my left side?  but that is ok.  really.  just'
  c4='leave me alone.  i will be fine.'
  message,c1,/info & message,c2,/info & message,c3,/info & message,c4,/info
  return
endif
pp=pp(1:*) & ww=ww(1:*) & hh=hh(1:*)
gg=indgen(nn)+1 & dp=0.*pp		;each component in a group by itself
if keyword_set(missing) then begin	;correct for missing width
  oo=where(ww eq 0) & if oo(0) ne -1 then ww(oo)=missing(1)
endif

;combine the interactively obtained list with the list read in from file
if keyword_set(p) then begin
  pp=[p,pp] & ww=[w,ww] & hh=[h,hh] & gg=[group,gg+max(group)] & dp=[delp,dp]
endif

;generate new model and overplot
m=mk_3model(xx,pp,ww,hh,gg,dp,missing=missing, _extra=e)
oplot,xx,m,linestyle=1

;figure out the output filenames and open for writing
pfile=root+'_pos.par' & if keyword_set(pardir) then pfile=pardir+'/'+pfile
wfile=root+'_wdt.par' & if keyword_set(pardir) then wfile=pardir+'/'+wfile
hfile=root+'_hgt.par' & if keyword_set(pardir) then hfile=pardir+'/'+hfile
openw,upf,pfile,/get_lun
openw,uwf,wfile,/get_lun
openw,uhf,hfile,/get_lun

;figure out how many unique groups
grp=gg(uniq(gg,sort(gg))) & ng=n_elements(grp)

;for each group, collect and write out
for i=0,ng-1 do begin
  oo=where(gg eq grp(i),nh)			;collect group
  pg=pp(oo) & wg=ww(oo) & hg=hh(oo) & dg=dp(oo)	;collect p,w,h,dp for group
  cp='' & cw='' & ch=''
  for j=0,nh-1 do begin
    if j gt 0 then pg(j)=dg(j)			;reset p's from dp's
    ;			output
    cp=cp+' '+strtrim(pg(j),2)
    cw=cw+' '+strtrim(wg(j),2)
    ch=ch+' '+strtrim(hg(j),2)
  endfor
  printf,upf,strtrim(cp,2) & printf,uwf,strtrim(cw,2) & printf,uhf,strtrim(ch,2)
endfor

;close output files
close,upf & free_lun,upf
close,uwf & free_lun,uwf
close,uhf,/all & free_lun,uhf

return
end

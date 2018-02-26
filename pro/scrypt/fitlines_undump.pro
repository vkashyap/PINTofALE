;+
;script	fitlines_undump
;	restores the variables saved to disk from FITLINES_EVENT
;	and calls FITLINES again.
;
;warning
;	this will restore all saved variables from within FITLINES_EVENT.
;	all of your existing variables with the same names will be
;	overwritten.  best to always start from a clean environment.
;
;	save files written in IDL 5.4 are incompatible with earlier versions,
;	but can be restored under certain circumstances using the CMSVLIB
;	package of Craig Markwardt
;	(http://cow.physics.wisc.edu/~craigm/idl/idl.html)
;	This program has a "restore" command right at the start that
;	may fail, but if the variables have already been loaded in,
;	just type .SKIP and .CON, and bob's your uncle.
;
;usage
;	fitsavfil='fitlines.save'
;	.run fitlines_undump
;
;input
;	FITSAVFIL	name of IDL save file containing the variables
;			dumped from within the FITLINES GUI
;
;output
;	FITSTR	structure containing the fit parameters
;
;history
;	vinay kashyap (MMJul)
;	cosmetic surgery (VK; DecMM)
;	included keyword HISTERIX (VK; FebMMIV)
;-

if not keyword_set(fitsavfil) then begin
  help,fitsavfil
  message,'fitsavfil: missing IDL savfile of FITLINES dump'
endif else restore,fitsavfil,verbose=verbose

c1=''
if MOO lt NX then begin
  read,prompt='restore only selected region [y] or full spectrum [n]? ',c1
  cc=strlowcase(strtrim(c1,2))
  if cc ne '' and cc ne 'y' then begin
    MOO=NX & OO=lindgen(MOO)
  endif
endif

read,prompt='zero out all curvature errors [y/n]? ',c1
cc=strlowcase(strtrim(c1,2))
if cc eq '' or cc eq 'y' then begin
  hp=where(perr eq perrm and perrc eq 0,mhp)
  hf=where(ferr eq ferrm and ferrc eq 0,mhf)
  hw=where(werr eq werrm and werrc eq 0,mhw)
  if mhp gt 0 then begin & perr[hp]=-1 & perrm[hp]=-1 & endif
  if mhf gt 0 then begin & ferr[hf]=-1 & ferrm[hf]=-1 & endif
  if mhw gt 0 then begin & werr[hw]=-1 & werrm[hw]=-1 & endif
endif

read,prompt='zero out all errors [n/y]? ',c1
cc=strlowcase(strtrim(c1,2))
if cc eq 'y' then begin
  perr[*]=-1 & perrm[*]=-1 & perrc[*]=0.0
  ferr[*]=-1 & ferrm[*]=-1 & ferrc[*]=0.0
  werr[*]=-1 & werrm[*]=-1 & werrc[*]=0.0
endif

if n_tags(histerix) eq 0 then histerix=0

fitstr=fitlines(x[oo],y[oo],ysig=sigy[oo],$
	funcs=funcs,intens=intens,dchi=delchi,$
	pos=pos,wdt=wdt,flx=flx,$
	perrp=perr,werrp=werr,ferrp=ferr,$
	perrm=perrm,werrm=werrm,ferrm=ferrm,$
	perrc=perrc,werrc=werrc,ferrc=ferrc,$
	thaw=thaw,type=type,ties=ties,epithet=epithet,$
	conlev=conlev[oo],consig=consig[oo],$
	comment=comment,histerix=histerix, _extra=e)

message,'The output is in the structure FITSTR',/info,/noname

end

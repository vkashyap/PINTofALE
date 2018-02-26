function squishem,idstr,abund=abund,mproxy=mproxy, _extra=e
;+
;function	squishem
;	returns a new ID structure from the old one with all the
;	multiple IDs (if any) squished into a single ID pointing
;	to the one with the highest peak emissivity.  If the IDs
;	include different elements, their emissivities are combined
;	in the ratio of their abundances.
;
;	a few more words about this combining:
;	o the >fluxes< are combined without regard to abundance, i.e.,
;	  are added up directly.  this is because the fluxes are assumed
;	  to already include all info about DEMs, abundances, ion balance,
;	  and so on. besides, UPDATID splits a measured flux value simply
;	  on the basis of the ion-balanced emissivities.
;	o the >emissivities< are combined by normalizing the weaker IDs to
;	  that of the strongest ID in the ratio of the abundances.  because
;	  the main ID itself only retains the strongest ID, the absolute
;	  abundance on the emissivity is not set, but relative abundances
;	  are fixed.
;
;syntax
;	sqID=squishem(idstr,abund=abund)
;
;parameters
;	idstr	[INPUT; required] an ID structure, e.g., output of LINEID
;
;keywords
;	abund	[INPUT] abundances relative to H (abund(0)=1)
;		* abund(Z-1) contains the abundance for element Z
;		* if array size is smaller than the largest present Z,
;		  the last element is replicated to fill the gap
;		* default: Anders & Grevesse (1989)
;       mproxy  [OUTPUT] if set, will output 2 dimensional MPROXY string 
;               as prescribed by MIXIE()
;	_extra	[JUNK] here only to avoid crashing the program
;
;subroutines
;	GETABUND
;	ZION2SYMB
;       INICON
;history
;	vinay kashyap (FebMM)
;	bug correction with flux err (VK; MarMM)
;	unknown IDs were crashing at abundances; added call to ZION2SYMB
;	  (VK; MMJul)
;       ADDED MPROXY keyword (Feb2004)
;-

;	usage
ok='ok' & np=n_params() & nid=n_tags(idstr)
if np eq 0 then ok='Insufficient parameters' else $
 if nid eq 0 then ok='IDs should be in a structure' else $
  if nid eq 1 then ok='input format incompatible' else begin
    nw=n_elements(idstr.(0))
    if nid ne nw+1L then ok='Input in unknown format'
  endelse
if ok ne 'ok' then begin
  print,'Usage: squished_ID=squishem(idstr,abund=abund)'
  print,"  convert multiply ID'd lines to simple single IDs"
  if np ne 0 then message,ok,/info
  if np eq 0 then return,0L
  return,idstr
endif

;	abundances
nabu=n_elements(abund)
defabu=getabund('anders & grevesse') & abu=defabu
if nabu eq 0 then abu=defabu else begin
  abu(*)=abund(nabu-1L)
  if nabu lt n_elements(defabu) then abu(0L:nabu-1L)=abund else $
   abu(*)=abund(0L:n_elements(defabu)-1L)
endelse

;	start copying the old ID structure into new, squished, one
tnames=tag_names(idstr)
sqid=create_struct(tnames(0),idstr.(0))
tmproxy = strarr(2,2)
inicon, atom=atom,roman=roman
for id=1L,nid-1L do begin			;{step through IDSTR
  tmp=idstr.(id) & mtmp=n_tags(tmp)
  ww=tmp.WVL & zz=tmp.Z & ii=tmp.ION & ll=tmp.LABL
  ff=tmp.FLUX & ffee=tmp.FLUXERR & tt=tmp.LOGT & ee=tmp.EMIS
  if mtmp gt 8 then ss=tmp.NOTES else ss=''
  ;
  zab=abu([zz-1]) & zion2symb,zz,ii,zymb,ziform='Z0'
  mw=n_elements(ww) & ml=n_elements(ll) & emx=dblarr(mw)
  for i=0L,mw-1L do emx(i)=max(reform(ee(*,i)))
  eemx=max(zab*emx,jmx)
  w0=ww(jmx) & z0=zz(jmx) & i0=ii(jmx) & f0=0*ff(jmx) & fe0=0*ffee(jmx)
  e0=0.D*tt
  if ml eq 1 then l0=ll
  if ml eq 2 then l0=ll
  if ml gt 2 then l0=ll(*,jmx)
  zab=abu([zz-1]) & if z0 gt 0 then zab=zab/abu(z0-1)
  ; 
;stop
  for i=0L,mw-1L do begin	;{add up multiple components
      if arg_present(mproxy) then begin;and i ne jmx then begin  
          bse = strcompress(string(id-1),/remove_all) 
          prxee = atom[tmp.z(i)-1]+roman[tmp.ion(i)-1]+' '+$
                  strjoin(tmp.labl(*,i),' ')
          tmproxy = [tmproxy,transpose([bse,prxee])]
      endif 
    emis=reform(ee(*,i))
    e0=e0+zab(i)*emis
    ;f0=f0+zab(i)*ff(i)
    ;fe0=fe0+(zab(i)*ffee(i))^2
    f0=f0+ff(i)
    fe0=fe0+ffee(i)^2
    ;if i eq 0 then s0=strtrim(zab(i),2)+'*'+strtrim(ff(i),2)
    ;if i gt 0 then s0=s0+' +'+strtrim(zab(i),2)+'*'+strtrim(ff(i),2)
    if i eq 0 then s0=strtrim(zab(i),2)+'* E('+zymb(i)+')'
    if i gt 0 then s0=s0+' +'+strtrim(zab(i),2)+'* E('+zymb(i)+')'
  endfor			;I=0,MW-1}
  fe0=sqrt(fe0)
  s0='EMIS = '+s0
  ;
  if arg_present(mproxy) then begin 
     if n_elements(tmproxy) ne 4 then mproxy = transpose(tmproxy[2:*,*]) else $
        mproxy='none'
  endif
  if mw gt 1 then ss=ss+string("12b)+s0 else ss=s0
  tstr=create_struct('WVL',[w0],'Z',[z0],'ION',[i0],$
	'LABL',l0,'FLUX',[f0],'FLUXERR',[fe0],$
	'LOGT',tt,'EMIS',e0,'NOTES',ss)
  sqid=create_struct(sqid,tnames(id),tstr)
endfor						;ID=1,NID}

return,sqid
end

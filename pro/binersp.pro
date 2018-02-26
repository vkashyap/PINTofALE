function binersp,nrg,binmn,binmx,bstr=bstr,nbin=nbin, _extra=e
;+
;function	binErsp
;	returns position indices of input energ(y/ies) appropriately
;	binned into the specified binning scheme
;
;	returns -1 where E is outside the bin range
;
;syntax
;	iE=binersp(nrg,binmn,maxbin,bstr=bstr,nbin=nbin)
;
;parameters
;	nrg	[INPUT; required] photon energies for which to find
;		the binning numbers
;	binmn	[INPUT] array of bin minimum values
;		* overrides all keywords
;		* if not given, see description below on how it is
;		  derived from BSTR.ELO, and NBIN & NRG.
;	binmx	[INPUT] maximum of the bin max values.
;		* if not a scalar, then the last element of array is used
;		* see description below for how it is determined if
;		  not given.
;
;keywords
;	bstr	[INPUT] structure containing the following pieces of
;		information, as obtained from an OGIP-compliant response
;		matrix (see RDRESP.PRO/RDARF.PRO): {NNRG, ELO, EHI}
;	nbin	[INPUT; default=101] number of equally spaced bins into
;		which to split the input array of photon energies
;		* if -ve, binning will be logarithmic!
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	if BINMN is not given, then BINMN=BSTR.ELO
;	  if BSTR is also not given (or is incorrect), then
;	  	BINMN=FINDGEN(ABS(NBIN))*DE+EMIN
;	  where DE=(EMAX-EMIN)/(ABS(NBIN)-1), EMIN=MIN(EE,MAX=EMAX),
;	  and EE=NRG if NBIN>0, EE=ALOG10(NRG) if NBIN<0
;	if BINMX < MAX(BINMN), BINMN is truncated
;	if BINMX < MIN(BINMN), BINMX is reset to MIN(BINMN)
;	if BINMX is not given, then BINMX=MAX(BINMN)+BINMX(LAST)-BINMX(LAST-1)
;	  if BINMN is not given, BINMX=MAX(BSTR.EHI)
;	    if BSTR is also not given (or is incorrect), then BINMX=EMAX
;
;history
;	vinay kashyap (Jul97)
;	converted to IDL5 notation (VK; OctMM)
;-

;	usage
nph=n_elements(nrg)
if nph eq 0 then begin
  print,'Usage: iE=binersp(E,binmn,maxbin,bstr=bstr,nbin=nbin)'
  print,'  returns bin numbers of given photon energies'
  return,-1L
endif

;	grab inputs
nelo=n_elements(binmn) & nehi=n_elements(binmx) & nbstr=n_tags(bstr)
ee=[nrg[*]]
if nelo gt 0 then elo=[binmn[*]]
if nehi gt 0 then ehi=binmx
if nbstr gt 0 then tbstr=strtrim(tag_names(bstr),2) else tbstr=['NONE']
if not keyword_set(nbin) then nbin=101L & nbn=abs(nbin)

;	define checks for bad inputs
if nelo gt 0 then begin					;(BINMN given
  emin=min(elo,max=emax)
  if nelo gt 1 then emax=emax+(elo[nelo-1]-elo[nelo-2])
  if nehi gt 0 then emax=ehi[nehi-1]
endif else begin		;)(BINMN not given
  ol=where(tbstr eq 'ELO',mol) & oh=where(tbstr eq 'EHI',moh)
  if (mol eq 0 or moh eq 0) then begin	;(BSTR not given
    emin=min(ee,max=emax)
    if nbin lt 0 then begin
      if emin le 0 then ee=ee-emin+1e-6*abs(emax)	;shift
      ee=alog10(abs(ee))
      emin=min(ee,max=emax)
    endif
    dE=emax-emin
    if nbn gt 1 then begin
      dE=dE/(nbn-1) & elo=findgen(nbn)*dE+emin
    endif else return,long(0*ee)
  endif else begin			;)(BSTR given
    elo=bstr.ELO & emin=min(elo) & emax=max(bstr.EHI)
  endelse				;BSTR)
endelse				;BINMN)

;	problems with min v/s max
if emax lt max(elo) then begin
  oo=where(elo lt emax,moo)
  if moo eq 0 then emax=emin else $	;BINMX < MIN(BINMN)
    elo=elo[oo]				;cut off binmin array
endif

;	quickie out
if (emin eq emax) then begin
  oo=where(ee ne emin,moo) & ie=lonarr(nph)
  if moo gt 0 then ie[oo]=ie[oo]-1
  return,ie
endif

;	and now for the actual binning!
;nelo=n_elements(elo) & jlo=lindgen(nelo) & jee=lonarr(nph)-1
;allE=[elo,ee] & allI=[jlo,jee] & oee=sort(allE)	;combine arrays & sort
;allE=allE(oee) & allI=allI(oee) & oee=where(allI lt 0)
;imn=(where(allI ge 0))[0]	;(first non-zero element)
;for i=imn+1,nelo+nph-1 do $
;	if allI[i] lt 0 then allI[i]=allI[i-1]
;ie=allI[oee]

nelo=n_elements(elo) & allE=[elo,ee]
ilo=lindgen(nelo) & iee=lindgen(nph)+nelo & allI=[ilo,iee]
jlo=lindgen(nelo) & jee=lonarr(nph)-1 & allJ=[jlo,jee]
oee=sort(allE)
allE=allE(oee) & allI=allI(oee) & allJ=allJ(oee)
imn=(where(allJ ge 0))[0] & oj=where(allJ lt 0)
for i=imn+1,nelo+nph-1 do $
	if allJ[i] lt 0 then allJ[i]=allJ[i-1]
ie=allJ[oj] & oee=allI[oj] & ie=ie[sort(oee)]

return,ie
end

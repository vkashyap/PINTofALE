;+
;EG_IDSTR.PRO
;	an example program that extracts information out of an ID structure
;
;	requires an input ID structure in the variable IN_IDSTR
;	extracts all the ID'd wavelengths, atomic numbers, ionic states,
;	and places them in arrays named OUT_<var>
;	look in OUT_IDX for a pointer to the original component
;
;	TBD: need to be cleverer about labels and notes.
;
;subroutines
;	ZION2SYMB
;
;vinay kashyap (AugMM)
;-

;	first check that the input exists
OK='ok' & NID=n_tags(IN_IDSTR)
if not keyword_set(IN_IDSTR) then message,'IN_IDSTR: missing input ID structure'
if nid eq 0 then message,'IN_IDSTR: not a structure'
IDNAM=tag_names(IN_IDSTR)
if nid eq 1 then ok='IN_IDSTR: in unknown format' else $
 if idnam[0] ne 'WVL' then ok='IN_IDSTR: missing wavelengths' else $
  if idnam[1] eq 'WVL_COMMENT' then ok='IN_IDSTR: contains no data'
if ok ne 'ok' then message,ok

;	warning
message,'WARNING: overwrites the following intermediate variables',/info
help,idnam,idtmp,k,mw,ncomp,nid,nlogT,nwvl,slabl,ok

;	initialize
OUT_OBSWVL=IN_IDSTR.WVL		;wavelengths of the observed features
NCOMP=n_elements(out_obswvl)	;number of features ID'd
NWVL=0L				;total number of wavelengths
for ICOMP=0L,ncomp-1L do nwvl=nwvl+n_elements(IN_IDSTR.(icomp+1L).WVL)
OUT_IDX=lonarr(nwvl)		;index pointing ID wavelength to ID component
OUT_WVL=fltarr(nwvl)		;ID wavelengths
OUT_Z=intarr(nwvl)		;atomic numbers of ID
OUT_ION=out_z			;ionic states of ID
OUT_LABL=strarr(2,nwvl)		;level designation/e-configuration of ID
				;*** NEED IMPROVEMENT IN HANDLING ***
OUT_ELEM=strarr(nwvl)		;atomic symbol for ID
OUT_NOTES=out_elem		;notes, if any *** NOT IMPLEMENTED ***
OUT_FLUX=fltarr(nwvl)		;fluxes (see UPDATID/SQUISHEM) of ID
OUT_FLUXERR=out_flux		;errors on OUT_FLUX
OUT_LOGT=IN_IDSTR.(1).LOGT & nlogT=n_elements(out_logt)	;temperature grid
				;OUT_LOGT *must* be identical for all
OUT_EMIS=dblarr(nlogT,nwvl)	;emissivities, including ion balance

;	now step through the structure and extract the contents into arrays
K=0L
for icomp=0L,ncomp-1L do begin
  IDTMP=IN_IDSTR.(icomp+1L)
  MW=n_elements(idtmp.WVL) & SLABL=size(idtmp.LABL)
  OUT_IDX[k:k+mw-1L]=icomp+1L
  OUT_WVL[k:k+mw-1L]=idtmp.WVL
  OUT_Z[k:k+mw-1L]=idtmp.Z
  OUT_ION[k:k+mw-1L]=idtmp.ION
  OUT_FLUX[k:k+mw-1L]=idtmp.FLUX
  OUT_FLUXERR[k:k+mw-1L]=idtmp.FLUXERR
  OUT_EMIS[*,k:k+mw-1L]=(idtmp.EMIS)
  if slabl[0] eq 1 then begin
    if slabl[1] eq 1 then OUT_LABL[*,k:k+mw-1L]=(idtmp.LABL)[0]
    if slabl[1] eq mw then OUT_LABL[0,k:k+mw-1L]=idtmp.LABL
  endif
  if slabl[0] eq 2 then begin
    if slabl[2] eq mw then begin
      OUT_LABL[0,k:k+mw-1L]=(idtmp.LABL)[0,*]
      OUT_LABL[1,k:k+mw-1L]=(idtmp.LABL)[1,*]
    endif
  endif
  k=k+mw
endfor
zion2symb,OUT_Z,OUT_ION,OUT_ELEM,ziform='Z ION'

message,'OUTPUTS are in the variables:',/info
help,OUT_OBSWVL,OUT_IDX,OUT_WVL,OUT_Z,OUT_ION,OUT_ELEM,OUT_LABL,$
	OUT_FLUX,OUT_FLUXERR,OUT_LOGT,OUT_EMIS

end

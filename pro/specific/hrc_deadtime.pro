pro hrc_deadtime,dtffil,dtf,edtf,gate1=gate1,gate2=gate2,$
	dtfus=dtfus,edtfus=edtfus,dtfs=dtfs,edtfs=edtfs,dtf1=dtf1,$
	forceI=forceI,verbose=verbose, _extra=e
;+
;procedure	hrc_deadtime
;	compute deadtime correction given the total, valid, and telemetered counts
;
;syntax
;	hrc_deadtime,dtffil,dtf,edtf,gate1=gate1,gate2=gate2,$
;	dtfus=dtfus,edtfus=edtfus,dtfs=dtfs,edtfs=edtfs,dtf1=dtf1,$
;	verbose=verbose
;
;parameters
;	dtffil	[INPUT; required] DTF file, must contain the columns
;		PROC_EVT_COUNT, TOTAL_EVT_COUNT, and VALID_EVT_COUNT
;	dtf	[OUTPUT; required] computed DTF
;	edtf	[OUTPUT; required] error on DTF
;
;keywords
;	gate1	[INPUT; default=19.5 musec] time for hardware checks
;	gate2	[INPUT; default=49.0 musec] time to complete event
;		information processing
;	dtfus	[OUTPUT] DTF for unsaturated case
;	edtfus	[OUTPUT] error on DTFUS
;	dtfs	[OUTPUT] DTF for saturated case
;	edtfs	[OUTPUT] error on DTFS
;	dtf1	[OUTPUT] structure of contents of DTFFIL
;	forceI	[INPUT] if set, ignores DETNAM and assumes that
;		only the HRC-I specific logic for determining
;		which bins are saturated and which are not apply
;	verbose	[INPUT] controls chatter
;		- stop if >1000
;		- plot if >100
;	_extra	[JUNK] here only to prevent crashing the program
;
;usage
;	.run hrc_deadtime
;		reads from dtf1 file in /data/hrc/ hierarchy, computes
;		new DTF and eDTF, and plots it along with eDTFUS (green)
;		eDTFS (red), and original (yellow) (requires PEASECOLR)
;
;history
;	vinay kashyap (2010nov)
;	now computes DTF by adding the deadtimes for saturated and
;	  unsaturated cases; included usage example (VK; 2010dec)
;-

np=n_params() & ok='ok' & nf=n_elements(dtffil) & szf=size(dtffil,/type)
if np lt 3 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='DTFFIL: filename not specified' else $
  if nf gt 1 then ok='DTFFIL must be a scalar' else $
   if szf[0] ne 7 then ok='DTFFIL must be a character string'
if ok ne 'ok' then begin
  print,'Usage: hrc_deadtime,dtffil,dtf,edtf,gate1=gate1,gate2=gate2,$'
  print,'	dtfus=dtfus,edtfus=edtfus,dtfs=dtfs,edtfs=edtfs,dtf1=dtf1,$'
  print,'       /forceI,verbose=verbose'
  print,'  compute deadtime correction factor and error on it'
  if np ne 0 then message,ok,/informational
  return
endif

vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	read in the DTFFIL
if vv gt 0 then message,'reading DTF from '+dtffil[0],/informational
dtf1=mrdfits(dtffil[0],1,hdtf) & tnam=tag_names(dtf1)
ip=where(strpos(tnam,'PROC_EVT_COUNT') ge 0)
it=where(strpos(tnam,'TOTAL_EVT_COUNT') ge 0)
iv=where(strpos(tnam,'VALID_EVT_COUNT') ge 0)
if ip lt 0 or it lt 0 or iv lt 0 then begin
  message,dtffil[0]+': not a standard DTF file; does not contain proper columns',/informational
  return
endif
time=dtf1.TIME & delT=sxpar(hdtf,'TIMEDEL') & detnam=strtrim(sxpar(hdtf,'DETNAM'),2)
tot=1.0*dtf1.TOTAL_EVT_COUNT & pp=1.0*dtf1.PROC_EVT_COUNT & vld=1.0*dtf1.VALID_EVT_COUNT

;	compute dtf
;	unsaturated case
;	dtf = g1*(TOTAL*PROC)/VALID + g2*PROC
g1=19.5e-6 & if keyword_set(gate1) then g1=gate1[0]*1e-6
g2=49e-6 & if keyword_set(gate2) then g2=gate2[0]*1e-6
gg=g1+g2
dtfus = 0.*pp+1.0
o0=where(vld gt 0,mo0)
if mo0 gt 0 then dtfus[o0] = g1 * ((tot[o0]*pp[o0])/vld[o0])
dtfus= dtfus + g2 * pp
dtfus = 1.-dtfus/delT

;	and if it is saturated..
;	dtf ~ PROC/VALID
dtfs = 0.*pp
if mo0 gt 0 then dtfs[o0] = 1.-(vld[o0]-pp[o0])/vld[o0]

;	acc. to FAP, in hrc_dtfstats, DTF is computed
;	based on unsaturated deadtime below VALID_EVT_COUNT=368
;	and on the sum of the deadtimes above that
;old code;old code;;osat=where(dtf1.DTF lt 0.9,mosat)
;old code;old code;;if mosat gt 0 then dtf[osat]=dtfs[osat]
;old code;dtf=1.-(2.-(dtfus+dtfs))
if keyword_set(forceI) then detnam='HRC-I'	;bypass the HRC-S "kludge"
if detnam eq 'HRC-S' then begin
  ;	by trial and error, it seems that pipeline also
  ;	includes an extra threshold of PROC_EVT_COUNT>366
  ;	for HRC-S.  this is the only way to get the computed
  ;	DTFs to match up with the pipeline DTFs. -VK
  os=where(vld gt 368 and pp gt 366,mos,complement=ous,ncomplement=mous)
endif else begin
  os=where(vld gt 368,mos,complement=ous,ncomplement=mous)
endelse
if vv gt 0 then print,mos,mous
dtf=0.*dtfus+1.
if mous gt 0 then dtf[ous]=dtfus[ous]
if mos gt 0 then dtf[os]=1.-(2.-(dtfus[os]+dtfs[os]))

;	compute dtf error
nn=vld-pp & mm=tot-vld
x11=gg*(pp+nn)^2 & x12=g1*nn*mm & x1=(x11+x12)^2*pp
x2=(g1*pp*mm)^2*nn
x3=(g1*pp)^2*(pp+nn)^2*mm
edtfus = 0.*dtfus+1.0
if mo0 gt 0 then edtfus[o0]=sqrt(x1[o0]+x2[o0]+x3[o0])/(pp[o0]+nn[o0])^2
edtfus=edtfus/delT

;	and errors when saturated
edtfs = 0.*dtfs+1.0
o0=where(nn gt 0 and pp gt 0,mo0)
if mo0 gt 0 then edtfs[o0]=((nn[o0]*pp[o0])/(nn[o0]+pp[o0])^2)*sqrt(1./nn[o0]+1./pp[o0])

;	and combine the errors
edtf=sqrt(edtfus^2+edtfs^2)
if mous gt 0 then edtf[ous]=edtfus[ous]
if mos gt 0 then edtf[os] = sqrt(edtfus[os]^2+edtfs[os]^2)

if vv gt 100 then begin
  plot,dtf1.DTF,/ynoz,_extra=e
  oplot,dtf,color=100
  oplot,dtf+edtf,col=150,psym=3
  oplot,dtf-edtf,col=150,psym=3
endif
if vv gt 1000 then begin
  if vv gt 10000 then begin
    oy=where(abs(dtf-dtf1.DTF) gt 0.01,moy) & help,moy
    plot,vld,pp,psym=3 & if moy gt 0 then oplot,vld[oy],pp[oy],psym=3,col=2
  endif
  stop,'HALTing; type .CON to continue'
endif

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	this is a usage example
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

peasecolr & loadct,3 & peasecolr	;requires PEASECOLR

if not keyword_set(obsid) then obsid=3698
sobsid=string(obsid[0],'(i5.5)')
dtffil=findfile('/data/hrc/?/'+sobsid+'/primary/*dtf1*',count=nfil)
if nfil eq 0 then message,'DTF file not found for ObsID: '+sobsid

;	call hrc_deadtime
hrc_deadtime,dtffil[0],dtf,edtf,gate1=gate1,gate2=gate2,$
	dtfus=dtfus,edtfus=edtfus,dtfs=dtfs,edtfs=edtfs,dtf1=dtf1,$
	forceI=forceI,verbose=verbose

;	make plots comparing computed DTFERR with saturated, unsaturated, and pipeline
if n_elements(xrange) ne 2 then xrange=[0,1]
if n_elements(yrange) ne 2 then yrange=[0,0.03]
plot,dtf,edtf,psym=3,xr=xrange,yr=yrange,xtitle='DTF',ytitle='DTF_ERR',title=dtffil[0]
oplot,dtf,edtfs,psym=3,col=2
oplot,dtf,edtfus,psym=3,col=3
oplot,dtf,dtf1.DTF_ERR,psym=3,col=4
print,'averaged DTF (old vs new)',mean(dtf1.DTF,/nan),mean(dtf,/nan)
print,'error weighted averaged DTF (old vs new)',$
	total(dtf1.DTF/(dtf1.DTF_ERR)^2,/nan)/total(1./(dtf1.DTF_ERR)^2,/nan),$
	total(dtf/edtf^2,/nan)/total(1./edtf^2,/nan)

end

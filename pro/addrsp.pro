function addrsp,rm1,rm2,ea1=ea1,ea2=ea2,eps=eps,morecol=morecol,$
	verbose=verbose,ea3=ea3,ea4=ea4,ea5=ea5,ea6=ea6, _extra=e
;+
;function	addrsp
;	combine two response matrices to make a composite RSP and
;	return a structure in the same format as RD_OGIP_RMF()
;
;	steps through each row of an RMF and tacks on the array
;	from a 2nd RMF, updating the N_GRP, N_CHAN, F_CHAN elements
;	as appropriate.
;
;syntax
;	rsp=addrsp(rm1,rm2,ea1=ea1,ea2=ea2,eps=eps,verbose=verbose)
;
;parameters
;	rm1	[INPUT; required] response matrix structure as
;		read in by RD_OGIP_RMF
;	rm2	[INPUT; required] response matrix structure to
;		be added to RM1
;		* RM1 and RM2 are expected to have the following fields:
;	        NNRG,ELO,EHI,NCHAN,EMN,EMX,N_GRP,F_CHAN,N_CHAN,MATRIX,FIRSTCHAN
;
;keywords
;	ea1	[INPUT] if set to a scalar or size matches RM1.ELO,
;		then multiplies RM1 by EA1 prior to addition
;	ea2	[INPUT] as EA1, but for RM2
;		* default for EA1 and EA2 is 1.0
;	eps	[INPUT] a small number, below which to set everything
;		to zero
;		* default is 1e-6
;	morecol	[INPUT] designed as a buffer so that a matrix won't
;		overflow the bounds automatically set.  by default,
;		the new matrix is defined at an early stage to be
;		max(RM1.N_CHAN)+max(RM2.N_CHAN)+100+MORECOL.
;		* note that MORECOL can even be negative.
;		* use wisely.
;	verbose	[INPUT] controls chatter
;	ea3	[JUNK] here just to prevent user error
;	ea4	[JUNK] here just to prevent user error
;	ea5	[JUNK] here just to prevent user error
;	ea6	[JUNK] here just to prevent user error
;		* if EA3, EA4, EA5, or EA6 are set, the program will
;		  complain and quit.  this is done to avoid the obvious
;		  user error of trying to specify EA3=EA3, etc. if
;		  adding a third RMF/ARF to a preexisting one.
;
;	_extra	[JUNK] here only to avoid crashing the program
;
;subroutines
;	KILROY
;	OGIPZIP
;
;history
;	vinay kashyap (Nov2001)
;	bug correction re max possible matrix elements when NGROUP>>1;
;	  added dummy keywords EA3, EA4, EA5, EA6 and extra checks
;	  against tricky inputs (VK; Aug02)
;	bug correction: wasn't handling N_GRP=1 (VK; Aug04)
;-

;	usage
ok='ok' & np=n_params()
n1=n_elements(rm1) & n2=n_elements(rm2) & nr1=n_tags(rm1) & nr2=n_tags(rm2)
if np lt 2 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='RM1 is undefined' else $
  if n2 eq 0 then ok='RM2 is undefined' else $
   if nr1 eq 0 then ok='RM1 is not a structure' else $
    if nr2 eq 0 then ok='RM2 is not a structure'
if ok ne 'ok' then begin
  print,'Usage: rsp=addrsp(rm1,rm2,ea1=ea1,ea2=ea2,eps=eps,verbose=verbose)
  print,'  combine two response matrices to make a composite RSP'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	verify RM inputs
ok='ok' & t1=tag_names(rm1) & t2=tag_names(rm2) & k1=0 & k2=0
for i=0,nr1-1 do begin
  if t1[i] eq 'ELO' then k1=k1+1
  if t1[i] eq 'EHI' then k1=k1+1
  if t1[i] eq 'EMN' then k1=k1+1
  if t1[i] eq 'EMX' then k1=k1+1
  if t1[i] eq 'N_GRP' then k1=k1+1
  if t1[i] eq 'F_CHAN' then k1=k1+1
  if t1[i] eq 'N_CHAN' then k1=k1+1
  if t1[i] eq 'MATRIX' then k1=k1+1
  if t1[i] eq 'FIRSTCHAN' then k1=k1+1
endfor
for i=0,nr2-1 do begin
  if t2[i] eq 'ELO' then k2=k2+1
  if t2[i] eq 'EHI' then k2=k2+1
  if t2[i] eq 'EMN' then k2=k2+1
  if t2[i] eq 'EMX' then k2=k2+1
  if t2[i] eq 'N_GRP' then k2=k2+1
  if t2[i] eq 'F_CHAN' then k2=k2+1
  if t2[i] eq 'N_CHAN' then k2=k2+1
  if t2[i] eq 'MATRIX' then k2=k2+1
  if t2[i] eq 'FIRSTCHAN' then k2=k2+1
endfor
if k1 lt 9 then ok='RM1 is not a standard response matrix structure' else $
 if k2 lt 9 then ok='RM2 is not a standard response matrix structure'
if ok ne 'ok' then begin
  message,ok,/info
  return,-1L
endif

;	check keywords
thr=1e-6 & if keyword_set(eps) then thr=eps[0]+0.0
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
if not keyword_set(morecol) then morecol=0L
elo1=rm1.ELO & elo2=rm2.ELO & ehi1=rm1.EHI & ehi2=rm2.EHI
na1=n_elements(elo1) & na2=n_elements(elo2)
ar1=fltarr(na1)+1. & ar2=fltarr(na2)+1.
ma1=n_elements(ea1) & ma2=n_elements(ea2)
if ma1 eq 0 then message,'EA1 is undefined, assuming EA1=1',/info
if ma2 eq 0 then message,'EA2 is undefined, assuming EA2=1',/info
if ma1 eq 1 then begin
  message,'EA1 is a scalar, assuming constant',/info
  ar1=ar1*ea1[0] ;& ma1=na1
endif
if ma2 eq 1 then begin
  message,'EA2 is a scalar, assuming constant',/info
  ar2=ar2*ea2[0] ;& ma2=na2
endif
if na1 eq ma1 then ar1=ea1*1.0	;(multiply by 1.0 to make sure ..
if na2 eq ma2 then ar2=ea2*1.0	;.. that it is at least float)
firstchan1=rm1.FIRSTCHAN & firstchan2=rm2.FIRSTCHAN
ok='ok'
if keyword_set(ea3) then ok='EA3: unrecognized keyword' else $
 if keyword_set(ea4) then ok='EA4: unrecognized keyword' else $
  if keyword_set(ea5) then ok='EA5: unrecognized keyword' else $
   if keyword_set(ea6) then ok='EA6: unrecognized keyword'
if ok ne 'ok' then begin
  message,ok,/info
  message,'Exiting without doing anything; please check documentation',/info
  return,-1L
endif

;	make arrays of appropriate sizes to hold everything
;	n_grp contains the number of different groups
;	f_chan contains the index of the first channel of each group
;	n_chan contains the number of channels in each group
n_grp1=rm1.N_GRP & n_grp2=rm2.N_GRP
f_chan1=rm1.F_CHAN & f_chan2=rm2.F_CHAN
n_chan1=rm1.N_CHAN & n_chan2=rm2.N_CHAN
matrix1=rm1.MATRIX & matrix2=rm2.MATRIX
;
;	complication: max number of groups (and hence final size of
;	the matrix) may be smaller than sum if there is overlap between
;	f_chan1:f_chan1+n_chan1 and f_chan2:f_chan2+n_chan2
;
nnrg=na1 & nchan=n_elements(rm1.EMN)
n_grp=n_grp1+n_grp2 & n_grpmax=max(n_grp) & maxn_grp=0L
n_grpmax1=max(n_grp1) & n_grpmax2=max(n_grp2)
if n_grpmax1 gt 1 then begin
  sumnchan=0.*n_grp
  for i=0L,n_elements(sumnchan)-1L do begin
    if n_grp1[i] eq 1 then sumnchan[i]=sumnchan[i]+n_chan1[i]
    if n_grp1[i] gt 1 then sumnchan[i]=sumnchan[i]+total(n_chan1[*,i])
    if n_grp2[i] eq 1 then sumnchan[i]=sumnchan[i]+n_chan2[i]
    if n_grp2[i] gt 1 then sumnchan[i]=sumnchan[i]+total(n_chan2[*,i])
    ;	sumnchan[i]=total(n_chan1[*,i])+total(n_chan2[*,i])
  endfor
  nmatmax=long(max(sumnchan))
endif else nmatmax=max(n_chan1)+max(n_chan2)
maxnmat=0L
nmatmax=nmatmax+100L+morecol	;just for the heck of it
f_chan=intarr(n_grpmax,nnrg) & n_chan=intarr(n_grpmax,nnrg)
matrix=fltarr(nmatmax,nnrg)

;	here's the deal
for i=0L,na1-1L do begin	;{step through RM1 rows
  if vv gt 1 and i eq 100*long(i/100.) then kilroy,dot=strtrim(i,2)+'..'
  e10=elo1[i] & e11=ehi1[i] & ng1=n_grp1[i]

  ;	response for this row, for RM1
  rsp1=fltarr(nchan) & jfc=0L
  for ig=0,ng1-1 do begin
    if n_grpmax1 gt 1 then begin
      ifc=f_chan1[ig,i] & inc=n_chan1[ig,i]
    endif else begin
      ifc=f_chan1[i] & inc=n_chan1[i]
    endelse
    if keyword_set(firstchan1) then ifc=ifc-firstchan1
    if inc gt 0 then rsp1[ifc:ifc+inc-1L]=matrix1[jfc:jfc+inc-1L,i]
    jfc=jfc+inc
  endfor
  resp=rsp1*ar1[i]*(e11-e10)

  o20=where(ar2 gt 0,mo20) & j0=0L & melo2=min(elo2) & mehi2=max(ehi2)
  if mo20 gt 0 then begin
    j0=o20[0] & melo2=elo2[o20[0]]
  endif
  if melo2 gt e11 or mehi2 lt e10 then goto,skip2
  if vv gt 2 then kilroy,dot=strtrim(i,2)+'..'

  for j=j0,na2-1L do begin	;{step through RM2 rows
    e20=elo2[j] & e21=ehi2[j] & ng2=n_grp2[j]
    if e20 lt e11 and e21 gt e10 and ar2[j] gt 0 then begin	;(ok, this bit is relevant

      ;	response for this row, for RM2
      rsp2=fltarr(nchan) & jfc=0L
      for ig=0,ng2-1 do begin
	if n_grpmax2 gt 1 then begin
	  ifc=f_chan2[ig,j] & inc=n_chan2[ig,j]
	endif else begin
	  ifc=f_chan2[j] & inc=n_chan2[j]
	endelse
	if keyword_set(firstchan2) then ifc=ifc-firstchan2
	if inc gt 0 then rsp2[ifc:ifc+inc-1L]=matrix2[jfc:jfc+inc-1L,j]
	jfc=jfc+inc
      endfor

      ;	add the two responses
      ;	now, in general, the 2nd RMF does not have exactly the same
      ;	energy binning as the 1st.  for e.g., 3 bins of the 2nd may
      ;	correspond to just 1 of the 1st.  in that case, if nothing is
      ;	done, the response due to the 2nd will be 3 times too large
      ;	when added to the 1st.  we must therefore weight them by the
      ;	fractional overlap of the corresponding bins.  that is, if
      ;	only a small fraction of the original bin is covered by the
      ;	bin of the 2nd RMF, then weight by the ratio of the overlap
      ;	to the original.  if the 2nd completely overlaps the 1st,
      ;	that automatically means include all of the 2nd (i.e., max of
      ;	ratio=1)
      f=[e10,e11,e20,e21] & f=f[sort(f)]
      resp=resp+rsp2*ar2[j]*(f[2]-f[1])

      if vv gt 5 then kilroy,dot=strtrim(j,2)+','+strtrim((f[2]-f[1])/(e11-e10),2)+' '

    endif					;E20<E11 and E21>E10)
  endfor			;J=0,NA2-1}
  if vv gt 5 then plot,rm1.EMN,resp/(e11-e10),/yl,/xl,yr=[thr,max(resp)/(e11-e10)],title=i,/xs
  skip2:	;yeah, a goto

  ;	compress the response
  resp=resp/(e11-e10)
  ogipzip,resp,mat,ng,fc,nc,eps=thr,chan0=firstchan1

  ;	I am betting that the number of groups won't exceed estimated maximum
  if ng gt 0 then begin
    n_grp[i]=ng & f_chan[0:ng-1,i]=fc & n_chan[0:ng-1,i]=nc
    matrix[0L:long(total(nc))-1L,i]=mat
  endif else begin
    n_grp[i]=ng & f_chan[*,i]=0+firstchan1 & n_chan[*,i]=0
    matrix[*,i]=0.
  endelse

endfor				;I=0,NA1-1}

;	shrink further if possible
maxn_grp=max(n_grp)
if maxn_grp lt n_grpmax then begin
  f_chan=f_chan[0L:maxn_grp-1L,*]
  n_chan=n_chan[0L:maxn_grp-1L,*]
endif

;	shrink matrix if possible
ok=lonarr(nmatmax)
for i=0L,nmatmax-1L do if total(matrix[i,*]) eq 0 then ok[i]=-1
ook=where(ok ge 0,mook)
if mook eq 0 then message,'nothing in matrix?!'
matrix=matrix[ook,*]

;	reconstruct the structure
rsp=create_struct('NNRG',nnrg,'ELO',elo1,'EHI',ehi1,$
	'NCHAN',nchan,'EMN',rm1.EMN,'EMX',rm1.EMX,$
	'N_GRP',n_grp,'F_CHAN',f_chan,'N_CHAN',n_chan,$
	'MATRIX',matrix,'FIRSTCHAN',firstchan1)

return,rsp
end

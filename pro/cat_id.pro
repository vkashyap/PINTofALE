function cat_id,AA,B,pick=pick,ask=ask,okev=okev,comm=comm,flst=flst,$
	mplet=mplet,dbdir=dbdir,fluxes=fluxes,fluxerr=fluxerr, _extra=e
;+
;function	cat_id
;	returns a concatenated structure of Line IDs
;
;syntax
;	C=cat_id(A,B,pick=pick,ask=ask,/okev,comm=comm,flst=flst,dbdir=dbdir,$
;	fluxes=fluxes,fluxerr=fluxerr)
;
;parameters
;	A	[INPUT; required] output of LINEID
;	B	[INPUT] output of LINEID to be merged with A
;		* if any wavelengths are common, the IDs in A take
;		  precedence (but see keyword ASK)
;		* if not given, then the program *prints* a summary
;		  of the IDs to the screen
;
;keywords
;	pick	[INPUT; default: all] position indices of elements in A
;		to be selected (or not).
;		* if scalar, returns all elements EXCEPT the specified one
;		* has no effect on B
;		* use these two to delete entries, for e.g.
;	ask	[INPUT] if set, asks for confirmation before throwing
;		away an ID from B that also exists in A.  allows overriding
;		the default action.
;		* ASK='n' will reverse the default by throwing away the
;		  ID from A instead of the one from B
;		* ASK='r' *reverses* the default, in that if there is a
;		  component in B that's also in A, the one in A is >replaced<
;		  by the one in B (as opposed to if ASK='n', then the one in
;		  A is deleted, and B is just appended at the end)
;		* ASK='k' overrides the default by keeping BOTH
;	okev	[INPUT] if set, converts wavelengths [Ang] to energy [keV]
;	comm	[INPUT] appends a comment to each line while printing
;		* takes the comments from A.LABL(1,i)+A.LABL(0,i)
;		* truncates the comment to a maximum of 30 characters
;		  unless COMM is set to a higher number
;		* for this to work when FLST is set, make sure DBDIR is
;		  also defined
;	flst	[INPUT] if set, prints out the line list in the form that
;		RD_LIST likes
;		* if FLST is a string, writes to file FLST
;       mplet   [OUTPUT] if set, will contain on output a 2 dimensional MPLET
;               string array as prescribed by MIXIE()
;	dbdir	[INPUT] if set and FLST is set, appends DBDIR to IDs
;		* 1 => $SCAR , 2 => $SPEX , 3 => $CHIANTI
;	fluxes	[INPUT] if set and matches either the number of components
;		or the number of matches, then renormalizes/overwrites the
;		line fluxes in the output.
;		* ONLY if IDs are being printed
;	fluxerr	[INPUT] if set and matches either the number of components
;		or the number of matches, then appends an error to output
;		* ignored if FLST is not set
;		* if FLUXERR is given for each match, square-adds the errors
;		  for the components
;		* if single element, but there are more than 1 component to
;		  ID structure, this is assumed to be the fractional error
;		  (if < 1) or S/N (if > 1) for each match
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	* works ONLY on the outputs of LINEID
;	* requires subroutines:
;	  -- INICON
;	  -- CREATE_STRUCT
;	  -- LINEID
;
;history
;	vinay kashyap (Feb97)
;	bug -- was breaking down in simplest case (VK; Mar97)
;	was breaking down if ID was unknown (VK; Jul97)
;	added keyword FLST (VK; Oct98)
;	added call to INITSTUFF; added Tmax to output; added keywords DBDIR,
;	  FLUXES, FLUXERR; corrected bug re undefined output when FLST was set;
;	  increased number of significant digits in FLST output; generalized
;	  ASK to "Yes/No/Replace/Keep" (VK; Nov98)
;	updated FLUXERR behavior to account for new field (VK; Dec98)
;	modified to allow for field NOTES to be printed out (VK; Mar99)
;	changed call to INITSTUFF to INICON (VK; 99May)
;	added FLUXERR to output (VK; MarMM)
;	converted to IDL5 array notation (VK; OctMM)
;	changed keyword KEV to OKEV (VK; JanMMI)
;	now /COMM also works when FLST (and DBDIR) are defined (VK;Nov01)
;       added MPLET keyword (LL; Feb04) 
;	adjusted print format a bit to avoid **** in FLUXERR (VK; Jul12)
;-

;	Usage
sza=size(AA) & szb=size(B) & nsa=n_elements(sza) & nsb=n_elements(szb)
if sza[nsa-2] ne 8 then begin
  print,'Usage: C=cat_id(A,B,pick=pick,ask=ask,/okev,comm=comm,flst=flst,$'
  print,'       dbdir=dbdir,fluxes=fluxes)'
  print,'  returns duplicate-less concatenation of LINE ID structures A and B'
  print,'--------------------------------------------------------------------'
  print,form="($,a7)",'LINEID ' & return,lineid()
endif
A=AA

;	catch input errors
c1='ok'
;
if szb[0] ne 0 and szb[nsa-2] ne 8 then begin
  message,'Input not a structure',/info & return,A
endif
;
tagA=tag_names(A)
if tagA[0] ne 'WVL' then c1='Structure not in right format'
if tagA[1] eq 'WVL_COMMENT' then c1='Structure contains no data'
if tagA[1] eq 'Z' then c1='Structure has no wrapper'
;
if szb[nsa-2] eq 8 then begin
  tagB=tag_names(B)
  if tagB[0] ne 'WVL' then c1='Structure in incorrect form'
  if tagB[1] eq 'WVL_COMMENT' then c1='Structure only has descriptions'
  if tagB[1] eq 'Z' then c1='Structure has been stripped'
endif
;
if c1 ne 'ok' then begin & message,c1,/info & return,A & endif

;{	a small digression to "reset" A
if n_elements(pick) ne 0 then begin
  szo=size(pick)
  wA=A.wvl & nwA=n_elements(wA) & iA=lindgen(nwA)
  if szo[0] eq 0 then begin		;{scalar -- delete this element!
    oA=where(iA ne pick,moA)
    if moA eq 0 then begin		;nothing to return?
      if n_tags(B) ne 0 then return,B else return,A
    endif else begin
      wA=wA[oA] & tmp=create_struct(tagA[1],A.(oA[0]+1))
      for i=1,moA-1 do tmp=create_struct(tmp,tagA[i+1],A.(oA[i]+1))
      A=create_struct(tagA[0],wA,tmp)
    endelse
  endif else begin			;}{vector -- pick only these!
    oA=pick & moA=n_elements(oA)
    wA=wA[oA] & tmp=create_struct(tagA[1],A.(oA[0]+1))
    for i=1,moA-1 do tmp=create_struct(tmp,tagA[i+1],A.(oA[i]+1))
    A=create_struct(tagA[0],wA,tmp)
  endelse				;}
endif
;}

wA=A.wvl & ww=[wA] & nwA=n_elements(wA)
if szb[0] ne 0 then begin
  wB=B.wvl & ww=[ww,wB] & nwB=n_elements(wB)
endif else nwB=0L
wC=ww(uniq(ww,sort(ww))) & nwC=n_elements(wC)
if nwC gt 126 then c1='Too many fields; cannot concatenate'
ncomm=30 & if keyword_set(comm) then ncomm=fix(comm)>ncomm

;	initialize
k=0L & outw=ww
atom=1 & rom=1 & inicon,atom=atom,roman=rom
atom=['X',atom]		;allow for "UNKNOWN"
rom=['',rom]		;allow for "missing"

;	merge
if nwB gt 0 then begin				;{there be sumpin' to merge
  noB=lonarr(nwB)-1				;keep track of common elements

  for i=0,nwA-1 do begin			;(step through A to fill C
    tmp=A.(i+1) & setask='y'
    oo=where(wA[i] eq wB,moo)			;any in common?
    if moo gt 0 then begin			;{yep
      print,string(wA[i],'(f10.4)')
      xw=tmp.wvl & xz=tmp.z & xi=tmp.ion & nx=n_elements(xw)
      for ix=0,nx-1 do print,'A:'+strtrim(ix+1,2)+' ',$
	string(xw[ix],'(f10.4)')+' '+atom(xz[ix])+rom(xi[ix])
      for j=0,moo-1 do begin
	xw=B.(oo[j]+1).wvl & xz=B.(oo[j]+1).z & xi=B.(oo[j]+1).ion
	nx=n_elements(xw)
        for ix=0,nx-1 do print,'B:'+strtrim(ix+1,2)+' ',$
	  string(xw[ix],'(f10.4)')+' '+atom(xz[ix])+rom(xi[ix])
        message,tagA[i+1]+' in A =?= '+tagB[oo[j]+1]+' in B',/info
      endfor
      if keyword_set(ask) then begin		;(ask to replace
	ca=strtrim(ask,2) & cc=strmid(strlowcase(ca),0,1)
	case cc of
	  'n': setask='n'	;reverse default
	  'r': setask='r'	;replace A's ID by B's
	  'k': setask='k'	;keep both
	  else: begin
	    if cc ne 'y' then begin
	      c1='A over B? (n:opposite, r:replace A by B, k:keep both)'+$
		' [y/n/r/k]' & print,c1 & c1=strlowcase(get_kbrd(1))
	      if c1 ne 'n' and c1 ne 'r' and c1 ne 'k' then setask='y' else $
		setask=c1
	    endif
	  end
	endcase
      endif					;asked)
    endif					;common IDs}
    case setask of				;{what to do with overlaps
      'n': outw[i]=0.			;keep B, not A
      'r': begin			;replace A by B
	noB[oo]=i
	for j=0,moo-1 do begin
	  tmpB=B.(oo[j]+1)
	  idname='ID'+strtrim(k+1,2)
	  if k eq 0 then C=create_struct(idname,tmpB) else $
		C=create_struct(C,idname,tmpB)
	  k=k+1L
	endfor
      end
      'k': begin			;keep both
	idname='ID'+strtrim(k+1,2)
	if k eq 0 then C=create_struct(idname,tmp) else $
		C=create_struct(C,idname,tmp)
	k=k+1L
      end
      else: begin			;keep A, not B
	if moo gt 0 then noB[oo]=i
	idname='ID'+strtrim(k+1,2)
	if k eq 0 then C=create_struct(idname,tmp) else $
		C=create_struct(C,idname,tmp)
	k=k+1L
      end
    endcase					;SETASK}
  endfor					;i=0,nwA-1)

  for i=0,nwB-1 do begin			;(step through B to fill C
    tmp=B.(i+1)
    if noB[i] lt 0 then begin
      idname='ID'+strtrim(k+1,2)
      if k eq 0 then C=create_struct(idname,tmp) else $
	C=create_struct(C,idname,tmp)
      k=k+1
    endif else outw[nwA+i]=0.
  endfor					;i=0,nwB-1)

  oo=where(outw ne 0,moo)
  if moo gt 0 then begin
    wC=outw[oo] & C=create_struct('WVL',wC,C)
  endif else message,'bug in program'

endif else begin				;}{just print, OK?

  C=A

  ;	extract fluxes and fluxerr
  flx=fltarr(nwA) & fle=flx
  for i=0,nwA-1 do begin
    tmp=A.(i+1) & tgnam=tag_names(tmp) & mW=n_elements(tmp.WVL)
    iflux=(where(tgnam eq 'FLUX'))[0]
    iflxe=(where(tgnam eq 'FLUXERR'))[0]
    if iflux ge 0 then flx[i]=total(tmp.FLUX)
    if iflxe ge 0 then fle[i]=sqrt(total((tmp.FLUXERR)^2))
  endfor

  ;	update fluxes or not?
  nfx=n_elements(fluxes) & updateflx=0
  if nfx gt 0 then begin
    if nfx eq nwA then updateflx=1
    if nfx gt nwA then begin
      ;	maybe size matches number of matches?
      mW=0L & for i=0,nwA-1 do mW=mW+n_elements(A.(i+1).WVL)
      if nfx eq mW then updateflx=2
    endif
    ;
    ;	so update 'em
    if updateflx gt 0 then begin
      j=0L
      for i=0,nwA-1 do begin
        tmp=A.(i+1) & tgnam=tag_names(tmp) & mW=n_elements(tmp.WVL)
        iflux=(where(tgnam eq 'FLUX'))[0]
        if iflux ge 0 then flx=tmp.FLUX else flux=fltarr(mW)
        tflx=total(flx) & if tflx eq 0 then tflx=1.
        if updateflx eq 1 then flx=flx*fluxes[i]/tflx
        if updateflx eq 2 then begin
	  flx=fluxes[j:j+mW-1] & j=j+mW
        endif
        if iflux ge 0 then C.(i+1).flux=flx
      endfor
    endif else message,'input fluxes are being ignored',/info
  endif

  ;	update fluxerr or not?
  nfe=n_elements(fluxerr) & updatefle=0
  if nfe gt 0 then begin
    ;	exact match?
    if nfe eq nwA then updatefle=1	;as many FLUXERR as components
    ;	maybe size matches number of matches?
    if nfe gt nwA then begin
      mW=0L & for i=0,nwA-1 do mW=mW+n_elements(A.(i+1).WVL)
      if nfe eq mW then updatefle=2	;as many FLUXERR as IDs
    endif
    ;	use FLUXERR[0] as fractional error
    if nfe eq 1 and nfe ne nwA then updatefle=3		;fixed error
    ;
    ;	so update 'em
    if updatefle gt 0 then begin
      j=0L & fle=fltarr(nwA)
      for i=0,nwA-1 do begin
        tmp=A.(i+1) & tgnam=tag_names(tmp) & mW=n_elements(tmp.WVL)
        iflxe=(where(tgnam eq 'FLUXERR'))[0]
        if iflxe ge 0 then flxe=tmp.FLUXERR else flxe=fltarr(mW)
	if updatefle eq 1 then begin
	  flxx=tmp.FLUX
	  tflxx=total(flxx) & if tflxx eq 0 then tflxx=1.
	  fle[i]=fluxerr[i]
	  flxe=fle[i]*flxx/tflxx
	  if iflxe ge 0 then C.(i+1).FLUXERR=flxe
	endif
	if updatefle eq 2 then begin
	  flxe=fluxerr[j:j+mW-1L] & j=j+mW
	  fle[i]=sqrt(total(flxe^2))
	  if iflxe ge 0 then C.(i+1).FLUXERR=flxe
	endif
	if updatefle eq 3 then begin
	  flxx=tmp.FLUX
	  if abs(fluxerr[0]) lt 1 then flxe=abs(fluxerr[0])*flxx else $
	    flxe=flxx/abs(fluxerr[0])
	  fle[i]=sqrt(total(flxe^2))
	  if iflxe ge 0 then C.(i+1).FLUXERR=flxe
	endif
      endfor
    endif else message,'input FLUXERR are being ignored',/info
  endif

  if not keyword_set(flst) then begin		;(normal output
    c1=string('WVL [Ang]','(a10)')
    if keyword_set(okev) then c1=string('[keV]','(a10)')
    c1=c1+'  '+string('MATCH','(a10)')+' '+string('Z ION','(a10)')+$
      string('FLUX','(a12)')+string('+-','(a12)')+string('Tmax','(a10)')
    print,c1
    for i=0,nwA-1 do begin
      tmp=C.(i+1) & tgnam=tag_names(tmp)
      ;
      iwvl=(where(tgnam eq 'WVL'))[0]
      iZ=(where(tgnam eq 'Z'))[0]
      iion=(where(tgnam eq 'ION'))[0]
      ilabl=(where(tgnam eq 'LABL'))[0]
      iflux=(where(tgnam eq 'FLUX'))[0]
      iflxe=(where(tgnam eq 'FLUXERR'))[0]
      ilogT=(where(tgnam eq 'LOGT'))[0]
      iemis=(where(tgnam eq 'EMIS'))[0]
      inote=(where(tgnam eq 'NOTES'))[0]
      if iwvl ge 0 then wvl=tmp.WVL else wvl=[0.] & mw=n_elements(wvl)
      if iZ ge 0 then Z=tmp.Z else Z=intarr(mw)
      if iion ge 0 then ion=tmp.ION else ion=intarr(mw)
      if ilabl ge 0 then labl=tmp.LABL else labl=strarr(2,mw)
      if iflux ge 0 then flux=tmp.FLUX else flux=fltarr(mw)
      if iflxe ge 0 then flxerr=tmp.FLUXERR else flxerr=fltarr(mw)
      if ilogT ge 0 then logT=tmp.logT else logT=[0.] & nt=n_elements(logT)
      if iemis ge 0 then emis=tmp.EMIS else emis=fltarr(nt,mw)
      if inote ge 0 then scratch=tmp.NOTES else scratch=''
      if labl[0] eq 'Unknown' then labl=['Un','known']
      ;
      print,strtrim(i,2)+$
	  '-----------------------------------------------------'+strtrim(i+1,2)
      for j=0,mw-1 do begin
        if keyword_set(okeV) then begin
	  wAk=12.3985/(abs(wA[i])>0.001) & wvlk=12.3985/(abs(wvl[j])>0.001)
        endif else begin
	  wAk=wA[i] & wvlk=wvl[j]
        endelse
        emismx=max(emis[*,j],imx) & tmax=logT(imx)
        c1=string(wAk,'(f10.4)')+': '+string(wvlk,'(f10.4)')+' '+$
	  string(atom(z[j])+' '+rom(ion[j]),'(a10)')+' '+$
	  string(flux[j],'(g11.5)')+' +- '+string(flxerr[j],'(g10.4)')
	if tmax ne 0 then c1=c1+string(tmax,'(g8.3)')
        c2='  '+labl[1,j]+' --> '+labl[0,j] & c2=strmid(c2,0,ncomm)
        if keyword_set(comm) then begin
	  c1=c1+c2
	  if j eq mw-1L and scratch ne '' then c1=c1+$
		string("12b)+string("15b)+scratch(0)
	endif
        print,c1
      endfor
    endfor
  endif else begin				;)(RD_LIST compatible output
    szl=size(flst) & nszl=n_elements(szl) & ulst=-1L
    if szl[nszl-2] eq 7 then openw,ulst,strtrim(flst[0],2),/get_lun
    kdb=0L & ndb=n_elements(dbdir)
    for i=0,nwA-1 do begin
      tmp=C.(i+1) & tgnam=tag_names(tmp)
      ;
      iwvl=(where(tgnam eq 'WVL'))[0]
      iZ=(where(tgnam eq 'Z'))[0]
      iion=(where(tgnam eq 'ION'))[0]
      ilabl=(where(tgnam eq 'LABL'))[0]
      iflux=(where(tgnam eq 'FLUX'))[0]
      iflxe=(where(tgnam eq 'FLUXERR'))[0]
      ilogT=(where(tgnam eq 'LOGT'))[0]
      iemis=(where(tgnam eq 'EMIS'))[0]
      inote=(where(tgnam eq 'NOTES'))[0]
      if iwvl ge 0 then wvl=tmp.WVL else wvl=[0.] & mw=n_elements(wvl)
      if iZ ge 0 then Z=tmp.Z else Z=intarr(mw)
      if iion ge 0 then ion=tmp.ION else ion=intarr(mw)
      if ilabl ge 0 then labl=tmp.LABL else labl=strarr(2,mw)
      if iflux ge 0 then flx=tmp.FLUX else flx=fltarr(mw)
      if iflxe ge 0 then flxerr=tmp.FLUXERR else flxerr=fle(i)
      if ilogT ge 0 then logT=tmp.logT else logT=[0.] & nt=n_elements(logT)
      if iemis ge 0 then emis=tmp.EMIS else emis=fltarr(nt,mw)
      if inote ge 0 then scratch=tmp.NOTES else scratch=''
      if labl[0] eq 'Unknown' then labl=['Un','known']
      ;
      mw=n_elements(wvl) & tflux=total(flx) & tflxe=sqrt(total(flxerr^2))
      c1='/{CAT_ID	' & c2=strarr(mw)
      wAk=abs(wA[i]) > 0.001
      if keyword_set(okeV) then wAk=12.3985/wAk
      c1=c1+strtrim(wAk,2)
      for j=0,mw-1 do begin
	wvlk=abs(wvl[j]) > 0.001
	if keyword_set(okeV) then wvlk=12.3985/wvlk
	c2[j]=atom(z[j])+' '+rom(ion[j])+'	'+$
		strtrim(string(wvlk,'(f10.4)'),2)
	if ndb ne 0 then begin		;(add DBDIR
	  cdb=strtrim(dbdir([kdb]),2)
	  case cdb(0) of
	   '0': ;nothing
	   '1': c2[j]=c2[j]+'	$SCAR'
	   '2': c2[j]=c2[j]+'	$SPEX'
	   '4': c2[j]=c2[j]+'	$CHIANTI'
	   else: c2[j]=c2[j]+'	'+cdb
	  endcase
	  if ndb eq nwA then begin
	    if j eq mw-1 then kdb=kdb+1		;one DBDIR per WVL
	  endif else begin
	    kdb=kdb+1		;one DBDIR per ID
	  endelse
          c3='	'+labl[0,j]+' '+labl[1,j] & c3=strmid(c3,0,ncomm)
          if keyword_set(comm) then begin
	    c2[j]=c2[j]+c3
	    if j eq mw-1L and scratch ne '' then c2[j]=c2[j]+'	'+scratch(0)
	  endif
	endif				;DBDIR)
      endfor
      c1=c1+'	'+strtrim(tflux,2)
      ;if updatefle gt 0 then c1=c1+'	'+strtrim(tflxe,2)
      c1=c1+'	'+strtrim(tflxe,2)
      printf,ulst,c1 & for j=0,mw-1 do printf,ulst,c2[j]
      if scratch[0] ne '' then printf,ulst,'//	'+scratch[0]
      printf,ulst,'/}CAT_ID'
    endfor
    if szl[nszl-2] eq 7 then begin & close,ulst,/all & free_lun,ulst & endif
  endelse					;FLST)

endelse						;nwB=0}


if arg_present(mplet) then begin 
    nid = n_tags(C) 
    tmplet = strarr(2,2)
    contr = -1L
    for i=1L,nid-1L do begin 
        tmp=C.(i) & zz = tmp.z 
        ii = tmp.ion & ll = tmp.labl
        nl  = n_elements(zz) 
        for j = 0,nl-1L do begin  
            contr = contr+1L
            for k = 0,nl-1L do begin 
                if k ne j then begin 
                    pckg=[strcompress(string(contr),/remove_all),atom(tmp.z(k)) $
                         + rom(tmp.ion(k)) + ' ' + strjoin(tmp.labl(k))]
                    tmplet = [tmplet,transpose(pckg)]
                endif
            endfor
        endfor  
    endfor 
     if n_elements(tmplet) ne 4L then mplet = transpose(tmplet[2:*,*]) else mplet='none'
endif 
return,C
end

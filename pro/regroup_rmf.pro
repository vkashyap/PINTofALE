function regroup_rmf, rmfstr, enrg, spec,nb, numatrix,noeffar=noeffar

;+
;function   regroup_rmf
;  
;        Converts RMF structure created with  RD_OGIP_RMF() to a 
;        one group per photon energy RMF structure. Done so to 
;        maximize convolution efficiency with SCRMF() e.g. for 
;        model fitting via RMFs in FITLINES(). Will return input
;        if RMF is already optimal. If effective area is folded 
;        into matrix and this is undesired, REGROUP_RMF is also 
;        capable of returning an unfolded RMF via the NOEFFAR 
;        keyword.
;
;        If ENRG paramater is fed,  will be used to define the photon
;        energy grid of the reconstructed RMF. 
; 
;        If a spectrum SPEC is fed together with ENRG, then 
;        the photon energy range will be constructed so that
;        bin-density is highest at energy ranges with largest 
;        signal. 
;
;        If NB is set, will be used as the total number of bins 
;        for resulting photon energy grid.
;
;syntax  nrmfstr = regroup_rmf(rmfstr,enrg,spec,nb,/noeffar)
;
;parameters 
;        rmfstr  [INPUT;required] structure containing info of 
;                response matrix, the output of RD_OGIP_RMF()      
;        enrg    [INPUT] desired energy grid
;        spec    [INPUT] spectrum of interest.  if given, should
;		 be on the same grid as ENRG, and is used to tweak
;		 the energy grid in the output such that the
;		 bin-density is correlated with SPEC.
;		 * overrides ENRG except that the energy range in
;		   the output is taken from ENRG.
;        nb      [INPUT] total number of photon energy bins desired.
;		 * is used only when SPEC is defined, and by default
;		   is set to N(ENRG) if not given
;
;keywords 
;        noeffar [INPUT] if set, the LSF corresponding to each 
;                photon energy will be normalized to 1 i.e. 
;                effective area is unfolded.
;		 WARNING: Because photons are _lost_ at the edges
;		 of applicability of the RMFs, the renormalization
;		 will not necessarily be correct in these regions.
;history 
;        liwei lin (Aug 03) 
;               ADDED option to regrid rmf in enrgy space via 
;                enrg, spec, and nb (LL Sep03) 
;               BUG FIX CF handling/interpolation re-worked to avoid 
;                shifting of grid (LL Sep03)  
;               BUG FIX now robust to zero vavlues for N_GRP. 
;                thanks Jan-Uwe (LL Jun06)
;		BUG FIX changed condition for when regroup is done
;		 thanks Jan-Uwe (VK Dec07)
;- 

ok = 'ok' & np=n_params() & mr = n_tags(rmfstr) 
if np lt 0 then ok = 'Insufficient parameters' else $ 
if mr eq 0 then ok = 'RMFSTR must be a structure' else begin 
 rname=tag_names(rmfstr)
       o1=(where(strpos(rname,'NNRG',0) ge 0))[0]
       o2=(where(strpos(rname,'ELO',0) ge 0))[0]
       o3=(where(strpos(rname,'EHI',0) ge 0))[0]
       o4=(where(strpos(rname,'NCHAN',0) ge 0))[0]
       o5=(where(strpos(rname,'EMN',0) ge 0))[0]
       o6=(where(strpos(rname,'EMX',0) ge 0))[0]
       o7=(where(strpos(rname,'N_GRP',0) ge 0))[0]
       o8=(where(strpos(rname,'F_CHAN',0) ge 0))[0]
       o9=(where(strpos(rname,'N_CHAN',0) ge 0))[0]
       o10=(where(strpos(rname,'MATRIX',0) ge 0))[0]
       o11=(where(strpos(rname,'FIRSTCHAN',0) ge 0))[0]
   if o1 lt 0 or o2 lt 0 or o3 lt 0 or o4 lt 0 or o5 lt 0 or o6 lt 0 $
       or o7 lt 0 or o8 lt 0 or o9 lt 0 or o10 lt 0 or o11 lt 0 then $
       ok='RMFSTR not in standard format, see RD_OGIP_RMF()'  
endelse

if ok ne 'ok' then begin 
   print, 'Usage: nrmfstr=regroup_rmf(rmfstr,enrg,spec,nb,/noeffar)'
   if np ne 0 then message, ok, /info 
   return, -1
endif

firstchan = rmfstr.firstchan ;	first channel index -- 0-based or 1-based
nnrg   = rmfstr.nnrg   ;number of nrg bins 
nchan  = rmfstr.nchan  ;number of total channels
n_grp  = rmfstr.n_grp  ;number of groups for each phtn nrg INT   [nnrg]
f_chan = rmfstr.f_chan ;first channel of group             INT   [ngrp, nnrg]
n_chan = rmfstr.n_chan ;number of channels in each grp     INT   [ngrp, nnrg]
matrix = rmfstr.matrix ;response matrix                    FLOAT [ !varies! , nnrg] 
elo    = rmfstr.elo    ;lower phtn nrg bin bounds
ehi    = rmfstr.ehi    ;upper phtn nrg bin bounds 
emn    = rmfstr.emn    ;lower channel nrg bin bounds 
emx    = rmfstr.emx    ;upper channel nrg bin bounds 
   
egrid = (elo+ehi)/2    ; matrix phtn nrg grid
changrd = (emn+emx)/2  ; matrix channel nrg grid 
nspc = n_elements(spec)
nenrg = n_elements(enrg)

if nenrg gt 0 then onrg = enrg            ; avoid returning alterred enrg array 
if nenrg eq nspc+1 then enrg = (enrg[0:nenrg-1]+enrg[1:*])/2 & nenrg = n_elements(enrg)
if nenrg ne nspc and np gt 2 then begin 
  message, 'enrg and spec array mismatch', /info 
  return, -1 
endif 

if nspc gt 0 then begin 
   if not n_elements(nb) gt 0 then  nb = nspc 
   if (min((enrg(1:*)-enrg)<0)) ne 0 then begin ;more general than  total((xx(1:*)-xx)<0) ne 0
   message, 'resorting enrg and spec to be in ascending order', /info
   esort = sort(enrg) & enrg = enrg(esort) & spec = spec(esort) 
   endif 
   CF = fltarr(nspc+1)
   for i = 1,nspc do CF[i] = CF[i-1L]+spec(i-1L) ; + 1 ?!! who put this in and why?
   CF = CF[1:*]/max(CF)
   enrg = interpol(enrg,CF,(1-min(cf))*findgen(nb)/nb) 
   nenrg = n_elements(enrg) 
endif

if nenrg gt 0 then begin
   jnk = min(abs(max(enrg)-egrid),bg) & jnk = min(abs(min(enrg)-egrid),sm)
   mnie = sm<bg  & mxie = sm>bg  
   nunnrg = nenrg 
   regrid = 1 
   fmatrix = fltarr(nchan,nnrg)
endif else begin 
   regrid =0
   mnie = 0 & mxie = nnrg-1 & nunnrg = nnrg        
endelse
 
  nuf_chan = intarr(nnrg)+1  ;new f_chan
  nun_chan = intarr(nnrg)+1  ;new n_chan
  dn =  mnie-1>0 & up = mxie<(nnrg-1) 
;{ORIGINAL if (total(n_grp) ne n_elements(n_grp)) then begin ; if regroup neccessary
;this is not general enough.  you could have a case where N_GRP is 0 for some
;energies, which will violate the equality condition, but there is still
;nothing that needs to be done.  so instead check to see that F_CHAN is 2D.
;VLK}
if (size(f_chan))(0) ne 1 then begin	; if regroup necessary
  regroup = 1 
  for bb = mnie, mxie do begin        ;
      tmp = max(f_chan(*,bb),mxf) ;
      nuf_chan(bb) = min(f_chan(0:n_grp(bb)-1>0,bb)) ;prelim runthrough to see  
      nun_chan(bb) = (tmp+n_chan(mxf,bb))-nuf_chan(bb) ;how small a matrix we can 
  endfor                                             ;afford. also calculate new 
  nulngth = max(nun_chan)       ;f_chan and n_chan
  numatrix = fltarr(nulngth, nnrg) 
  for cc = dn, up do begin   ;loop through phtn nrgs 
      profc = fltarr(nchan) & ipcw = 0 
      for pp = 0, n_grp(cc)-1 do begin ;loop through the grps for this nrg bin
          cbeg = f_chan(pp,cc) & cw = n_chan(pp,cc)
          if firstchan eq 1L then cbeg=cbeg-1L>0 ;IDL index correction
          if pp eq 0 then pcw=0 else pcw=n_chan((pp-1>0),cc)
          ipcw = ipcw+pcw
          profc[cbeg:cbeg+(cw-1>0)] = matrix[ipcw:ipcw+(cw-1>0),cc]
      endfor                    ;end loop through groups 
      if regrid eq 1 then fmatrix[*,cc] = profc
      if firstchan eq 1L then corig=nuf_chan(cc)-1>0 else corig=nuf_chan(cc) 
      profc = profc[corig:corig+nun_chan(cc)-1>0] ;smush so that profile support
      if n_elements(profc) gt nulngth then profc = profc[0:nulngth-1>0]; overcompensation corction
      crctn = nulngth-n_elements(profc)         ;begins at index of 0 and is 
      if crctn gt 0 then profc = [profc, fltarr(crctn)] ;of adequate size for new matrix
      numatrix[*,cc] = profc
  endfor; end loop through phtn nrgs
  matrix=numatrix[*,mnie:mxie]  & f_chan=nuf_chan[mnie:mxie] & n_chan = nun_chan[mnie:mxie]
  elo = elo[mnie:mxie] & ehi=ehi[mnie:mxie] & nnrg = nunnrg & n_grp=intarr(nnrg)+1L  
endif else regroup =0                           ; end if regroup needed

if regroup eq 0 then message, 'No re-grouping necessary...', /info
if not arg_present(enrg) and regroup eq 0 and not keyword_set(noeffar) then begin 
  message, 'RMF already optimal. Nothing to be done', /info 
return, rmfstr
endif 

if regrid eq 1 then begin ;if regrid required then 
  if regroup ne 1 then begin ;if regouping was not done thus fmatrix not populated then populate
    for cc = dn, up do begin ; loop through phtn nrgs  
       cbeg = f_chan(cc) & cw = n_chan(cc)          
       profc = fltarr(nchan)
       if firstchan eq 1L then cbeg=cbeg-1L ;IDL index correction
       profc[cbeg:cbeg+cw-1] = matrix[0:cw-1,cc]               
       fmatrix[*,cc] = profc
    endfor ; end loop through phtn nrgs
  endif    ; end if fmatrix was not populated during regrouping
  ; adjust channel energy grids to accomodate interpolation 
    elo = rmfstr.elo[dn:up] & ehi=rmfstr.ehi[dn:up] ;& nnrg = nunnrg
    fmatrix = fmatrix[*,dn:up] ; chop to necessary size
    egrid =(elo+ehi)/2   ; matrix phtn nrg grid
  ; get index mapping of phtn nrg in spec grd index space
    row_grid = interpol(findgen(n_elements(egrid)), egrid, enrg) 
  ; use linear interpolation to create new response matrix 
    numatrix = interpolate(fmatrix, row_grid,/GRID)  
    nnrg = nenrg   
     f_chan = intarr(nnrg) & n_chan = intarr(nnrg) 
    for cc = 0, nnrg-1 do begin ;loop through phtn nrgs 
        gt0 = where(numatrix[*,cc] gt 0) 
        f_chan(cc) = min(gt0) & n_chan(cc) = max(gt0)-min(gt0)+1
    endfor  ; end loop through phtn nrgs
    mn_chan = max(n_chan)
  ; smush matrix 
    finmatrix = fltarr(mn_chan, nnrg)
    for cc = 0, nnrg-1 do finmatrix[0:n_chan(cc)-1,cc] = $
    numatrix[f_chan(cc):f_chan(cc)+n_chan(cc)-1,cc]
    matrix = finmatrix
  ; redifine elo, ehi, emn, emx, n_grp 
    n_grp = intarr(nnrg) + 1L
    new_grd = mid2bound(enrg, _extra = e)  
  if firstchan eq 1L then f_chan = f_chan +1
  if new_grd(1) gt new_grd(0) then begin 
  elo = new_grd[0:nnrg-1] & ehi = new_grd[1:*] 
  endif else begin  
  elo = new_grd[1:*] & ehi = new_grd[0:nnrg-1]
  endelse
   ; emn = new_grd[0:nchan-1] & emx = new_grd[1:*] 
endif; if regrid required if

if nenrg gt 0 then enrg = onrg ; avoid returning alterred enrg array 

if keyword_set(noeffar) then begin ; unfold effective area if asked to
  for dd = 0, nnrg-1 do matrix[*,dd] = matrix[*,dd]/total(matrix[*,dd])
endif; end if unfold needed 

  nrmfstr=create_struct('NNRG',nnrg,'ELO',elo,'EHI',ehi,$
                        'NCHAN',nchan,'EMN',emn,'EMX',emx,$
                        'N_GRP',n_grp,'F_CHAN',f_chan,'N_CHAN',n_chan,$
                        'MATRIX',matrix,'FIRSTCHAN',firstchan)  
return, nrmfstr
end 

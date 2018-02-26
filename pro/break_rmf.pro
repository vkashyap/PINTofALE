function break_rmf, rmfstr, erng,effar=effar, _extra=e 

;+ 
;procedure    break_rmf 
;       extract a sub-matrix from response matrix and  associated 
;       axes in OGIP compliant form. both input and output are given 
;       in a structure of the form: 
;
;	{NNRG, ELO, EHI, NCHAN, EMN, EMX, N_GRP, F_CHAN, N_CHAN, MATRIX, FIRSTCHAN}
;       
;	where ELO and EHI refer to range of photon energies at which the
;	response is valid, EMN and EMX refer to bin boundaries for each
;	channel, N_GRP refers to number of groups of non-zero data, F_CHAN
;	refer to beginning indices of channels, N_CHAN are the number of
;	channels in each group, and the output response MATRIX excludes
;	zeros to save space.  The value of FIRSTCHAN indicates whether F_CHAN
;	indices are 0-based or 1-based. SEE RD_OGIP_RMF().
;
;	http://legacy.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html
;	for a description of the format
;
;syntax
;	break_rmf, rmfstr, erng, _extra=e 
;
;parameters
;	rmfstr	[INPUT; required] response matrix structure of the form : 
;               {NNRG, ELO, EHI, NCHAN, EMN, EMX, N_GRP, F_CHAN, N_CHAN, MATRIX,
;               FIRSTCHAN} as read in by e.g. RD_OGIP_RMF() 
;       erng    [INPUT; required] two element array designating the channel 
;               energy range to which to limit the output rmf 
;
;keywords 
;       effar	[OUTPUT] returns the effective area as a function of ELO
;		* if an RMF file is read instead of an RSP file, EFFAR are
;		uniformly 1, unless it is ASCA in which case it is the QEs
;      
;history 
;       liwei lin (Dec 2004)  

;       usage 
ok = 'ok' & np=n_params() & nrng = n_elements(erng)
if np lt 2 then ok='Insufficient parameters' else $ 
  if nrng ne 2 then ok='ERNG must be 2-element array'
if ok ne 'ok' then begin 
  print, 'Usage: break_rmf, rmfstr, erng' 
  if np ne 0 then message,ok,/info 
  return, -1L 
endif 
if n_tags(rmfstr) eq 0 then return, -1L 

;       deconstruct rmf structure 
NNRG=rmfstr.NNRG & ELO=rmfstr.ELO & EHI=rmfstr.EHI
N_GRP=rmfstr.N_GRP & F_CHAN=rmfstr.F_CHAN 
N_CHAN=rmfstr.N_CHAN & FIRSTCHAN=rmfstr.FIRSTCHAN
MATRIX=rmfstr.MATRIX & NCHAN=rmfstr.NCHAN 
EMN=rmfstr.EMN & EMX=rmfstr.EMX
szn=size(N_CHAN)

;       set to zero all elements ouside erng  
chan = [emn+emx]/2  
oo   = where(chan lt erng(0) or chan gt erng(1),n)  
ooi  = where(chan ge erng(0) and chan lt erng(1),ni)  
matrix[oo,*] = -1.0 
nmatrix = fltarr(nchan-n,nnrg)  & nngrp = n_grp
nf_chan = f_chan  & nn_chan = n_chan 

;       compress response matrix 
for j = 0, nnrg-1 do begin  
    ogipzip, matrix[*,j],carr,ngrp,fchan,nchano,eps=-0.01,chan0=firstchan,cchan=cchan
    nmatrix[*,j]=carr  & nngrp(j) = ngrp 
    nf_chan(j) = fchan & nn_chan(j) = nchano      
endfor

nrmfstr=create_struct('NNRG',nnrg,'ELO',elo,'EHI',ehi,$
                      'NCHAN',ni,'EMN',emn(ooi),'EMX',emx(ooi),$
                      'N_GRP',nngrp,'F_CHAN',nf_chan,'N_CHAN',nn_chan,$
                      'MATRIX',nmatrix,'FIRSTCHAN',firstchan)  

;      effar  
    effar = fltarr(nnrg)  
    for j = 0, nnrg-1 do effar(j) = total(nmatrix[*,j]) 

return, nrmfstr

end 

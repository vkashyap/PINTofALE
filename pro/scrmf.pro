pro scrmf, nrg,flx,chan,spec,rmf,effar=effar,nrgar=nrgar,rmfstr=rmfstr,$
	_extra=e
;+
;procedure	scrmf 
;	Convolves input energy spectrum with the response matrix
;       
;       SCRMF is a version of CONV_RMF optimized for speed. 
;       It maximized speed at the expense of memory. 
;       A faster version of CONV_RMF was necessary for  
;       use with e.g., FITLINES() because the convolution procedure is
;       called literally tens of thousands of times. 
; 
;       On input, the RMF is required to be an RMF structure created
;       by RD_OGIP_RMF() from an OGIP-compliant FITS file AND
;       optimized with REGROUP_RMF(). 
;
;       NOTE: If the RMF is not optimized prior to
;       the SCRMF call, REGROUP_RMF() will be called 
;       internally, and SCRMF will be slowed to a crawl. 
;                      
;syntax
;	scrmf,nrg,flx,chan,spec,rmf,effar=effar,nrgar=nrgar,$
;	rmfstr=rmfstr,fchcol=fchcol,/shift1
;
;parameters
;	nrg	[INPUT; required] mid-bin values at which input FLX is defined.
;               * preferably NRG is identical to RMF channel energy grid 
;		* depends on the RMF, but is invariably [keV]
;		* if size = N(FLX)+1, assumed to be bin-boundaries
;	flx	[INPUT; required] fluxes in [ph/bin/(s/sr/...)]
;		* note that this assumes that the effective areas are
;		  already folded in.  this assumption is made in order
;		  to speed up the calculation some more
;		* the input should also be integrated over the bin width,
;		  i.e., do not give [ph/keV/...] or [ph/Ang/...] as input
;		* also note, the units should not be [erg/...] -- be sure
;		  to divide the input energy spectrum by the energy of
;		  the photons corresponding to that bin before calling
;		  the program
;	chan	[OUTPUT; required] bin boundary values at which output SPEC
;		is defined. 
;               * same energy range as NRG (so size will not necesarily match
;		  EMN,EMX of RMF)
;		* same units as RMF.EMN, RMF.EMX
;	spec	[OUTPUT; required] output, convolved spectrum
;	rmf	[I/O; required] rmf response structure in RD_OGIP_RMF() format. 
;		WARNING: Unlike CONV_RMF, this _cannot_ be the name of a file.
;		Also, the input will be overwritten if necessary.
;               NOTE: For most efficient use with FITLINES(), RMFs which
;               are defined with only one group per photon energy bin are
;               preferred. If RMF contains more than one group per
;               photon energy, the RMF will be optimized for FITLINES()
;               use with REGROUP_RMF(). 
; 
;keywords
;       verbose [INPUT] set verbosity level
;	effar	[INPUT] (not implemented yet)
;	nrgar	[INPUT] (not implemented yet)
;	rmfstr	[OUTPUT] structure containing info of response matrix,
;		the output of RD_OGIP_RMF()
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;
;subroutines
;	REGROUP_RMF
;
;these are all the smoothing tools in PINTofALE
;	ALALOESS() makes a loess curve
;	CLTSMOOTH() removes -ves from background subtracted spectra or light curves
;	CONV_RMF convolves with specified RMF
;	HAARTRAN() smoothes by filtering on Haar wavelet coefficients
;	LINEREM() removes lines from spectrum by iterative S/N filtering
;	NOISMOOTH() does boxcar accumulation a la Ebeling's asmooth
;	REGROUP() accumulates from one end
;	SCRMF does fast convolution using optimized RMF
;	SMOOTHIE does peak-descent accumulation to build up S/N
;	SPLAC computes a segmented piecewise linear approximations
;	UNKINK() removes sharp kinks from a curve
;	VARSMOOTH() does local smoothing at specified scales
;	VOORSMOOTH() does fast local smoothing at specified scales
;
;history
;      liwei lin  (Aug 03)    
;            BUGFIX crashed with no warning when input spectrum 
;            not fully covered. Thanks Jan-Uwe. LL (Jun 06) 
;            BUGFIX crashed with zero values for N_CHAN. 
;            Thanks Jan-Uwe. LL (Jun 06) 
;	BUGFIX crashed when FIRSTCHAN=1 for IBEG=0 (thanks Jan-Uwe; Vinay K; May 07)
;-

;	usage
ok='ok'
np=n_params() & mnrg=n_elements(nrg) & mflx=n_elements(flx)
szr=size(rmf) & nszr=n_elements(szr)
nrmf=n_elements(rmf) & mrmf=n_tags(rmf)
if np lt 5 then ok='missing parameters' else $
 if mnrg eq 0 then ok='NRG is undefined: require spectrum grid' else $
  if mflx eq 0 then ok='FLX is undefined: require spectrum' else $
   if mnrg ne mflx and mnrg ne (mflx+1L) then ok='input spectrum garbled' else $
    if arg_present(rmf) and szr[nszr-2] ne 7 and mrmf eq 0 then $
     ok='RMF must be either a filename or a response matrix structure ' 
if ok ne 'ok' then begin
  print,'Usage: scrmf, nrg,flx,chan,spec,rmf,rmfstr=rmfstr'
  print,'convolve spectrum with response'
  if np gt 0 then message,ok,/info
  return
endif 

if not keyword_set(verbose) then verbose = 2 

ee=[nrg[*]] & fe=[flx[*]] 
if mnrg eq mflx+1 then begin
  ee=[0.5*(nrg[1:*]+nrg)]
  mnrg=mnrg-1L
  if mnrg ne mflx then begin
  if verbose gt 4 then message,'Energy grid and flux grid do not match',/info 
  return
  endif
endif
NNRG=rmf.nnrg 

if total(rmf.n_grp) ne nnrg then begin 
  if verbose gt 4 then message,'RMF is not optimized for SCRMF(). Optimizing with REGROUP_RMF()...',/informational
  rmf = regroup_rmf(rmf)
  rmfstr = rmf
endif 

 EHI       = rmf.EHI & ELO = rmf.ELO 
 FIRSTCHAN = rmf.FIRSTCHAN  & MATRIX = rmf.MATRIX
 EMN = rmf.EMN & EMX=rmf.EMX & nchan = rmf.nchan
 f_chan    = rmf.f_chan & n_chan = rmf.n_chan 
 chan = [EMN,EMX[nchan-1L]]
 chanmd = [chan[1:*]+chan]/2 
 spec = fltarr(nchan) 
 ie   = binersp(ee,elo,ehi)     
if min(ie) eq -1 then $
 message, 'WARNING: RMF does not cover full range of input spectrum.',/info

if (mflx gt nnrg) then begin  
 mnie = min(ie) & mxie = max(ie) 
  for inrg=mnie,mxie do begin        
      oi=where(ie eq inrg,moi)
      ibeg = f_chan[inrg]
      iw   = n_chan[inrg] & rsp = fltarr(nchan)
      if keyword_set(firstchan) then ibeg=ibeg-1	;IDL index correction 
      if moi gt 0 then rsp[ibeg:ibeg+iw-1L>0]=matrix[0:iw-1>0,inrg]*total(fe(oi))
      spec = spec+rsp
  endfor 
endif else begin 
    for i=0L,mflx-1L do begin
      inrg=ie[i] & rsp = fltarr(nchan)
      ;this produces a crash if FIRSTCHAN=1 (thx Jan-Uwe): if inrg ne -1 then begin
      if inrg gt firstchan-1 then begin
          ibeg =f_chan[inrg] & iw = n_chan[inrg] 
          if keyword_set(firstchan) then ibeg=ibeg-1 ;IDL index correction
	  ;	JWN also suggested if statement and removal of ">0" checks:
          ;	rsp[ibeg:ibeg+iw-1L>0]=matrix[0:iw-1>0,inrg]*fe(i)
          if iw gt 0 then rsp[ibeg:ibeg+iw-1L]=matrix[0:iw-1,inrg]*fe(i)
          spec=spec+rsp
      endif
    endfor  
endelse

jnk1 = min( abs(max(nrg)-chanmd), mxieo) 
jnk2 = min( abs(min(nrg)-chanmd), mnieo) 
spec = spec(mnieo<mxieo:mnieo>mxieo) 
chan = chan(mnieo<mxieo:mnieo>mxieo)
message, 'done convolving...',/info
return
end


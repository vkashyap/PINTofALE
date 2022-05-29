function colbehr,Ssrc,Msrc,Hsrc, Sbkg=Sbkg,Mbkg=Mbkg,Hbkg=Hbkg,$
	Sarea=Sarea,Marea=Marea,Harea=Harea, Seff=Seff,Meff=Meff,Heff=Heff,$
	Sidx=Sidx,Midx=Midx,Hidx=Hidx, Sscl=Sscl,Mscl=Mscl,Hscl=Hscl,$
	Stable=Stable,Mtable=Mtable,Htable=Htable,$
	nsim=nsim,nburnin=nburnin,post=post,level=level,$
	outputF=outputF,BEHRdir=BEHRdir,verbose=verbose, _extra=e
;+
;function	colbehr
;	computes hardness ratios based on 3-band data by invoking
;	BEHR (Bayesian Estimate of Hardness Ratios) twice, and
;	returns the inputs as well as the following quantities
;	and their credible ranges in a structure:
;	  {S, M, H, lS=log10(S), lM=log10(M), lH=log10(H),
;	  SpM=S+M, MpH=M+H, SpH=S+H, T=S+M+H,
;	  SmM=S-M, MmH=M-H, SmH=S-H,
;	  R1=S/M, R2=M/H, R3=S/H,
;	  C1=C_SM, C2=C_MH, C3=C_SH,
;	  HR1=(S-M)/(S+M), HR2=(M-H)/(M+H), HR3=(S-H)/(S+H),
;	  HRA=(S-M)/(S+M+H), HRB=(M-H)/(S+M+H), HRC=(S-H)/(S+M+H)}
;
;	Warning: Because of the necessity of combining two separate
;	runs of BEHR, only the MCMC option is used.  Thus, if the
;	length of the chain is small, the computed values may be
;	subject to computational instability.
;
;	Warning: BEHR assumes that the passbands in question
;	do not overlap, and that the counts input to the program
;	are statistically independent.  It is up to the users to
;	ensure the validity of this assumption.  No checks are
;	made to verify it.
;
;	Reference:
;	"Bayesian Estimation of Hardness Ratios: Modeling and Computations",
;	  Park, T., Kashyap, V.L., Siemiginowska, A., van Dyk, D., Zezas, A.,
;	  Heinke, C., and Wargelin, B., 2006, ApJ, 652, 610
;	http://hea-www.harvard.edu/AstroStat/BEHR/
;
;syntax
;	behr3=colbehr(Ssrc,Msrc,Hsrc,Sbkg=Sbkg,Mbkg=Mbkg,Hbkg=Hbkg,$
;	Sarea=Sarea,Marea=Marea,Harea=Harea,$
;	Seff=Seff,Meff=Meff,Heff=Heff,$
;	Sidx=Sidx,Midx=Midx,Hidx=Hidx,$
;	Sscl=Sscl,Mscl=Mscl,Hscl=Hscl,$
;	Stable=Stable,Mtable=Mtable,Htable=Htable,$
;	/post,level=level,nsim=nsim,nburnin=nburnin,hpd=hpd,$
;	outputf=outputf,BEHRdir=BEHRdir, verbose=verbose)
;
;parameters
;	Ssrc	[INPUT; required] counts in source region in the soft (S) band
;	Msrc	[INPUT; required] counts in source region in the medium (M) band
;	Hsrc	[INPUT; required] counts in source region in the hard (H) band
;		* can be arrays; if so, the array with the most elements
;		  determines the size of the output and the shortfalls in the
;		  others, if any, are made up by replicating the first elements
;
;keywords
;	Sbkg	[INPUT] counts in background region in the S band
;	Mbkg	[INPUT] counts in background region in the M band
;	Hbkg	[INPUT] counts in background region in the H band
;		* if not given, assumed to be 0
;		* if size smaller than Xsrc, first element gets replicated
;	Sarea	[INPUT] background scaling factor in the S band
;	Marea	[INPUT] background scaling factor in the M band
;	Harea	[INPUT] background scaling factor in the H band
;		* (background region area)/(source region area)
;		* if not given, assumed to be 1
;		* can also include differences in exposure time into
;		  the ratio, in the same manner as geometric area
;		* if size smaller than Xsrc, first element gets replicated
;	Seff	[INPUT] effective area in S band
;	Meff	[INPUT] effective area in M band
;	Heff	[INPUT] effective area in H band
;		* if none are set, all are assumed to be 1,
;		  else if one is set, all are assumed to be equal to that one,
;		  else if two are set and unequal, third is assumed to be 1
;		* can also be the effective area relative to some
;		  special point on the detector (e.g., aimpoint)
;		  or even some specific detector (e.g., ACIS-I v/s ACIS-S)
;		* if size smaller than Xsrc, first element gets replicated
;	Sidx	[INPUT] index of prior on S (range = 0+)
;	Midx	[INPUT] index of prior on M (range = 0+)
;	Hidx	[INPUT] index of prior on H (range = 0+)
;		* if none are set, all are assumed to be 0.5,
;		  else if one is set, all are assumed to be equal to that one
;		  else if two are set and are unequal, third is assumed to be 0.5
;		* if size smaller than Xsrc, first element gets replicated
;		* similar to AGAMMA of PPD_SRC()
;	Sscl	[INPUT] scaling index of prior on Ssrc
;	Mscl	[INPUT] scaling index of prior on Msrc
;	Hscl	[INPUT] scaling index of prior on Hsrc
;		* if none are set, all are assumed to be 0
;		  else if one is set, all are assumed to be equal to that one
;		  else if two are set and are unequal, third is assumed to be 0
;		* if size smaller than Xsrc, first element gets replicated
;		* similar to BGAMMA of PPD_SRC()
;	Stable	[INPUT] filename containing a tabulated prior for Ssrc
;	Mtable	[INPUT] filename containing a tabulated prior for Msrc
;	Htable	[INPUT] filename containing a tabulated prior for Hsrc
;		* the table prior must be an ascii file with the following format:
;		  line 1: number of entries, say NLIN
;		  line 2: labels for the columns, ignored
;		  lines 3..NLIN+2: two whitespace separated columns of numbers,
;	            with each row containing the source intensity and the posterior
;	            density, in that order
;		* the default filenames are "./tblprior_{soft|med|hard}.txt"
;		* the default filenames are used iff Stable, Mtable, and Htable are set
;		  but are not found
;		* WARNING: if regex is used in the filename specification, only the
;		  first file from the list will be used.  furthermore, if specified,
;		  the table priors are applied to _all_ SSRC, MSRC, and HSRC
;	post	[INPUT] if set, suggests the values of (Sidx,Sscl), (Midx,Mscl),
;		and (Hidx,Hscl) going forward, i.e., what you should set the
;		priors to in your next calculation for the same source -- the
;		suggested values are stored in the output structure
;	level	[INPUT] percentage confidence level at which to report error
;		(default = 68)
;	details	[INPUT] compute various ratios (true/false)?
;		(default = true)
;	nsim	[INPUT] number of draws if algo=gibbs (default=10000)
;	nburnin	[INPUT] number of burn-in draws if algo=gibbs
;		(default=5000 or NSIM/2, whichever is smaller)
;	outputF	[INPUT] root of filename in which to place output
;		(default = 'none')
;		* output will be placed in the files OUTPUTF.txt and OUTPUTF_draws.txt
;		* NOTE: if OUTPUTF='none', then MC draws will be in BEHR_draws.txt
;	BEHRdir	[INPUT] full path to directory where BEHR executable resides
;		(default = '/fubar/kashyap/AstroStat/BEHR')
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;requirements
;	uses subroutines HIPD_INTERVAL() and MODALPOINT()
;	BEHR executable must be installed in BEHRDIR
;	BEHR should be executable under the shell via SPAWN
;	BEHR output assumed to be compatible with 12-19-2005 version
;
;side-effects
;	potentially creates numerous ascii files in $cwd or `basedir OUTPUTF`
;
;history
;	vinay kashyap (Mar07; based on behr_hug.pro)
;	bug correction with NaNs not being caught in some cases
;	  (VK; Mar07)
;	added keywords Stable,Mtable,Htable (VK; Feb08)
;	bugfix: wasn't reading HTABLE (VK; Oct16)
;
;etymology
;	getting color-color diagrams using BEHR
;	(that's my story and I'm sticking to it)
;-

vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	where is it?
if not keyword_set(BEHRdir) then begin
  stdout='BEHR: Command not found.'
  if !version.os_family eq 'unix' then spawn,'which BEHR',stdout
  if strpos(strlowcase(stdout[0]),'command not found',0) ge 0 then $
	BEHRdir='/fubar/kashyap/AstroStat/BEHR' else $
	BEHRdir=stdout[0]
endif
BEHRexe=filepath('BEHR',root_dir=BEHRdir)
if vv ge 5 then message,'Using : '+BEHRexe,/informational

;	usage
ok='ok'
np=n_params() & nss=n_elements(Ssrc) & nms=n_elements(Msrc) & nhs=n_elements(Hsrc)
if np lt 3 then ok='Insufficient parameters'
;
spawn,BEHRexe,usage
if strpos(strlowcase(usage[0]),'command not found',0) ge 0 then $
	ok='BEHR: command not found'
;
if nss gt 0 then cmd=BEHRexe+' Ssrc='+strtrim(Ssrc[0],2) else $
  cmd=BEHRexe+' Ssrc=0'
spawn,cmd,help

if ok ne 'ok' then begin
  print,'Usage: behr3=colbehr(Ssrc,Msrc,Hsrc, Sbkg=Sbkg,Mbkg=Mbkg,Hbkg=Hbkg,$
  print,'       Sarea=Sarea,Marea=Marea,Harea=Harea, Seff=Seff,Meff=Meff,Heff=Heff,$
  print,'       Sidx=Sidx,Midx=Midx,Hidx=Hidx, Sscl=Sscl,Mscl=Mscl,Hscl=Hscl,$
  print,'       Stable=Stable,Mtable=Mtable,Htable=Htable, /post,level=level,$'
  print,'       nbin=nbins,outputf=outputf,BEHRdir=BEHRdir,verbose=verbose)'
  print,''
  print,'  an IDL wrapper to Bayesian Estimate of Hardness Ratios (BEHR)'
  print,'  returns a structure containing the relevant outputs neatly'
  print,'  summarized into fields for IDL consumption'
  print,''
  print,'  REFERENCES:'
  print,'  Park, T., Kashyap, V.L., Siemiginowska, A., van Dyk, D., Zezas, A.,
  print,'    Heinke, C., and Wargelin, B., 2006, ApJ, 652, 610'
  print,'    http://hea-www.harvard.edu/AstroStat/BEHR/'
  print,''
  if np ne 0 then message,ok,/informational
  print,''
  if vv ne 0 then begin
    if mp eq 0 then for i=0,n_elements(usage)-1 do print,usage[i] else $
  	for i=0,n_elements(help)-1 do print,help[i]
  endif
  return,-1L
endif

;	how many are input?
ndat=nss > nms > nhs

;	set up the input parameters
;	source counts
if nss eq 0 then xSsrc=0 else xSsrc=Ssrc
if nms eq 0 then xMsrc=0 else xMsrc=Msrc
if nhs eq 0 then xHsrc=0 else xHsrc=Hsrc
ss=lonarr(ndat)+xSsrc[0] & if nss gt 0 then ss[0L:nss-1L]=xSsrc[*]
ms=lonarr(ndat)+xMsrc[0] & if nms gt 0 then ms[0L:nms-1L]=xMsrc[*]
hs=lonarr(ndat)+xHsrc[0] & if nhs gt 0 then hs[0L:nhs-1L]=xHsrc[*]
;
;	background counts
nsb=n_elements(Sbkg) & if nsb eq 0 then xSbkg=0 else xSbkg=Sbkg
nmb=n_elements(Mbkg) & if nmb eq 0 then xMbkg=0 else xMbkg=Mbkg
nhb=n_elements(Hbkg) & if nhb eq 0 then xHbkg=0 else xHbkg=Hbkg
sb=lonarr(ndat)+xSbkg[0] & if nsb gt 0 then sb[0L:nsb-1L]=xSbkg[*]
mb=lonarr(ndat)+xMbkg[0] & if nmb gt 0 then mb[0L:nmb-1L]=xMbkg[*]
hb=lonarr(ndat)+xHbkg[0] & if nhb gt 0 then hb[0L:nhb-1L]=xHbkg[*]
;
;	background area scaling factors
nsa=n_elements(Sarea) & if nsa eq 0 then xSarea=1. else xSarea=Sarea
nma=n_elements(Marea) & if nma eq 0 then xMarea=1. else xMarea=Marea
nha=n_elements(Harea) & if nha eq 0 then xHarea=1. else xHarea=Harea
sa=dblarr(ndat)+xSarea[0] & if nsa gt 0 then sa[0L:nsa-1L]=xSarea[*]
ma=dblarr(ndat)+xMarea[0] & if nma gt 0 then ma[0L:nma-1L]=xMarea[*]
ha=dblarr(ndat)+xHarea[0] & if nha gt 0 then ha[0L:nha-1L]=xHarea[*]
;
;	effective area
nse=n_elements(Seff) & if nse eq 0 then xSeff=1.0 else xSeff=Seff
nme=n_elements(Meff) & if nme eq 0 then xMeff=1.0 else xMeff=Meff
nhe=n_elements(Heff) & if nhe eq 0 then xHeff=1.0 else xHeff=Heff
if nse eq 0 and nme eq 0 and nhe ne 0 then begin & xSeff=xHeff & xMeff=xHeff & endif
if nse eq 0 and nme ne 0 and nhe eq 0 then begin & xSeff=xMeff & xHeff=xMeff & endif
if nse ne 0 and nme eq 0 and nhe eq 0 then begin & xMeff=xSeff & xHeff=xSeff & endif
se=dblarr(ndat)+xSeff[0] & if nse gt 0 then se[0L:nse-1L]=xSeff[*]
me=dblarr(ndat)+xMeff[0] & if nme gt 0 then me[0L:nme-1L]=xMeff[*]
he=dblarr(ndat)+xHeff[0] & if nhe gt 0 then he[0L:nhe-1L]=xHeff[*]
;
;	index of prior
nsi=n_elements(Sidx) & if nsi eq 0 then xSidx=0.5 else xSidx=Sidx
nmi=n_elements(Midx) & if nmi eq 0 then xMidx=0.5 else xMidx=Midx
nhi=n_elements(Hidx) & if nhi eq 0 then xHidx=0.5 else xHidx=Hidx
if nsi eq 0 and nmi eq 0 and nhi ne 0 then begin & xSidx=xHidx & xMidx=xHidx & endif
if nsi eq 0 and nmi ne 0 and nhi eq 0 then begin & xSidx=xMidx & xHidx=xMidx & endif
if nsi ne 0 and nmi eq 0 and nhi eq 0 then begin & xMidx=xSidx & xHidx=xSidx & endif
si=fltarr(ndat)+xSidx[0] & if nsi gt 0 then si[0L:nsi-1L]=xSidx[*]
mi=fltarr(ndat)+xMidx[0] & if nmi gt 0 then mi[0L:nmi-1L]=xMidx[*]
hi=fltarr(ndat)+xHidx[0] & if nhi gt 0 then hi[0L:nhi-1L]=xHidx[*]
;
;	scale index of prior
nsc=n_elements(Sscl) & if nsc eq 0 then xSscl=0. else xSscl=Sscl
nmc=n_elements(Mscl) & if nmc eq 0 then xMscl=0. else xMscl=Mscl
nhc=n_elements(Hscl) & if nhc eq 0 then xHscl=0. else xHscl=Hscl
if nsc eq 0 and nmc eq 0 and nhc ne 0 then begin & xSscl=xHscl & xMscl=xHscl & endif
if nsc eq 0 and nmc ne 0 and nhc eq 0 then begin & xSscl=xMscl & xHscl=xMscl & endif
if nsc ne 0 and nmc eq 0 and nhc eq 0 then begin & xMscl=xSscl & xHscl=xSscl & endif
sc=fltarr(ndat)+xSscl[0] & if nsc gt 0 then sc[0L:nsc-1L]=xSscl[*]
mc=fltarr(ndat)+xMscl[0] & if nmc gt 0 then mc[0L:nmc-1L]=xMscl[*]
hc=fltarr(ndat)+xHscl[0] & if nhc gt 0 then hc[0L:nhc-1L]=xHscl[*]
;
nlev=n_elements(level) & if nlev eq 0 then xlevel=68. else xlevel=level
clev=fltarr(ndat)+xlevel[0] & if nlev gt 0 then clev[0L:nlev-1L]=xlevel[*]
;
msim=n_elements(nsim) & if msim eq 0 then nnsim=10000L else nnsim=nsim
numsim=lonarr(ndat)+nnsim[0] & if msim gt 0 then numsim[0L:msim-1L]=nnsim[*]
;
mburn=n_elements(nburnin) & if mburn eq 0 then nnburnin=(5000L < numsim/2) else nnburnin=nburnin
numburn=lonarr(ndat)+nnburnin[0] & if mburn gt 0 then numburn[0L:mburn-1L]=nnburnin[*]
;
nfmt=strtrim(fix(alog10(ndat)+1),2) & fmt='(i'+nfmt+'.'+nfmt+')'
nout=n_elements(outputf) & if nout eq 0 then outputf='none'
if ndat gt 1 then filroot=outputf[0]+'_'+string(lindgen(ndat)+1L,fmt) else $
	filroot=outputf[0]
if ndat eq 1 then if nout gt 1 then filroot[0L:nout-1L]=outputf[*]
;
mPost=n_elements(post) & if mPost eq 0 then iPost=0 else iPost=Post
outPost=bytarr(ndat)+iPost[0] & if mPost gt 0 then outPost[0L:mPost-1L]=iPost[*]
soutPost=strarr(ndat)+'false' & oo=where(outPost ne 0,moo) & if moo gt 0 then soutPost[oo]='true'
;
sinthpd=strarr(ndat)+'true'
;
;	table priors
if keyword_set(Stable) then begin
  stbl='tblprior_soft.txt'
  sfil=findfile(string(Stable),count=nsfil)
  if nsfil ne 0 then stbl=sfil[0]
endif
if keyword_set(Mtable) then begin
  Mtbl='tblprior_med.txt'
  Mfil=findfile(string(Mtable),count=nMfil)
  if nMfil ne 0 then Mtbl=Mfil[0]
endif
if keyword_set(Htable) then begin
  Htbl='tblprior_hard.txt'
  Hfil=findfile(string(Htable),count=nHfil)
  if nHfil ne 0 then Htbl=Hfil[0]
endif

;	define the output
dat=fltarr(ndat) & nat=dat+!values.F_NAN
Sstr=create_struct('mode',nat,'mean',nat,'median',nat,'lowerbound',nat,'upperbound',nat)
Mstr=Sstr & Hstr=Sstr					;S, M, H,
lSstr=Sstr & lMstr=Sstr & lHstr=Sstr			;logS, logM, logH,
SmMstr=Sstr & MmHstr=Sstr & SmHstr=Sstr			;SmM=S-M, MmH=M-H, SmH=S-H,
SpMstr=Sstr & MpHstr=Sstr & SpHstr=Sstr & Tstr=Sstr	;SpM=S+M, MpH=M+H, SpH=S+H, T=S+M+H
R1str=Sstr & R2str=Sstr & R3str=Sstr			;R1=S/M, R2=M/H, R3=S/H,
C1str=Sstr & C2str=Sstr & C3str=Sstr			;C1=C_SM, C2=C_MH, C3=C_SH,
HR1str=Sstr & HR2str=Sstr & HR3str=Sstr			;HR1=(S-M)/(S+M), HR2=(M-H)/(M+H), HR3=(S-H)/(S+H),
HRAstr=Sstr & HRBstr=Sstr & HRCstr=Sstr			;HRA=(S-M)/(S+M+H), HRB=(M-H)/(S+M+H), HRC=(S-H)/(S+M+H)}
;
postlamSidx=si & postlamMidx=mi & postlamHidx=hi
postlamSscl=dat & postlamMscl=dat & postlamHscl=dat

;	call BEHR for each data point
for i=0L,ndat-1L do begin
  if vv ne 0 and vv le 4 then print,strtrim(i,2),format='($,a10)'
  algo_method='gibbs'

  if ss[i]-sb[i]/sa[i] lt 15 or ms[i]-mb[i]/ma[i] lt 15 or hs[i]-hb[i]/ha[i] lt 15 then begin
    if ss[i]-sb[i]/sa[i] gt 100 or ms[i]-mb[i]/ma[i] gt 100 or hs[i]-hb[i]/ha[i] gt 100 then begin
      if numsim[i] lt 50000 then message,'WARNING: '+$
	'consider raising NSIM; too many counts in one band and too few in the other',/informational
    endif
  endif

  cmd1=BEHRexe+$
	' softsrc='+strtrim(ss[i],2)+$
	' hardsrc='+strtrim(ms[i],2)+$
	' softbkg='+strtrim(sb[i],2)+$
	' hardbkg='+strtrim(mb[i],2)+$
	' softarea='+strtrim(sa[i],2)+$
	' hardarea='+strtrim(ma[i],2)+$
	' softeff='+strtrim(se[i],2)+$
	' hardeff='+strtrim(me[i],2)+$
	' softidx='+strtrim(si[i],2)+$
	' hardidx='+strtrim(mi[i],2)+$
	' softscl='+strtrim(sc[i],2)+$
	' hardscl='+strtrim(mc[i],2)
  if keyword_set(Stable) then cmd1=cmd1+' softtbl='+stbl
  if keyword_set(Mtable) then cmd1=cmd1+' hardtbl='+mtbl
  cmd1=cmd1+$
	' post='+soutPost[i]+$
	' level='+strtrim(clev[i],2)+$
	' algo=gibbs'+$
	' details=true'+$
	' nsim='+strtrim(numsim[i],2)+$
	' nburnin='+strtrim(numburn[i],2)+$
	' HPD=true'+$
	' output='+filroot[i]+$
	' outputMC=true'
  if vv gt 4 then print,cmd1
  ;
  if numburn[i] ge numsim[i] then begin
    message,'NBURNIN must be less than NSIM; skipping this command:',$
    /informational
    print,cmd1
  endif else spawn,cmd1,stdout

  MCdrawsfil='BEHR_draws.txt'
  if strpos(strlowcase(filroot[i]),'none') lt 0 then MCdrawsfil=filroot[i]+'_draws.txt'
  var=fltarr(2,numsim[i]-numburn[i])
  openr,umc,MCdrawsfil,/get_lun
  readf,umc,var & Sdraw=reform(var[0,*]) & Mdraw=reform(var[1,*])
  close,umc & free_lun,umc

  nn=n_elements(stdout) & j=0L & k=-1
  while k lt 0 do begin
    if strpos(stdout[j],'Mode',0) ge 0 and $
       strpos(stdout[j],'Mean',0) ge 0 and $
       strpos(stdout[j],'Median',0) ge 0 and $
       strpos(stdout[j],'Lower Bound',0) ge 0 and $
       strpos(stdout[j],'Upper Bound',0) ge 0 then k=j
    if k gt 0 and vv ge 10 then print,stdout[j]
    j=j+1L
    if j eq nn then message,'BUG?  BEHR returns nothing'
  endwhile
  for j=k,nn-1L do begin
    if strpos(stdout[j],'p(lamS|data)',0) ge 0 then begin
      if vv ge 10 then print,stdout[j]
      i0=strpos(stdout[j],'(',strpos(stdout[j],'p(lamS|data)',0)+strlen('p(lamS|data)'))
      i1=strpos(stdout[j],')',i0+1)
      cc=strmid(stdout[j],i0+1,i1-i0-1)
      ccc=strsplit(cc,',',/extract)
      postlamSidx[i]=float(ccc[0])
      postlamSscl[i]=float(ccc[1])
    endif
    if strpos(stdout[j],'p(lamH|data)',0) ge 0 then begin
      if vv ge 10 then print,stdout[j]
      i0=strpos(stdout[j],'(',strpos(stdout[j],'p(lamH|data)',0)+strlen('p(lamH|data)'))
      i1=strpos(stdout[j],')',i0+1)
      cc=strmid(stdout[j],i0+1,i1-i0-1)
      ccc=strsplit(cc,',',/extract)
      postlamMidx[i]=float(ccc[0])
      postlamMscl[i]=float(ccc[1])
    endif
  endfor

  cmd2=BEHRexe+$
	' softsrc='+strtrim(ss[i],2)+$
	' hardsrc='+strtrim(hs[i],2)+$
	' softbkg='+strtrim(sb[i],2)+$
	' hardbkg='+strtrim(hb[i],2)+$
	' softarea='+strtrim(sa[i],2)+$
	' hardarea='+strtrim(ha[i],2)+$
	' softeff='+strtrim(se[i],2)+$
	' hardeff='+strtrim(he[i],2)+$
	' softidx='+strtrim(si[i],2)+$
	' hardidx='+strtrim(hi[i],2)+$
	' softscl='+strtrim(sc[i],2)+$
	' hardscl='+strtrim(hc[i],2)
  if keyword_set(Stable) then cmd2=cmd2+' softtbl='+mtbl
  if keyword_set(Htable) then cmd2=cmd2+' hardtbl='+htbl
  cmd2=cmd2+$
	' post='+soutPost[i]+$
	' level='+strtrim(clev[i],2)+$
	' algo=gibbs'+$
	' details=true'+$
	' nsim='+strtrim(numsim[i],2)+$
	' nburnin='+strtrim(numburn[i],2)+$
	' HPD=true'+$
	' output='+filroot[i]+$
	' outputMC=true'
  if vv gt 4 then print,cmd2
  ;
  if numburn[i] ge numsim[i] then begin
    message,'NBURNIN must be less than NSIM; skipping this command:',$
    /informational
    print,cmd2
  endif else spawn,cmd2,stdout

  MCdrawsfil='BEHR_draws.txt'
  if strpos(strlowcase(filroot[i]),'none') lt 0 then MCdrawsfil=filroot[i]+'_draws.txt'
  var=fltarr(2,numsim[i]-numburn[i])
  openr,umc,MCdrawsfil,/get_lun
  readf,umc,var & Sdraw=reform(var[0,*]) & Hdraw=reform(var[1,*])
  close,umc & free_lun,umc

  nn=n_elements(stdout) & j=0L & k=-1
  while k lt 0 do begin
    if strpos(stdout[j],'Mode',0) ge 0 and $
       strpos(stdout[j],'Mean',0) ge 0 and $
       strpos(stdout[j],'Median',0) ge 0 and $
       strpos(stdout[j],'Lower Bound',0) ge 0 and $
       strpos(stdout[j],'Upper Bound',0) ge 0 then k=j
    if k gt 0 and vv ge 10 then print,stdout[j]
    j=j+1L
    if j eq nn then message,'BUG?  BEHR returns nothing'
  endwhile
  for j=k,nn-1L do begin
    if strpos(stdout[j],'p(lamH|data)',0) ge 0 then begin
      if vv ge 10 then print,stdout[j]
      i0=strpos(stdout[j],'(',strpos(stdout[j],'p(lamH|data)',0)+strlen('p(lamH|data)'))
      i1=strpos(stdout[j],')',i0+1)
      cc=strmid(stdout[j],i0+1,i1-i0-1)
      ccc=strsplit(cc,',',/extract)
      postlamHidx[i]=float(ccc[0])
      postlamHscl[i]=float(ccc[1])
    endif
  endfor

  y=Sdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  Sstr.mode[i]=fmode & Sstr.mean[i]=mean(y) & Sstr.median[i]=median(y) & Sstr.lowerbound[i]=z[0] & Sstr.upperbound[i]=z[1]
  y=Mdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  Mstr.mode[i]=fmode & Mstr.mean[i]=mean(y) & Mstr.median[i]=median(y) & Mstr.lowerbound[i]=z[0] & Mstr.upperbound[i]=z[1]
  y=Hdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  Hstr.mode[i]=fmode & Hstr.mean[i]=mean(y) & Hstr.median[i]=median(y) & Hstr.lowerbound[i]=z[0] & Hstr.upperbound[i]=z[1]

  y=alog10(Sdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    lSstr.mode[i]=fmode & lSstr.mean[i]=mean(y[os]) & lSstr.median[i]=median(y[os]) & lSstr.lowerbound[i]=z[0] & lSstr.upperbound[i]=z[1]
  endif
  y=alog10(Mdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    lMstr.mode[i]=fmode & lMstr.mean[i]=mean(y[os]) & lMstr.median[i]=median(y[os]) & lMstr.lowerbound[i]=z[0] & lMstr.upperbound[i]=z[1]
  endif
  y=alog10(Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    lHstr.mode[i]=fmode & lHstr.mean[i]=mean(y[os]) & lHstr.median[i]=median(y[os]) & lHstr.lowerbound[i]=z[0] & lHstr.upperbound[i]=z[1]
  endif

  y=Sdraw-Mdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  SmMstr.mode[i]=fmode & SmMstr.mean[i]=mean(y) & SmMstr.median[i]=median(y) & SmMstr.lowerbound[i]=z[0] & SmMstr.upperbound[i]=z[1]
  y=Mdraw-Hdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  MmHstr.mode[i]=fmode & MmHstr.mean[i]=mean(y) & MmHstr.median[i]=median(y) & MmHstr.lowerbound[i]=z[0] & MmHstr.upperbound[i]=z[1]
  y=Sdraw-Hdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  SmHstr.mode[i]=fmode & SmHstr.mean[i]=mean(y) & SmHstr.median[i]=median(y) & SmHstr.lowerbound[i]=z[0] & SmHstr.upperbound[i]=z[1]

  y=Sdraw+Mdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  SpMstr.mode[i]=fmode & SpMstr.mean[i]=mean(y) & SpMstr.median[i]=median(y) & SpMstr.lowerbound[i]=z[0] & SpMstr.upperbound[i]=z[1]
  y=Mdraw+Hdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  MpHstr.mode[i]=fmode & MpHstr.mean[i]=mean(y) & MpHstr.median[i]=median(y) & MpHstr.lowerbound[i]=z[0] & MpHstr.upperbound[i]=z[1]
  y=Sdraw+Hdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  SpHstr.mode[i]=fmode & SpHstr.mean[i]=mean(y) & SpHstr.median[i]=median(y) & SpHstr.lowerbound[i]=z[0] & SpHstr.upperbound[i]=z[1]
  y=Sdraw+Mdraw+Hdraw & z=hipd_interval(y,/fsample,fmode=fmode,clev=clev[i])
  Tstr.mode[i]=fmode & Tstr.mean[i]=mean(y) & Tstr.median[i]=median(y) & Tstr.lowerbound[i]=z[0] & Tstr.upperbound[i]=z[1]

  y=float(Sdraw)/float(Mdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    R1str.mode[i]=fmode & R1str.mean[i]=mean(y[os]) & R1str.median[i]=median(y[os]) & R1str.lowerbound[i]=z[0] & R1str.upperbound[i]=z[1]
  endif
  y=float(Mdraw)/float(Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    R2str.mode[i]=fmode & R2str.mean[i]=mean(y[os]) & R2str.median[i]=median(y[os]) & R2str.lowerbound[i]=z[0] & R2str.upperbound[i]=z[1]
  endif
  y=float(Sdraw)/float(Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    R3str.mode[i]=fmode & R3str.mean[i]=mean(y[os]) & R3str.median[i]=median(y[os]) & R3str.lowerbound[i]=z[0] & R3str.upperbound[i]=z[1]
  endif

  y=alog10(Sdraw)-alog10(Mdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    C1str.mode[i]=fmode & C1str.mean[i]=mean(y[os]) & C1str.median[i]=median(y[os]) & C1str.lowerbound[i]=z[0] & C1str.upperbound[i]=z[1]
  endif
  y=alog10(Mdraw)-alog10(Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    C2str.mode[i]=fmode & C2str.mean[i]=mean(y[os]) & C2str.median[i]=median(y[os]) & C2str.lowerbound[i]=z[0] & C2str.upperbound[i]=z[1]
  endif
  y=alog10(Sdraw)-alog10(Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    C3str.mode[i]=fmode & C3str.mean[i]=mean(y[os]) & C3str.median[i]=median(y[os]) & C3str.lowerbound[i]=z[0] & C3str.upperbound[i]=z[1]
  endif

  y=(Sdraw-Mdraw)/(Sdraw+Mdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    HR1str.mode[i]=fmode & HR1str.mean[i]=mean(y[os]) & HR1str.median[i]=median(y[os]) & HR1str.lowerbound[i]=z[0] & HR1str.upperbound[i]=z[1]
  endif
  y=(Mdraw-Hdraw)/(Mdraw+Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    HR2str.mode[i]=fmode & HR2str.mean[i]=mean(y[os]) & HR2str.median[i]=median(y[os]) & HR2str.lowerbound[i]=z[0] & HR2str.upperbound[i]=z[1]
  endif
  y=(Sdraw-Hdraw)/(Sdraw+Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    HR3str.mode[i]=fmode & HR3str.mean[i]=mean(y[os]) & HR3str.median[i]=median(y[os]) & HR3str.lowerbound[i]=z[0] & HR3str.upperbound[i]=z[1]
  endif

  y=(Sdraw-Mdraw)/(Sdraw+Mdraw+Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    HRAstr.mode[i]=fmode & HRAstr.mean[i]=mean(y[os]) & HRAstr.median[i]=median(y[os]) & HRAstr.lowerbound[i]=z[0] & HRAstr.upperbound[i]=z[1]
  endif
  y=(Mdraw-Hdraw)/(Sdraw+Mdraw+Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    HRBstr.mode[i]=fmode & HRBstr.mean[i]=mean(y[os]) & HRBstr.median[i]=median(y[os]) & HRBstr.lowerbound[i]=z[0] & HRBstr.upperbound[i]=z[1]
  endif
  y=(Sdraw-Hdraw)/(Sdraw+Mdraw+Hdraw) & os=where(finite(y) ne 0,mos)
  if mos gt 2 then begin
    z=hipd_interval(y[os],/fsample,fmode=fmode,clev=clev[i])
    HRCstr.mode[i]=fmode & HRCstr.mean[i]=mean(y[os]) & HRCstr.median[i]=median(y[os]) & HRCstr.lowerbound[i]=z[0] & HRCstr.upperbound[i]=z[1]
  endif

endfor

;	define the output
BEHR3=create_struct('softsrc',ss,'mediumsrc',ms,'hardsrc',hs,$
	'softbkg',sb,'mediumbkg',mb,'hardbkg',hb,$
	'softarea',sa,'mediumarea',ma,'hardarea',ha,$
	'softeff',se,'mediumeff',me,'hardeff',he,$
	'softidx',si,'mediumidx',mi,'hardidx',hi,$
	'softscl',sc,'mediumscl',mc,'hardscl',hc,$
	'level',clev,'algo','gibbs','nsim',numsim,'nburnin',numburn,$
	'nbins',0,'HPD',sinthpd,'filroot',filroot,$
	'S',Sstr,'M',Mstr,'H',Hstr,'lS',lSstr,'lM',lMstr,'lH',lHstr,$
	'SmM',SmMstr,'MmH',MmHstr,'SmH',SmHstr,'SpM',SpMstr,'MpH',MpHstr,'SpH',SpHstr,'T',Tstr,$
	'R1',R1str,'R2',R2str,'R3',R3str,'C1',C1str,'C2',C2str,'C3',C3str,$
	'HR1',HR1str,'HR2',HR2str,'HR3',HR3str,'HRA',HRAstr,'HRB',HRBstr,'HRC',HRCstr,$
	'postlamS_idx',postlamSidx,'postlamS_scl',postlamSscl,$
	'postlamM_idx',postlamMidx,'postlamM_scl',postlamMscl,$
	'postlamH_idx',postlamHidx,'postlamH_scl',postlamHscl,$
	'help', 'S, M, H ; lS=log10(S), lM=log10(M), lH=log10(H) ; '+$
		'SpM=S+M, MpH=M+H, SpH=S+H, T=S+M+H ; '+$
		'SmM=S-M, MmH=M-H, SmH=S-H ; '+$
		'R1=S/M, R2=M/H, R3=S/H ; C1=C_SM, C2=C_MH, C3=C_SH ; '+$
		'HR1=(S-M)/(S+M), HR2=(M-H)/(M+H), HR3=(S-H)/(S+H) ; '+$
		'HRA=(S-M)/(S+M+H), HRB=(M-H)/(S+M+H), HRC=(S-H)/(S+M+H)')

return,BEHR3
end

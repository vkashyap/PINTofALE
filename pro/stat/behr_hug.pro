function behr_hug,softsrc,hardsrc,softbkg,hardbkg,softarea,hardarea,$
	softeff=softeff,hardeff=hardeff,softidx=softidx,hardidx=hardidx,$
	softscl=softscl,hardscl=hardscl,softtbl=softtbl,hardtbl=hardtbl,$
	post=post,level=level,algo=algo,details=details,nsim=nsim,nburnin=nburnin,$
	hpd=hpd,nbins=nbins,outputf=outputf,outputR=outputR,outputHR=outputHR,$
	outputC=outputC,outputMC=outputMC,outputPr=outputPr,BEHRdir=BEHRdir,$
	verbose=verbose
;+
;function	BEHR_hug
;	IDL wrapper to Bayesian Estimate of Hardness Ratios (BEHR)
;	returns a structure containing the relevant outputs neatly
;	summarized into fields for IDL consumption
;
;	references:
;
;	"Bayesian Estimation of Hardness Ratios: Modeling and Computations",
;	  Park, T., Kashyap, V.L., Siemiginowska, A., van Dyk, D., Zezas, A.,
;	  Heinke, C., and Wargelin, B., 2006, ApJ, 652, 610
;
;	"BEHR: Bayesian Estimation of Hardness Ratios", Park, T., Kashyap, V.,
;	  Zezas, A., Siemiginowska, A., van Dyk, D., Connors, A., and Heinke, C.,
;	  2005, at Six Years of Science with Chandra Symposium, Nov 2-5, 2005, #5.6
;
;	"Computing Hardness Ratios with Poissonian Errors",
;	  van Dyk, D.A., Park, T., Kashyap, V.L., \& Zezas, A.,
;	  2004, AAS-HEAD #8, 16.27
;	  http://www.aas.org/publications/baas/v36n3/head2004/137.htm
;
;syntax
;	behr=behr_hug(softsrc,hardsrc,softbkg,hardbkg,softarea,hardarea,$
;	softeff=softeff,hardeff=hardeff,softidx=softidx,hardidx=hardidx,$
;	softscl=softscl,hardscl=hardscl,softtbl=softtbl,hardtbl=hardtbl,/post,$
;	level=level,algo=algo,details=details,nsim=nsim,nburnin=nburnin,$
;	/hpd,nbins=nbins,outputf=outputf,outputR=outputR,outputHR=outputHR,$
;	outputC=outputC,outputMC=outputMC,outputPr=outputPr,BEHRdir=BEHRdir,$
;	verbose=verbose)
;
;parameters
;	softsrc	[INPUT; required] counts in source region in soft (S) band
;	hardsrc	[INPUT; required] counts in source region in hard (H) band
;	softbkg	[INPUT; required] counts in background region in S band
;	hardbkg	[INPUT; required] counts in background region in H band
;	softarea[INPUT; required] background scaling factor in S band
;	hardarea[INPUT; required] background scaling factor in H band
;		* (background region area)/(source region area)
;		* can also include differences in exposure time into
;		  the ratio, in the same manner as geometric area
;
;keywords
;	softeff	[INPUT] effective area in soft (S) band
;		* if not set, then set equal to HARDEFF if given, or
;		  1 otherwise
;	hardeff	[INPUT] effective area in hard (H) band
;		* if not set, then set equal to SOFTEFF if given, or
;		  1 otherwise
;		* can also be the effective area relative to some
;		  special point on the detector (e.g., aimpoint)
;		  or even some specific detector (e.g., ACIS-I v/s ACIS-S)
;	softidx	[INPUT] index of prior on S (range = 0+)
;		* if not set, then set equal to HARDIDX if given, or
;		  0.5 otherwise
;	hardidx	[INPUT] index of prior on H (range = 0+)
;		* if not set, then set equal to SOFTIDX if given, or
;		  0.5 otherwise
;		* similar to AGAMMA of PPD_SRC()
;	softscl	[INPUT] scaling index of prior on S
;		* if not set, then set to HARDSCL if given, or
;		  0 otherwise
;	hardscl	[INPUT] scaling index of prior on H
;		* if not set, then set to SOFTSCL if given, or
;		  0 otherwise
;		* similar to BGAMMA of PPD_SRC()
;	softtbl	[INPUT] filename containing a tabulated prior for S
;	hardtbl	[INPUT] filename containing a tabulated prior for H
;		* the table prior must be an ascii file with the following format:
;		  line 1: number of entries, say NLIN
;		  line 2: labels for the columns, ignored
;		  lines 3..NLIN+2: two whitespace separated columns of numbers,
;	            with each row containing the source intensity and the posterior
;	            density, in that order
;		* the default filenames are "./tblprior_{soft|hard}.txt"
;		* the default filenames are used iff SOFTTBL and HARDTBL are set
;		  but do not appear to exist
;		* WARNING: if regex is used in the filename specification, only the
;		  first file from the list will be used.  furthermore, if specified,
;		  the table priors are applied to _all_ SOFTSRC and HARDSRC
;	post	[INPUT] if set, suggests the values of (SOFTIDX,SOFTSCL)
;		and (HARDIDX,HARDSCL) going forward, i.e., what you should
;		set the priors to in your next calculation for the same
;		source
;	level	[INPUT] percentage confidence level at which to report error
;		(default = 68)
;	details	[INPUT] compute various ratios (true/false)?
;		(default = true)
;	algo	[INPUT] calculation method, GIBBS (default) or QUAD
;		* if set to AUTO or AUTO=N, then uses GIBBS for any
;		  case where {SOFT|HARD}SRC > N, and QUAD below that
;		  unless one of {SOFT|HARD}SRC > 99, in which case
;		  GIBBS is set automatically
;		  - if N is not set, assumed to be 15
;		  - N can be set to any integer less than 100
;	nsim	[INPUT] number of draws if algo=gibbs (default=10000)
;	nburnin	[INPUT] number of burn-in draws if algo=gibbs
;		(default=5000 or NSIM/2, whichever is smaller)
;	nbins	[INPUT] number of bins in integration if algo=quad (default=500)
;	hpd	[INPUT] if set, computes the highest posterior density interval,
;		which also is the smallest interval that includes the mode.
;		* this is set by default
;		* only has an effect if ALGO=QUAD
;	outputF	[INPUT] root of filename in which to place output
;		(default = 'none')
;		* output will be placed in the file, OUTPUTF.txt
;		* unless OUTPUTF='none', in which case nothing is written out,
;		    unless one of OUTPUT{R|C|HR|MC|PR} are set, in which case
;		    the corresponding output is put in file BEHR_{R|C|HR|MC|PR}.txt
;	outputR	[INPUT] if set, writes output for R in OUTPUTF_R.txt
;	outputC	[INPUT] if set, writes output for C in OUTPUTF_C.txt
;	outputHR[INPUT] if set, writes output for HR in OUTPUTF_HR.txt
;	outputMC[INPUT] if set, writes Monte Carlo draws for lamS and lamH
;		to OUTPUTF_draws.txt when algo=gibbs
;	outputPr[INPUT] if set, writes the probability distributions for
;		R, HR, C, lamS, and lamH to OUTPUTF_prob.txt when algo=quad
;	BEHRdir	[INPUT] full path to directory where BEHR executable resides
;		(default = '/fubar/kashyap/AstroStat/BEHR')
;	verbose	[INPUT] controls chatter
;
;requirements
;	BEHR executable must be installed in BEHRDIR
;	BEHR should be executable under the shell via SPAWN
;	BEHR output assumed to be compatible with 12-19-2005 version
;
;side-effects
;	potentially creates numerous ascii files in $cwd
;
;history
;	vinay kashyap (SepMMV)
;	bugfix: output file format was overflowing; force NBURNIN < NSIM;
;	  outputF for 1st iteration; added more references (VK; NovMMV)
;	added keywords SOFTSCL,HARDSCL,POST,HPD; changed DETAILS default
;	  to true, OUTPUTF default to none; added extra fields to output
;	  structure; added ALGO option AUTO (VK; DecMMV)
;	altered behavior of ALGO=AUTO to use GIBBS if counts in one of the
;	  bands exceeds 100 (VK; OctMMVI)
;	added keywords SOFTTBL and HARDTBL (VK; FebMMVIII)
;	updated default BEHRdir (VK; OctMMXV)
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
np=n_params() & if np lt 6 then ok='Insufficient parameters'
spawn,BEHRexe,usage
if strpos(strlowcase(usage[0]),'command not found',0) ge 0 then $
	ok='BEHR: command not found'
cmd=BEHRexe
nss=n_elements(softsrc) & if nss eq 0 then xsoftsrc=0 else xsoftsrc=softsrc
nhs=n_elements(hardsrc) & if nhs eq 0 then xhardsrc=0 else xhardsrc=hardsrc
nsb=n_elements(softbkg) & if nsb eq 0 then xsoftbkg=0 else xsoftbkg=softbkg
nhb=n_elements(hardbkg) & if nhb eq 0 then xhardbkg=0 else xhardbkg=hardbkg
nsa=n_elements(softarea) & if nsa eq 0 then xsoftarea=0 else xsoftarea=softarea
nha=n_elements(hardarea) & if nha eq 0 then xhardarea=0 else xhardarea=hardarea
mp=0
if nss gt 0 then mp=mp+1 & if nhs gt 0 then mp=mp+1
if nsb gt 0 then mp=mp+1 & if nhb gt 0 then mp=mp+1
if nsa gt 0 then mp=mp+1 & if nha gt 0 then mp=mp+1
if mp lt 6 then begin
  if nss ne 0 then cmd=cmd+' softsrc='+strtrim(xsoftsrc[0],2)
  if nhs ne 0 then cmd=cmd+' hardsrc='+strtrim(xhardsrc[0],2)
  if nsb ne 0 then cmd=cmd+' softbkg='+strtrim(xsoftbkg[0],2)
  if nhb ne 0 then cmd=cmd+' hardbkg='+strtrim(xhardbkg[0],2)
  if nsa ne 0 then cmd=cmd+' softarea='+strtrim(xsoftarea[0],2)
  if nha ne 0 then cmd=cmd+' hardarea='+strtrim(xhardarea[0],2)
  spawn,cmd,help
endif
;
if ok ne 'ok' then begin
  print,'Usage: behr=behr_hug(softsrc,hardsrc,softbkg,hardbkg,softarea,hardarea,$'
  print,'       softeff=softeff,hardeff=hardeff,softidx=softidx,hardidx=hardidx,$'
  print,'       softscl=softscl,hardscl=hardscl,softtbl=softtbl,hardtbl=hardtbl,$
  print,'       /post,level=level,algo=algo,details=details,nsim=nsim,nburnin=nburnin,$'
  print,'       /hpd,nbin=nbins,outputf=outputf,outputR=outputR,outputHR=outputHR,$'
  print,'       outputC=outputC,outputMC=outputMC,outputPr=outputPr,BEHRdir=BEHRdir,$'
  print,'       verbose=verbose)'
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
nse=n_elements(softeff) & if nse eq 0 then xsofteff=1.0 else xsofteff=softeff
nhe=n_elements(hardeff) & if nhe eq 0 then xhardeff=1.0 else xhardeff=hardeff
nsi=n_elements(softidx) & if nsi eq 0 then xsoftidx=0.5 else xsoftidx=softidx
nhi=n_elements(hardidx) & if nhi eq 0 then xhardidx=0.5 else xhardidx=hardidx
nsc=n_elements(softscl) & if nsi eq 0 then xsoftscl=0. else xsoftscl=softscl
nhc=n_elements(hardscl) & if nhi eq 0 then xhardscl=0. else xhardscl=hardscl
ndat=nss > nhs > nsb > nhb > nsa > nha > nse > nhe > nsi > nhi > nsc > nhc

;	set up the input parameters
ss=lonarr(ndat)+xsoftsrc[0] & if nss gt 0 then ss[0L:nss-1L]=xsoftsrc[*]
hs=lonarr(ndat)+xhardsrc[0] & if nhs gt 0 then hs[0L:nhs-1L]=xhardsrc[*]
sb=lonarr(ndat)+xsoftbkg[0] & if nsb gt 0 then sb[0L:nsb-1L]=xsoftbkg[*]
hb=lonarr(ndat)+xhardbkg[0] & if nhb gt 0 then hb[0L:nhb-1L]=xhardbkg[*]
sa=dblarr(ndat)+xsoftarea[0] & if nsa gt 0 then sa[0L:nsa-1L]=xsoftarea[*]
ha=dblarr(ndat)+xhardarea[0] & if nha gt 0 then ha[0L:nha-1L]=xhardarea[*]
;	set up the input keywords
se=dblarr(ndat)+xsofteff[0] & if nse gt 0 then se[0L:nse-1L]=xsofteff[*]
he=dblarr(ndat)+xhardeff[0] & if nhe gt 0 then he[0L:nhe-1L]=xhardeff[*]
si=fltarr(ndat)+xsoftidx[0] & if nsi gt 0 then si[0L:nsi-1L]=xsoftidx[*]
hi=fltarr(ndat)+xhardidx[0] & if nhi gt 0 then hi[0L:nhi-1L]=xhardidx[*]
sc=fltarr(ndat)+xsoftscl[0] & if nsc gt 0 then sc[0L:nsc-1L]=xsoftscl[*]
hc=fltarr(ndat)+xhardscl[0] & if nhc gt 0 then hc[0L:nhc-1L]=xhardscl[*]
;
nlev=n_elements(level) & if nlev eq 0 then xlevel=68. else xlevel=level
clev=fltarr(ndat)+xlevel[0] & if nlev gt 0 then clev[0L:nlev-1L]=xlevel[*]
;
nalg=n_elements(algo) & if nalg eq 0 then xalgo='gibbs' else xalgo=algo
alg=strarr(ndat)+xalgo[0] & if nalg gt 0 then alg[0L:nalg-1L]=xalgo[*]
alg=strlowcase(alg)
;
msim=n_elements(nsim) & if msim eq 0 then nnsim=10000L else nnsim=nsim
numsim=lonarr(ndat)+nnsim[0] & if msim gt 0 then numsim[0L:msim-1L]=nnsim[*]
;
mburn=n_elements(nburnin) & if mburn eq 0 then nnburnin=(5000L < numsim/2) else nnburnin=nburnin
numburn=lonarr(ndat)+nnburnin[0] & if mburn gt 0 then numburn[0L:mburn-1L]=nnburnin[*]
;
mbin=n_elements(nbins) & if mbin eq 0 then nnbins=500L else nnbins=nbin
numbin=lonarr(ndat)+nnbins[0] & if mbin gt 0 then numbin[0L:mbin-1L]=nnbins[*]
;
nfmt=strtrim(fix(alog10(ndat)+1),2) & fmt='(i'+nfmt+'.'+nfmt+')'
nout=n_elements(outputf) & if nout eq 0 then outputf='none'
if ndat gt 1 then filroot=outputf[0]+'_'+string(lindgen(ndat)+1L,fmt) else $
	filroot=outputf[0]
if ndat eq 1 then if nout gt 1 then filroot[0L:nout-1L]=outputf[*]
;
mR=n_elements(outputR) & if mR eq 0 then xoutputR=0 else xoutputR=outputR
outR=bytarr(ndat)+xoutputR[0] & if mR gt 0 then outR[0L:mR-1L]=xoutputR[*]
soutR=strarr(ndat)+'false' & oo=where(outR ne 0,moo) & if moo gt 0 then soutR[oo]='true'
;
mC=n_elements(outputC) & if mC eq 0 then xoutputC=0 else xoutputC=outputC
outC=bytarr(ndat)+xoutputC[0] & if mC gt 0 then outC[0L:mC-1L]=xoutputC[*]
soutC=strarr(ndat)+'false' & oo=where(outC ne 0,moo) & if moo gt 0 then soutC[oo]='true'
;
mHR=n_elements(outputHR) & if mHR eq 0 then xoutputHR=0 else xoutputHR=outputHR
outHR=bytarr(ndat)+xoutputHR[0] & if mHR gt 0 then outHR[0L:mHR-1L]=xoutputHR[*]
soutHR=strarr(ndat)+'false' & oo=where(outHR ne 0,moo) & if moo gt 0 then soutHR[oo]='true'
;
mMC=n_elements(outputMC) & if mMC eq 0 then xoutputMC=0 else xoutputMC=outputMC
outMC=bytarr(ndat)+xoutputMC[0] & if mMC gt 0 then outMC[0L:mMC-1L]=xoutputMC[*]
soutMC=strarr(ndat)+'false' & oo=where(outMC ne 0,moo) & if moo gt 0 then soutMC[oo]='true'
;
mPr=n_elements(outputPr) & if mPr eq 0 then xoutputPr=0 else xoutputPr=outputPr
outPr=bytarr(ndat)+xoutputPr[0] & if mPr gt 0 then outPr[0L:mPr-1L]=xoutputPr[*]
soutPr=strarr(ndat)+'false' & oo=where(outPr ne 0,moo) & if moo gt 0 then soutPr[oo]='true'
;
mPost=n_elements(post) & if mPost eq 0 then iPost=0 else iPost=Post
outPost=bytarr(ndat)+iPost[0] & if mPost gt 0 then outPost[0L:mPost-1L]=iPost[*]
soutPost=strarr(ndat)+'false' & oo=where(outPost ne 0,moo) & if moo gt 0 then soutPost[oo]='true'
;
mdet=n_elements(details) & if mdet eq 0 then idetails=1 else idetails=details
dotell=bytarr(ndat)+idetails[0] & if mdet gt 0 then dotell[0L:mdet-1L]=idetails[*]
sdotell=strarr(ndat)+'true' & oo=where(dotell eq 0,moo) & if moo gt 0 then sdotell[oo]='false'
;
mhpd=n_elements(hpd) & if mhpd eq 0 then ihpd=1 else ihpd=hpd
inthpd=bytarr(ndat)+ihpd[0] & if mhpd gt 0 then inthpd[0L:mhpd-1L]=ihpd[*]
sinthpd=strarr(ndat)+'true' & oo=where(inthpd eq 0,moo) & if moo gt 0 then sinthpd[oo]='false'
;
if keyword_set(softtbl) then begin
  stbl='tblprior_soft.txt'
  sfil=findfile(string(softtbl),count=nsfil)
  if nsfil ne 0 then stbl=sfil[0]
endif
if keyword_set(hardtbl) then begin
  htbl='tblprior_hard.txt'
  hfil=findfile(string(hardtbl),count=nhfil)
  if nhfil ne 0 then htbl=hfil[0]
endif

;	define the output
dat=fltarr(ndat) & nat=dat+!values.F_NAN & warn=strarr(ndat)
Sstr=create_struct('mode',nat,'mean',nat,'median',nat,'lowerbound',nat,'upperbound',nat)
Hstr=Sstr & HmSstr=Sstr & HpSstr=Sstr & lSstr=Sstr & lHstr=Hstr
Rstr=create_struct('mode',dat,'mean',dat,'median',dat,'lowerbound',dat,'upperbound',dat)
Cstr=Rstr & HRstr=Rstr
postlamSidx=si & postlamHidx=hi & postlamSscl=dat & postlamHscl=dat

;	call BEHR for each data point
for i=0L,ndat-1L do begin
  if vv ne 0 and vv le 4 then print,strtrim(i,2),format='($,a10)'
  algo_method='gibbs'
  if strpos(strlowcase(alg[i]),'gibbs') lt 0 then begin
    if strpos(strlowcase(alg[i]),'auto') lt 0 then algo_method='quad' else begin
      cc=strsplit(alg[i],'=',/extract) & ncc=n_elements(cc)
      if ncc eq 1 then ctthresh=15L else ctthresh=long(cc[1]) < 99
      if ss[i]-sb[i]/sa[i] lt ctthresh or hs[i]-hb[i]/ha[i] lt ctthresh then begin
	algo_method='quad'
        if ss[i]-sb[i]/sa[i] gt 100 or hs[i]-hb[i]/ha[i] gt 100 then begin
	  algo_method='gibbs'
	  if numsim[i] lt 50000 then message,'WARNING: '+$
	    'Choosing GIBBS when there are too many counts in one band and too few in the other; consider raising NSIM',/informational
        endif
      endif
      if vv gt 0 then print,alg[i]+'::',ss[i],hs[i],ctthresh,' :: Choosing algo='+algo_method
    endelse
  endif
  cmd=BEHRexe+$
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
  if keyword_set(softtbl) then cmd=cmd+' softtbl='+stbl
  if keyword_set(hardtbl) then cmd=cmd+' hardtbl='+htbl
  cmd=cmd+$
	' post='+soutPost[i]+$
	' level='+strtrim(clev[i],2)+$
	' algo='+algo_method+$
	' details='+sdotell[i]+$
	' nsim='+strtrim(numsim[i],2)+$
	' nburnin='+strtrim(numburn[i],2)+$
	' HPD='+sinthpd[i]+$
	' nbins='+strtrim(numbin[i],2)+$
	' output='+filroot[i]+$
	' outputR='+soutR[i]+$
	' outputC='+soutC[i]+$
	' outputHR='+soutHR[i]+$
	' outputMC='+soutMC[i]+$
	' outputPr='+soutPr[i]
  if vv gt 4 then print,cmd
  ;
  if numburn[i] ge numsim[i] then begin
    message,'NBURNIN must be less than NSIM; skipping this command:',$
    /informational
    print,cmd
  endif else spawn,cmd,stdout
  ;
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
    if strpos(stdout[j],'S/H',0) ge 0 and $
       strpos(stdout[j],'log10',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      Rstr.mode[i]=float(cc[1])
      Rstr.mean[i]=float(cc[2])
      Rstr.median[i]=float(cc[3])
      Rstr.lowerbound[i]=float(cc[4])
      Rstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' R '
    endif
    if strpos(stdout[j],'H-S',0) ge 0 and $
       strpos(stdout[j],'H+S',0) ge 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      HRstr.mode[i]=float(cc[1])
      HRstr.mean[i]=float(cc[2])
      HRstr.median[i]=float(cc[3])
      HRstr.lowerbound[i]=float(cc[4])
      HRstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' HR '
    endif
    if strpos(stdout[j],'log10',0) ge 0 and $
       strpos(stdout[j],'S/H',0) ge 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      Cstr.mode[i]=float(cc[1])
      Cstr.mean[i]=float(cc[2])
      Cstr.median[i]=float(cc[3])
      Cstr.lowerbound[i]=float(cc[4])
      Cstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' C '
    endif
    if strpos(stdout[j],'S',0) ge 0 and $
       strpos(stdout[j],'H',0) lt 0 and strpos(stdout[j],'log10',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      Sstr.mode[i]=float(cc[1])
      Sstr.mean[i]=float(cc[2])
      Sstr.median[i]=float(cc[3])
      Sstr.lowerbound[i]=float(cc[4])
      Sstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' S '
    endif
    if strpos(stdout[j],'H',0) ge 0 and $
       strpos(stdout[j],'S',0) lt 0 and strpos(stdout[j],'log10',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      Hstr.mode[i]=float(cc[1])
      Hstr.mean[i]=float(cc[2])
      Hstr.median[i]=float(cc[3])
      Hstr.lowerbound[i]=float(cc[4])
      Hstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' H '
    endif
    if strpos(stdout[j],'H-S',0) ge 0 and $
       strpos(stdout[j],'H+S',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      HmSstr.mode[i]=float(cc[1])
      HmSstr.mean[i]=float(cc[2])
      HmSstr.median[i]=float(cc[3])
      HmSstr.lowerbound[i]=float(cc[4])
      HmSstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' H-S '
    endif
    if strpos(stdout[j],'H+S',0) ge 0 and $
       strpos(stdout[j],'H-S',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      HpSstr.mode[i]=float(cc[1])
      HpSstr.mean[i]=float(cc[2])
      HpSstr.median[i]=float(cc[3])
      HpSstr.lowerbound[i]=float(cc[4])
      HpSstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' H+S '
    endif
    if strpos(stdout[j],'S',0) ge 0 and strpos(stdout[j],'log10',0) ge 0 and $
       strpos(stdout[j],'H',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      lSstr.mode[i]=float(cc[1])
      lSstr.mean[i]=float(cc[2])
      lSstr.median[i]=float(cc[3])
      lSstr.lowerbound[i]=float(cc[4])
      lSstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' logS '
    endif
    if strpos(stdout[j],'H',0) ge 0 and strpos(stdout[j],'log10',0) ge 0 and $
       strpos(stdout[j],'S',0) lt 0 then begin
      if vv ge 10 then print,stdout[j]
      cc=strsplit(stdout[j],/extract)
      lHstr.mode[i]=float(cc[1])
      lHstr.mean[i]=float(cc[2])
      lHstr.median[i]=float(cc[3])
      lHstr.lowerbound[i]=float(cc[4])
      lHstr.upperbound[i]=float(cc[5])
      if strpos(stdout[j],'*',0) ge 0 then warn[i]=warn[i]+' logH '
    endif
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
      postlamHidx[i]=float(ccc[0])
      postlamHscl[i]=float(ccc[1])
    endif
  endfor
endfor

;	define the output
BEHR=create_struct('softsrc',ss,'hardsrc',hs,'softbkg',sb,'hardbkg',hb,$
	'softarea',sa,'hardarea',ha,'softeff',se,'hardeff',he,$
	'softidx',si,'hardidx',hi,'softscl',sc,'hardscl',hc,'level',clev,$
	'algo',alg,'nsim',numsim,'nburnin',numburn,'nbins',numbin,'HPD',sinthpd,$
	'filroot',filroot,'R',Rstr,'C',Cstr,'HR',HRstr,$
	'S',Sstr,'H',Hstr,'HmS',HmSstr,'HpS',HpSstr,'lS',lSstr,'lH',lHstr,$
	'warning',warn,$
	'postlamS_idx',postlamSidx,'postlamS_scl',postlamSscl,$
	'postlamH_idx',postlamHidx,'postlamH_scl',postlamHscl)

return,BEHR
end

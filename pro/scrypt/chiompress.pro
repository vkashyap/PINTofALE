;+
;chiompress.pro
;	a script to compress the CHIANTI databases into a smaller dataset
;	to make it easy to transport and easy to read.  The max emissivity
;	of each line is tested against the max emissivities of nearby lines
;	and the line is retained iff it is above the locally defined threshold.
;
;usage
;	chifulldir='/data/fubar/SCAR/emissivity/chi713'	;full database dir
;	outdir='/data/fubar/SCAR/emissivity/chiompress713'	;compressed dir
;	dynrng=1e2	;dynamic range
;	hwdth=1.0	;half-width of window in AA
;	.run chiompress
;
;side effects
;	creates files in OUTDIR
;	creates a postscript file in $CWD
;
;vinay kashyap (2010jul)
;	bug fix: density sensitive lines were being ignored in the
;	write (VK; 2010sep)
;-

;	initialize
;chifulldir='/data/fumtu/kashyap/chiompress/chi521'
;ioneqf='ioneq/mazzotta_etal.ioneq'
if not keyword_set(chifulldir) then chifulldir='/data/fubar/SCAR/emissivity/chi713
if not keyword_set(outdir) then outdir='/data/fubar/SCAR/emissivity/chiompress713'
outroot='chiom'
if not keyword_set(hwdth) then hwdth=1.0	;moving window half-width over which to check for useful lines
if not keyword_set(dynrng) then dynrng=1e2	;keep all lines down to 1/dynrng of the max line in moving window
peasecolr & loadct,3 & peasecolr

;	first figure out how many elements there are
afils=findfile(chifulldir+'D'+'/*_wvl',count=natoms)
atoms=strarr(natoms)
for i=0,natoms-1 do begin
  cc=afils[i] & l1=strlen(chifulldir+'D') & i1=strpos(cc,'_',l1) & atoms[i]=strmid(cc,l1+1,i1-l1-1)
endfor
help,atoms

;	figure out how many densities there are
efils=findfile(chifulldir+'D'+'/'+atoms[0]+'_??.?',count=nedens)
edens=fltarr(nedens)
for i=0,nedens-1 do begin & $
  cc=efils[i] & l1=strlen(chifulldir+'D') & i1=strpos(cc,'_',l1) & edens[i]=float(strmid(cc,i1+1)) & $
endfor
n_e=10.^(edens)
help,edens

set_plot,'ps' & device,file='chiompress.ps',/landscape,/color
;	read in all the emissivities for all lines at all the densities, one element at a time
for i=0L,natoms-1L do begin

  ;	read it in
  for j=0L,nedens-1L do begin
    emis=rd_line(atoms[i],n_e=n_e[j],wrange=[-1,1e5],dbdir=chifulldir,fstr=fstr,verbose=verbose)
    if j eq 0 then begin
      nw=n_elements(fstr.WVL) & nT=n_elements(fstr.LOGT)
      allemis=dblarr(nT,nw,nedens)
      allemis[*,*,0]=emis
      allwvl=abs(fstr.wvl)
    endif else allemis[*,*,j]=emis
  endfor

  ;	now figure out how many lines to keep
  print,nw
  maxemis=dblarr(nw) & for j=0L,nw-1L do maxemis[j]=max(allemis[*,j,*])
  if not keyword_set(hwdth) then hwdth=0.5
  wmin=min(allwvl,max=wmax) & iw=lonarr(nw)
  for j=0L,nw-1L do begin
    w0=allwvl[j] & w0m=w0-hwdth & w0p=w0+hwdth & ow=where(allwvl ge w0m and allwvl le w0p,mow)
    maxmaxemis=max(maxemis[ow])
    ;ok=where(maxemis ge maxmaxemis/dynrng,mok)
    if maxemis[j] ge maxmaxemis/dynrng then iw[j]=1
    if j eq 1000*long(j/1000L) then kilroy
  endfor

  ow=where(iw ne 0,mw)
  plot,allwvl,maxemis,psym=1,/xl,/yl,title=atoms[i]
  oplot,allwvl[ow],maxemis[ow],psym=1,col=2

  ;	now read back in the emissivities for all densities
  for j=0L,nedens-1L do begin
    emis=rd_line(atoms[i],n_e=n_e[j],wrange=[-1,1e5],dbdir=chifulldir,fstr=fstr,verbose=verbose)
    linstr=cat_ln(fstr,pick=ow)
    wrt_ln_generic,linstr,outdir+'D',n_e[j],verbose=10
  endfor
  print,'compression factor for '+atoms[i]+': '+strtrim(nw,2)+' -> '+strtrim(mw,2)

endfor
device,/close & set_plot,'x' & print,'' & print,'$gv chiompress.ps &' & print,''

end

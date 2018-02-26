pro xspecplot,x,y1,y2,y1err=y1err,y2err=y2err,xerr=xerr,$
	y1over=y1over,y2over=y2over,$
	fraclow=fraclow,verbose=verbose,outfile=outfile,$
	xtitle=xtitle,title=title,$
	Ulstyl=Ulstyl,Upsym=Upsym,Ucol=Ucol,Uylog=Uylog,Uytitle=Uytitle,Uyr=Uyr,$
	Llstyl=Llstyl,Lpsym=Lpsym,Lcol=Lcol,Lylog=Lylog,Lytitle=Lytitle,Lyr=Lyr,$
	oUlstyl=oUlstyl,oUpsym=oUpsym,oUcol=oUcol,oLlstyl=oLlstyl,oLpsym=oLpsym,oLcol=oLcol,$
	linestyle=linestyle,psym=psym,color=color,ylog=ylog,ytitle=ytitle,yrange=yrange,$
	_extra=e
;+
;procedure	xspecplot
;	plot two curves in two panels the way XSPEC does spectra and
;	residuals
;
;syntax
;	xspecplot,x,y1,y2,y1err=y1err,y2err=y2err,xerr=xerr,$
;	y1over=y1over,y2over=y2over,$
;	fraclow=fraclow,verbose=verbose,outfile=outfile,$
;	Ulstyl=Ulstyl,Upsym=Upsym,Ucol=Ucol,Uylog=Uylog,Uytitle=Uytitle,Uyr=Uyr,$
;	Llstyl=Llstyl,Lpsym=Lpsym,Lcol=Lcol,Lylog=Lylog,Lytitle=Lytitle,Lyr=Lyr,$
;	oUlstyl=oUlstyl,oUpsym=oUpsym,oUcol=oUcol,oLlstyl=oLlstyl,oLpsym=oLpsym,oLcol=oLcol,$
;	[PLOT keywords]
;	
;
;parameters
;	x	[INPUT; required] abscissa array
;	y1	[INPUT; required] array that goes in the top panel
;	y2	[INPUT; required] array that goes in the bottom panel
;		* sizes of all arrays _must_ match
;
;keywords
;	y1err	[INPUT] error bar on Y1
;	y2err	[INPUT] error bar on Y2
;	xerr	[INPUT] error bar on Y2
;		* if scalar and +ve, assumed to be a fractional error
;		  if between 0 and 1, and %error if greater than 1
;		* if scalar and -ve, abs value is assumed to be a constant
;		  error
;		* if array of same size as X, assumed to be symmetric
;		* if array of size [2,NX] or [NX,2], assumed to be upper
;		  and lower bounds
;		* ignored if array size doesn't match X
;	y1over	[INPUT] additional arrays to overplot in upper panel
;	y2over	[INPUT] additional arrays to overplot in lower panel
;		* size must be [m,NX] or [NX,m]
;		* if not set, both assumed to be 0*NX
;	fraclow	[INPUT] what fraction of the plot should the lower panel
;		take up?
;		* default is 0.35
;		* hard minimum of 0.05 and hard maximum of 0.95
;	verbose	[INPUT] controls chatter
;	outfile	[INPUT] if set, makes a postscript plot
;		* if set to a string, assumes that that is the name of
;		  the file, otherwise writes to ./xspecplot.ps
;	Ulstyl	[INPUT] line style for upper panel
;	Llstyl	[INPUT] line style for lower panel
;	linestyle	[INPUT] caught and discarded
;	Upsym	[INPUT] PSYM for upper panel
;	Lpsym	[INPUT] PSYM for lower panel
;	psym	[INPUT] caught and discarded
;	Ucol	[INPUT] COLOR for upper panel
;	Lcol	[INPUT] COLOR for lower panel
;	color	[INPUT] caught and discarded
;	Uylog	[INPUT] YLOG for upper panel
;	Lylog	[INPUT] YLOG for lower panel
;	ylog	[INPUT] caught and discarded
;	Uytitle	[INPUT] YTITLE for upper panel
;	Lytitle	[INPUT] YTITLE for lower panel
;	ytitle	[INPUT] caught and discarded
;	Uyr	[INPUT] YRANGE for upper panel
;	Lyr	[INPUT] YRANGE for lower panel
;	oUlstyl	[INPUT] line styles of Y1OVER (must be array of size M)
;	oLlstyl	[INPUT] line styles of Y2OVER (must be array of size M)
;	oUpsym	[INPUT] PSYM for Y1OVER (must be array of size M)
;	oLpsym	[INPUT] PSYM for Y2OVER (must be array of size M)
;	oUcol	[INPUT] COLOR for Y1OVER (must be array of size M)
;	oLcol	[INPUT] COLOR for Y2OVER (must be array of size M)
;	yrange	[INPUT] caught and discarded
;	
;	_extra	[INPUT ONLY] pass defined keywords to PLOT
;
;history
;	vinay kashyap (2012-jun)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & n1=n_elements(y1) & n2=n_elements(y2)
if np lt 3 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if n1 eq 0 then ok='Y1 is not defined' else $
   if n2 eq 0 then ok='Y2 is not defined' else $
    if nx ne n1 then ok='X and Y1 do not match' else $
     if nx ne n2 then ok='X and Y2 do not match'
if ok ne 'ok' then begin
  print,'Usage: xspecplot,x,y1,y2,y1err=y1err,y2err=y2err,xerr=xerr,$'
  print,'       y1over=y1over,y2over=y2over,$'
  print,'       fraclow=fraclow,verbose=verbose,outfile=outfile,$'
  print,'       Ulstyl=Ulstyl,Upsym=Upsym,Ucol=Ucol,Uylog=Uylog,Uytitle=Uytitle,Uyr=Uyr,$'
  print,'       Llstyl=Llstyl,Lpsym=Lpsym,Lcol=Lcol,Lylog=Lylog,Lytitle=Lytitle,Lyr=Lyr,$'
  print,'       oUlstyl=oUlstyl,oUpsym=oUpsym,oUcol=oUcol,oLlstyl=oLlstyl,oLpsym=oLpsym,oLcol=oLcol,$'
  print,'       [PLOT keywords]'
  print,'  make XSPEC-like 2-panel plots'
  if np ne 0 then message,ok,/informational
  return
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=(long(verbose[0]))>1L

if keyword_set(outfile) then begin
  dname=!d.NAME
  if size(outfile,/type) eq 7 then psfile=strtrim(outfile[0],2) else $
  	psfile='xspecplot.ps'
  if vv gt 1 then message,'writing to '+psfile,/informational
  set_plot,'ps' & device,file=psfile,/inches,ysize=10,yoffset=1,xsize=8,xoffset=1,/color
endif
;
n1e=n_elements(y1err) & n2e=n_elements(y2err) & nxe=n_elements(xerr)
y1u=0.*y1 & y1l=y1u & y2u=y1u & y2l=y1u & xu=y1u & xl=y1u
;
if n1e eq 1 then begin	;(scalar
  y1e=y1err[0]
  if y1e lt 0 then begin
    y1u[*]=abs(y1e) & y1l[*]=abs(y1e)
  endif
  if y1e gt 0 and y1e lt 1 then begin
    y1u=abs(y1e*y1) & y1l=abs(y1e*y1)
  endif
endif else begin	;scalar)(array
  if n1e eq nx then begin	;(size matches NX
    y1u=y1err & y1l=y1err
  endif else begin		;N=NX)(does not match NX
    if n1e eq 2*nx then begin	;(N=2*NX
      sz=size(y1err) & m1=sz[1] & m2=sz[2]
      if m1 eq 2 then begin & y1u=reform(y1err[0,*]) & y1l=reform(y1err[1,*]) & endif
      if m2 eq 2 then begin & y1u=reform(y1err[*,0]) & y1l=reform(y1err[*,1]) & endif
    endif			;N=NX*2)
  endelse			;N.ne.NX)
endelse			;array)
;
if n2e eq 1 then begin	;(scalar
  y2e=y2err[0]
  if y2e lt 0 then begin
    y2u[*]=abs(y2e) & y2l[*]=abs(y2e)
  endif
  if y2e gt 0 and y2e lt 1 then begin
    y2u=abs(y2e*y2) & y2l=abs(y2e*y2)
  endif
endif else begin	;scalar)(array
  if n2e eq nx then begin	;(size matches NX
    y2u=y2err & y2l=y2err
  endif else begin		;N=NX)(does not match NX
    if n2e eq 2*nx then begin	;(N=2*NX
      sz=size(y2err) & m1=sz[1] & m2=sz[2]
      if m1 eq 2 then begin & y2u=reform(y2err[0,*]) & y2l=reform(y2err[1,*]) & endif
      if m2 eq 2 then begin & y2u=reform(y2err[*,0]) & y2l=reform(y2err[*,1]) & endif
    endif			;N=NX*2)
  endelse			;N.ne.NX)
endelse			;array)
;
if nxe eq 1 then begin	;(scalar
  xe=xerr[0]
  if xe lt 0 then begin
    xu[*]=abs(xe) & xl[*]=abs(xe)
  endif
  if xe gt 0 and xe lt 1 then begin
    xu=abs(xe*x) & xl=abs(xe*x)
  endif
endif else begin	;scalar)(array
  if nxe eq nx then begin	;(size matches NX
    xu=xerr & xl=xerr
  endif else begin		;N=NX)(does not match NX
    if nxe eq 2*nx then begin	;(N=2*NX
      sz=size(xerr) & m1=sz[1] & m2=sz[2]
      if m1 eq 2 then begin & xu=reform(xerr[0,*]) & xl=reform(xerr[1,*]) & endif
      if m2 eq 2 then begin & xu=reform(xerr[*,0]) & xl=reform(xerr[*,1]) & endif
    endif			;N=NX*2)
  endelse			;N.ne.NX)
endelse			;array)
;
flo=0.35 & if keyword_set(fraclow) then flo=(abs(fraclow[0])>0.05)<0.95
;
if not keyword_set(Upsym) then Upsym=10 & if not keyword_set(Lpsym) then Lpsym=4
if not keyword_set(Ucol) then Ucol=255-!p.background & if not keyword_set(Lcol) then Lcol=255-!p.background
if not keyword_set(Uylog) then Uylog=0 & if not keyword_set(Lylog) then Lylog=0
if not keyword_set(Uytitle) then Uytitle='Y1' & if not keyword_set(Lytitle) then Lytitle='Y2'
if not keyword_set(Uyr) then Uyr=minmax(y1) & if not keyword_set(Lyr) then Lyr=minmax(y2)
if n_elements(Uyr) ne 2 then begin
  if Uyr[0] gt min(y1) then Uyr=[min(y1),Uyr[0]] else Uyr=[Uyr[0],max(y1)]
endif
if n_elements(Lyr) ne 2 then begin
  if Lyr[0] gt min(y2) then Lyr=[min(y2),Lyr[0]] else Lyr=[Lyr[0],max(y2)]
endif
;
if not keyword_set(xtitle) then xtitle='X'
if not keyword_set(title) then title=''
;
y1o=0*y1 & y2o=0*y2 & n1o=n_elements(y1over) & n2o=n_elements(y2over)
if n1o gt 0 then begin
  if n1o eq 1 then y1o[*]=y1over[0] else begin
    if n1o eq nx then y1o=y1over else begin
      sz=size(y1over)
      if sz[0] eq 2 then begin
        if sz[1] eq nx or sz[2] eq nx then begin
	  if sz[1] eq nx then y1o=transpose(y1over) else y1o=y1over
	endif else begin
          if vv gt 1 then message,'cannot understand format of Y1OVER; ignoring',/informational
	endelse
      endif else begin
        if vv gt 1 then message,'cannot understand format of Y1OVER; ignoring',/informational
      endelse
    endelse
  endelse
endif
if n2o gt 0 then begin
  if n2o eq 1 then y2o[*]=y2over[0] else begin
    if n2o eq nx then y2o=y2over else begin
      sz=size(y2over)
      if sz[0] eq 2 then begin
        if sz[1] eq nx or sz[2] eq nx then begin
	  if sz[1] eq nx then y2o=transpose(y2over) else y2o=y2over
	endif else begin
          if vv gt 1 then message,'cannot understand format of Y2OVER; ignoring',/informational
	endelse
      endif else begin
        if vv gt 1 then message,'cannot understand format of Y2OVER; ignoring',/informational
      endelse
    endelse
  endelse
endif
sz=size(y1o) & if sz[0] eq 2 then mn1=sz[1] else mn1=1
sz=size(y2o) & if sz[0] eq 2 then mn2=sz[1] else mn2=1
;
zoUlstyl=intarr(mn1) & zoLlstyl=intarr(mn2)
if n_elements(oUlstyl) eq 1 then zoUlstyl[*]=oUlstyl[0] else $
 if n_elements(oUlstyl) eq mn1 then zoUlstyl=oUlstyl
if n_elements(oLlstyl) eq 1 then zoLlstyl[*]=oLlstyl[0] else $
 if n_elements(oLlstyl) eq mn1 then zoLlstyl=oLlstyl
;
zoUpsym=intarr(mn1)-3 & zoLpsym=intarr(mn2)-3
if n_elements(oUpsym) eq 1 then zoUpsym[*]=oUpsym[0] else $
 if n_elements(oUpsym) eq mn1 then zoUpsym=oUpsym
if n_elements(oLpsym) eq 1 then zoLpsym[*]=oLpsym[0] else $
 if n_elements(oLpsym) eq mn1 then zoLpsym=oLpsym
;
zoUcol=intarr(mn1)+(255-!P.BACKGROUND) & zoLcol=intarr(mn2)+(255-!P.BACKGROUND)
if n_elements(oUcol) eq 1 then zoUcol[*]=oUcol[0] else $
 if n_elements(oUcol) eq mn1 then zoUcol=oUcol
if n_elements(oLcol) eq 1 then zoLcol[*]=oLcol[0] else $
 if n_elements(oLcol) eq mn1 then zoLcol=oLcol

;	initialize
pmulti=!p.multi & pposition=!p.position

!p.multi=[0,1,2]

;!p.position=[0.12,flo*0.84,0.97,0.94]
plot,/nodata,x,y1,ylog=Uylog,ytitle=Uytitle,yrange=Uyr,title=title,xticks=1,xtickname=replicate(' ',2),$
	position=[0.12,flo*0.84,0.97,0.94], _extra=e
oplot,x,y1,psym=Upsym,color=Ucol
for i=0,nx-1 do if y1u[i]+y1l[i] gt 0 then oplot,x[i]*[1,1],y1[i]+[-y1l[i],y1u[i]],color=Ucol
;
if mn1 eq 1 then oplot,x,y1o,linestyle=zoUlstyl[0],psym=zoUpsym[0],color=zoUcol[0], _extra=e else $
  for i=0,mn1-1 do oplot,x,y1o[i,*],linestyle=zoUlstyl[i],psym=zoUpsym[i],color=zoUcol[i], _extra=e

;!p.position=[0.12,0.1,0.97,flo*0.84]
plot,/nodata,x,y2,ylog=Lylog,ytitle=Lytitle,yrange=Lyr,xtitle=xtitle,$
	position=[0.12,0.1,0.97,flo*0.84], _extra=e
for i=0,nx-1 do if y2u[i]+y2l[i] gt 0 then oplot,x[i]*[1,1],y2[i]+[-y2l[i],y2u[i]],color=Lcol
;
oplot,x,y2,psym=Lpsym,color=Lcol
if mn2 eq 1 then oplot,x,y2o,linestyle=zoLlstyl[0],psym=zoLpsym[0],color=zoLcol[0], _extra=e else $
  for i=0,mn2-1 do oplot,x,y2o[i,*],linestyle=zoLlstyl[i],psym=zoLpsym[i],color=zoLcol[i], _extra=e

if vv gt 1000 then stop,'halting; type .CON to continue'

if keyword_set(outfile) then begin
  device,/close & set_plot,dname
  if vv gt 0 then spawn,'ls -l '+psfile
endif

!p.multi=pmulti & !p.position=pposition

return
end

peasecolr & loadct,3 & peasecolr

nx=361 & x=findgen(nx)*!pi/180. & y0=sin(x) & ysig=randomn(seed,nx)*0.1
y1=y0+ysig & y2=ysig
xspecplot,x,y1,y2,y1err=0.1,y2err=y2err,xerr=xerr,$
	y1over=y0,y2over=0,$
	fraclow=0.35,verbose=101,outfile=1,$
	xtitle=xtitle,title='X vs Y1,Y2',$
	Ulstyl=Ulstyl,Upsym=Upsym,Ucol=Ucol,Uylog=Uylog,Uytitle=Uytitle,Uyr=Uyr,$
	Llstyl=Llstyl,Lpsym=Lpsym,Lcol=Lcol,Lylog=Lylog,Lytitle=Lytitle,Lyr=Lyr,$
	oUlstyl=oUlstyl,oUpsym=oUpsym,oUcol=1,oLlstyl=1,oLpsym=oLpsym,oLcol=2,$
	/xs,/ys,xthick=2,ythick=2,thick=2,charthick=2,charsize=1.2

end

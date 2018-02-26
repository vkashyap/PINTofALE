function ti_degti,tarr,tstart,tstop,tdelt,igti=igti,startt=startt,stopt=stopt,$
	verbose=verbose, _extra=e
;+
;function	ti_degti
;	given the event times and GTI start and stop times, removes the
;	gaps from the input and returns a set of points that run continuously.
;	useful for when there are large gaps in the data or large numbers of
;	breaks, and the source does not have any periodicities
;
;syntax
;	tx=ti_degti(tarr,tstart,tstop,tdelt,igti=igti,startt=startt,stopt=stopt,$
;	verbose=verbose, leap=leap)
;
;parameters
;	tarr	[INPUT; required] input time array
;	tstart	[INPUT; required] array of GTI start times
;	tstop	[INPUT; required] array of GTI stop times
;		* TSTART and TSTOP must match in size
;	tdelt	[OUTPUT] add this to the output to get back the input
;
;keywords
;	igti	[OUTPUT] index array that says which interval each element of TARR is in
;		* NOTE: this does not necessarily match TSTART and TSTOP indices,
;		  because this program first cleans them up, and could shrink the
;		  number of distinct intervals.
;	startt	[OUTPUT] cleaned TSTART
;	stopt	[OUTPUT] cleaned TSTOP
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] pass defined keywords to subroutines
;		TI_CLEAN: LEAP
;
;subroutines
;	TI_CLEAN
;
;example
;	.run ti_degti
;
;history
;	Vinay Kashyap (2015sep)
;-

;	usage
ok='ok' & np=n_params() & nt=n_elements(tarr) & n0=n_elements(tstart) & n1=n_elements(tstop)
sz0=size(tstart,/type) & sz1=size(tstop,/type)
if np lt 3 then ok='Insufficient parameters' else $
 if nt eq 0 then ok='TARR is not defined' else $
  if n0 eq 0 then ok='TSTART is not defined' else $
   if n1 eq 0 then ok='TSTOP is not defined' else $
    if n0 ne n1 then ok='TSTART and TSTOP do not match' else $
     if sz0 ne sz1 then ok='TSTART and TSTOP are of different types' else $
      if sz0 eq 6 then ok='TSTART cannot be complex' else $
       if sz0 eq 7 then ok='TSTART cannot be char'
if ok ne 'ok' then begin
  print,'Usage: tx=ti_degti(tarr,tstart,tstop,tdelt,igti=igti,startt=startt,stopt=stopt,$
  print,'       verbose=verbose, leap=leap)'
  print,'  removes data gaps from list of times'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L

;	make sure GTIs are in proper order
startt=tstart & stopt=tstop & ti_clean,startt,stopt, _extra=e	;accepts LEAP on input
ngti=n_elements(startt)	;in case some intervals got merged
tmin=min(startt) & dgti=dblarr(ngti)
;
tx=tarr-tmin & tdelt=0.D*tarr+tmin & igti=lonarr(nt)	;output arrays
for i=1L,ngti-1L do begin
  dgti[i]=startt[i]-stopt[i-1L]
  ok=where(tarr ge startt[i],mok)
  if mok gt 0 then tdelt[ok]=tdelt[ok]+dgti[i]
  if mok gt 0 then igti[ok]=i
endfor
tx=tarr-tdelt

if vv gt 1000 then stop,'halting; type .CON to continue'

return,tx
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	made up example
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	let's make 3 intervals with different rates
tarr=[randomu(seed,100)*10.+20., randomu(seed,200)*10.+50., randomu(seed,100)*20.+70.]
tstart=[20.,50.,70.] & tstop=[30.,60.,90.]
pmulti=!p.multi & !p.multi=[0,1,2]
h0=histogram(tarr,min=0.,max=100.,bin=1.) & xh0=findgen(n_elements(h0))*1.
plot,xh0,h0,psym=10,title='example light curve with data gaps',xtitle='time',ytitle='rate'

;	print call ti_degti()
jnk=ti_degti()
tx=ti_degti(tarr,tstart,tstop,tdelt,igti=igti,verbose=1)
hx=histogram(tx,min=0.,max=total(tstop-tstart),bin=1.) & xhx=findgen(n_elements(hx))*1.

ugti=igti[uniq(igti,sort(igti))] & nugti=n_elements(ugti)

plot,xhx,hx,psym=10,title='light curve with data gaps removed',xtitle='collapsed time',ytitle='rate'
for i=0L,nugti-1L do begin
  oo=where(igti eq ugti[i],moo)
  if moo gt 0 then begin
    xt=tx[oo] & oplot,min(xt)*[1,1],!y.crange,line=2
  endif
endfor
!p.multi=pmulti

end

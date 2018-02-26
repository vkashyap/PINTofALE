;+
;EG_TRACESPEC.PRO
;
;example program for calling LINESPEC for TRACE data
;							vinay kashyap
;-

restore,'trace.save'	;read in TRACE data
nbin=1000		;number of bins

;	for 173 A passband
w0=min(w173,max=w1)
sp173=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww173,effar=a173,wvlar=w173)

plot,w173,sp173,psym=10,title='TRACE 173A passband',xtitle='[Ang]',$
	ytitle='[ph/s/Ang]'

;	for 195 A passband
w0=min(w195,max=w1)
sp195=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww195,effar=a195,wvlar=w195)

plot,w195,sp195,psym=10,title='TRACE 195A passband',xtitle='[Ang]',$
	ytitle='[ph/s/Ang]'

;	for 284 A passband
w0=min(w284,max=w1)
sp284=linespec(atom,wmin=w0,wmax=w1,nbin=nbin,ww=ww284,effar=a284,wvlar=w284)

plot,w284,sp284,psym=10,title='TRACE 284A passband',xtitle='[Ang]',$
	ytitle='[ph/s/Ang]'

end

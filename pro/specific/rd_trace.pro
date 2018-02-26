;+
;RD_TRACE.PRO
;	read in effective areas for the TRACE mirrors
;-

dir='/data/sdata1/'
m173='mirrors_173.save' & m195='mirrors_195.save' & m284='mirrors_284.save'

restore,dir+m173
w173=pri_lam			;wavelength [Ang]
a173=current_ea/100.		;effective area [cm^2]

restore,dir+m195
w195=pri_lam			;wavelength [Ang]
a195=current_ea/100.		;effective area [cm^2]

restore,dir+m284
w284=pri_lam			;wavelength [Ang]
a284=current_ea/100.		;effective area [cm^2]

save,file='trace.save',w173,a173,w195,a195,w284,a284

end

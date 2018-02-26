;+
;EG_FIDGIT.PRO
;
;example program for calling FIDGIT
;							vinay kashyap
;-

;	initialize
;atom=['FeXX','Ni18']	;uncomment to restrict to specified species
logP=15							;pressure = nT [cm^-3 K]
dbdir='$SCAR'				;line DB directory
dem=0.5*(findgen(81)*0.05+4)+11 & dem=10.D^(dem)        ;DEM
abund=getabund('anders')				;abundances
euve_sw,swea,sw						;EUVE SW effective area
exptim=1e5						;exposure time [s]
loadct,3						;set color table
fwhm=10							;line width in pixels

;	read in EUVE SW data
!path = !path + ':/data/drake7/vinay/euve_fit/pro'
@/data/drake7/vinay/euve_fit/read_new2d
infil='/data/drake7/vinay/HR1099/swsp'
read_new2d,infil,oswv,scts,sbkg,simhdr,dir=imdir,genx=usegen
x=oswv & y=scts/exptim

;	pick a small range
oo=pickrange(x,y) & x=x(oo) & y=y(oo)
plot,x,y,xtitle='[Ang]',ytitle='ct s!u-1!n'

;	digitally fit DEM to spectrum
fidgit,x,y,atom,dem=dem,logt=logt,fwhm=fwhm,$
	logP=logP,dbdir=dbdir,abund=abund,effar=swea,wvlar=sw

end

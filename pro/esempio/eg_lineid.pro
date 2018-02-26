;+
;EG_LINEID.PRO
;
;example program to call LINEID
;				-vinay kashyap
;-

;	initialize
if not keyword_set(stor) then stor=1	;to avoid calling RD_LINE many times
wrng=1.1	;search for IDs within 10% of observed line
best=1		;find "best" match
logP=15		;pressure = nT [cm^-3 K]
n_e=1e10	;density = n [cm^-3]
dbdir='$CHIANTI'			;line DB directory
dem=0.5*(findgen(81)*0.05+4)+11 & dem=10.D^(dem)	;DEM
abund=getabund('anders')				;abundances
euve_sw,swea,sw						;EUVE SW effective area

;	mark some features
markw=[132.,140.,180.]
feature=[' Big Fe line',' No lines from here..',' ..till here']

;	read in EUVE SW data
!path = !path + ':/data/drake7/vinay/euve_fit/pro'
;@/data/drake7/vinay/euve_fit/pro/read_new2d
infil='/data/drake7/vinay/HR1099/swsp'
read_new2d,infil,oswv,scts,sbkg,simhdr,dir=imdir,genx=usegen
x=oswv & y=scts

;	get LINE IDs
lid=lineid(x,y,wrng=wrng,stor=stor,best=best,logP=logP,dbdir=dbdir,$
	n_e=n_e,dem=dem,abund=abund,effar=swea,wvlar=sw,$
	markw=markw,feature=feature)

;	save for future reference
lamda=x & spec=y
save,file='tmp/eg_lineid.sav',lamda,spec,wrng,stor,best,logP,dbdir,dem,$
	swea,sw,lid

end

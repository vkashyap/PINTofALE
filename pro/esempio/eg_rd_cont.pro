eden=1e10 & fcstr=1 & wrange=[1.,500.] & nbin=-500L
t0=6.5
wmn=(wrange(0) < (wrange(1))) > (0.1)
wmx=(wrange(1) > (wrange(0))) > (0.1)
defabu=getabund('anders & grevesse')
!p.multi=[0,2,2] & loadct,3

;	only H+He
atm='He' & symb2zion,atm,zz,ii
elem=[atm] & Z=[1,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,fcstr=fcstr,abund=abnd)
f02=reform(fc(50,*)) & print,''
wc=fcstr.wvl & nw=n_elements(wc) & wc=wc(0:nw-2)
;
g02=f02
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
cd,'../emissivity/'
c02=reform(cf(0,*))
;
d02=c02
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c02-f02,/xl,title=tt,psym=7
plot,wc,c02,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f02,col=150
plot,wc,d02-g02,/xl,psym=7,title='H+He'
plot,wc,d02,/xl & oplot,wc,g02,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	you can rearrange the blocks below in any order
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	add O
atm='O' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f08=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g08=reform(tmp(50,*)) & print,''
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c08=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d08=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c08-f08,/xl,title=tt,psym=7
plot,wc,c08,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f08,col=150
plot,wc,d08-g08,/xl,psym=7,title='H+He+'+atm
plot,wc,d08,/xl & oplot,wc,g08,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add C
atm='C' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f06=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g06=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c06=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d06=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c06-f06,/xl,title=tt,psym=7
plot,wc,c06,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f06,col=150
plot,wc,d06-g06,/xl,psym=7,title='H+He+'+atm
plot,wc,d06,/xl & oplot,wc,g06,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add N
atm='N' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f07=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g07=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c07=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d07=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c07-f07,/xl,title=tt,psym=7
plot,wc,c07,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f07,col=150
plot,wc,d07-g07,/xl,psym=7,title='H+He+'+atm
plot,wc,d07,/xl & oplot,wc,g07,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Na
atm='Na' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f11=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g11=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c11=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d11=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c11-f11,/xl,title=tt,psym=7
plot,wc,c11,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f11,col=150
plot,wc,d11-g11,/xl,psym=7,title='H+He+'+atm
plot,wc,d11,/xl & oplot,wc,g11,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Ne
atm='Ne' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f10=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g10=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c10=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d10=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c10-f10,/xl,title=tt,psym=7
plot,wc,c10,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f10,col=150
plot,wc,d10-g10,/xl,psym=7,title='H+He+'+atm
plot,wc,d10,/xl & oplot,wc,g10,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Al
atm='Al' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f13=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g13=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c13=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d13=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c13-f13,/xl,title=tt,psym=7
plot,wc,c13,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f13,col=150
plot,wc,d13-g13,/xl,psym=7,title='H+He+'+atm
plot,wc,d13,/xl & oplot,wc,g13,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Si
atm='Si' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f14=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g14=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c14=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d14=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c14-f14,/xl,title=tt,psym=7
plot,wc,c14,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f14,col=150
plot,wc,d14-g14,/xl,psym=7,title='H+He+'+atm
plot,wc,d14,/xl & oplot,wc,g14,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add S
atm='S' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f16=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g16=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c16=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d16=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c16-f16,/xl,title=tt,psym=7
plot,wc,c16,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f16,col=150
plot,wc,d16-g16,/xl,psym=7,title='H+He+'+atm
plot,wc,d16,/xl & oplot,wc,g16,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Ar
atm='Ar' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f18=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g18=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c18=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d18=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c18-f18,/xl,title=tt,psym=7
plot,wc,c18,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f18,col=150
plot,wc,d18-g18,/xl,psym=7,title='H+He+'+atm
plot,wc,d18,/xl & oplot,wc,g18,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Ca
atm='Ca' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f20=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g20=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c20=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d20=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c20-f20,/xl,title=tt,psym=7
plot,wc,c20,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f20,col=150
plot,wc,d20-g20,/xl,psym=7,title='H+He+'+atm
plot,wc,d20,/xl & oplot,wc,g20,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Fe
atm='Fe' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f26=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g26=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c26=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d26=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c26-f26,/xl,title=tt,psym=7
plot,wc,c26,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f26,col=150
plot,wc,d26-g26,/xl,psym=7,title='H+He+'+atm
plot,wc,d26,/xl & oplot,wc,g26,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Ni
atm='Ni' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f28=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g28=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c28=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d28=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c28-f28,/xl,title=tt,psym=7
plot,wc,c28,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f28,col=150
plot,wc,d28-g28,/xl,psym=7,title='H+He+'+atm
plot,wc,d28,/xl & oplot,wc,g28,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

;	add Mg
atm='Mg' & symb2zion,atm,zz,ii
elem=[elem,atm] & Z=[Z,zz]
;
abnd=defabu & abnd(Z-1)=-abnd(Z-1) & abnd=(-abnd)>0
fc=rd_cont('cie',n_e=eden,abund=abnd)
f12=reform(fc(50,*)) & print,''
;
Z2=[1,2,zz]
abnd=defabu & abnd(Z2-1)=-abnd(Z2-1) & abnd=(-abnd)>0
tmp=rd_cont('cie',n_e=eden,abund=abnd)
g12=reform(tmp(50,*))
;
cd,'../CIE/'
cf=cont_cie(elem,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
c12=reform(cf(0,*))
;
elem2=['He',atm]
tmp=cont_cie(elem2,pr,t0,cw,wmn=wmn,wmx=wmx,nbin=nbin,n_e=eden)
d12=reform(tmp(0,*))
cd,'../emissivity/'
;
tt='H' & nelm=n_elements(elem) & for i=0,nelm-1 do tt=tt+'+'+elem(i)
plot,wc,c12-f12,/xl,title=tt,psym=7
plot,wc,c12,/xl,title='CIE in white, RD_CONT in color' & oplot,wc,f12,col=150
plot,wc,d12-g12,/xl,psym=7,title='H+He+'+atm
plot,wc,d12,/xl & oplot,wc,g12,col=120
c1='any key to continue, q to stop' & print,c1 & c1=get_kbrd(1)
if c1 eq 'q' then stop

end

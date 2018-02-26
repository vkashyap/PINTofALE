;pro exgu,filename,linstr
;+
;-

if not keyword_set(filename) then message,'FILENAME is undefined'

tmp=read_ascii(filename)
gudata=tmp.(0)
nn=reform(fix(gudata[0,*]))
ll=reform(fix(gudata[1,*]))
ul=reform(fix(gudata[2,*]))
tlog=reform(float(gudata[3,*]))
kev=reform(float(gudata[4,*]))
ang=reform(float(gudata[5,*]))
de=reform(double(gudata[6,*]))
re=reform(double(gudata[7,*]))
rec=reform(double(gudata[8,*]))
ion=reform(double(gudata[9,*]))
numlin=n_elements(nn)

;	how many unique temperature bins?
oT=uniq(tlog,sort(tlog)) & mT=n_elements(oT) & utlog=tlog[oT]
;	expand to default
nT=81L & logT=findgen(nT)*0.05+4.

;	how many unique lines?
nwvl=long(float(numlin)/float(mT))
if nwvl ne n_elements(uniq(kev,sort(kev))) then message,'BUG!'
uwvl=fltarr(nwvl) & ull=intarr(nwvl) & uul=ull & k=0L
uwvl[k]=ang[0] & ull[k]=ll[0] & uul[k]=ul[0]
for i=1L,numlin-1L do begin
  if ang[i] ne uwvl[k] or ll[i] ne ull[k] or ul[i] ne uul[k] then begin
    k=k+1L
    uwvl[k]=ang[i] & ull[k]=ll[i] & uul[k]=ul[i]
  endif
endfor
if k ne nwvl-1L then message,'BUG!'

;	emissivity matrix
emis_dere=dblarr(mT,nwvl) & emis_rec=dblarr(mT,nwvl) & emis_ion=dblarr(mT,nwvl)
emisdere=dblarr(nT,nwvl) & emisrec=dblarr(nT,nwvl) & emision=dblarr(nT,nwvl)

;	put the data into the arrays
for i=0L,nwvl-1L do begin
  ok=where(ang eq uwvl[i] and ll eq ull[i] and ul eq uul[i],mok)
  if mok eq 0 then message,'BUG!'
  emis_dere[*,i]=alog10(de[ok]+re[ok] > 1e-30)
  emis_rec[*,i]=alog10(rec[ok] > 1e-30)
  emis_ion[*,i]=alog10(ion[ok] > 1e-30)
  emisdere[*,i]=(interpol(emis_dere[*,i],utlog,logt)>(-30)) < max(emis_dere[*,i])
  emisrec[*,i]=(interpol(emis_rec[*,i],utlog,logt)>(-30)) < max(emis_rec[*,i])
  emision[*,i]=(interpol(emis_ion[*,i],utlog,logt)>(-30)) < max(emis_ion[*,i])
endfor

;	construct LINSTR
emis=dblarr(nT,nwvl*3)
wvl=fltarr(nwvl*3)
Z=intarr(nwvl*3)+26
ion=26-[nn,nn,nn]+1
jon=ion
desig=strarr(2,nwvl*3)
config=strarr(2,nwvl*3)
src=intarr(nwvl*3)+5
for i=0,nwvl-1L do begin
  emis[*,3*i]=10.D^(emisdere[*,i])
  emis[*,3*i+1L]=10.D^(emisrec[*,i]) & jon[3*i+1L]=jon[3*i+1L]+1L
  emis[*,3*i+2L]=10.D^(emision[*,i]) & jon[3*i+2L]=jon[3*i+2L]-1L
  wvl[3*i:3*i+2L]=uwvl[i]
  desig[0,3*i:3*i+2L]=strtrim(ull[i],2)
  desig[1,3*i:3*i+2L]=strtrim(uul[i],2)
endfor
linstr=create_struct('LINE_INT',emis,'LOGT',logt,'WVL',wvl,'Z',Z,$
	'ION',ion,'DESIG',desig,'CONFIG',config,'SRC',src,'JON',jon)

end

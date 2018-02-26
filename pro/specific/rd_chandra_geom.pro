function rd_chandra_geom,help=help,caldb=caldb,version=version,$
	verbose=verbose, _extra=e
;+
;function	rd_chandra_geom
;	reads in the entire geometry file from Chandra's CALDB
;	into an IDL structure and spits it out
;
;syntax
;	geom=rd_chandra_geom(/help,caldb=caldb,version=version,verbose=verbose)
;
;parameters	NONE
;
;keywords
;	help	[INPUT] if set, prints out the usage help
;	caldb	[INPUT] top level directory containing the Chandra CALDB
;		* default is to look for environment variable CALDB
;		* if $CALDB is not defined, then assumes /soft/ciao/CALDB/
;	version	[INPUT] version number to use
;		* default is 5
;		* the geometry file is assumed to be at
;		  $CALDB/data/chandra/tel/bcf/geom/telD1999-07-23geomN####.fits
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Oct04)
;	changed default VERSION to 5 (VK; Jul05)
;-

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
if not keyword_set(caldb) then begin
  caldb=getenv('CALDB')
  if not keyword_set(caldb) then caldb='/soft/ciao/CALDB/'
endif
;
nver=5 & if keyword_set(version) then nver=fix(version[0])>2

;	usage
if keyword_set(help) or vv gt 10 then begin
  print,'Usage: geom=rd_chandra_geom(/help,caldb=caldb,version=version,verbose=verbose)'
  print,'  read the Chandra geometry file into an IDL structure'
endif

;	define the input file
dir=caldb+'/data/chandra/tel/bcf/geom/'
geomfil='telD1999-07-23geomN'+string(nver,'(i4.4)')+'.fits'

;	there are 7 extensions.  read them all in
if vv gt 0 then message,'Reading from: '+dir+geomfil,/informational
ext1=mrdfits(dir+geomfil,1,hext1)
ext2=mrdfits(dir+geomfil,2,hext2)
ext3=mrdfits(dir+geomfil,3,hext3)
ext4=mrdfits(dir+geomfil,4,hext4)
ext5=mrdfits(dir+geomfil,5,hext5)
ext6=mrdfits(dir+geomfil,6,hext6)
ext7=mrdfits(dir+geomfil,7,hext7)

;	figure out their names
extnam=['INSTRUMENTS','GRATINGS',$
	'INSTGEOM1','INSTGEOM2','INSTGEOM3','INSTGEOM4','INSTGEOM5']
iext=intarr(n_elements(extnam))
o1=where(strtrim(sxpar(hext1,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext1,'HDUNAME'),2) eq extnam,mo1) & iext[0]=o1[0]
o2=where(strtrim(sxpar(hext2,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext2,'HDUNAME'),2) eq extnam,mo2) & iext[1]=o2[0]
o3=where(strtrim(sxpar(hext3,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext3,'HDUNAME'),2) eq extnam,mo3) & iext[2]=o3[0]
o4=where(strtrim(sxpar(hext4,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext4,'HDUNAME'),2) eq extnam,mo4) & iext[3]=o4[0]
o5=where(strtrim(sxpar(hext5,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext5,'HDUNAME'),2) eq extnam,mo5) & iext[4]=o5[0]
o6=where(strtrim(sxpar(hext6,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext6,'HDUNAME'),2) eq extnam,mo6) & iext[5]=o6[0]
o7=where(strtrim(sxpar(hext7,'EXTNAME'),2) eq extnam or $
	strtrim(sxpar(hext7,'HDUNAME'),2) eq extnam,mo7) & iext[6]=o7[0]
ok='ok'
if mo1 eq 0 then ok='Extension '+extnam[1-1]+' is missing' else $
 if mo2 eq 0 then ok='Extension '+extnam[2-1]+' is missing' else $
  if mo3 eq 0 then ok='Extension '+extnam[3-1]+' is missing' else $
   if mo4 eq 0 then ok='Extension '+extnam[4-1]+' is missing' else $
    if mo5 eq 0 then ok='Extension '+extnam[5-1]+' is missing' else $
     if mo6 eq 0 then ok='Extension '+extnam[6-1]+' is missing' else $
      if mo7 eq 0 then ok='Extension '+extnam[7-1]+' is missing'
if ok ne 'ok' then begin
  message,ok,/informational
  if vv gt 100 then stop
  return,-1L
endif

;	make the output
if iext[0] ne 0 then j=execute('ht0=hext'+strtrim(iext[0]+1,2)) else ht0=hext1
if iext[1] ne 1 then j=execute('ht1=hext'+strtrim(iext[1]+1,2)) else ht1=hext2
if iext[2] ne 2 then j=execute('ht2=hext'+strtrim(iext[2]+1,2)) else ht2=hext3
if iext[3] ne 3 then j=execute('ht3=hext'+strtrim(iext[3]+1,2)) else ht3=hext4
if iext[4] ne 4 then j=execute('ht4=hext'+strtrim(iext[4]+1,2)) else ht4=hext5
if iext[5] ne 5 then j=execute('ht5=hext'+strtrim(iext[5]+1,2)) else ht5=hext6
if iext[6] ne 6 then j=execute('ht6=hext'+strtrim(iext[6]+1,2)) else ht6=hext7
;
if iext[0] ne 0 then j=execute('xt0=ext'+strtrim(iext[0]+1,2)) else xt0=ext1
if iext[1] ne 1 then j=execute('xt1=ext'+strtrim(iext[1]+1,2)) else xt1=ext2
if iext[2] ne 2 then j=execute('xt2=ext'+strtrim(iext[2]+1,2)) else xt2=ext3
if iext[3] ne 3 then j=execute('xt3=ext'+strtrim(iext[3]+1,2)) else xt3=ext4
if iext[4] ne 4 then j=execute('xt4=ext'+strtrim(iext[4]+1,2)) else xt4=ext5
if iext[5] ne 5 then j=execute('xt5=ext'+strtrim(iext[5]+1,2)) else xt5=ext6
if iext[6] ne 6 then j=execute('xt6=ext'+strtrim(iext[6]+1,2)) else xt6=ext7

geom=create_struct(	'INSTRUMENTS',xt0,$
			'GRATINGS',xt1,$
			'INSTGEOM1',xt2,$
			'INSTGEOM2',xt3,$
			'INSTGEOM3',xt4,$
			'INSTGEOM4',xt5,$
			'INSTGEOM5',xt6,$
			'hINSTRUMENTS',ht0,$
			'hGRATINGS',ht1,$
			'hINSTGEOM1',ht2,$
			'hINSTGEOM2',ht3,$
			'hINSTGEOM3',ht4,$
			'hINSTGEOM4',ht5,$
			'hINSTGEOM5',ht6 )

return,geom
end

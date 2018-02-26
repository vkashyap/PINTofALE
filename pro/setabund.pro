;+
;function	setabund
;	widget-based editor to set abundances
;
;syntax
;	abund=setabund(init=init,norm=norm,/help)
;
;parameters	NONE
;
;keywords
;	init	[INPUT] initial set of abundances -- default is Anders &
;		Gervasse (1989)
;	help	[INPUT] as it says
;	norm	[INPUT] same as in GETABUND, to decide whether the
;		abundances are in normal or log format.
;		* abs(NORM) > 1 ==> log
;	_extra	[INPUT] junk -- here only to prevent program from crashing
;
;restrictions
;	widget-based
;
;subroutines
;	GETABUND
;	GETEMALL, SETEMALL (appended)
;
;history
;	vinay kashyap (Dec96)
;	bug correction (MULT not being set) (VK; Jan96)
;-

;----------------------------------------------------------
pro getemall,wid,abund,frac,nelem
;
;get all the values.
;	WID	[INPUT] all the widget IDs
;	ABUND	[OUTPUT] all the abundances
;	FRAC	[OUTPUT] all the deviations
;	NELEM	[optional INPUT] number of elements (this is really useless)
;

;initialize
if n_elements(nelem) eq 0 then nelem=30
abund=fltarr(nelem) & frac=abund

;unpack WID
ii=0
hx_abu=wid(ii) & ii=ii+1
he_abu=wid(ii) & ii=ii+1
li_abu=wid(ii) & ii=ii+1
be_abu=wid(ii) & ii=ii+1
bx_abu=wid(ii) & ii=ii+1
cx_abu=wid(ii) & ii=ii+1
nx_abu=wid(ii) & ii=ii+1
ox_abu=wid(ii) & ii=ii+1
fx_abu=wid(ii) & ii=ii+1
ne_abu=wid(ii) & ii=ii+1
na_abu=wid(ii) & ii=ii+1
mg_abu=wid(ii) & ii=ii+1
al_abu=wid(ii) & ii=ii+1
si_abu=wid(ii) & ii=ii+1
px_abu=wid(ii) & ii=ii+1
sx_abu=wid(ii) & ii=ii+1
cl_abu=wid(ii) & ii=ii+1
ar_abu=wid(ii) & ii=ii+1
kx_abu=wid(ii) & ii=ii+1
ca_abu=wid(ii) & ii=ii+1
sc_abu=wid(ii) & ii=ii+1
ti_abu=wid(ii) & ii=ii+1
vx_abu=wid(ii) & ii=ii+1
cr_abu=wid(ii) & ii=ii+1
mn_abu=wid(ii) & ii=ii+1
fe_abu=wid(ii) & ii=ii+1
co_abu=wid(ii) & ii=ii+1
ni_abu=wid(ii) & ii=ii+1
cu_abu=wid(ii) & ii=ii+1
zn_abu=wid(ii) & ii=ii+1
hx_twk=wid(ii) & ii=ii+1
he_twk=wid(ii) & ii=ii+1
li_twk=wid(ii) & ii=ii+1
be_twk=wid(ii) & ii=ii+1
bx_twk=wid(ii) & ii=ii+1
cx_twk=wid(ii) & ii=ii+1
nx_twk=wid(ii) & ii=ii+1
ox_twk=wid(ii) & ii=ii+1
fx_twk=wid(ii) & ii=ii+1
ne_twk=wid(ii) & ii=ii+1
na_twk=wid(ii) & ii=ii+1
mg_twk=wid(ii) & ii=ii+1
al_twk=wid(ii) & ii=ii+1
si_twk=wid(ii) & ii=ii+1
px_twk=wid(ii) & ii=ii+1
sx_twk=wid(ii) & ii=ii+1
cl_twk=wid(ii) & ii=ii+1
ar_twk=wid(ii) & ii=ii+1
kx_twk=wid(ii) & ii=ii+1
ca_twk=wid(ii) & ii=ii+1
sc_twk=wid(ii) & ii=ii+1
ti_twk=wid(ii) & ii=ii+1
vx_twk=wid(ii) & ii=ii+1
cr_twk=wid(ii) & ii=ii+1
mn_twk=wid(ii) & ii=ii+1
fe_twk=wid(ii) & ii=ii+1
co_twk=wid(ii) & ii=ii+1
ni_twk=wid(ii) & ii=ii+1
cu_twk=wid(ii) & ii=ii+1
zn_twk=wid(ii) & ii=ii+1

;get abund
ii=0
widget_control,hx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,he_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,li_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,be_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,bx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,cx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,nx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,ox_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,fx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,ne_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,na_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,mg_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,al_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,si_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,px_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,sx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,cl_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,ar_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,kx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,ca_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,sc_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,ti_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,vx_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,cr_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,mn_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,fe_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,co_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,ni_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,cu_abu,get_value=ab & abund(ii)=ab & ii=ii+1
widget_control,zn_abu,get_value=ab & abund(ii)=ab & ii=ii+1

;get frac
ii=0
widget_control,hx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,he_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,li_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,be_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,bx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,cx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,nx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,ox_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,fx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,ne_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,na_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,mg_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,al_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,si_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,px_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,sx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,cl_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,ar_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,kx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,ca_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,sc_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,ti_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,vx_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,cr_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,mn_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,fe_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,co_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,ni_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,cu_twk,get_value=fr & frac(ii)=fr & ii=ii+1
widget_control,zn_twk,get_value=fr & frac(ii)=fr & ii=ii+1

return
end
;----------------------------------------------------------
pro setemall,wid,abund,frac
;
;set all the values.
;	WID	[INPUT] all the widget IDs
;	ABUND	[INPUT] all the abundances
;	FRAC	[INPUT] all the deviations
;

;unpack WID
ii=0
hx_abu=wid(ii) & ii=ii+1
he_abu=wid(ii) & ii=ii+1
li_abu=wid(ii) & ii=ii+1
be_abu=wid(ii) & ii=ii+1
bx_abu=wid(ii) & ii=ii+1
cx_abu=wid(ii) & ii=ii+1
nx_abu=wid(ii) & ii=ii+1
ox_abu=wid(ii) & ii=ii+1
fx_abu=wid(ii) & ii=ii+1
ne_abu=wid(ii) & ii=ii+1
na_abu=wid(ii) & ii=ii+1
mg_abu=wid(ii) & ii=ii+1
al_abu=wid(ii) & ii=ii+1
si_abu=wid(ii) & ii=ii+1
px_abu=wid(ii) & ii=ii+1
sx_abu=wid(ii) & ii=ii+1
cl_abu=wid(ii) & ii=ii+1
ar_abu=wid(ii) & ii=ii+1
kx_abu=wid(ii) & ii=ii+1
ca_abu=wid(ii) & ii=ii+1
sc_abu=wid(ii) & ii=ii+1
ti_abu=wid(ii) & ii=ii+1
vx_abu=wid(ii) & ii=ii+1
cr_abu=wid(ii) & ii=ii+1
mn_abu=wid(ii) & ii=ii+1
fe_abu=wid(ii) & ii=ii+1
co_abu=wid(ii) & ii=ii+1
ni_abu=wid(ii) & ii=ii+1
cu_abu=wid(ii) & ii=ii+1
zn_abu=wid(ii) & ii=ii+1
hx_twk=wid(ii) & ii=ii+1
he_twk=wid(ii) & ii=ii+1
li_twk=wid(ii) & ii=ii+1
be_twk=wid(ii) & ii=ii+1
bx_twk=wid(ii) & ii=ii+1
cx_twk=wid(ii) & ii=ii+1
nx_twk=wid(ii) & ii=ii+1
ox_twk=wid(ii) & ii=ii+1
fx_twk=wid(ii) & ii=ii+1
ne_twk=wid(ii) & ii=ii+1
na_twk=wid(ii) & ii=ii+1
mg_twk=wid(ii) & ii=ii+1
al_twk=wid(ii) & ii=ii+1
si_twk=wid(ii) & ii=ii+1
px_twk=wid(ii) & ii=ii+1
sx_twk=wid(ii) & ii=ii+1
cl_twk=wid(ii) & ii=ii+1
ar_twk=wid(ii) & ii=ii+1
kx_twk=wid(ii) & ii=ii+1
ca_twk=wid(ii) & ii=ii+1
sc_twk=wid(ii) & ii=ii+1
ti_twk=wid(ii) & ii=ii+1
vx_twk=wid(ii) & ii=ii+1
cr_twk=wid(ii) & ii=ii+1
mn_twk=wid(ii) & ii=ii+1
fe_twk=wid(ii) & ii=ii+1
co_twk=wid(ii) & ii=ii+1
ni_twk=wid(ii) & ii=ii+1
cu_twk=wid(ii) & ii=ii+1
zn_twk=wid(ii) & ii=ii+1

;set abund
ii=0
widget_control,hx_abu,set_value=abund(ii) & ii=ii+1
widget_control,he_abu,set_value=abund(ii) & ii=ii+1
widget_control,li_abu,set_value=abund(ii) & ii=ii+1
widget_control,be_abu,set_value=abund(ii) & ii=ii+1
widget_control,bx_abu,set_value=abund(ii) & ii=ii+1
widget_control,cx_abu,set_value=abund(ii) & ii=ii+1
widget_control,nx_abu,set_value=abund(ii) & ii=ii+1
widget_control,ox_abu,set_value=abund(ii) & ii=ii+1
widget_control,fx_abu,set_value=abund(ii) & ii=ii+1
widget_control,ne_abu,set_value=abund(ii) & ii=ii+1
widget_control,na_abu,set_value=abund(ii) & ii=ii+1
widget_control,mg_abu,set_value=abund(ii) & ii=ii+1
widget_control,al_abu,set_value=abund(ii) & ii=ii+1
widget_control,si_abu,set_value=abund(ii) & ii=ii+1
widget_control,px_abu,set_value=abund(ii) & ii=ii+1
widget_control,sx_abu,set_value=abund(ii) & ii=ii+1
widget_control,cl_abu,set_value=abund(ii) & ii=ii+1
widget_control,ar_abu,set_value=abund(ii) & ii=ii+1
widget_control,kx_abu,set_value=abund(ii) & ii=ii+1
widget_control,ca_abu,set_value=abund(ii) & ii=ii+1
widget_control,sc_abu,set_value=abund(ii) & ii=ii+1
widget_control,ti_abu,set_value=abund(ii) & ii=ii+1
widget_control,vx_abu,set_value=abund(ii) & ii=ii+1
widget_control,cr_abu,set_value=abund(ii) & ii=ii+1
widget_control,mn_abu,set_value=abund(ii) & ii=ii+1
widget_control,fe_abu,set_value=abund(ii) & ii=ii+1
widget_control,co_abu,set_value=abund(ii) & ii=ii+1
widget_control,ni_abu,set_value=abund(ii) & ii=ii+1
widget_control,cu_abu,set_value=abund(ii) & ii=ii+1
widget_control,zn_abu,set_value=abund(ii)

;set frac
ii=0
widget_control,hx_twk,set_value=frac(ii) & ii=ii+1
widget_control,he_twk,set_value=frac(ii) & ii=ii+1
widget_control,li_twk,set_value=frac(ii) & ii=ii+1
widget_control,be_twk,set_value=frac(ii) & ii=ii+1
widget_control,bx_twk,set_value=frac(ii) & ii=ii+1
widget_control,cx_twk,set_value=frac(ii) & ii=ii+1
widget_control,nx_twk,set_value=frac(ii) & ii=ii+1
widget_control,ox_twk,set_value=frac(ii) & ii=ii+1
widget_control,fx_twk,set_value=frac(ii) & ii=ii+1
widget_control,ne_twk,set_value=frac(ii) & ii=ii+1
widget_control,na_twk,set_value=frac(ii) & ii=ii+1
widget_control,mg_twk,set_value=frac(ii) & ii=ii+1
widget_control,al_twk,set_value=frac(ii) & ii=ii+1
widget_control,si_twk,set_value=frac(ii) & ii=ii+1
widget_control,px_twk,set_value=frac(ii) & ii=ii+1
widget_control,sx_twk,set_value=frac(ii) & ii=ii+1
widget_control,cl_twk,set_value=frac(ii) & ii=ii+1
widget_control,ar_twk,set_value=frac(ii) & ii=ii+1
widget_control,kx_twk,set_value=frac(ii) & ii=ii+1
widget_control,ca_twk,set_value=frac(ii) & ii=ii+1
widget_control,sc_twk,set_value=frac(ii) & ii=ii+1
widget_control,ti_twk,set_value=frac(ii) & ii=ii+1
widget_control,vx_twk,set_value=frac(ii) & ii=ii+1
widget_control,cr_twk,set_value=frac(ii) & ii=ii+1
widget_control,mn_twk,set_value=frac(ii) & ii=ii+1
widget_control,fe_twk,set_value=frac(ii) & ii=ii+1
widget_control,co_twk,set_value=frac(ii) & ii=ii+1
widget_control,ni_twk,set_value=frac(ii) & ii=ii+1
widget_control,cu_twk,set_value=frac(ii) & ii=ii+1
widget_control,zn_twk,set_value=frac(ii)

return
end

;----------------------------------------------------------
;	the main subroutine starts here
;----------------------------------------------------------
function setabund,init=init,norm=norm,help=help,_extra=e

;	usage
if keyword_set(help) then begin
  print,'Usage: abund=setabund(init=init,norm=norm,/help)'
  print,'  widget-based editor to set abundances'
endif

;	initialize
defabu=getabund('anders & gervasse',norm=norm,_extra=e)
if not keyword_set(init) then ab0=defabu else ab0=init
elem=[	'H','He','Li','Be','B', 'C','N','O','F','Ne','Na','Mg','Al',$
	'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',$
	'Ni','Cu','Zn']
nab=n_elements(ab0) & nelem=n_elements(elem)
if nab lt nelem then ab0=[ab0,defabu(nab:*)]	;fill out with defaults
if nab gt nelem then ab0=ab0(0:nelem-1)		;ignore the rest
if keyword_set(norm) then begin			;add or multiply deltas?
  if abs(norm) gt 1 then mult='+' else mult='X'
endif else mult='X'
frac=fltarr(nelem) & if mult eq 'X' then frac(*)=1. & fr0=frac

;	set up widget layout
;1st row: DONE, RESET, TWEAK
;all other rows: ELEMENT: FIELD, RESET, ELEMENT: DEVIATION

;	base widget
base=widget_base(frame=2,/column,/map,title='SET ABUNDANCES',uvalue='base')

;	control menu
base_menu=widget_base(base,frame=2,/row,uvalue='base_menu')
;
menu_done=widget_button(base_menu,value='DONE',uvalue='menu_done')
menu_rset=widget_button(base_menu,value='RESET',uvalue='menu_rset')
menu_nudg=widget_button(base_menu,value='TWEAK',uvalue='menu_nudg')

;	set abundances
base_abnd=widget_base(base,frame=2,/scroll,scr_xsize=600,scr_ysize=600,$
	/column,uvalue='base_abnd') & ii=0
;
abnd_hx=widget_base(base_abnd,/row,uvalue='abnd_hx') & el='H'
  hx_abu=cw_field(abnd_hx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='hx_abu')
  hx_def=widget_button(abnd_hx,value='RESET',uvalue='hx_def') & ii=ii+1
  hx_twk=cw_field(abnd_hx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='hx_twk')
abnd_he=widget_base(base_abnd,/row,uvalue='abnd_he') & el='He'
  he_abu=cw_field(abnd_he,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='he_abu')
  he_def=widget_button(abnd_he,value='RESET',uvalue='he_def') & ii=ii+1
  he_twk=cw_field(abnd_he,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='he_twk')
abnd_li=widget_base(base_abnd,/row,uvalue='abnd_li') & el='Li'
  li_abu=cw_field(abnd_li,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='li_abu')
  li_def=widget_button(abnd_li,value='RESET',uvalue='li_def') & ii=ii+1
  li_twk=cw_field(abnd_li,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='li_twk')
abnd_be=widget_base(base_abnd,/row,uvalue='abnd_be') & el='Be'
  be_abu=cw_field(abnd_be,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='be_abu')
  be_def=widget_button(abnd_be,value='RESET',uvalue='be_def') & ii=ii+1
  be_twk=cw_field(abnd_be,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='be_twk')
abnd_bx=widget_base(base_abnd,/row,uvalue='abnd_bx') & el='B'
  bx_abu=cw_field(abnd_bx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='bx_abu')
  bx_def=widget_button(abnd_bx,value='RESET',uvalue='bx_def') & ii=ii+1
  bx_twk=cw_field(abnd_bx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='bx_twk')
abnd_cx=widget_base(base_abnd,/row,uvalue='abnd_cx') & el='C'
  cx_abu=cw_field(abnd_cx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='cx_abu')
  cx_def=widget_button(abnd_cx,value='RESET',uvalue='cx_def') & ii=ii+1
  cx_twk=cw_field(abnd_cx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='cx_twk')
abnd_nx=widget_base(base_abnd,/row,uvalue='abnd_nx') & el='N'
  nx_abu=cw_field(abnd_nx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='nx_abu')
  nx_def=widget_button(abnd_nx,value='RESET',uvalue='nx_def') & ii=ii+1
  nx_twk=cw_field(abnd_nx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='nx_twk')
abnd_ox=widget_base(base_abnd,/row,uvalue='abnd_ox') & el='O'
  ox_abu=cw_field(abnd_ox,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='ox_abu')
  ox_def=widget_button(abnd_ox,value='RESET',uvalue='ox_def') & ii=ii+1
  ox_twk=cw_field(abnd_ox,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='ox_twk')
abnd_fx=widget_base(base_abnd,/row,uvalue='abnd_fx') & el='F'
  fx_abu=cw_field(abnd_fx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='fx_abu')
  fx_def=widget_button(abnd_fx,value='RESET',uvalue='fx_def') & ii=ii+1
  fx_twk=cw_field(abnd_fx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='fx_twk')
abnd_ne=widget_base(base_abnd,/row,uvalue='abnd_ne') & el='Ne'
  ne_abu=cw_field(abnd_ne,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='ne_abu')
  ne_def=widget_button(abnd_ne,value='RESET',uvalue='ne_def') & ii=ii+1
  ne_twk=cw_field(abnd_ne,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='ne_twk')
abnd_na=widget_base(base_abnd,/row,uvalue='abnd_na') & el='Na'
  na_abu=cw_field(abnd_na,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='na_abu')
  na_def=widget_button(abnd_na,value='RESET',uvalue='na_def') & ii=ii+1
  na_twk=cw_field(abnd_na,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='na_twk')
abnd_mg=widget_base(base_abnd,/row,uvalue='abnd_mg') & el='Mg'
  mg_abu=cw_field(abnd_mg,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='mg_abu')
  mg_def=widget_button(abnd_mg,value='RESET',uvalue='mg_def') & ii=ii+1
  mg_twk=cw_field(abnd_mg,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='mg_twk')
abnd_al=widget_base(base_abnd,/row,uvalue='abnd_al') & el='Al'
  al_abu=cw_field(abnd_al,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='al_abu')
  al_def=widget_button(abnd_al,value='RESET',uvalue='al_def') & ii=ii+1
  al_twk=cw_field(abnd_al,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='al_twk')
abnd_si=widget_base(base_abnd,/row,uvalue='abnd_si') & el='Si'
  si_abu=cw_field(abnd_si,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='si_abu')
  si_def=widget_button(abnd_si,value='RESET',uvalue='si_def') & ii=ii+1
  si_twk=cw_field(abnd_si,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='si_twk')
abnd_px=widget_base(base_abnd,/row,uvalue='abnd_px') & el='P'
  px_abu=cw_field(abnd_px,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='px_abu')
  px_def=widget_button(abnd_px,value='RESET',uvalue='px_def') & ii=ii+1
  px_twk=cw_field(abnd_px,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='px_twk')
abnd_sx=widget_base(base_abnd,/row,uvalue='abnd_sx') & el='S'
  sx_abu=cw_field(abnd_sx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='sx_abu')
  sx_def=widget_button(abnd_sx,value='RESET',uvalue='sx_def') & ii=ii+1
  sx_twk=cw_field(abnd_sx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='sx_twk')
abnd_cl=widget_base(base_abnd,/row,uvalue='abnd_cl') & el='Cl'
  cl_abu=cw_field(abnd_cl,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='cl_abu')
  cl_def=widget_button(abnd_cl,value='RESET',uvalue='cl_def') & ii=ii+1
  cl_twk=cw_field(abnd_cl,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='cl_twk')
abnd_ar=widget_base(base_abnd,/row,uvalue='abnd_ar') & el='Ar'
  ar_abu=cw_field(abnd_ar,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='ar_abu')
  ar_def=widget_button(abnd_ar,value='RESET',uvalue='ar_def') & ii=ii+1
  ar_twk=cw_field(abnd_ar,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='ar_twk')
abnd_kx=widget_base(base_abnd,/row,uvalue='abnd_kx') & el='K'
  kx_abu=cw_field(abnd_kx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='kx_abu')
  kx_def=widget_button(abnd_kx,value='RESET',uvalue='kx_def') & ii=ii+1
  kx_twk=cw_field(abnd_kx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='kx_twk')
abnd_ca=widget_base(base_abnd,/row,uvalue='abnd_ca') & el='Ca'
  ca_abu=cw_field(abnd_ca,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='ca_abu')
  ca_def=widget_button(abnd_ca,value='RESET',uvalue='ca_def') & ii=ii+1
  ca_twk=cw_field(abnd_ca,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='ca_twk')
abnd_sc=widget_base(base_abnd,/row,uvalue='abnd_sc') & el='Sc'
  sc_abu=cw_field(abnd_sc,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='sc_abu')
  sc_def=widget_button(abnd_sc,value='RESET',uvalue='sc_def') & ii=ii+1
  sc_twk=cw_field(abnd_sc,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='sc_twk')
abnd_ti=widget_base(base_abnd,/row,uvalue='abnd_ti') & el='Ti'
  ti_abu=cw_field(abnd_ti,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='ti_abu')
  ti_def=widget_button(abnd_ti,value='RESET',uvalue='ti_def') & ii=ii+1
  ti_twk=cw_field(abnd_ti,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='ti_twk')
abnd_vx=widget_base(base_abnd,/row,uvalue='abnd_vx') & el='V'
  vx_abu=cw_field(abnd_vx,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='vx_abu')
  vx_def=widget_button(abnd_vx,value='RESET',uvalue='vx_def') & ii=ii+1
  vx_twk=cw_field(abnd_vx,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='vx_twk')
abnd_cr=widget_base(base_abnd,/row,uvalue='abnd_cr') & el='Cr'
  cr_abu=cw_field(abnd_cr,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='cr_abu')
  cr_def=widget_button(abnd_cr,value='RESET',uvalue='cr_def') & ii=ii+1
  cr_twk=cw_field(abnd_cr,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='cr_twk')
abnd_mn=widget_base(base_abnd,/row,uvalue='abnd_mn') & el='Mn'
  mn_abu=cw_field(abnd_mn,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='mn_abu')
  mn_def=widget_button(abnd_mn,value='RESET',uvalue='mn_def') & ii=ii+1
  mn_twk=cw_field(abnd_mn,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='mn_twk')
abnd_fe=widget_base(base_abnd,/row,uvalue='abnd_fe') & el='Fe'
  fe_abu=cw_field(abnd_fe,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='fe_abu')
  fe_def=widget_button(abnd_fe,value='RESET',uvalue='fe_def') & ii=ii+1
  fe_twk=cw_field(abnd_fe,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='fe_twk')
abnd_co=widget_base(base_abnd,/row,uvalue='abnd_co') & el='Co'
  co_abu=cw_field(abnd_co,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='co_abu')
  co_def=widget_button(abnd_co,value='RESET',uvalue='co_def') & ii=ii+1
  co_twk=cw_field(abnd_co,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='co_twk')
abnd_ni=widget_base(base_abnd,/row,uvalue='abnd_ni') & el='Ni'
  ni_abu=cw_field(abnd_ni,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='ni_abu')
  ni_def=widget_button(abnd_ni,value='RESET',uvalue='ni_def') & ii=ii+1
  ni_twk=cw_field(abnd_ni,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='ni_twk')
abnd_cu=widget_base(base_abnd,/row,uvalue='abnd_cu') & el='Cu'
  cu_abu=cw_field(abnd_cu,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='cu_abu')
  cu_def=widget_button(abnd_cu,value='RESET',uvalue='cu_def') & ii=ii+1
  cu_twk=cw_field(abnd_cu,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='cu_twk')
abnd_zn=widget_base(base_abnd,/row,uvalue='abnd_zn') & el='Zn'
  zn_abu=cw_field(abnd_zn,/float,/return_events,/row,title=el,$
	value=ab0(ii),uvalue='zn_abu')
  zn_def=widget_button(abnd_zn,value='RESET',uvalue='zn_def')
  zn_twk=cw_field(abnd_zn,/float,/return_events,/row,title=mult,$
	value=frac(ii),uvalue='zn_twk')

;	start up widgets
widget_control,base,/realize,/hourglass

;	pack widget IDs of interest
wid=[ hx_abu, he_abu, li_abu, be_abu, bx_abu, cx_abu, nx_abu, ox_abu,$
      fx_abu, ne_abu, na_abu, mg_abu, al_abu, si_abu, px_abu, sx_abu,$
      cl_abu, ar_abu, kx_abu, ca_abu, sc_abu, ti_abu, vx_abu, cr_abu,$
      mn_abu, fe_abu, co_abu, ni_abu, cu_abu, zn_abu,$
      hx_twk, he_twk, li_twk, be_twk, bx_twk, cx_twk, nx_twk, ox_twk,$
      fx_twk, ne_twk, na_twk, mg_twk, al_twk, si_twk, px_twk, sx_twk,$
      cl_twk, ar_twk, kx_twk, ca_twk, sc_twk, ti_twk, vx_twk, cr_twk,$
      mn_twk, fe_twk, co_twk, ni_twk, cu_twk, zn_twk ]

ok=''
while ok ne 'quit' do begin		;{I/O
  event=widget_event(base) & widget_control,event.id,get_uvalue=ev
  case ev of
  ;
  ;	menu bar
  ;
  'menu_done': begin
    getemall,wid,abund,frac,nelem
    ok='quit'
  end
  'menu_rset': begin
    setemall,wid,ab0,fr0
  end
  'menu_nudg': begin
    getemall,wid,abund,frac,nelem
    if mult eq '+' then abund=abund+frac else abund=abund*frac
    setemall,wid,abund,frac
  end
  ;
  ;	reset individual abundances
  ;
  'hx_def': widget_control,hx_abu,set_value=ab0(0)
  'he_def': widget_control,he_abu,set_value=ab0(1)
  'li_def': widget_control,li_abu,set_value=ab0(2)
  'be_def': widget_control,be_abu,set_value=ab0(3)
  'bx_def': widget_control,bx_abu,set_value=ab0(4)
  'cx_def': widget_control,cx_abu,set_value=ab0(5)
  'nx_def': widget_control,nx_abu,set_value=ab0(6)
  'ox_def': widget_control,ox_abu,set_value=ab0(7)
  'fx_def': widget_control,fx_abu,set_value=ab0(8)
  'ne_def': widget_control,ne_abu,set_value=ab0(9)
  'na_def': widget_control,na_abu,set_value=ab0(10)
  'mg_def': widget_control,mg_abu,set_value=ab0(11)
  'al_def': widget_control,al_abu,set_value=ab0(12)
  'si_def': widget_control,si_abu,set_value=ab0(13)
  'px_def': widget_control,px_abu,set_value=ab0(15)
  'sx_def': widget_control,sx_abu,set_value=ab0(16)
  'cl_def': widget_control,cl_abu,set_value=ab0(17)
  'ar_def': widget_control,ar_abu,set_value=ab0(18)
  'kx_def': widget_control,kx_abu,set_value=ab0(19)
  'ca_def': widget_control,ca_abu,set_value=ab0(20)
  'sc_def': widget_control,sc_abu,set_value=ab0(21)
  'ti_def': widget_control,ti_abu,set_value=ab0(22)
  'vx_def': widget_control,vx_abu,set_value=ab0(23)
  'cr_def': widget_control,cr_abu,set_value=ab0(24)
  'mn_def': widget_control,mn_abu,set_value=ab0(25)
  'fe_def': widget_control,fe_abu,set_value=ab0(26)
  'co_def': widget_control,co_abu,set_value=ab0(27)
  'ni_def': widget_control,ni_abu,set_value=ab0(28)
  'cu_def': widget_control,cu_abu,set_value=ab0(29)
  'zn_def': widget_control,zn_abu,set_value=ab0(30)
  else: begin
    getemall,wid,abund,frac,nelem
  end
  endcase
endwhile				;ok='quit'}

;	clean up
widget_control,base,/destroy

return,abund
end

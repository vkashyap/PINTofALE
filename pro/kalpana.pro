pro kalpana,x,y,ststr,idstr,oststr=oststr,wshow=wshow,wwork=wwork, _extra=e
;+
;procedure	kalpana
;	widget-based procedure to make an annotated plot, e.g., of a spectrum.
;
;	"kalpana" means "imagination".  or "delirium", if you prefer.
;	pronounced "cull-pun-ah"
;
;syntax
;	kalpana,x,y,ststr,idstr,oststr=oststr,wshow=wshow,wwork=wwork
;
;parameters
;	x	[INPUT; required] points where the curve is defined
;	y	[INPUT; required] the curve Y[X] (usually a spectrum)
;		* size of Y must match that of X
;	ststr	[I/O; required] the state structure containing all
;		the necessary information on what labels to put where.
;		* see >>description<< below
;	idstr	[INPUT] the output of LINEID, containing the wavelengths,
;		IDs, labels, etc. of identified features in the spectrum.
;		* if given, STSTR is suitably modified.
;
;keywords
;	oststr	[OUTPUT] returns the original STSTR, in case modifications
;		are not acceptable
;	wshow	[INPUT] display final result in this window number (def: 1)
;	wwork	[INPUT] display working plot in this window number (def: 0)
;	_extra	[INPUT] pass defined keywords to subroutines
;
;description
;	The state structure contains the following substructures
;	WINDOW:
;		MULTI: !P.MULTI
;		TITLE: overall title
;		SUBTITLE: array of subtitles for each plot if >1,
;			the subtitle for one if one.
;		XTITLE: string array of xtitles, for each of the plots
;		YTITLE: string array of ytitles
;		XSTYLE: integer aarray of xstyles
;		YSTYLE: integer array of ystyles
;		XLOG: integer array for x-axis in log style (1) or linear (0)
;		YLOG: integer array for y-axis in log style (1) or linear (0)
;		CHARS: character size for labels
;		XMIN: array of minimum values of XRANGE
;		XMAX: array of minimum values of XRANGE
;		YMIN: array of minimum values of YRANGE
;		YMAX: array of minimum values of YRANGE
;	LOC:
;		X: x-values to which the labels apply
;		Y: y-values to which the labels apply
;		GROUP: index in X from which to take attributes.
;		  e.g., if STSTR.LOC.GROUP(I+1)=J+1, then STSTR.LI.(ALIGN,
;		  ORIENT,SIZE,THICK,LABCOLOR,LINCOLOR,ARRANGE,UNDERLINE,
;		  SIDELINE) are taken from STSTR.LJ
;	L1..N:	one structure for each of LOC, containing
;		POS: [x,y] at which label is to be written
;		LABEL: the label text (could be a string array)
;		XPATH: array of intermediate x-positions of line connecting
;			POS to LOC
;		YPATH: as XPATH, for intermediate y-positions
;		ALIGN: the alignment of the label wrt POS (LEFT/RIGHT/UP/DOWN)
;		ORIENT: the orientation of the label text
;		SIZE: size of characters in LABEL
;		THICK: thickness of line in PATH
;		LABCOLOR: a number from 1..255
;		LINCOLOR: a number from 1..255
;		ARRANGE: if multiple fields in LABEL, arrange them in
;			a ROW or a COLUMN
;		UNDERLINE: underline (+ve), overline (-ve), no line (0)
;			with values >1 indicating extent of pincer tip
;		SIDELINE: same as UNDERLINE, but sideways
;
;restrictions
;	requires graphics capability
;	subroutines:
;		KALPANA_EVENT [PICKLABL, MOVELABL]
;		DARSHANA
;		MUDRA
;		STR_2_ARR
;	hardcopies are not working, needs a few more features to be
;	  really useful, and still somewhat klunky.
;
;history
;	vinay kashyap (SepMIM)
;-

;	usage
ok='ok'
nx=n_elements(x) & ny=n_elements(y) & np=n_params()
if nx eq 0 then ok='Missing X' else $
 if ny eq 0 then ok='Missing Y(X)' else $
  if nx ne ny then ok='Y(X) and X are incompatible' else $
   if np lt 3 then ok='insufficiant parameters'
if ok ne 'ok' then begin
  print,'Usage: kalpana,x,y,ststr,idstr,oststr=oststr'
  print,'  GUI to make annotated spectral plots'
  message,ok,/info & return
endif

;	insert ID structure data into state structure, if either exist
mudra,ststr,idstr,oststr=oststr, _extra=e

;	initialize displays
if n_elements(wshow) eq 0 then wshow=1
if n_elements(wwork) eq 0 then wwork=0
darshana,x,y,ststr,wid=wshow, _extra=e		;all
darshana,x,y,ststr,wid=wwork,pid=0		;just the first

;	define the widget base heirarchy to control the procedure
base=widget_base(frame=2,/row,/map,title='KALPANA',uvalue='base')
;organization:
;	1st column			2nd column
;	1st row -- menu buttons		tweak
;	2nd row -- status display
;	3rd row -- window controls
;	4th row -- plot controls
;	5th row -- label controls

;{first column
base_col1=widget_base(base,frame=2,/column,uvalue='base_col1')

;{first row
;	buttons for SHOW, QUIT, SAVE, RESTORE, HELP,
;	menu for PRINT,
;	entry field for OUTROOT,
;	button for DEVICE
col1_menu=widget_base(base_col1,frame=2,/row,uvalue='col1_menu')
;
menu_show=widget_button(col1_menu,value='SHOW',uvalue='menu_show')
menu_quit=widget_button(col1_menu,value='QUIT',uvalue='menu_quit')
menu_save=widget_button(col1_menu,value='SAVE',uvalue='menu_save')
menu_bkup=widget_button(col1_menu,value='RESTORE',uvalue='menu_bkup')
menu_help=widget_button(col1_menu,value='HELP',uvalue='menu_help')
menu_prnt=widget_button(col1_menu,/menu,value='PRINT',uvalue='menu_prnt')
prnt_post=widget_button(menu_prnt,value='POSTSCRIPT',uvalue='prnt_post')
prnt_gifs=widget_button(menu_prnt,value='GIF',uvalue='prnt_gifs')
menu_file=cw_field(col1_menu,/string,/return_events,xsize=10,$
	title=' ',value='root',uvalue='menu_file')
menu_dpar=widget_button(col1_menu,value='DEVICE',uvalue='menu_dpar')
;
menu=[col1_menu,menu_show,menu_quit,menu_save,menu_bkup,menu_help,$
	menu_prnt,prnt_post,prnt_gifs,menu_file,menu_dpar]
;first row}

;{second row
;	text status message
col1_txts=widget_base(base_col1,frame=2,/row,uvalue='col1_txts',/align_left)
;
txts_stts=widget_text(col1_txts,font='10x20',frame=2,/scroll,/wrap,$
	xsize=55,ysize=2,value='Widget-based graph annotation program',$
	uvalue='txts_stts')
;
txts=[col1_txts,txts_stts]
;second row}

;{third row
;	window controls
;	entry fields for NCOL, NROW, TITLE
;	entry fields for plot color, frame color, background color
col1_wndo=widget_base(base_col1,frame=2,/column,uvalue='col1_wndo')
;
ncol=(ststr.window.multi)(1)
nrow=(ststr.window.multi)(2)
;
wndo_wind=widget_base(col1_wndo,/row,uvalue='wndo_wind',/align_left)
wind_ncol=cw_field(wndo_wind,/integer,/return_events,title='# cols',$
	xsize=2,value=nrow,uvalue='wind_ncol')
wind_nrow=cw_field(wndo_wind,/integer,/return_events,title='# rows',$
	xsize=2,value=ncol,uvalue='wind_nrow')
wind_titl=cw_field(wndo_wind,/string,/return_events,title='Plot Title',$
	xsize=30,value=(ststr.window.title)(0),uvalue='wind_titl')
;
wndo_cols=widget_base(col1_wndo,/row,uvalue='wndo_cols',/align_left)
cols_pltn=cw_field(wndo_cols,/string,/return_events,title='Plot Color',$
	xsize=7,value='',uvalue='cols_pltn')
cols_frmn=cw_field(wndo_cols,/string,/return_events,title='Frame Color',$
	xsize=7,value='',uvalue='cols_frmn')
cols_bkgn=cw_field(wndo_cols,/string,/return_events,title='Background Color',$
	xsize=7,value='',uvalue='cols_bkgn')
;
wind=[col1_wndo,$
	wndo_wind,wind_ncol,wind_nrow,wind_titl,$
	wndo_cols,cols_pltn,cols_frmn,cols_bkgn]
;third row}

;{fourth row
;	plot controls
;	entry field for PLOT#, title, and character size, menu for DESTROY
;	X: entry fields for title and range, droplists for axis type and style 
;	Y: entry fields for title and range, droplists for axis type and style 
cxr=strtrim(min(x),2)+','+strtrim(max(x),2)
cyr=strtrim(min(y),2)+','+strtrim(max(y),2)
axistype=['LINEAR','LOG']
axistyle=['Default Axis','Exact Axis','Extend','Suppress','Lower Only']
col1_plot=widget_base(base_col1,frame=2,/column,uvalue='col1_plot',/align_left)
;
plot_row1=widget_base(col1_plot,/row,uvalue='plot_row1',/align_left)
row1_pnum=cw_field(plot_row1,/long,/return_events,title='Plot #',$
	xsize=2,value=-1L,uvalue='row1_pnum')
row1_subt=cw_field(plot_row1,/string,/return_events,title='Title',$
	xsize=20,value='',uvalue='row1_subt')
row1_csiz=cw_field(plot_row1,/float,/return_events,title='CharSize',$
	xsize=3,value=1,uvalue='row1_csiz')
row1_dstr=widget_button(plot_row1,value='DESTROY',/menu,uvalue='row1_dstr')
dstr_plot=widget_button(row1_dstr,value='Plot',uvalue='dstr_plot')
;
plot_rowx=widget_base(col1_plot,/row,uvalue='plot_rowx',/align_left)
rowx_xtit=cw_field(plot_rowx,/string,/return_events,title='Xtitle',$
	xsize=10,value='X',uvalue='rowx_xtit')
rowx_xrng=cw_field(plot_rowx,/string,/return_events,title='range',$
	xsize=15,value=cxr,uvalue='rowx_xrng')
rowx_xsty=cw_field(plot_rowx,/integer,/return_events,title='Axis:',$
	xsize=2,value=0,uvalue='rowx_xsty')
rowx_xlog=widget_droplist(plot_rowx,value=axistype,uvalue='rowx_xlog')
;
plot_rowy=widget_base(col1_plot,/row,uvalue='plot_rowy',/align_left)
rowy_ytit=cw_field(plot_rowy,/string,/return_events,title='Ytitle',$
	xsize=10,value='y',uvalue='rowy_ytit')
rowy_yrng=cw_field(plot_rowy,/string,/return_events,title='range',$
	xsize=15,value=cyr,uvalue='rowy_yrng')
rowy_ysty=cw_field(plot_rowy,/integer,/return_events,title='Axis:',$
	xsize=2,value=0,uvalue='rowy_ysty')
rowy_ylog=widget_droplist(plot_rowy,value=axistype,uvalue='rowy_ylog')
;
plot=[col1_plot,$
	plot_row1,row1_pnum,row1_subt,row1_csiz,row1_dstr,dstr_plot,$
	plot_rowx,rowx_xtit,rowx_xrng,rowx_xsty,rowx_xlog,$
	plot_rowy,rowy_ytit,rowy_yrng,rowy_ysty,rowy_ylog]
;fourth row}

;{fifth row
;	label controls
;	buttons for SELECT/NEXT/PREVIOUS LABEL, entry field for SELECT LABEL
;	  droplist for HIDE, menu for DESTROY
;	entry field for LABEL TEXT
;	droplist for ALIGN (l/r/u/d), ARRANGE (c/r)
;	entry fields for ORIENT, SIZE, THICK,
;	entry fields for LABEL COLOR/NAME, LINE COLOR/NAME
;	entry fields for UNDERLINE, SIDELINE
;;;	LABCOLOR, LINCOLOR, UNDERLINE, SIDELINE
col1_labl=widget_base(base_col1,frame=2,/column,uvalue='col1_labl')
;
hide=['SHOW','HIDE']
labl_row0=widget_base(col1_labl,/row,uvalue='labl_row0')
row0_sele=widget_button(labl_row0,value='SELECT',uvalue='row0_sele')
row0_next=widget_button(labl_row0,value='>>',uvalue='row0_next')
row0_prev=widget_button(labl_row0,value='<<',uvalue='row0_prev')
row0_labl=cw_field(labl_row0,/long,/return_events,title='Label #',$
	xsize=3,value=-1L,uvalue='row0_labl')
row0_hide=widget_droplist(labl_row0,value=hide,uvalue='row0_hide')
row0_newl=widget_button(labl_row0,value='ADD NEW',uvalue='row0_newl')
row0_dele=widget_button(labl_row0,value='DESTROY',/menu,uvalue='row0_dele')
dele_labl=widget_button(row0_dele,value='Label',uvalue='dele_labl')
;
labl_row1=widget_base(col1_labl,/row,uvalue='labl_row1')
row1_ltxt=cw_field(labl_row1,/string,/return_events,title='Label',$
	xsize=50,value='',uvalue='row1_ltxt')
;
qstyl=[	'STYLE',$
	'hang/top','hang/top/slide',$
	'vert/top','vert/top/slide',$
	'flat/top','flat/top/slide',$
	'straight line','half-T','inverted half-T','half-Y',$
	'inverted half-Y','drop+half-Y']
align=['LEFT','RIGHT','UP','DOWN']
arrange=['COLUMN','ROW']
labl_row2=widget_base(col1_labl,/row,uvalue='labl_row2')
row2_styl=widget_droplist(labl_row2,value=qstyl,title='Label: ',$
	uvalue='row2_styl')
row2_algn=widget_droplist(labl_row2,value=align,title='Align',$
	uvalue='row2_algn')
row2_arng=widget_droplist(labl_row2,value=arrange,title='Arrange',$
	uvalue='row2_arng')
;
labl_row3=widget_base(col1_labl,/row,uvalue='labl_row3')
row3_ornt=cw_field(labl_row3,/float,/return_events,title='Orientation',$
	xsize=3,value=90,uvalue='row3_ornt')
row3_size=cw_field(labl_row3,/float,/return_events,title='Charcter Size',$
	xsize=2,value=1,uvalue='row3_size')
row3_thck=cw_field(labl_row3,/float,/return_events,title='Line Thickness',$
	xsize=2,value=1,uvalue='row3_thck')
;
labl_row4=widget_base(col1_labl,/row,uvalue='labl_row4')
row4_labc=cw_field(labl_row4,/integer,/return_events,title='Label Color: #',$
	xsize=3,value=100,uvalue='row4_labc')
row4_labn=cw_field(labl_row4,/string,/return_events,title='Name',$
	xsize=7,value='',uvalue='row4_labn')
row4_linc=cw_field(labl_row4,/integer,/return_events,title='Line Color: #',$
	xsize=3,value=150,uvalue='row4_linc')
row4_linn=cw_field(labl_row4,/string,/return_events,title='Name',$
	xsize=7,value='',uvalue='row4_linn')
;
labl_row5=widget_base(col1_labl,/row,uvalue='labl_row5')
row5_ulin=cw_field(labl_row5,/float,/return_events,title='Underline',$
	xsize=4,value=1,uvalue='row5_ulin')
row5_slin=cw_field(labl_row5,/float,/return_events,title='Sideline',$
	xsize=4,value=0,uvalue='row5_slin')
;
labl=[col1_labl,$
	labl_row0,row0_sele,row0_next,row0_prev,row0_labl,$
	row0_hide,row0_newl,row0_dele,dele_labl, labl_row1,row1_ltxt,$
	labl_row2,row2_styl,row2_algn,row2_arng,$
	labl_row3,row3_ornt,row3_size,row3_thck,$
	labl_row4,row4_labc,row4_labn,row4_linc,row4_linn,$
	labl_row5,row5_ulin,row5_slin]
;fifth row}

;first column}{second column

base_col2=widget_base(base,frame=2,/column,uvalue='base_col2')
;	tweak attributes
;	button to ADJUST
;	button to UNDO
;	droplist to make a GROUP (ALL/SELECT/DISCARD/DISPERSE)
;	menu for ALIGN-LABEL (ALL/GROUPED/SELECTED)
;	menu for ORIENT (ALL/GROUPED/SELECTED)
;	menu for LABEL SIZE (ALL/GROUPED/SELECTED)
;	menu for LINE THICKNESS (ALL/GROUPED/SELECTED)
;	menu for ARRANGEMENT (ALL/GROUPED/SELECTED)
;	menu for UNDERLINE (ALL/GROUPED/SELECTED)
;	menu for SIDELINE (ALL/GROUPED/SELECTED)
;	droplist for COARSE/EXACT setting of XPATH/YPATH, and
;	  droplist for SPRINGY/HINGY/STICKY setting or adjusting XPATH/YPATH
;	entry field for COARSENESS
;	entry field for YTOP
;	entry field for YBOT
coarse=['COARSE','EXACT'] & tieup=['SPRINGY','HINGY','STICKY']
col2_move=widget_button(base_col2,value='ADJUST',uvalue='col2_move')
col2_undo=widget_button(base_col2,value='UNDO',uvalue='col2_undo')
;
col2_grps=widget_button(base_col2,value='GROUP',/menu,uvalue='col2_grps')
grps_alle=widget_button(col2_grps,value='ALL',uvalue='grps_alle')
grps_sele=widget_button(col2_grps,value='SELECT',uvalue='grps_sele')
grps_dele=widget_button(col2_grps,value='DISCARD',uvalue='grps_dele')
grps_undo=widget_button(col2_grps,value='DISPERSE',uvalue='grps_undo')
;
col2_alny=widget_button(base_col2,value='Align Labels',/menu,$
	uvalue='col2_alny')
alny_alle=widget_button(col2_alny,value='ALL',uvalue='alny_alle')
alny_grps=widget_button(col2_alny,value='GROUPED',uvalue='alny_grps')
alny_sele=widget_button(col2_alny,value='SELECTED',uvalue='alny_sele')
;
col2_ornt=widget_button(base_col2,value='Match Orientations',/menu,$
	uvalue='col2_ornt')
ornt_alle=widget_button(col2_ornt,value='ALL',uvalue='ornt_alle')
ornt_grps=widget_button(col2_ornt,value='GROUPED',uvalue='ornt_grps')
ornt_sele=widget_button(col2_ornt,value='SELECTED',uvalue='ornt_sele')
;
col2_lsiz=widget_button(base_col2,value='Match Sizes',/menu,$
	uvalue='col2_lsiz')
lsiz_alle=widget_button(col2_lsiz,value='ALL',uvalue='lsiz_alle')
lsiz_grps=widget_button(col2_lsiz,value='GROUPED',uvalue='lsiz_grps')
lsiz_sele=widget_button(col2_lsiz,value='SELECTED',uvalue='lsiz_sele')
;
col2_thck=widget_button(base_col2,value='Match Thicknesses',/menu,$
	uvalue='col2_thck')
thck_alle=widget_button(col2_thck,value='ALL',uvalue='thck_alle')
thck_grps=widget_button(col2_thck,value='GROUPED',uvalue='thck_grps')
thck_sele=widget_button(col2_thck,value='SELECTED',uvalue='thck_sele')
;
col2_arng=widget_button(base_col2,value='Match Arrangements',/menu,$
	uvalue='col2_arng')
arng_alle=widget_button(col2_arng,value='ALL',uvalue='arng_alle')
arng_grps=widget_button(col2_arng,value='GROUPED',uvalue='arng_grps')
arng_sele=widget_button(col2_arng,value='SELECTED',uvalue='arng_sele')
;
col2_ulin=widget_button(base_col2,value='Match Underlines',/menu,$
	uvalue='col2_ulin')
ulin_alle=widget_button(col2_ulin,value='ALL',uvalue='ulin_alle')
ulin_grps=widget_button(col2_ulin,value='GROUPED',uvalue='ulin_grps')
ulin_sele=widget_button(col2_ulin,value='SELECTED',uvalue='ulin_sele')
;
col2_slin=widget_button(base_col2,value='Match Sidelines',/menu,$
	uvalue='col2_slin')
slin_alle=widget_button(col2_slin,value='ALL',uvalue='slin_alle')
slin_grps=widget_button(col2_slin,value='GROUPED',uvalue='slin_grps')
slin_sele=widget_button(col2_slin,value='SELECTED',uvalue='slin_sele')
;
col2_posn=widget_droplist(base_col2,title='Positioning',value=coarse,$
	uvalue='col2_posn')
col2_ties=widget_droplist(base_col2,title='Rigidity',value=tieup,$
	uvalue='col2_ties')
col2_coar=cw_field(base_col2,/float,/return_events,title='Coarseness',$
	xsize=3,value=0.2,uvalue='col2_coar')
col2_ytop=cw_field(base_col2,/string,/return_events,title='Y-TOP',$
	xsize=10,value='0,0.1',uvalue='col2_ytop')
col2_ybot=cw_field(base_col2,/string,/return_events,title='Y-BOT',$
	xsize=10,value='0,1.1',uvalue='col2_ybot')
col2_yfrc=cw_field(base_col2,/float,/return_events,title='Y-Fraction',$
	xsize=3,value='0.1',uvalue='col2_yfrc')

tweak=[col2_move,col2_undo,col2_grps,grps_alle,grps_sele,grps_dele,grps_undo,$
	col2_alny,alny_alle,alny_grps,alny_sele,$
	col2_ornt,ornt_alle,ornt_grps,ornt_sele,$
	col2_lsiz,lsiz_alle,lsiz_grps,lsiz_sele,$
	col2_thck,thck_alle,thck_grps,thck_sele,$
	col2_arng,arng_alle,arng_grps,arng_sele,$
	col2_ulin,ulin_alle,ulin_grps,ulin_sele,$
	col2_slin,slin_alle,slin_grps,slin_sele,$
	col2_posn,col2_ties,col2_coar,col2_ytop,col2_ybot,col2_yfrc]
;second column}

;	the help section
hilfe=['',$
'To interactively tweak graph annotations',$
'']

;	this is how the widget IDs are transferred to KALPANA_EVENT
widg={base:base, col1:base_col1, col2:base_col2, menu:menu,$
	txts:txts, wind:wind, plot:plot, labl:labl, tweak:tweak,$
	help:hilfe}

;	start up the widgets
widget_control,base,/realize,/hourglass

;register with XMANAGER
;	this is the officially recommended way, but we're not going to do it
;	because there are too many parameters and keywords (not to mention
;	the _EXTRA keywords) to pass to the event handling routine.
;		xmanager,'fitlines',base
;therefore,
;	we handle the widget events in a separate subourtine, partly also
;	to keep things easy to read

kalpana_event,x,y,ststr,widg,oststr=oststr,wshow=wshow,wwork=wwork, _extra=e

;	clean up
widget_control,base,/destroy

return
end

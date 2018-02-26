function lineid_menu,menu,w0,w,f,wshift,topX,wrange,cord,code,scratch,$
	wvls=wvls,idwvls=idwvls,gwvl=gwvl,renorm=renorm
;+
;function	lineid_menu
;	function to select matches to wavelengths and return position
;	indices of selected matches.
;
;syntax
;	oo=lineid_menu(menu,w0,w,f,wshift,topX,wrange,cord,code,scratch,$
;	wvls=wvls,idwvls=idwvls,gwvl=gwvl,renorm=renorm)
;
;parameters
;	menu	[INPUT] list of items to choose from; string array
;	w0	[INPUT] wavelength that needs match
;	w	[INPUT] wvls of matched lines (must be of same size as MENU)
;	f	[INPUT] fluxes of all matched lines (same size as MENU)
;	wshift	[I/O] wavelength shift
;	topX	[I/O] maximum number of matches
;	wrange	[I/O] interval in which to search for matches
;	cord	[I/O] comma-separated list of grating orders
;	code	[OUTPUT] 0: NEXT WAVELENGTH
;			 1: QUIT LINEID
;			 2: REPEAT FOR NEW PARAMETERS
;	scratch	[OUTPUT] human readable notation
;
;keywords
;	wvls	[INPUT] wavelengths which are previously identified
;	idwvls	[INPUT] structure describing IDs 
;	gwvl	[INPUT] grating-order corrected wavelengths
;	renorm	[INPUT] renormalization factor
;
;restrictions
;	* widget based
;	* dedicated subroutine to LINEID (written simply to avoid
;	  making LINEID too large to comprehend)
;	* no error checking is done.  it is assumed that this function
;	  is called from LINEID, and naught else.
;	* requires subroutines
;	  -- CAT_ID
;	  -- CREATE_STRUCT
;
;warning
;	* here be hippogrif
;	* don't even <<think>> about understanding this program
;
;history
;	author of this hack prefers anonymity, thank you very much (Dec96)
;	some cosmetic surgery, added button STATUS and call to CAT_ID (Feb97)
;	corrected bug with reordering (Jul98)
;	added CORD, GWVL, changed TOP field from /LONG to /FLOAT (Nov98)
;	changed behavior of UNDO LAST; renamed NEXT as HELP; renamed DONE
;	  as NEXT; added keyword RENORM and entry field PARS_NORM; changed
;	  MBAR_NEXT to MBAR_HELP (Dec98)
;	added output SCRATCH (Mar99)
;	bug correction re SCRATCH (FebMM)
;-

if n_params(0) lt 8 then begin
  print,'Usage: oo=lineid_menu(menu,w0,wvl,flx,wshift,topX,wrange,code,scratch,$'
  print,'  wvls=wvls,idwvls=idwvls,gwvl=gwvl,renorm=renorm)
  return,-1L
endif

;	save inputs
cc=string(menu)					;input items
nn=n_elements(menu) & oh=lonarr(nn) & oldoh=oh
oo=lindgen(nn)					;initial order
if not keyword_set(gwvl) then gwvl=w		;error check
if not keyword_set(renorm) then renorm=0.	;error check
ow=sort(abs(gwvl-w0)) & kw=0			;sort by wavelength difference
oa=sort(gwvl) & ka=0				;sort by wavelength
oi=sort(abs(f)) & oi=reverse(oi) & ki=0		;sort by intensity
ws=wshift & top=topX & wrl=wrange(0) & wrr=wrange(1)
code=0 & kpars=0 & last1=0 & last2=nn-1

;	initialize scribbles
scratch=''

;	help array
lineid_help=[ $
'Pick IDs by clicking on the buttons next to each candidate line',$
'The columns are:',$
'        Wavelength, Atom, Ionic state, Tmax, Trange',$
'        optionally: level designation, e-configuration, matched wavelength',$
'',$
'NEXT: finish with IDing this feature and go to next one',$
'',$
'SELECT: select/unselect large numbers of candidates at a time',$
'-> ALL: Select all candidates as IDs',$
'-> Between Last 2: (not very robust)',$
'-> Undo Last: Undo previous SELECT action',$
'-> Reset All: Unselect all candidates',$
'',$
'SORT: Sort listed candidates in various orders',$
'-> By Wavelength difference: in increasing distance from feature',$
'-> By Intensity: in decreasing order of predicted flux',$
'-> By Wavelength: in increasing order of wavelength',$
'-> Default: restore original order',$
'',$
'UNKNOWN: Mark this feature as having an "Unknown" ID',$
'  *WARNING* This will result in discarding all other IDs for this feature',$
'',$
'HELP: as it says',$
'',$
'STATUS: print list of features that have been IDd',$
'',$
'QUIT: quit LINEID',$
'',$
'APPLY: Apply changes in parameters to generate new candidate list',$
'  *WARNING* All existing IDs for this feature will be lost!',$
'',$
'PARAMETERS: Click to hide/display global control parameters',$
'== Wavelength Shift: shift wavelength scale',$
'== Maximum matches: sets max number of candidates',$
'   (same as keyword TOPX of LINEID)',$
'   Some tricks: set to 0<val<1 for it to act as a relative threshold,',$
'                    i.e., include all lines with fluxes > val*max(flux)',$
'                set to -NN.FF for it to return at most NN matches with',$
'                    a threshold of 0.FF',$
'== Grating Orders: to include higher-order lines as candidates',$
'   (comma-separated string of integers)',$
'== Delta[Wvl] (Left): left-side range in which to search for candidates',$
'== Delta[Wvl] (Right): right-side range in which to search for candidates',$
'   (same as keyword WRNG of LINEID)',$
'== RENORM to: set to +ve number to normalize the strongest line among the',$
'   candidates to this number; to -1 to normalize said line to value of the',$
'   spectrum at this wavelength; to -(#.NE.1) to set multiplicative factor',$
'   for renormalization.',$
'   If -1 is set, "-1*multiplicative_factor" is echoed, so to recover the',$
'   theoretical fluxes, use -1/echoed_value',$
'' ]

;	widget layout
;ID selection commands window
;ID selections list
;parameters (unmapped)
;menu bar of global commands

c1='ID LINE @ '+strtrim(w0,2)+' A'
base=widget_base(frame=2,/column,/map,title=c1,uvalue='base')

;	menu bar for IDs
base_idnt=widget_base(base,frame=2,/row,uvalue='base_idnt')
;
idnt_done=widget_button(base_idnt,value='NEXT',uvalue='idnt_done')
;
idnt_slct=widget_button(base_idnt,value='SELECT',/menu,uvalue='idnt_slct')
slct_opt1=widget_button(idnt_slct,value='All',uvalue='slct_opt1')
slct_opt2=widget_button(idnt_slct,value='Between Last 2',uvalue='slct_opt2')
slct_opt3=widget_button(idnt_slct,value='Undo Last',uvalue='slct_opt3')
slct_opt4=widget_button(idnt_slct,value='Reset All',uvalue='slct_opt4')
;
idnt_sort=widget_button(base_idnt,value='SORT',/menu,uvalue='idnt_sort')
sort_wvls=widget_button(idnt_sort,value='by Wavelength difference',$
	uvalue='sort_wvls')
sort_ints=widget_button(idnt_sort,value='by Intensity',uvalue='sort_ints')
sort_angs=widget_button(idnt_sort,value='by Wavelength',uvalue='sort_angs')
sort_defs=widget_button(idnt_sort,value='Default',uvalue='sort_defs')
;
idnt_unkn=widget_button(base_idnt,value='UNKNOWN',uvalue='idnt_unkn')
;
idnt_note=widget_button(base_idnt,value='Add COMMENT',uvalue='idnt_note')

;	all the matches
base_mtch=widget_base(base,frame=2,/column,uvalue='base_mtch')
label_top='{?id?} : [*] WVL : Z, ION : Fx : Tmax <Trange> : [description] : [match]'
mtch_item=cw_bgroup(base_mtch,cc,/column,/nonexclusive,uvalue='mtch_item',$
	/scroll,y_scroll_size=300,x_scroll_size=700,set_value=oh,$
	label_top=label_top)

;	menu bar
base_mbar=widget_base(base,frame=2,/row,/map,uvalue='base_mbar')
mbar_help=widget_button(base_mbar,value='HELP',uvalue='mbar_help')
mbar_pars=widget_button(base_mbar,value='PARAMETERS',uvalue='mbar_pars')
mbar_stus=widget_button(base_mbar,value='STATUS',uvalue='mbar_stus')
mbar_quit=widget_button(base_mbar,value='QUIT',uvalue='mbar_quit')

;	global parameters
base_parm=widget_base(base,frame=2,map=0,/row,uvalue='base_parm')
;base_parm=widget_base(base,frame=2,/row,uvalue='base_parm')
par1_pars=widget_base(base_parm,/column,uvalue='par1_pars')
pars_wsft=cw_field(par1_pars,/float,/return_events,title='Wavelength Shift',$
	/row,value=ws,uvalue='pars_wsft')
pars_maxm=cw_field(par1_pars,/float,/return_events,title='Maximum Matches',$
	/row,value=top,uvalue='pars_maxm')
pars_ordr=cw_field(par1_pars,/return_events,title='Grating Orders',$
	/row,value=cord,uvalue='pars_ordr')
par2_pars=widget_base(base_parm,/column,uvalue='par2_pars')
pars_wr00=cw_field(par2_pars,/float,/return_events,$
	/row,title='Delta[Wvl] (Left)',value=wrl,uvalue='pars_wr00')
pars_wr01=cw_field(par2_pars,/float,/return_events,$
	/row,title='Delta[Wvl] (Right)',value=wrr,uvalue='pars_wr01')
pars_norm=cw_field(par2_pars,/float,/return_events,$
	/row,title='RENORM to',value=renorm,uvalue='pars_norm')
pars_appl=widget_button(base_parm,value='APPLY',uvalue='pars_appl')

;	start up widgets
widget_control,base,/realize,/hourglass

ok=''
while ok ne 'quit' do begin
  event=widget_event(base) & widget_control,event.id,get_uvalue=ev
  case ev of
    ;
    ;	menu bar
    ;
    'mbar_help': begin
      nhelp=n_elements(lineid_help)
      if fix(strmid(!version.release,0,1)) lt 5 then $
	for i=0,nhelp-1 do print,lineid_help(i) else xdisplayfile,$
	'help',text=lineid_help,title='LINEID HELP',width=75
      ;widget_control,mtch_item,get_value=oh		;read selections
      ;widget_control,pars_wsft,get_value=ws		;read wvl shift
      ;widget_control,pars_maxm,get_value=top		;read max matches
      ;widget_control,pars_ordr,get_value=cord		;read grating orders
      ;widget_control,pars_wr00,get_value=wrl		;read wvl range left
      ;widget_control,pars_wr01,get_value=wrr		;read wvl range right
      ;widget_control,pars_norm,get_value=renorm	;renormalize fluxes
      ;ok='quit'					;out
    end
    'mbar_pars': begin
      if kpars eq 0 then begin				;open parameter input
	widget_control,base_parm,/map & kpars=1
      endif else begin					;close parameter input
	widget_control,base_parm,map=0 & kpars=0
      endelse
    end
    'mbar_stus': begin					;show current IDs
      stus='ok'
      if not keyword_set(wvls) then stus='nope'		;no matched wvls
      if not keyword_set(idwvls) then stus='nope'	;no IDs
      if n_tags(idwvls) ne n_elements(wvls) then stus='nope'	;incorrect
      if stus eq 'ok' then begin			;all A-OK
	tmp=cat_id(create_struct('WVL',wvls,idwvls))	;print out
      endif else print,'None identified yet'
    end
    'mbar_quit': begin
      widget_control,mtch_item,get_value=oh		;read selections
      widget_control,pars_wsft,get_value=ws		;read wvl shift
      widget_control,pars_maxm,get_value=top		;read max matches
      widget_control,pars_ordr,get_value=cord		;read grating orders
      widget_control,pars_wr00,get_value=wrl		;read wvl range left
      widget_control,pars_wr01,get_value=wrr		;read wvl range right
      widget_control,pars_norm,get_value=renorm		;renormalize fluxes
      ok='quit' & code=1				;out
    end
    'idnt_note': begin					;add comment
      print,string("7b)
      print,'type in comment.'
      print,'  type "SHOW" or "?" to print'
      print,'  end with BLANK LINE, "EOF", "QUIT", or "DONE":'
      print,'  type "CLEAR" to start over'
      c='!EOF' & c1='' & k=0L
      while c ne 'EOF' do begin
	read,c1
	if strtrim(c1,2) eq '' then c='EOF'
	c2=strlowcase(strtrim(c1,2))
	if c2 eq 'eof' or c2 eq 'quit' or c2 eq 'done' then c='EOF'
	if c2 eq 'clear' then scratch=''
	if c2 eq 'show' or c2 eq '?' then begin
	  c='EOF' & if keyword_set(scratch) then print,scratch
	endif
	if c ne 'EOF' then begin
	  if keyword_set(scratch) then scratch=scratch+$
		string("12b)+c1 else scratch=c1
	endif
      endwhile
      widget_control,mtch_item,get_value=oh	;read selections
      oo=where(oh gt 0,moo)
      if moo eq 0 and keyword_set(scratch) then begin
	print,scratch
	print,string("7b)+'WARNING:'
	print,'comment will be ignored unless this feature is IDd'
      endif
    end
    ;
    ;	parameters
    ;
    'pars_appl': begin
      widget_control,pars_wsft,get_value=ws		;read wvl shift
      widget_control,pars_maxm,get_value=top		;read max matches
      widget_control,pars_ordr,get_value=cord		;read grating orders
      widget_control,pars_wr00,get_value=wrl		;read wvl range left
      widget_control,pars_wr01,get_value=wrr		;read wvl range right
      widget_control,pars_norm,get_value=renorm		;read renormalization
      ok='quit' & code=2				;out
      widget_control,base_parm,map=0 & kpars=0
    end
    'pars_wsft': widget_control,pars_wsft,get_value=ws	;read wvl shift
    'pars_maxm': widget_control,pars_maxm,get_value=top	;read max matches
    'pars_ordr': widget_control,pars_ordr,get_value=cord ;read grating orders
    'pars_wr00': widget_control,pars_wr00,get_value=wrl	;read wvl range left
    'pars_wr01': widget_control,pars_wr01,get_value=wrr	;read wvl range right
    'pars_norm': widget_control,pars_norm,get_value=renorm ;renormalize fluxes
    ;
    ;	ID parameters
    ;
    'idnt_done': begin					;done!
      widget_control,mtch_item,get_value=oh		;read selections
      ok='quit' & code=0
    end
    'idnt_unkn': begin
      widget_control,mtch_item,set_value=lonarr(nn)
      ok='quit' & code=3				;unknown ID
      wait,1
    end
    'slct_opt1': begin					;select ALL
      oldoh=oh & oh=oh>1
      widget_control,mtch_item,set_value=oh
    end
    'slct_opt2': begin					;all between last 2
      oldoh=oh
      if last2 gt last1 then oh(last1:last2)=1
      if last2 lt last1 then oh(last2:last1)=1
      widget_control,mtch_item,set_value=oh
    end
    'slct_opt3': begin					;undo last
      ;oh=oldoh
      widget_control,mtch_item,set_value=oldoh
      oldoh=oh
    end
    'slct_opt4': begin					;reset all to 0
      oh(*)=0 & widget_control,mtch_item,set_value=oh
    end
    'sort_wvls': begin					;sort by wvl. diffs
      if kw eq 0 then begin
        widget_control,mtch_item,get_value=oh		;read selections
        message,'	please wait...',/info
        if ki eq 1 then begin				;reorder
	  oh=oh(sort(oi)) & oldoh=oldoh(sort(oi))
        endif
        if ka eq 1 then begin				;reorder
          oh=oh(sort(oa)) & oldoh=oldoh(sort(oa))
        endif
        kw=1 & ki=0 & ka=0
        oldoh=oldoh(ow) & oh=oh(ow)			;rearrange
        widget_control,mtch_item,/destroy		;destroy
	widget_control,event.top,update=0
        mtch_item=cw_bgroup(base_mtch,cc(ow),/column,/nonexclusive,$
	  uvalue='mtch_item',/scroll,y_scroll_size=300,x_scroll_size=700,$
	  map=0,set_value=oh,label_top=label_top)
	widget_control,event.top,/update
        widget_control,mtch_item,/realize,/map,/hour	;recreate
      endif
    end
    'sort_ints': begin					;sort by intensity
      if ki eq 0 then begin
        message,'	please wait...',/info
        widget_control,mtch_item,get_value=oh		;read selections
        if kw eq 1 then begin				;reorder
	  oh=oh(sort(ow)) & oldoh=oldoh(sort(ow))
        endif
        if ka eq 1 then begin				;reorder
          oh=oh(sort(oa)) & oldoh=oldoh(sort(oa))
        endif
        kw=0 & ki=1 & ka=0
        oldoh=oldoh(oi) & oh=oh(oi)			;rearrange
	widget_control,event.top,update=0
        widget_control,mtch_item,/destroy		;destroy
        mtch_item=cw_bgroup(base_mtch,cc(oi),/column,/nonexclusive,$
	  uvalue='mtch_item',/scroll,y_scroll_size=300,x_scroll_size=700,$
	  map=0,set_value=oh,label_top=label_top)
        widget_control,mtch_item,/realize,/map,/hour	;recreate
	widget_control,event.top,/update
      endif
    end
    'sort_angs': begin					;sort by wavelength
      if ka eq 0 then begin
        message,'	please wait...',/info
        widget_control,mtch_item,get_value=oh		;read selections
        if kw eq 1 then begin				;reorder
	  oh=oh(sort(ow)) & oldoh=oldoh(sort(ow))
        endif
        if ki eq 1 then begin				;reorder
          oh=oh(sort(oi)) & oldoh=oldoh(sort(oi))
        endif
        kw=0 & ki=0 & ka=1
        oldoh=oldoh(oa) & oh=oh(oa)			;rearrange
	widget_control,event.top,update=0
        widget_control,mtch_item,/destroy		;destroy
        mtch_item=cw_bgroup(base_mtch,cc(oa),/column,/nonexclusive,$
	  uvalue='mtch_item',/scroll,y_scroll_size=300,x_scroll_size=700,$
	  map=0,set_value=oh,label_top=label_top)
        widget_control,mtch_item,/realize,/map,/hour	;recreate
	widget_control,event.top,/update
      endif
    end
    'sort_defs': begin
      if ki eq 1 or kw eq 1 or ka eq 1 then begin
        message,'	please wait...',/info
        widget_control,mtch_item,get_value=oh		;read selections
        if ki eq 1 then begin				;reorder
	  oh=oh(sort(oi)) & oldoh=oldoh(sort(oi))
        endif
        if kw eq 1 then begin				;reorder
          oh=oh(sort(ow)) & oldoh=oldoh(sort(ow))
        endif
        if ka eq 1 then begin				;reorder
          oh=oh(sort(oa)) & oldoh=oldoh(sort(oa))
        endif
	widget_control,event.top,update=0
	widget_control,mtch_item,/destroy		;destroy
	mtch_item=cw_bgroup(base_mtch,cc,/column,/nonexclusive,$
          uvalue='mtch_item',/scroll,y_scroll_size=300,x_scroll_size=700,$
          map=0,set_value=oh,label_top=label_top)
	widget_control,mtch_item,/realize,/map,/hour	;recreate
	widget_control,event.top,/update
      endif
      ki=0 & kw=0 & ka=0
    end
    'mtch_item': begin					;update
      last2=last1 & last1=event.value
    end
    else: ;nothing
  endcase
endwhile

;	ok, lesse what we got
if kw eq 1 then oh=oh(sort(ow))			;sorted by dwvl -- reorder
if ki eq 1 then oh=oh(sort(oi))			;sorted by flx -- reorder
if ka eq 1 then oh=oh(sort(oa))			;sorted by wvl -- reorder
wshift=ws & topX=top & wrange=[wrl,wrr]
widget_control,base,/destroy

oo=where(oh gt 0)

return,oo
end

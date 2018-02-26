restore,'plotid.sav'
x=lamda & y=spec & lid=line & loadct,3
plotid,lid,x,y,xpos=xp,ypos=yp,hylyt=1,style=sty
xp2=xp+randomn(seed,4)*10 & yp2=yp+randomn(seed,4)*20

plotid,lid,x,y,xpos=xp2,ypos=yp2,style=sty,/finger

end

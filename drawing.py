import math
from tkinter import *
from Biologic import get_bonds

# Startposition für das Zeichnen
sx = 600
sy = 500

# zeichne eine Base
def drawBase(canvas,x,y,b):
    canvas.create_oval(x-7,y-7,x+7,y+7, width = 1)
    canvas.create_text(x,y,text = str(b))

# zeichne eine Linie
def drawLine(canvas,lx, ly, rx, ry, n):
	canvas.create_line(0.8*lx+0.2*rx,0.8*ly+0.2*ry,0.2*lx+0.8*rx,0.2*ly+0.8*ry, width=n)

# zeichne eine gegebene RNA nach gegebener Struktur (rekursive Funktion)
def drawRna(canvas, rna, struct, l, r, lx, ly, rx, ry):
	level = 0
	count = 2.2

	for i in range(l+1,r):
		if struct[i] == ')':
				level = level-1
		if level == 0:
				count = count+1
		if struct[i] == '(':
				level = level+1

	f = 25
	theta = 2 * 3.14159 / count
	rad = f / (2*  math.sin(theta/2))
	h = rad * math.cos(theta/2)
	le = math.sqrt((ly-ry)*(ly-ry) + (lx-rx)*(lx-rx))
	cx = int((lx+rx)/2 + h*((ly-ry)/le))
	cy = int((ly+ry)/2) + int(h*((rx-lx)/le))
	alpha = math.atan2(ly-cy, lx-cx)

	ii = l
	xx = lx
	yy = ly

	for i in range(l+1,r+1):
		if struct[i] == ')':
				level = level-1
		if level == 0:
				alpha = alpha-theta
				x = int(cx + rad*math.cos(alpha))
				y = int(cy + rad*math.sin(alpha))
				if struct[i] == ')':
						drawRna(canvas, rna, struct, ii, i, xx, yy, x, y)
				else:
						canvas.create_line(0.8*xx + 0.2*x, 0.8*yy+0.2*y, 0.2*xx+0.8*x, 0.2*yy + 0.8*y)
				drawBase(canvas, xx,yy,rna[ii])
				xx = x
				yy = y
				ii = i
		if struct[i] == '(':
			level = level+1

	canvas.create_line(0.8*xx+0.2*rx,0.8*yy+0.2*ry,0.2*xx+0.8*rx,0.2*yy+0.8*ry)
	drawBase(canvas,xx,yy,rna[r-1])
	drawBase(canvas,rx, ry, rna[r])

	#if not (lx == sx or ly == sy):
	drawLine(canvas, lx, ly, rx, ry, get_bonds().get(str(rna[l])+str(rna[r])))

# zeichne eine gegebene RNA nach gegebener Struktur (Einstiegsfunktion)
def drawRnaStructure(rna, strukt, canvas):

    rna = '5'+rna+'3'
    struct = '('+strukt+')'
    canvas.delete(ALL)
    print(struct)
    drawRna(canvas, rna, struct, 0, len(struct)-1, sx, sy, sx+20,sy)
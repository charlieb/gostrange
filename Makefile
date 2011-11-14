include $(GOROOT)/src/Make.inc

TARG=strange
GOFILES= main.go matrix.go

INCLUDES=-I../Go-SDL/sdl/_obj \
				 -I../Go-SDL/mixer/_obj \
				 -I../Go-SDL/ttf/_obj \
				 -I../Go-SDL/gfx/_obj \
				 -I../Go-OpenGL/gl/_obj \
				 -I../Go-OpenGL/glu/_obj

LIBS=-L../Go-SDL/sdl/_obj \
		 -L../Go-SDL/mixer/_obj \
		 -L../Go-SDL/ttf/_obj \
		 -L../Go-SDL/gfx/_obj \
		 -L../Go-OpenGL/gl/_obj \
		 -L../Go-OpenGL/glu/_obj

GC=${O}g $(INCLUDES)
LD=${O}l $(LIBS)

include $(GOROOT)/src/Make.cmd
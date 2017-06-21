# glitch-finder
John Arthur's glitch finder program

Here is the fortran source code for my glitch finder program.  It may look
a little overly complicated, because I made it from pieces of code that I
already had, and wanted to use a standard input data format for all my
crystal programs, etc, etc…

But I have simplified some of the subroutines and attached them, so that
this is a standalone version.

Attached, find:
 glitch_SLAC.f   The main program, with some subroutines attached at the
end
 mtrx_SLAC.f     Some additional subroutines that are needed
 glitch.dat      The input data file (a text file)

Not included: 
 The output graphing routine.  I use the PGPLOT graphics library — if you
have that then just link to it and this will run.  Otherwise, you can use
you own graphics package, but you will need to rewrite the subroutine
GRAPH that is attached at the end of GLITCH_SLAC.f

This is written in very plain vanilla f77 fortran.  It should compile with
just about any fortran compiler.  But I guess you will want to port it to
Python.  In that case, refer to the comments within the source files.  And
feel free to contact me if you have any questions.

Best regards,
John

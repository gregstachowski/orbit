# orbit
Orbit determination from  three  sets of observations (t, R.A., dec.)
assuming heliocentric orbit and geocentric observer (i.e. requires R)
and using Laplacian methods. See "Astrodynamics: Orbit Determination,
Space Navigation,  Celestial Mechanics,  Volume 1" by Samuel Herrick,
Van Nostrand Reinhold, 1971;  esp.  Ch. 10 and 12 (page numbers refer
to the appropriate pages in this edition).

Requires an ASCII datafile, eight (8) lines:

   comment
   header
   t(1)    R.A.(1)   dec(1)  Rx(1)   Ry(1)   Rz(1)
   t(1)    R.A.(1)   dec(1)  Rx(2)   Ry(2)   Rz(2)
   t(1)    R.A.(1)   dec(1)  Rx(3)   Ry(3)   Rz(3)
   k*
   alpha
   epsilon

where t is in JD and R.A., dec. are both in decimal degrees,

and k* = k = 1/(unit of canonical time) :

     k = 0.017 202 09895 if canonical time in days
     k = 0.001 239 444   if canonical time in seconds

note: tau = k*(t - to).

and alpha is the light-time constant:

     alpha = 0.0057755 days/au or alpha = 0.021275 sec/gu

epsilon is the obliquity (inclination) of the ecliptic plane,
in decimal degrees.

GSS, 12 Mar 1998, 3 Jun 1998, 7-9 Oct 1998 


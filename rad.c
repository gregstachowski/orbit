/* <html>
   <title>rad.c - (c) G.S.Stachowski, 1998</title>
   <body bgcolor="white" fgcolor="black"><p>
   <b><font color="darkred">rad.c</font> - (c) <a href="mailto:greg@froggie.freeservers.com">G.S.Stachowski</a>, 1998</b><br>
   Functions to do <font color="darkred">degree -> radian</font> and <font color="darkred">radian -> degree</font> conversion.<p>
   <i>All HTML in this program is contained within </i> C <i> comment markers.  
   If you either copy-paste from the screen, &quot;save as text&quot; or even 
   save the HTML source, you will have a working </i>C<i> program source code.</i>
   <xmp>*/


#include <math.h>

double rad(double val) {	/* degree -> radian conversion */
	return(val*M_PI/180.0);
}

double deg(double val) {        /* radian -> degree conversion */
        return(val*180.0/M_PI);
}


/*</xmp><! ------------------------------------------------------------------>
<hr><img align=right border=0 src="../../graphics/fish-s.jpg"></body></html>*/


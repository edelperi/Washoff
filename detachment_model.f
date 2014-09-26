* ====================================================================
*                          DETACHMENT MODEL
* ====================================================================
*  This is the model to calculate the pesticide detachment by rainfall
*  based on the water interception by canopy (Calder, 1986)
*  Here is implemented a two-phase loss
*  The first stage accounts for the loss of loose attached pesticide
*  The second stage accounts for the loss of tight attached pesticide
*  Note lun jul 29 18:38:53 CEST 2013
*  This version is for fiting multiple sets namely particulate and soluble
*  losses by using the integer tag data_type(i_data)  
*  See formatting in Pages 7-8in the LM-OPT manual.
*  Note:
*  mar abr 22 13:24:43 CEST 2014
*  Aditional first order detachment process included.
* 
*  This subroutine must be implemented with the lm_opt optimizacion
*  algorithm  (Clausnitzer, V.  and Hopmans, J.W. 1995)
* ====================================================================
*  José Eugenio López Periago
*  University of Vigo 2013
* GNU GENERAL PUBLIC LICENSE
 *                      Version 2, June 1991

*Copyright (C) 1989, 1991 Free Software Foundation, Inc., <http://fsf.org/>
*51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*Everyone is permitted to copy and distribute verbatim copies
*of this license document, but changing it is not allowed.
* ====================================================================
      INTEGER FUNCTION model(n_par,par,
     &                       n_data,abscissa,y_calc,data_type)

********************************************************************
* Definition of variables
* these should not be modified:

      REAL*8 par(*),abscissa(*),y_calc(*),parm,temporal
      REAL factorial,effn
      INTEGER n_par,n_data,data_type(*),N,R,ISET
      REAL*8 tmp
*
********************************************************************
*     ===========================================
*              Key of model parameters 
*     ===========================================
*     par(1) number of elements in the surface parameter (L)
*     par(2) number of effective drop impacts parameter (q) (stage 1)
*     par(3) detachment efficiency per impact (stage 1)
*     par(4) number of effective drop impacts (stage 2)
*     par(5) detachment efficiency per impact (stage 2)
*     par(6) first order detachment factor (stage 1)
*     par(7) first order detachment factor (stage 2)
*     ===========================================
*
*     Read the experimental data
      DO 1 i_data=1,n_data

      ISET = data_type(i_data)ISET = data_type(i_data)

      temporal = 0.0
      factorial = 1.0
       parm= abscissa(i_data)/par(1)
      if (parm .gt. (par(2*ISET)+1)) then
       R= INT(par(2*ISET))
       else
       R = INT(parm)
       endif
*     Using 16 terms would be enough for a good approximation for the
*     Asymptotic term as follows in the DO 10  LOOP
      if ( R .gt. 1) THEN
      Do 10 N=1,R
      factorial = factorial*real(N)
        temporal = temporal
     +   +(real(N) - par(2*ISET)) * (parm**real(N))/factorial
   10 CONTINUE
      ENDIF
        effn  =
     +  (par(2*ISET) * ( 1.0 - exp(-parm) ) + exp(-parm)* temporal)
*
        y_calc(i_data) = par(1) * effn * par(2*ISET+1) *
     +  exp(-effn*par(5+ISET))
    1 CONTINUE
      model=0
      RETURN
      END
* =================================================================
*                     END OF DETACHMENT MODEL
* =================================================================

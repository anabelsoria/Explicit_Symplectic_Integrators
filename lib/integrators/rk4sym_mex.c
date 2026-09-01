/* Name:     rk4sym_mex()
Purpose:     Fixed-step 4th-order, symplectic  Runge-Kutta initial value 
             ODE solver.
             Converts MATLAB inputs into C inputs for the rk4sym.c
             integrator.  Calls the rk4sym() integrator and returns
             outputs to MATLAB.
             TurboProp Version 4.0.

Call:        [T,Y] = rk4sym_mex(odefun,tspan,y0,options,extras)

Inputs:      odefun  = string with the name of the derivative function 
                       (e.g. 'TwoBody')
             tspan   = vector of times for which an integrated 
                       state will be output.  The intervals must be 
                       evenly divisible by
                       dt, except for the last interval.
             y0      = the initial state in a vector.
             options = Array of two doubles.
                       options(1) contains dt, the time stepsize.
                       options(2) contains a tolerance used to verify
                       the fixed point iteration required for the implicit
                       Runge-Kutta scheme.
             extras  = a structure data type sent on to the derivative
                       function.  See the derivative function help files
                       for more information.
 
Outputs:     T       = a column vector of times for which the integrated
                       state is output.
             Y       = a matrix where row i is the state at the time 
                       given in T(i).    

Notes:       

References:  

Subroutines: User derivative function specified by the string odefun.
             For a list of compiled derivative functions, type
                help odefun
             FillExtras()

Examples:    

Author:      Brandon Jones

Revisions History:
7/20/05:     Initial Submission, see svn log for further changes
*/

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <matrix.h>
#include "TurboTypes.h"
#include "findptr.h"
#include "FillExtras.h"
#include "initExtras.h"
 
/* -------------------------------------------------------- */   
/* #define the name of the C integrator function: */ 
#define SOLVER rk4sym
#define SOLVERSTR "rk4sym"

/* -------------------------------------------------------- */   
/* define the number of required options for the integrator: */ 
/* If options aren't used by the integrator, put optsize.    */
#define NUMOPT 2

/* The C integrator function prototype: */ 
int SOLVER(double *T, double *Y, derivFcn_T funcptr,
      double *tspan,int tspansize, double *y0,
      int statesize, double *options, int optsize,
      extrastype *extras);
      
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
   double *y0, *Y, *tspan, *options, *T;
   double  dt, itnum, frac;
   char *odefun;
   int tspansize, statesize, optsize, buflen, status, i;
   derivFcn_T funcptr = 0;
   char solvername[] = SOLVERSTR;
   char errmsg[120];
   extrastype extras;
   int fieldnum;
   mxArray *fieldptr;
   double *arrayptr;

/* Check for proper number of arguments. */
   if (nrhs != 5) {
      sprintf(errmsg,"\n%s error: Inproper number of inputs.\n",solvername);
      mexErrMsgTxt(errmsg);
   }  else if (nlhs != 2) {
      sprintf(errmsg,"\n%s error: Inproper number of outputs.\n",solvername);
      mexErrMsgTxt(errmsg);
   }
  
/* Check the size, type of the odefun input and allocate memory */
   if (mxIsChar(prhs[0]) != 1 || mxGetM(prhs[0]) != 1)
   {
      sprintf(errmsg,"\n%s error: odefun must be a string.\n",solvername);
      mexErrMsgTxt(errmsg);
   }
   buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
   odefun = mxCalloc(buflen, sizeof(char));
   status = mxGetString(prhs[0], odefun, buflen);
   if (status != 0) 
   {
      sprintf(errmsg,"\n%s error: odefun string did not load correctly.\n",solvername);
      mexWarnMsgTxt(errmsg);
   }
  
/* Check size, type of tspan */
   if (mxGetNumberOfDimensions(prhs[1])>2 || 
     (mxGetN(prhs[1])>1 && mxGetM(prhs[1])>1) 
     || !(mxGetNumberOfElements(prhs[1])>1)
     || mxIsComplex(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      sprintf(errmsg,"\n%s error: tspan should be a vector of real doubles.\n",solvername);
      mexErrMsgTxt(errmsg);
   }
   tspansize = mxGetNumberOfElements(prhs[1]);

/* Check size, type of y0 */
   if (mxGetNumberOfDimensions(prhs[2])>2 || 
      (mxGetN(prhs[2])>1 && mxGetM(prhs[2])>1) 
     || mxIsComplex(prhs[2]) || !mxIsDouble(prhs[2]))
   {
       sprintf(errmsg,"\n%s error: y0 should be a vector of real doubles.\n",solvername);
       mexErrMsgTxt(errmsg);
   }
   statesize = mxGetNumberOfElements(prhs[2]);
  
/* Check size, type of options */
   if (mxGetNumberOfDimensions(prhs[3])>2 || 
      (mxGetN(prhs[3])>1 && mxGetM(prhs[3])>1)
     || mxIsComplex(prhs[3]) || !mxIsDouble(prhs[3]))
   {
       sprintf(errmsg,"\n%s error: options should be a vector of real doubles.\n",solvername);
       mexErrMsgTxt(errmsg);
   }
   optsize = mxGetNumberOfElements(prhs[3]);
   if (optsize != NUMOPT)
   {
       sprintf(errmsg,"\n%s error: options is not the correct size.\n",solvername);
       mexErrMsgTxt(errmsg);
   }
   
/* Check size, type of extras structure*/  
   if (!mxIsStruct(prhs[4]) || mxGetNumberOfElements(prhs[4])>1)
   {
       sprintf(errmsg,"\n%s error: extras should be a scalar structure.\n",solvername);
       mexErrMsgTxt(errmsg);
   }   

/* Call initExtras to initialize the extras structure */
   initExtras(&extras);   
   
/* Call FillExtras to load the extras structure */
   FillExtras(&extras,prhs[4]);

/* Assign pointers to each input. */
   tspan = (double *)mxGetPr(prhs[1]);
   y0 = (double *)mxGetPr(prhs[2]);
   options = (double *)mxGetPr(prhs[3]);

/* Check the tspan for strict increase or decrease: */   
   dt = (tspan[1] - tspan[0]);
   
   for (i=2; i<tspansize; i++)
   {
      if ((tspan[i] - tspan[i-1])/dt<=0.0)
      {
          sprintf(errmsg,"\n%s error: tspan should strictly increase or decrease.\n",solvername);
          mexErrMsgTxt(errmsg);
      }
   }
   
/*check tspan to make sure dt evenly divides the intervals:*/   
   dt = options[0];
   for (i=0; i<tspansize-2; i++)
   {
      /*find number of dt's in the interval:*/
      frac = modf(((tspan[i+1] - tspan[i])/dt),&itnum);
      /*check to see if dt evenly divides the interval:*/
      /*round off number of dt's in the interval:*/
      if (frac>0.5)
         itnum++;
      if (fabs(itnum-(tspan[i+1] - tspan[i])/dt) > options[1])
      {
         sprintf(errmsg,"\n%s error: %s\n%s (err = %.2e, tol = %.2e)\n",solvername,
         "An interval in tspan, which is not the last, is not",
         "evenly divisible by dt.",fabs(itnum-(tspan[i+1] - tspan[i])/dt),options[1]);
         mexErrMsgTxt(errmsg);
      }
   }
   
/* Create the output matrices. */
   plhs[0] = mxCreateDoubleMatrix(tspansize, 1, mxREAL);
   plhs[1] = mxCreateDoubleMatrix(tspansize, statesize, mxREAL);

/* Assign pointers to each output. */
   T = mxGetPr(plhs[0]);
   Y = mxGetPr(plhs[1]);
  
/* Find the function pointer corresponding to odefun */
   funcptr = findptr(odefun,tspan,tspansize,y0,statesize,solvername,1,&extras);

/* Call the integrate subroutine. */
   status = SOLVER(T, Y, funcptr, tspan, tspansize, y0, statesize, options, 
       optsize, &extras);
   
/* Use findptr for any termination instructions: */
   funcptr = findptr(odefun,tspan,tspansize,y0,statesize,solvername,0,&extras);

   if (status==TP_FAILURE)
   {
      sprintf(errmsg,"\n%s error: integration did not complete properly.\n",solvername);
      mexErrMsgTxt(errmsg);
   }
   
   return;
   
} /* mexFunction() */

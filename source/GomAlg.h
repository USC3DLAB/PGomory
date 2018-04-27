//
//  LGomFixedRecourse.h
//  LGomFixedRecourse
//
//  Created by Dinakar Gade on 9/8/12.
//  Copyright (c) 2012 Dinakar Gade. All rights reserved.
//

#ifndef LGomFixedRecourse_LGomFixedRecourse_h
#define LGomFixedRecourse_LGomFixedRecourse_h

#include <stdio.h>
#include "cplex.h"

#define TRUE 1
#define FALSE 0

#define NAMELEN 50
#define LENGTH 400
#define ROWSIZE_W2 65
#define ROWSIZE_W  1500              // Initial size of row (constraint) arrays for W
#define FIELDLEN 500
#define ROWSIZE_A 1000              // Initial size of row (constraint) arrays for A
#define INIT_LOWER_BOUND -1000000     //Added by Dinakar Gade. Lower bound on the optimal objective function value of SIP. Recourse can be now negative
#define PENALTY 100000
#define INT_PRECISION 1e-6
#define EPSILON_TOLERANCE 1e-6
#define NONZERO_LB 1e-4
#define INIT_ROWSIZE_FACTOR 100     //Row factor fo for initial change_sense, change_rhs, pi, etc
#define MAX_CUTS_PER_SCENARIO 5000    // Maximum number of cuts allowed per scenario - Before memory is realloced
// SMH: MAX_CUTS_PER_SCENARIO is increased from 500 to 5000 to avoid memory corruption. I believe no new memory is allocated when this bound is exceeded.
#define CUTS_PER_ROUND 1
#define TIME_LIMIT 3600
#define PERCENT_GAP 0.0025
#define MAX_ITER 10000
#define ZERO_TOL 1e-6
#define UB_ITER   200                      //Compute Upper bound after every UB_ITER iterations
#define GAP_UBOUNDING 25
// SMH: GAP_UBOUNDING has been changed to 25 from 30000





/**
 * Structure to hold subproblem info
 */
typedef struct  {
	int nrowsW;      	// an integer indicating the number of rows in matrix W (constant)
	int nrowsWmip;      // an integer indicating the number of rows in matrix W for mip (constant)
	int nrows_mip;      	// an integer indicating the number of rows in matrix W^k for mip
	int nrows_total;     // an integer indicating the number of rows in matrix W^k for C^3 lp
	int nrows;      	// an integer indicating the number of rows in matrix W^k
	int ncols;	   	// an integer indicating the length of the arrays cmatbeg_W, cmatind_W
	double *obj;    	// objective vector for the subproblem
	char *ctype;    	// Stores ctypes for for subproblem lp for determining integer variables
	char *ctype_ctns;	// Array for storing ctypes 'C'=`- continuous for subproblem lp
	char *sense;    	// Temp array for storing the contraint senses read from the core file
	char *sense_geq;    	// Temp array for storing the contraint senses >= as req'd by the D2 algorithm
	double *rhs;   	// RHS vector for subproblem r as read from the core file
	double **rhsRho;   	// RHS vector for subproblem scenario omega rh0(w) = r(w) - T(w)x
	char **colnames;	// Column or variable name array
	char **rownames;	// Row or constraint name array
	char *colnamestore;	// column name array
	char *rownamestore;	// row name array
	double *lb;		// Array of lower bounds on the cols
	double *ub;		// Array of upper bounds on the cols
	int *indices_ctype;	// Array of length ncols containing the indices for the ctypes
	int *indices_row;	// Array of length nrows containing the numerical indices
	// of the rows corresponding to the constraints for which
	// the rhs coefs are to be changed for the subproblem lp
	int *indices_rowmip;	// Array of length nrows containing the numerical indices
	// of the rows corresponding to the constraints for which
	// the rhs coefs are to be changed for the subproblem mip
	double *duals;	// Array of dual multipliers
	
	/* Constraint matrix W storage variables : Initialy it is a column-by-column sparse
	 matrix amd then it is changed to row-by-row to allow for addition of D^2 cut pi
	 coefs (this is more convenient)
	 */
	int nzcnt_W;		// number of nonzeros in matrix W
	int *cmatbeg_W;     // an array containing indices of where each col begins in the array cmatval and cmatind
	int *cmatcnt_W;     // an array containing the number of entries in each column of W
	int *cmatind_W;     // an array containing the row indices associated with the elements of cmatval
	double *cmatval_W;	// an array containing nonzero coefficients of the specified columns
	int cmatspace_W;	// an integer indicating the size of matrix W
	int spacesize_W;	// Current allocated space for W
	
	int *disjVarWindex; // For each row of pi coefs this array contains indices indicating
	// indicating the disjunction variable used to generate that row
	// Need this when forming the C3-LP. For a given disjunction var
	int nd2cuts;        // Number of D^2 cuts added thus far
	
	// Constant technology matrix T(w) = T storage variables
	int nrows_T;	   	// an integer indicating the number of rows in the matrix T
	int ncols_T;        //Added by Dinakar Gade
	int nzcnt_T;	// Number of nonzeros in matrix T for scenario w
	int *cmatbeg_T;	// An array containing indices of where each col begins in the array cmatval
	// and cmatind
	int *cmatcnt_T;     	// Array containing the number of entries in each column of T
	int *cmatind_T;	// Array containing the row indices associated with the elements of cmatval
	double*cmatval_T;  	// Array containing nonzero coefficients of the specified columns
	int cmatspace_T;    // An integers indicating the size of matrix T for each scenario w
	int spacesize_T;	// Current allocated space for T
	/**************************************************************
	 Added by Dinakar Gade - Row sparse matrices for T matrix
	 **************************************************************/
	int *rmatbeg, *rmatcnt;
	double *rmatval;
	int *rmatind;
	int rmatspace;
	
} subproblem_t;

typedef struct  {
	char *probname;       	// The name of the problem
	char *content;        	// File content: PERIODS
	char *distn;		// Support distrn type: DSCRETE
	char **scenName;		// Scenario name array
	double *coef;		// coef for variable 'col in contraint row
	double **rhs;		// rhs elements for each scenario
	// corresponding to the exact location in the subprobPtr->rhs
	double *scenProb;		// Array of probabilities for each scenario
	double *scenCondProb;	// Conditional probability for each scenario
	int nscens;			// total number of scenarios read from stoch file
	int nrows;           	// Number of rows in the T(w) matrix
	int ncols;	   	        // an integer indicating the length of the arrays cmatbeg_T(w), cmatind_T(w)
	
	// Random Objective: Added Dec 26, 2003
	int obj_cnt;	// Number of scenario subproblem random obj coefs
	double **obj;	// Random obj elements for each scenario subproblem
	int *obj_index;	// Array of indices of random obj elements for scenario subproblem
	
	
	
	// Technology matrix T(w) storage variables read from the stoch file
	// This is a col by col sparse matrix
	int *cnzcnt_T;	 	// vector of number of nonzeros in matrix T for scenario w
	int **cmatbeg_T;		// an array containing indices of where each col begins in the array cmatval
	// and cmatind for each scenario
	int **cmatcnt_T;     	// an array containing the number of entries in each column of T
	int **cmatind_T;		// an array containing the row indices associated with the elements of cmatval
	// for each scenario w
	double **cmatval_T;  	// an array containing nonzero coefficients of the specified columns for each
	// scenario w
	int *cmatspace_T;        	// an array of integers indicating the size of matrix T
	// for each scenario w
	int nscenspace_T;		// Current allocated space for this T in terms of scenarios
	
	// Technology matrix T(w) for newly added rows
	// This is a row by row sparse matrix
	int rnrows;           	// Number of rows in the rowT(w) matrix
	// cmatbeg_T(w), cmatind_T(w)
	int *rmatspace_T;
	int *rnzcnt_T;	 	// vector of number of nonzeros in matrix T for scenario w
	int **rmatbeg_T;		// an array containing indices of where each col begins in the array cmatval
	// and cmatind for each scenario
	int **rmatcnt_T;     	// an array containing the number of entries in each row of T
	int **rmatind_T;		// an array containing the row indices associated with the elements of cmatval
	// for each scenario w
	double **rmatval_T;  	// an array containing nonzero coefficients of the specified columns for each
	// scenario w
} stochfile_info;

typedef struct  {
	int nrows;      	// an integer indicating the number of rows in matrix A
	int ncols;	   	// an integer indicating the length of the arrays cmatbeg_A, cmatind_A
	double *obj;    	// objective vector for the subproblem
	char *ctype;    	// Stores ctypes for for subproblem lp for determining integer variables
	char *sense;    	// Temp array for storing the contraint senses read from the core file
	double *rhs;   	// RHS vector for subproblem
	char **colnames;	// Array of pointers to column names in colnamestore array
	char *colnamestore;	// Column name array
	//char **rownames;	// Row or constraint name array
	
	/* Constraint matrix A storage variables */
	// Row sparse matrix format
	int nzcnt_A;		// number of nonzeros in matrix W
	int *rmatbeg_A;	// an array containing indices of where each col begins in the array cmatval
	// and cmatind
	int *rmatcnt_A;     // an array containing the number of entries in each column of W
	int *rmatind_A;	// an array containing the row indices associated with the elements of cmatval
	double *rmatval_A;	// an array containing nonzero coefficients of the specified columns
	int rmatspace_A;	// an integer indicating the size of matrix A
	int spacesize_A;	// Current allocated space for A
	
	double rhsCoef;	// Rhs coefs in optimality cut
	double *cutCoefs;	// Array to store optimality cut coefs
	
} masterproblem_t;


/**
 Added by Dinakar Gade
 This structure will hold cut coefficients for the (FR:=) fixed recourse matrix
 Will hold cut coefficients for the Fixed Technology Matrix
 SCENARIO_RHS = Scenario Righ-Hand-Sides
 */
typedef struct {
	int ncuts;                              //Number of cuts
	int numnz_W, numnz_T;                   //Number of Nonzeros of cuts fixed recourse W and Technology Matrices (T)
	int *rmatbeg_W, *rmatind_W, *rmatcnt_W;  //Sparse index matrices for cuts W
	int *rmatbeg_T, *rmatind_T, *rmatcnt_T;  //Sparse index matrices for T
	double *rmatval_W;                       //Spare matrix values for W
	double *rmatval_T;                       //Sparse matrix values for T
	double **cutrhs;                         // Cut RHS by scenario
} CUTS_FR;


int loadTimeFile(char *rowname_start, char *colname_start, char *filename, FILE *fpout);
int getnumscenarios(char *filename);
int memAllocStochFileStruct(stochfile_info *stochdataPtr, int nrows_sub, int ncols_master,int ncols_sub);
int memAllocMasterProblemStructs(masterproblem_t *masterprobPtr, int nrows_master, int ncols_master);
int memAllocSubProblemStruct(subproblem_t *subprobPtr, int nrows_sub, int ncols_sub, int ncols_master,int numscens);
int loadMasterProblemData(CPXENVptr env, CPXLPptr lp_core, masterproblem_t *masterprobPtr, int nrows_master,int ncols_master,
						  int nrows_core, int ncols_core, FILE *fpout);
int setupSubProbMip(CPXENVptr env, CPXLPptr lp_submip, subproblem_t *subprobPtr, int nrows_master,
					int ncols_master, int nrows_sub, int ncols_sub,  FILE *fpout);
int addOptColToMasterProbLP(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr);
int loadSubProblemData(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
					   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr,
					   FILE *fpout);
int changeSlacksToIntegers(CPXENVptr env, CPXLPptr lp_submip);
int addObjConstraint2(CPXENVptr env, CPXLPptr lp);
int convert_col_to_row(int *cmatbeg, int *cmatcnt, int *cmatind, double *cmatval, int *rmatbeg, int *rmatcnt, int *rmatind, double *rmatval, int nrows, int ncols, int numnz);
int
loadStochFile(stochfile_info *stochdataPtr, subproblem_t *subprobPtr, char *filename, int *random_T, int *random_obj, int ncols_master, FILE *fpout);
int memAllocFixedRecourseStruct(CUTS_FR *frcuts, int ncols_W, int ncols_T, int nscens);
void freeFixedRecourseStructs(CUTS_FR *frcuts, int nscens);

int compute_rhs_rows(subproblem_t *subprobPtr, stochfile_info *stochdataPtr, CUTS_FR *frcuts, double *rhs, int rhssize, int omega, double *solnX, int x_size);
int is_element_int(double x);
int check_sol_frac(double *x, int n);
int generate_fixed_recourse_cuts(CPXENVptr env, CPXLPptr lp, subproblem_t *subprobPtr, stochfile_info *stochdataPtr, CUTS_FR *frcuts, int ncuts, double *solnX, int scen);

void print_CPX_sparse_matrix(int *rmatbeg, int *rmatcnt, int *rmatind, double *rmatval, int nrows, int ncols, int numnz);






void freeSubProbStruct(subproblem_t *subprobPtr, int nscens);
void free_stochDataPtr(stochfile_info *stochdataPtr);


double myfabs(double x);



#include <stdio.h>
#include <cplex.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h> 		/* Computation time functions */

class ParametricGomoryAlg
{

	stochfile_info    struct_s;             // STOCH file data structure
	subproblem_t      sub_prob_t ;
	masterproblem_t   master_prob_t ;       // Master problem lp data structure

public:
	ParametricGomoryAlg ();
	
	bool solve (int argc, const char * argv[]);
	
private:
	
	bool solve_subproblems();
	
};




#endif

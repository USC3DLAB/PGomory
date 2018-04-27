//
//  LGomFixedRecourseFuncs.c
//  LGomFixedRecourse
//
//  Created by Dinakar Gade on 9/8/12.
//  Copyright (c) 2012 Dinakar Gade. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cplex.h>
#include <math.h>

#include "GomAlg.h"

int resetconstraints(CPXENVptr env, CPXLPptr lp);
int addFeasColToSubProbLP(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr);

double myfabs(double x) {
	if (x < 0.0) {
		return -x;
	}
	else {
		return x;
	}
}

double myceil(double x) {
	double ceiling = -CPX_INFBOUND;
	if (is_element_int(x)) {
		ceiling = x;
	}
	else {
		ceiling = (double) ceil(x);
	}
	return ceiling;
}

double myfloor(double x) {
	double myfloornum = -CPX_INFBOUND;
	if (is_element_int(x)) {
		myfloornum = x;
	}
	else {
		myfloornum = (double)floor(x);
	}
	return myfloornum;
}

int loadTimeFile(char *rowname_start, char *colname_start, char *filename, FILE *fpout)
/**
 * Reads the TIME file data and puts the data in the structure timedataPtr
 * @param rowname_start subproblem first row name
 * @param colname_start subproblem first column name
 * @param filename  TIME file name
 * @param fpout output file pointer
 * @return 0 if TIME file is successfully read, otherwise return a nonzero integer
 */
{
	FILE *time;
	char field1[NAMELEN], field2[NAMELEN];
	char buffer[LENGTH];
	
	//int j;
	
	// Open the TIME file
	time = fopen(filename,"r" );
	if(time == NULL) {
		fprintf (stderr, "\n loadTimeFile():\n");
		fprintf(stderr, "Could not open the TIME file %s for reading!\n", filename);
		return(1);
	}
	
	// Process the TIME file
	if (fgets (buffer, LENGTH, time) != NULL) {
		//fprintf(stdout, "%s \n", buffer);
		sscanf(buffer, "%s %s", field1, field2);
		//fprintf(stdout, "field1: %s \n", field1);
		//fprintf(stdout, "field2: %s \n", field2);
		if ( strcmp(field1, "TIME") != 0) {
			fprintf (stderr, "\n loadTimeFile():\n");
			fprintf(stderr, "The first line of file %s must start with ""TIME"" and not %s \n!", filename, field1);
			fprintf(stderr, "Exiting ...\n");
			return(1);
		}
	}
	
	if (fgets (buffer, LENGTH, time) != NULL){
		sscanf(buffer, "%s%s", field1, field2);
	}
	
	if (fgets (buffer, LENGTH, time) == NULL){
		fprintf (stderr, "\n loadTimeFile():\n");
		fprintf(stderr, "The time file %s must have a description line for the first colname and row names", filename);
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	if (fgets (buffer, LENGTH, time) == NULL){
		fprintf (stderr, "\n loadTimeFile():\n");
		fprintf(stderr, "The time file %s must have a description line for the first subproblem colname and row names", filename);
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	sscanf(buffer, "%s%s", colname_start, rowname_start);
	//fprintf(fpout, "%s  %s \n", colname_start, rowname_start);;
	
	
	if (fgets (buffer, LENGTH, time) != NULL){
		sscanf(buffer, "%s", field1);
		if ( strcmp(field1, "ENDATA") != 0) {
			fprintf (stderr, "\n loadTimeFile():\n");
			fprintf(stderr, "The last line of file %s must start with ""ENDATA"" \n!", filename);
			fprintf(stderr, "Exiting ...\n");
			return(0);
		}
	} else {
		fprintf (stderr, "\n loadTimeFile():\n");
		fprintf(stderr, "The last line of file %s must start with ""ENDATA"" \n!", filename);
		fprintf(stderr, "Exiting ...\n");
		return(1);
	}
	
	// close the TIME file
	fclose(time);
	return (0);
	
} /************************** End loadTimeFile() function *****************************/

int getnumscenarios(char *filename)
/**
 * Reads the STOCH file data and determines the total number of scenarios
 * @param filename  STOCH file name
 * @return 0 number of scenarios
 */
{
	
	FILE *stoch;
	char field1[NAMELEN], field2[NAMELEN];
	char buffer[LENGTH];
	
	int numscenarios = 0;  // counters
	
	// Open the TIME file
	stoch = fopen(filename,"r" );
	if(stoch == NULL) {
		fprintf(stderr, "\ngetnumscenarios(...): \n");
		fprintf(stderr, " Could not open the STOCH file %s for reading!\n", filename);
		return(1);
	}
	
	
	//************************* Read the TIME file ************************************
	// Read first line: e.g. STOCH	example
	if (fgets (buffer, LENGTH, stoch) != NULL) {
		sscanf(buffer, "%s", field1);
		if (strcmp(field1, "STOCH") != 0) {
			fprintf(stderr, "\ngetnumscenarios(...): \n");
			fprintf(stderr, " The first line of the STOCH file %s must start with ""STOCH"" \n!", filename);
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	}
	// Read the second line: e.g SCENARIOS	DISCRETE
	if (fgets (buffer, LENGTH, stoch) != NULL){
		sscanf(buffer, "%s%s", field1, field2);
		if ( strcmp(field1, "SCENARIOS") != 0){
			fprintf(stderr, "\nloadStochFile(...): \n");
			fprintf(stderr, " The second line of file %s must start with ""SCENARIOS""!\n", filename);
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
		if ( strcmp(field2, "DISCRETE") != 0){
			fprintf(stderr, "\ngetnumscenarios(...): \n");
			fprintf(stderr, " The second word in the second line of file %s must be ""DISCRETE""!\n", filename);
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	}
	
	// Read the rest of the file and count scenarios
	while (fgets (buffer, LENGTH, stoch) != NULL) {
		//fprintf(stderr, "BUFFER: %s\n", buffer);
		sscanf(buffer, "%s", field1);
		if (strcmp(field1, "SC") == 0)  // Read in new scenario data
			numscenarios++; 	// count scenario
	}
	
	// close the STOCH file
	fclose(stoch);
	
	return numscenarios;
	
} //************************** End getnumscenarios() function *****************************/

int memAllocStochFileStruct(stochfile_info *stochdataPtr, int nrows_sub, int ncols_master,int ncols_sub) {
	/**
	 * Allocates memory to SMPS Format STOCH file data structure variables
	 * @param stochdataPtr  pointer to the STOCH file data structure
	 * @param nrows_sub  integer indicating number of rows in subproblem lp
	 * @param ncols_master  integer indicating number of columns in master lp
	 * @param ncols_sub  integer indicating number of columns in subproblem lp
	 * @return 0 if memory allocation is success, otherwise return a nonzero integer
	 */
	
	int j;
	
	// STOCH file data structure variables
	//stochdataPtr->nrows  = nrows_sub;
	stochdataPtr->ncols  = ncols_master;
	int rmatspace_T;
	int space;
	int numscens = stochdataPtr->nscens;
	//printf("nrows_sub = %d\n", nrows_sub);
	//printf("numscens = %d\n", numscens);
	//printf("ncols_sub = %d\n", ncols_sub);
	
	stochdataPtr->probname = (char*)malloc(NAMELEN*sizeof(char));
	stochdataPtr->content  = (char*)malloc(FIELDLEN*sizeof(char));
	stochdataPtr->distn = (char*)malloc(FIELDLEN*sizeof(char));
	stochdataPtr->scenName  = (char**)malloc(numscens*sizeof(char*));
	stochdataPtr->rhs      = (double**)malloc(numscens*sizeof(double*));
	stochdataPtr->scenProb = (double*)malloc(numscens*sizeof(double));
	stochdataPtr->scenCondProb = (double*)malloc(numscens*sizeof(double));
	
	//Added DEC 26, 2003
	stochdataPtr->obj = (double**)malloc(numscens*sizeof(double*));
	stochdataPtr->obj_index = (int*)malloc(ncols_sub*sizeof(int));
	
	if (stochdataPtr->probname == NULL || stochdataPtr->content   == NULL ||
		stochdataPtr->distn    == NULL || stochdataPtr->scenName  == NULL ||
		stochdataPtr->rhs       == NULL || stochdataPtr->scenProb == NULL ||
		stochdataPtr->scenCondProb == NULL || stochdataPtr->obj == NULL   ||
		stochdataPtr->obj_index == NULL)
	{
		fprintf (stderr, "memAllocStochFileStruct(...):\n");
		fprintf(stderr, "Failure to allocate memory to stochfile pointers\n");
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	
	for (j = 0; j < numscens; j++) {
		stochdataPtr->scenName[j] = (char*)malloc(FIELDLEN*sizeof(char));
		stochdataPtr->rhs[j]      = (double*)malloc((ROWSIZE_W)*sizeof(double));
		stochdataPtr->obj[j]      = (double*)malloc(ncols_sub*sizeof(double));
		if(stochdataPtr->scenName[j] == NULL || stochdataPtr->rhs[j] == NULL ||
		   stochdataPtr->obj[j] == NULL) {
			fprintf (stderr, "memAllocStochFileStructs(...):\n");
			fprintf(stderr, "\nFailure to allocate memory to stochfile scenName[]  and rhs[] \n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	}// End for loop
	
	// printf("Number 9b\n");
	
	// Random Technology Matrix as read from the stoch file if available
	// This is a col by col sparse matrix
	
	stochdataPtr->rnrows  = 0;	// Initialize to zero
	
	stochdataPtr->cnzcnt_T = (int*)malloc(numscens*sizeof(int));
	stochdataPtr->cmatbeg_T = (int**)malloc(numscens*sizeof(int*));
	stochdataPtr->cmatcnt_T = (int**)malloc(numscens*sizeof(int*));
	stochdataPtr->cmatind_T = (int**)malloc(numscens*sizeof(int*));
	stochdataPtr->cmatval_T = (double**)malloc(numscens*sizeof(double*));
	
	if (stochdataPtr->cnzcnt_T == NULL || stochdataPtr->rhs == NULL ||
		stochdataPtr->cmatbeg_T == NULL || stochdataPtr->cmatbeg_T == NULL ||
		stochdataPtr->cmatind_T == NULL || stochdataPtr->cmatval_T == NULL) //||       stochdataPtr->cmatspace_T == NULL)
	{
		fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
		fprintf(stderr, "Failure to allocate memory to subproblem T(w) data pointers\n");
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	
	space = ncols_master*ROWSIZE_W2;
	for (j = 0; j < numscens; j++){
		stochdataPtr->cmatbeg_T[j]  = (int*)malloc(stochdataPtr->ncols*sizeof(int));
		stochdataPtr->cmatcnt_T[j]  = (int*)malloc(stochdataPtr->ncols*sizeof(int));
		stochdataPtr->cmatind_T[j]  = (int*)malloc(space*sizeof(int));
		stochdataPtr->cmatval_T[j]  = (double*)malloc(space*sizeof(double));
		
		if (stochdataPtr->cmatbeg_T[j] == NULL || stochdataPtr->cmatcnt_T[j] == NULL ||
			stochdataPtr->cmatind_T[j] == NULL || stochdataPtr->cmatval_T[j] == NULL)
		{
			fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
			fprintf(stderr, "Failure to allocate memory to subproblem stochastic T(w)[] data arrays\n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	} // end for loop
	// Random Technology Matrix for update by the D2 algorithm
	// This is a row by row sparse matrix
	stochdataPtr->nscenspace_T = numscens;
	
	stochdataPtr->rnzcnt_T = (int*)calloc(numscens, sizeof(int));
	stochdataPtr->rmatbeg_T = (int**)malloc(numscens*sizeof(int*));
	stochdataPtr->rmatind_T = (int**)calloc(numscens, sizeof(int*));
	stochdataPtr->rmatval_T = (double**)malloc(numscens*sizeof(double*));
	
	if (stochdataPtr->rnzcnt_T  == NULL || stochdataPtr->rmatbeg_T == NULL ||
		stochdataPtr->rmatbeg_T == NULL || stochdataPtr->rmatind_T == NULL ||
		stochdataPtr->rmatval_T == NULL ) // stochdataPtr->rmatspace_T == NULL) //|| stochdataPtr->rmatcnt_T == NULL)
	{
		fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
		fprintf(stderr, "Failure to allocate memory to subproblem T(w) row-by-row data arrays\n");
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	
	// printf("Number 10\n");
	rmatspace_T  = ncols_master*ROWSIZE_W2;      // Initial assumed size
	for (j = 0; j < numscens; j++){
		stochdataPtr->rmatbeg_T[j]  = (int*)malloc(ROWSIZE_W2*sizeof(int));
		stochdataPtr->rmatind_T[j]  = (int*)malloc(rmatspace_T*sizeof(int));
		stochdataPtr->rmatval_T[j]  = (double*)malloc(rmatspace_T*sizeof(double));
		
		if (stochdataPtr->rmatbeg_T[j] == NULL || stochdataPtr->rmatind_T[j] == NULL
			|| stochdataPtr->rmatval_T[j] == NULL)
		{
			fprintf(stderr, "\n memAllocStochFileStruct(...): \n");
			fprintf(stderr, "Failure to allocate memory to subproblem stochastic T(w)[] row-by-row data arrays\n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	} // end for loop
	
	
	return(0); // successfull return
	
}//************************** End memAllocStochFileStruct function *****************************/

int memAllocMasterProblemStructs(masterproblem_t *masterprobPtr, int nrows_master, int ncols_master)
/**
 * Allocates memory to subproblem data structure variables
 * @param masterprobPtr  pointer to the master problem data structure
 * @param nrows_master  integer indicating number of stage 1 rows (constraints)
 * @param ncols_master  integer indicating number of stage 1 columns (vars)
 * @return 0 if memory allocation is success, otherwise return a nonzero integer
 */
{
	// initializations
	//int status;
	int j;
	
	// Initialize subproblem data structure vars
	// Will add a column for the theta variable, i.e, optimality variable
	masterprobPtr->nrows  = nrows_master;
	masterprobPtr->ncols  = ncols_master;
	masterprobPtr->rhsCoef = 0;
	masterprobPtr->rmatspace_A  = (ncols_master+1)*ROWSIZE_A; // assumed max number of nonzeros in A
	
	//printf("masterprobPtr->nrows = %d\n masterprobPtr->ncols = %d\n", masterprobPtr->nrows, masterprobPtr->ncols);
	
	masterprobPtr->obj = (double*)malloc((ncols_master+1)*sizeof(double));
	masterprobPtr->ctype = (char*)malloc((ncols_master+1)*sizeof(char));
	masterprobPtr->sense = (char*)malloc(ROWSIZE_A*sizeof(char));
	masterprobPtr->rhs = (double*)malloc(ROWSIZE_A*sizeof(double));
	masterprobPtr->colnames = (char**)malloc((ncols_master+1)*sizeof(char*));
	masterprobPtr->colnamestore = (char*)malloc((ncols_master+1)*FIELDLEN*sizeof(char));
	masterprobPtr->cutCoefs = (double*)malloc((ncols_master+1)*sizeof(double));
	
	masterprobPtr->rmatbeg_A = (int*)malloc((ncols_master+1)*sizeof(int));
	masterprobPtr->rmatcnt_A = (int*)malloc((ncols_master+1)*sizeof(int));
	masterprobPtr->rmatind_A = (int*)malloc(masterprobPtr->rmatspace_A*sizeof(int));
	masterprobPtr->rmatval_A = (double*)malloc(masterprobPtr->rmatspace_A*sizeof(double));
	
	if (masterprobPtr->obj   == NULL || masterprobPtr->ctype == NULL || masterprobPtr->sense == NULL ||
		masterprobPtr->rhs   == NULL || masterprobPtr->rmatbeg_A  == NULL ||
		masterprobPtr->rmatind_A == NULL || masterprobPtr->rmatind_A == NULL ||
		masterprobPtr->rmatval_A == NULL || masterprobPtr->colnames == NULL ||
		masterprobPtr->cutCoefs == NULL || masterprobPtr->colnamestore == NULL)
	{
		fprintf(stderr, "memAllocMasterProblemStructs(...): \n");
		fprintf(stderr, "Failure to allocate memory to masterproblem data arrays\n");
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	for (j = 0; j < (ncols_master+1); j++){
		masterprobPtr->colnames[j] = (char*)malloc(NAMELEN*sizeof(char));
		if (masterprobPtr->colnames[j] == NULL)
		{
			fprintf(stderr, "\n memAllocMasterProblemStruct(...): \n");
			fprintf(stderr, "Failure to allocate memory to master problem colnames[].\n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	} // end for loop
	
	return(0); // successfull return
	
}//************************** End memAllocMasterProblemStruct function *****************************

void free_stochDataPtr(stochfile_info *stochdataPtr) {
	/**
	 * Free up the stoch data pointer memory that was allocated
	 */
	
	int i;
	if(stochdataPtr->probname != NULL ) {
		free(stochdataPtr->probname);
		stochdataPtr->probname = NULL;
	}
	
	
	if(stochdataPtr->content != NULL ) {
		free(stochdataPtr->content);
		stochdataPtr->content = NULL;
	}
	
	
	if(stochdataPtr->distn != NULL) {
		free(stochdataPtr->distn);
		stochdataPtr->distn = NULL;
	}
	if(stochdataPtr->scenName != NULL) {
		free(stochdataPtr->scenName);
	}
	
	if (stochdataPtr->scenProb != NULL) {
		free(stochdataPtr->scenProb);
		stochdataPtr->scenProb = NULL;
	}
	
	if (stochdataPtr->scenCondProb != NULL) {
		free(stochdataPtr->scenCondProb);
		stochdataPtr->scenCondProb = NULL;
	}
	
	if (stochdataPtr->obj_index != NULL) {
		free(stochdataPtr->obj_index);
		stochdataPtr->obj_index = NULL;
	}
	
	if (stochdataPtr->obj != NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->obj[i] != NULL) {
				free(stochdataPtr->obj[i]);
				stochdataPtr->obj[i] = NULL;
			}
		}
		free(stochdataPtr->obj);
		stochdataPtr->obj = NULL;
	}
	if(stochdataPtr->cnzcnt_T != NULL) {
		free(stochdataPtr->cnzcnt_T);
	}
	
	if(stochdataPtr->cmatbeg_T != NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->cmatbeg_T[i] != NULL) {
				free(stochdataPtr->cmatbeg_T[i]);
				stochdataPtr->cmatbeg_T[i] = NULL;
			}
		}
		free(stochdataPtr->cmatbeg_T);
		stochdataPtr->cmatbeg_T = NULL;
	}
	
	if(stochdataPtr->cmatcnt_T != NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->cmatcnt_T[i] != NULL) {
				free(stochdataPtr->cmatcnt_T[i]);
				stochdataPtr->cmatcnt_T[i] = NULL;
			}
		}
		free(stochdataPtr->cmatcnt_T);
		stochdataPtr->cmatcnt_T = NULL;
	}
	
	if(stochdataPtr->cmatind_T!= NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->cmatind_T[i] != NULL) {
				free(stochdataPtr->cmatind_T[i]);
				stochdataPtr->cmatind_T[i] = NULL;
			}
		}
		free(stochdataPtr->cmatind_T);
		stochdataPtr->cmatind_T = NULL;
	}
	
	if(stochdataPtr->cmatval_T!= NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->cmatval_T[i] != NULL) {
				free(stochdataPtr->cmatval_T[i]);
				stochdataPtr->cmatval_T[i] = NULL;
			}
		}
		free(stochdataPtr->cmatval_T);
		stochdataPtr->cmatval_T = NULL;
	}
	
	if(stochdataPtr->rnzcnt_T != NULL) {
		free(stochdataPtr->rnzcnt_T);
		stochdataPtr->rnzcnt_T = NULL;
	}
	
	if(stochdataPtr->rmatbeg_T != NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->rmatbeg_T[i] != NULL ) {
				free(stochdataPtr->rmatbeg_T[i]);
				stochdataPtr->rmatbeg_T[i] = NULL;
			}
		}
		free(stochdataPtr->rmatbeg_T);
		stochdataPtr->rmatbeg_T = NULL;
	}
	
	if(stochdataPtr->rmatind_T != NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if(stochdataPtr->rmatind_T[i] != NULL ) {
				free(stochdataPtr->rmatind_T[i]);
				stochdataPtr->rmatind_T[i] = NULL;
			}
		}
		free(stochdataPtr->rmatind_T);
		stochdataPtr->rmatind_T = NULL;
	}
	
	if(stochdataPtr->rmatval_T != NULL) {
		for (i = 0; i < stochdataPtr->nscens; i++) {
			if (stochdataPtr->rmatval_T[i] != NULL) {
				free(stochdataPtr->rmatval_T[i]);
				stochdataPtr->rmatval_T[i] = NULL;
			}
		}
		free(stochdataPtr->rmatval_T);
		stochdataPtr->rmatval_T = NULL;
	}
}

//************************** End free_stochDataPtr function *****************************

int memAllocSubProblemStruct(subproblem_t *subprobPtr, int nrows_sub, int ncols_sub, int ncols_master,int numscens)
/**
 * NEW VERSION FOR GDD ALGORITHM - DINAKAR GADE
 * Allocates memory to subproblem data structure variables
 * @param subprobPtr  pointer to the subproblem structure
 * @param nrows_sub   	integer indicating number of stage 2 rows (constraints)
 * @param ncols_sub   	integer indicating number of stage 2 columns (vars)
 * @param ncols_master   integer indicating number of stage 1 columns (vars)
 * @param numscens   number of scenario subproblems
 * @return 0 if memory allocation is successful, otherwise return a nonzero integer
 */
{
	// initializations
	//int status;
	int j;
	
	// Initialize the subproblem data structure
	//nrows_sub += ncols_sub; // Will add binary constraints for each var - UNCOMMENTED BY DINAKAR GADE: No need to add binary constraints for each
	//subprobPtr->nrowsW   = nrows_sub;
	//subprobPtr->nrows   = nrows_sub;
	//subprobPtr->nrows_T = nrows_sub;
	subprobPtr->ncols   = ncols_sub;
	subprobPtr->cmatspace_W  = (ncols_sub+1)*ROWSIZE_W; // assumed max number of nonzeros in W
	//subprobPtr->nd2cuts   = 0;
	
	subprobPtr->cmatspace_T = (nrows_sub+1)*ncols_master;
	subprobPtr->nzcnt_T = nrows_sub*ncols_master;
	
	//printf("subprobPtr->ncols = %d\n", subprobPtr->ncols);
	//fprintf(stdout, "In method memAllocSubProblemStruct_new(): subprobPtr->nrows = %d\n", nrows_sub);
	
	
	subprobPtr->obj = (double*)malloc((ncols_sub+1)*sizeof(double));
	subprobPtr->ctype = (char*)malloc((ncols_sub+1)*sizeof(char));
	subprobPtr->ctype_ctns = (char*)malloc((ncols_sub+1)*sizeof(char));
	subprobPtr->sense = (char*)malloc(ROWSIZE_W*sizeof(char));
	subprobPtr->sense_geq = (char*)malloc(ROWSIZE_W*sizeof(char));
	subprobPtr->rhs = (double*)malloc(nrows_sub*sizeof(double));
	subprobPtr->rhsRho = (double**)malloc(numscens*sizeof(double*));
	subprobPtr->rownames      = (char**)malloc(nrows_sub*sizeof(char*));
	subprobPtr->colnames      = (char**)malloc((subprobPtr->ncols+1)*sizeof(char*));
	subprobPtr->rownamestore  = (char*)malloc(nrows_sub*FIELDLEN*sizeof(char));
	subprobPtr->colnamestore  = (char*)malloc((subprobPtr->ncols+1)*FIELDLEN*sizeof(char));
	subprobPtr->lb = (double*)malloc((ncols_sub+1)*sizeof(double));
	subprobPtr->ub = (double*)malloc((ncols_sub+1)*sizeof(double));
	subprobPtr->indices_row = (int*)malloc(ROWSIZE_W*sizeof(int));
	subprobPtr->indices_rowmip = (int*)malloc(ROWSIZE_W*sizeof(int));
	subprobPtr->indices_ctype = (int*)malloc((ncols_sub+1)*sizeof(int));
	//subprobPtr->disjVarWindex = (int*)malloc(ROWSIZE_W*sizeof(int));//Dinakar Gade comment
	subprobPtr->duals = (double*)malloc(ROWSIZE_W*sizeof(double));
	
	subprobPtr->cmatbeg_W = (int*)malloc((ncols_sub+1)*nrows_sub*sizeof(int));
	subprobPtr->cmatcnt_W = (int*)malloc((ncols_sub+1)*nrows_sub*sizeof(int));
	subprobPtr->cmatind_W = (int*)malloc(subprobPtr->cmatspace_W*sizeof(int));
	subprobPtr->cmatval_W = (double*)malloc(subprobPtr->cmatspace_W*sizeof(double));
	
	
	subprobPtr->cmatbeg_T = (int*)malloc((ncols_master)*sizeof(int));
	subprobPtr->cmatcnt_T = (int*)malloc((ncols_master)*sizeof(int));
	subprobPtr->cmatind_T = (int*)malloc(subprobPtr->cmatspace_T*sizeof(int));
	subprobPtr->cmatval_T = (double*)malloc(subprobPtr->cmatspace_T*sizeof(double));
	
	if (subprobPtr->obj   == NULL || subprobPtr->ctype == NULL || subprobPtr->ctype_ctns == NULL ||
		subprobPtr->sense == NULL || subprobPtr->rhsRho   == NULL || subprobPtr->rownames   == NULL ||
		subprobPtr->colnames   == NULL || subprobPtr->cmatbeg_W  == NULL ||
		subprobPtr->cmatind_W == NULL || subprobPtr->cmatind_W == NULL ||
		subprobPtr->cmatval_W == NULL || subprobPtr->sense_geq == NULL ||
		subprobPtr->lb == NULL || subprobPtr->ub == NULL || subprobPtr->indices_row == NULL ||
		subprobPtr->indices_ctype == NULL || subprobPtr->rhs   == NULL ||
		subprobPtr->indices_rowmip == NULL)
	{
		fprintf(stderr, "Function memAllocSubProblemStruct(...): \n");
		fprintf(stderr, "Failure to allocate memory to subproblem data arrays\n");
		fprintf(stderr, "Exiting...\n");
		return(1);
	}
	//First initialize these to null
	subprobPtr->rmatbeg = NULL;
	subprobPtr->rmatcnt = NULL;
	subprobPtr->rmatind = NULL;
	subprobPtr->rmatval = NULL;
	//printf("DONE\n");
	//Added by Dinakar Gade
	//subprobPtr->rmatbeg = (int *)malloc(nrows_sub*sizeof(int));
	//subprobPtr->rmatind = (int *)malloc(subprobPtr->cmatspace_T*sizeof(int));
	//subprobPtr->rmatval = (double *)malloc(subprobPtr->cmatspace_T*sizeof(double));
	//if(subprobPtr->rmatbeg == NULL || subprobPtr->rmatind == NULL || subprobPtr->rmatval == NULL) {
	//    fprintf(stderr, "In Function memAllocSubProblemStruct():\n");
	//    fprintf(stderr, "\tFailure to allocate memory to subproblem row sparse arrays\n");
	//    return (1);
	//}
	
	for (j = 0; j < nrows_sub; j++){
		subprobPtr->rownames[j] = (char*)malloc(NAMELEN*sizeof(char));
		if (subprobPtr->rownames[j] == NULL)
		{
			fprintf(stderr, "\n memAllocSubProblemStruct(...): \n");
			fprintf(stderr, "Failure to allocate memory to sub problem rownames[].\n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	} // end for loop
	
	for (j = 0; j < subprobPtr->ncols+1; j++){
		subprobPtr->colnames[j] = (char*)malloc(NAMELEN*sizeof(char));
		if (subprobPtr->colnames[j] == NULL)
		{
			fprintf(stderr, "\n memAllocSubProblemStruct(...): \n");
			fprintf(stderr, "Failure to allocate memory to sub problem colnames[].\n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	} // end for loop
	
	for (j = 0; j < numscens; j++){
		subprobPtr->rhsRho[j] = (double*)malloc(ROWSIZE_W*sizeof(double));
		if (subprobPtr->rhsRho[j] == NULL)
		{
			fprintf(stderr, "\n memAllocSubProblemStruct(...): \n");
			fprintf(stderr, "Failure to allocate memory to sub problem subprobPtr->rhsRho[j].\n");
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	} // end for loop
	
	return(0); // successfull return
	
}//************************** End memAllocSubProblemStruct function *****************************/

void freeSubProbStruct(subproblem_t *subprobPtr, int nscens) {
	int i;
	if(subprobPtr->obj != NULL) {
		free(subprobPtr->obj);
		subprobPtr->obj = NULL;
	}
	if(subprobPtr->ctype != NULL) {
		free(subprobPtr->ctype);
		subprobPtr->ctype = NULL;
	}
	if(subprobPtr->ctype_ctns != NULL) {
		free(subprobPtr->ctype_ctns);
		subprobPtr->ctype_ctns = NULL;
	}
	if(subprobPtr->sense != NULL) {
		free(subprobPtr->sense);
		subprobPtr->sense = NULL;
	}
	if(subprobPtr->sense_geq != NULL) {
		free(subprobPtr->sense_geq);
		subprobPtr->sense_geq = NULL;
	}
	if(subprobPtr->rhs != NULL) {
		free(subprobPtr->rhs);
		subprobPtr->rhs = NULL;
	}
	if (subprobPtr->lb != NULL) {
		free(subprobPtr->lb);
		subprobPtr->lb = NULL;
	}
	if(subprobPtr->ub != NULL) {
		free(subprobPtr->ub);
		subprobPtr->ub = NULL;
	}
	if (subprobPtr->indices_row != NULL) {
		free(subprobPtr->indices_row);
		subprobPtr->indices_row = NULL;
	}
	if(subprobPtr->indices_rowmip != NULL) {
		free(subprobPtr->indices_rowmip);
	}
	if(subprobPtr->indices_ctype != NULL) {
		free(subprobPtr->indices_ctype);
		subprobPtr->indices_ctype = NULL;
	}
	
	if(subprobPtr->cmatbeg_W != NULL ){
		free(subprobPtr->cmatbeg_W);
		subprobPtr->cmatbeg_W = NULL;
	}
	if(subprobPtr->cmatind_W != NULL){
		free(subprobPtr->cmatind_W);
		subprobPtr->cmatind_W = NULL;
	}
	if(subprobPtr->cmatcnt_W != NULL) {
		free(subprobPtr->cmatcnt_W);
		subprobPtr->cmatcnt_W = NULL;
	}
	if(subprobPtr->cmatval_W != NULL) {
		free(subprobPtr->cmatval_W);
		subprobPtr->cmatval_W = NULL;
	}
	if(subprobPtr->cmatbeg_T != NULL ){
		free(subprobPtr->cmatbeg_T);
		subprobPtr->cmatbeg_T = NULL;
	}
	if(subprobPtr->cmatind_T != NULL){
		free(subprobPtr->cmatind_T);
		subprobPtr->cmatind_T = NULL;
	}
	if(subprobPtr->cmatcnt_T != NULL) {
		free(subprobPtr->cmatcnt_T);
		subprobPtr->cmatcnt_T = NULL;
	}
	if(subprobPtr->cmatval_T != NULL) {
		free(subprobPtr->cmatval_T);
		subprobPtr->cmatval_T = NULL;
	}
	
	if(subprobPtr->rmatbeg != NULL) {
		free(subprobPtr->rmatbeg);
		subprobPtr->rmatbeg = NULL;
	}
	if(subprobPtr->rmatind != NULL) {
		free(subprobPtr->rmatind);
		subprobPtr->rmatind = NULL;
	}
	if(subprobPtr->rmatval != NULL) {
		free(subprobPtr->rmatval);
		subprobPtr->rmatval = NULL;
	}
	if(subprobPtr->rhsRho != NULL) {
		for (i = 0; i < nscens; i++) {
			if(subprobPtr->rhsRho[i] != NULL) {
				free(subprobPtr->rhsRho[i]);
				subprobPtr->rhsRho[i] = NULL;
			}
		}
		free(subprobPtr->rhsRho);
		subprobPtr->rhsRho = NULL;
	}
	
}
//************************** End freeSubProbStruct() function *****************************/


int loadMasterProblemData(CPXENVptr env, CPXLPptr lp_core, masterproblem_t *masterprobPtr, int nrows_master,int ncols_master,
						  int nrows_core, int ncols_core, FILE *fpout)
/**
 * Extracts  master problem data from the lp read from the core file and puts it in a data structure
 * It is assumed that the master problem data structure has already been allocated memory
 * @param lp_core pointer to the lp model read from the core file
 * @param master prob_t  pointer to the master problem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_master   integer indicating number of master problem columns (vars)
 * @param nrows_core   integer indicating number of core file lp problem rows (constraints)
 * @param ncols_core   integer indicating number of core file lp problem columns (vars)
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */
{
	
	int status;
	int surplus; // To check for sufficiency of the aloocated array size
	//int i, j;
	
#ifdef DEBUG_FUNC
	fprintf (stderr, "loadMasterProblemData(...):\n");
	fprintf(stderr, "masterprobPtr->nrows = %d\n masterprobPtr->ncols = %d\n", masterprobPtr->nrows, masterprobPtr->ncols);
#endif
	
	
	//***** Create master problem lp by deleting stage_2 cols and rows from core file lp  *******
	status = CPXdelcols(env, lp_core, ncols_master, ncols_core-1);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to delete stage-2 cols from lp_core, error %d.\n", status);
		return status;
	}
	
	/* Delete stage-2 rows (constraints) */
	status = CPXdelrows(env, lp_core, nrows_master, nrows_core-1);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to delete stage-2 rows from lp_core rows, error %d.\n", status);
		return status;
	}
	//            Reset the <= contraints to >=            //
	status = resetconstraints(env, lp_core);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
		return status;
	}
	// Write master to file in lp format
#ifdef DEBUG_FUNC
	fprintf (stdout, "\nloadMasterProblemData():\n");
	status = CPXwriteprob(env, lp_core, "master.lp", NULL);
	if ( status ) {
		fprintf (stderr, "d2algMain():\n");
		fprintf (stderr, "Failure to write lp_core problem lp to file, error %d.\n", status);
		return status;
	}
#endif
	
	
	//************* Get and store the ctype (var type) array for the master problem lp ********
	status = CPXgetctype(env, lp_core, masterprobPtr->ctype, 0, ncols_master-1);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to get ctype from lp_core for the master lp, error %d.\n", status);
		return (status);
	}
#ifdef DEBUG_FUNC
	fprintf(fpout,"loadMasterProblemData(...): \n Master problem ctype: \n");
	for (j = 0; j < ncols_master; j++)
		fprintf(fpout, "col %d: %c\n", j, masterprobPtr->ctype[j]);
#endif
	
	//************* Extract the obj from the lp_core from the core file ****************
	status = CPXgetobj(env, lp_core, masterprobPtr->obj, 0, ncols_master-1);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to write get obj from lp_core for the master lp, error %d.\n", status);
		return(status);
	}
	
	//******************* Extract the A matrix from lp_core from the core file ****************
	// in row sparse format: more convenient for the formation of the RHS LP later
	if (nrows_master > 0) {
		status = CPXgetrows(env, lp_core, &masterprobPtr->nzcnt_A, masterprobPtr->rmatbeg_A, masterprobPtr->rmatind_A,
							masterprobPtr->rmatval_A, masterprobPtr->rmatspace_A, &surplus, 0, nrows_master-1);
	}
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to get rows from lp_core for A matrix, error %d.\n", status);
		fprintf (stderr, "Surplus value is: %d.\n", surplus);
		return(status);
	}
	
	//******************** Populate the cmatcnt_A matrix: **************************************
	// an array containing the number of entries in each column of A
	//masterprobPtr->rmatcnt_A[0] = masterprobPtr->rmatbeg_A[1];
	// for (i = 1; i < ncols_master-1; i++){
	//    masterprobPtr->rmatcnt_A[i] = masterprobPtr->rmatbeg_A[i+1]-masterprobPtr->rmatbeg_A[i];
	// }
	// masterprobPtr->rmatcnt_A[ncols_master-1] = masterprobPtr->nzcnt_A - masterprobPtr->rmatbeg_A[ncols_master-1];
	
	//****** Get the sense vector for the range of contraints for the subproblem ******/
	status = CPXgetsense(env, lp_core, masterprobPtr->sense, 0, nrows_master-1);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to get the senses from lp_core problem, error %d.\n", status);
		return(status);
	}
	
	//****** Get the rhs (b) vector for the master problem *******************************
	status = CPXgetrhs(env, lp_core, masterprobPtr->rhs, 0, nrows_master-1);
	if ( status ) {
		fprintf (stderr, "loadMasterProblemData(...):\n");
		fprintf (stderr, "Failure to get the rhs vector (b) from lp_core problem, error %d.\n", status);
		return(status);
	}
	
	
	//**************************** Print the A matrix, senses, and rhs ********************
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadMasterProblemData(...):\n");
	fprintf(stdout, "masterprobPtr->rmatspace_A = %d\n", masterprobPtr->rmatspace_A);
	fprintf(stdout, "masterprobPtr->nzcnt_A = %d\n", masterprobPtr->nzcnt_A);
	fprintf(stdout, "The A matrix is: \n");
	printSparseMatrix(nrows_master, masterprobPtr->nzcnt_A, masterprobPtr->rmatbeg_A, masterprobPtr->rmatcnt_A, masterprobPtr->rmatind_A,
					  masterprobPtr->rmatval_A, stdout);
	
	fprintf(stdout, "master problem senses > \n");
	for (j = 0; j < masterprobPtr->nrows; j++)
		fprintf(stdout, "Row %d: %c\n", j, masterprobPtr->sense[j]);
	
	fprintf(stdout, "master problem rhs (b) > \n");
	for (j = 0; j < masterprobPtr->nrows; j++)
		fprintf(stdout, "Row %d: %f\n", j, masterprobPtr->rhs[j]);
#endif
	
	return (0); // successful return
	
} //************************** End loadMasterProblemData function *****************************/

int resetconstraints(CPXENVptr env, CPXLPptr lp)
/**
 * Resets <= to <= constraints in the CORE FILE lp. This is required by the D2 Algorithm.
 * This is done row by row by multiplying both sides of the constraint by -1.
 * The scenario RHS must also be multiplied by -1 in the function loadStochFile function
 * @param env CPLEX environment
 * @param lp pointer to the lp model
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */
{
	int i, j;
	int cur_numcols;
	int cur_numrows;
	int status =  0;
	int nzcnt;
	int  *rmatbeg;
	int  *rmatind;
	double *rmatval;
	char *sense;
	
	int surplus;
	
	int *rowindex;
	char *newsense;
	double *rhs;
	
	int *rows;
	
	cur_numcols = CPXgetnumcols (env, lp);
	cur_numrows = CPXgetnumrows (env, lp);
	
	rmatbeg = (int *) malloc (cur_numrows*sizeof(int));
	rmatind = (int *) malloc (cur_numcols*sizeof(int));
	rmatval = (double *) malloc (cur_numcols*sizeof(double));
	rows    = (int *) malloc (cur_numcols*sizeof(int));
	sense   = (char *) malloc (cur_numrows*sizeof(char));
	
	if (rmatbeg == NULL || rmatind == NULL || rmatval == NULL || rows == NULL ||
		sense == NULL) {
		fprintf(stderr, "\nresetconstraints: ");
		fprintf(stderr, "\nMemory allocation failure for arrays.");
		return 1;
	}
	rowindex   = (int *) malloc (sizeof(int));
	rhs        = (double *) malloc (sizeof(double));
	newsense   = (char *) malloc (sizeof(char));
	if (rowindex == NULL || rhs == NULL || newsense == NULL) {
		fprintf(stderr, "\nresetconstraints: ");
		fprintf(stderr, "\n**Memory allocation failure for arrays.");
		return 1;
	}
	
	newsense[0] = 'G';
	
	//printf("\n  cur_numrows = %d\n", cur_numrows);
	//printf("   cur_numcols = %d\n", cur_numcols);
	
	// Get sense
	status = CPXgetsense (env, lp, sense, 0, cur_numrows-1);
	if (status) {
		fprintf(stderr, "\nresetconstraints: ");
		fprintf(stderr, "\nCPXgetsense failure at row. status = %d\n", status);
		return status;
	}
	
	// Loop thru all the contraints
	for (i = 0; i < cur_numrows; i++) {
		
		//printf("   cur_row = %d\n", i);
		// Change sense
		//printf("   sense[%d] = %c\n", i, sense[i]);
		if (sense[i] == 'L') {
			
			rowindex[0] = i;
			status =  CPXchgsense(env, lp, 1, rowindex, newsense);
			if (status) {
				fprintf(stderr, "\nresetconstraints: ");
				fprintf(stderr, "\nCPXchgsense failure at row %d. status = %d\n", i, status);
				return status;
			}
			
			// Change RHS
			status = CPXgetrhs (env, lp, rhs, i, i);
			if (status) {
				fprintf(stderr, "\nresetconstraints: ");
				fprintf(stderr, "\nCPXgetrhs failure at row %d. status = %d\n", i, status);
				return status;
			}
			rhs[0] *= -1;
			status = CPXchgrhs (env,lp,1,rowindex, rhs);
			if (status) {
				fprintf(stderr, "\nresetconstraints: ");
				fprintf(stderr, "\nCPXchgrhs failure at row %d. status = %d\n", i, status);
				return status;
			}
			
			// Change constraint
			status =  CPXgetrows (env,lp, &nzcnt, rmatbeg, rmatind, rmatval, cur_numcols,
								  &surplus, i, i);
			if (status) {
				fprintf(stderr, "\nresetconstraints: ");
				fprintf(stderr, "\nCPXgetrows failure at row %d. status = %d\n", i, status);
				return status;
			}
			for (j = 0; j < nzcnt; j++){
				rmatval[j] *= -1;
				rows[j] = i;
			}
			
			status = CPXchgcoeflist (env, lp, nzcnt, rows, rmatind, rmatval);
			if (status) {
				fprintf(stderr, "\nresetconstraints: ");
				fprintf(stderr, "\nCPXchgcoeflist failure at row %d. status = %d\n", i, status);
				return status;
			}
			
		} // End if
		
		
	} // End for loop
	
	
	
#ifdef DEBUG_FUNC
	//Write sub prob to file in lp format
	status = CPXwriteprob(env, lp, "prob.lp", NULL);
	if ( status ) {
		fprintf (stderr, "resetconstraints:\n");
		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
		return status;
	}
#endif
	
	
	free(rmatbeg);
	free(rmatind);
	free(rmatval);
	free(rows);
	free(sense);
	
	free(rowindex);
	free(rhs);
	free(newsense);
	
	
	// printf("\n Done resetting constraints...\n");
	
	return 0;
	
}

//************************** End resetconstraints function *****************************/

int
addOptColToMasterProbLP(CPXENVptr env, CPXLPptr lp_master, masterproblem_t *masterprobPtr)
/**
 * This function adds an optimality (theta) column to master problem CPLEX lp object
 * for adding optimality cuts as required by the D^2 algorithm
 * @param env  a pointer to the CPLEX environment
 * @param lp_master  a pointer to the CPLEX LP master problem object
 * @param masterprobPtr a pointer to the master problem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
//******************************************************************************
// ******************** Add theta (opt) variable to master problem ***************
{
	double obj[1];
	char **colname_theta;
	obj[0] = 1;
	int i;
	int status;
	char mastername_lp[NAMELEN];
	double *theta_lb;
	double *theta_ub;
	
	strcpy( mastername_lp, "mastermip.lp");
	theta_lb = (double *)malloc(sizeof(double));
	theta_ub = (double *)malloc(sizeof(double));
	
	if(theta_lb == NULL || theta_ub == NULL) {
		fprintf(stderr, "\n addOptColToMasterProbLP:\n");
		fprintf(stderr, " Failed to allocate memory to lower or upper bound for theta column");
		return 1;
	}
	//Make the lower bound for theta -infinity and upper bound +infty so the variable is unrestricted in sign
	theta_lb[0] = INIT_LOWER_BOUND;
	theta_ub[0] = CPX_INFBOUND;
	
	colname_theta = (char**)malloc(sizeof(char*));
	if ( colname_theta == NULL ) {
		fprintf (stderr, "\n addOptColToMasterProbLP:\n");
		fprintf (stderr, " Failed to allocate memory to master col name ""theta"".\n");
		return 1;
	}
	
	colname_theta[0] = (char*)malloc(NAMELEN*sizeof(char*));
	if ( colname_theta[0] == NULL ) {
		fprintf (stderr, "\n addOptColToMasterProbLP:\n");
		fprintf (stderr, " Failed to allocate memory to master problem optimality col name ""theta"".\n");
		return 1;
	}
	strcpy(colname_theta[0], "theta");
	
	status = CPXnewcols(env, lp_master, 1, obj, theta_lb, theta_ub, NULL, colname_theta);
	if ( status ) {
		fprintf (stderr, "\n addOptColToMasterProbLP():\n");
		fprintf (stderr, " Failed to add new theta column (var) to master LP.\n");
		return(status);
	}
	
	
	i = masterprobPtr->ncols;
	strcpy(masterprobPtr->colnames[i], "theta");
	masterprobPtr->ncols +=1;
	
	//********************** Write prob to file in lp format **************************
#ifdef DEBUG_FUNC
	status = CPXwriteprob(env, lp_master, mastername_lp, NULL);
	if ( status ) {
		fprintf (stderr, " addOptColToMasterProbLP:\n");
		fprintf (stderr, "Failure to write lp_master to file, error %d.\n", status);
		return(status);
	}
#endif
	
	// Free memory
	if(colname_theta[0] != NULL) {
		free(colname_theta[0]);
		colname_theta[0]= NULL;
	}
	if ( colname_theta != NULL ){
		free (colname_theta);
		colname_theta = NULL;
	}
	if( theta_lb != NULL) {
		free(theta_lb);
		theta_lb = NULL;
	}
	if(theta_ub != NULL) {
		free(theta_ub);
		theta_ub = NULL;
	}
	
	return (0);
	
}
//************************** // End addOptColToMasterProbLP function *****************************/


int setupSubProbMip(CPXENVptr env, CPXLPptr lp_submip, subproblem_t *subprobPtr, int nrows_master,
					int ncols_master, int nrows_sub, int ncols_sub,  FILE *fpout)
/**
 * Dinakar Gade - modified from Yang Yuan/Lewis Ntaimo's code
 * Sets up the subproblem mip CPLEX LP obj by deleting the first stage rows and cols
 * from the core file LP. This LP object will be used for solving subproblem mips for
 * upper bounding. This function also updates the subproblem data structure number
 * of rows in the mip obj. D^2 cuts will sequentially be added to this model in main.
 * @param lp_submip pointer to the subproblem mip model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_mp   integer indicating number of master problem columns (vars)
 * @param nrows_sub  integer indicating number of subproblem rows (constraints)
 * @param ncols_sub   integer indicating number of subproblem columns (vars)
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise returns a nonzero integer
 */
{
	
	int status;
	//int surplus; // To check for sufficiency of the aloocated array size
	//int i,
	int j;
	char mipname[20];
	//int curr_row;
	
	strcpy(mipname, "subprobmip1.lp");
	
	//Create subproblem lp by deleting stage_1 cols and rows from core file lp
	
	//***************** Delete stage-1 rows (constraints) *********************************
	if (nrows_master > 0) {
		status = CPXdelrows(env, lp_submip, 0, nrows_master-1);
		if ( status ) {
			fprintf(stderr,"\nFunction setupSubProbMip(...): \n");
			fprintf (stderr,"Failure to delete suproblem lp rows, error %d.\n", status);
			return(status);
		}
	}
	
	//********************** Delete stage-1 columns (vars) **********************************
	//fprintf(stdout, "Number of original sub cols = %d\n", CPXgetnumcols(env, lp_submip));
	status = CPXdelcols(env, lp_submip, 0, ncols_master-1);
	if ( status ) {
		fprintf(fpout, "\nFunction setupSubProbMip(...): \n");
		fprintf (stderr, "Failure to delete lp_core cols, error %d.\n", status);
		return(status);
	}
	//fprintf(stdout, "Number of sub cols after first stage cols deletion = %d\n", CPXgetnumcols(env, lp_submip));
	
	subprobPtr->nrows_mip = CPXgetnumrows (env, lp_submip);
	subprobPtr->nrowsWmip = subprobPtr->nrows_mip;
	subprobPtr->nrows = subprobPtr->nrows_mip;
	
	//////////////////// ADDED Jan 7, 2002 //////////////////
	//            Reset the <= contraints to >=            //
	/////////////////////////////////////////////////////////
	status = resetconstraints(env, lp_submip);
	if ( status ) {
		fprintf (stderr, "setupSubProbMip(...):\n");
		fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
		return status;
	}
	////////////////////////////////////////////////////////
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Add explicit feasibility column to the model to enable complete recourse //
	// as require by the D^2 algorithm                                          //
	//////////////////////////////////////////////////////////////////////////////
	
	status = addFeasColToSubProbLP(env, lp_submip, subprobPtr);
	if ( status ) {
		fprintf (stderr, "setupSubProbMip:\n");
		fprintf (stderr, "Failure to add feasibility column to lp_submip, error %d.\n", status);
		return(status);
	}
	//********************** Write prob to file in lp format **************************
#ifdef DEBUG_FUNC
	status = CPXwriteprob(env, lp_submip, mipname, NULL);
	if ( status ) {
		fprintf (stderr, "Function setupSubProbMip():\n");
		fprintf (stderr, "Failure to write lp_submip to file, error %d.\n", status);
		return(status);
	}
#endif
	
	//****** Initialize the indices array corresponding to the constraints *********//
	//************* for which the rhs coefs are to be changed *********************//
	for (j = 0; j < subprobPtr->nrows_mip; j++) {
		subprobPtr->indices_rowmip[j] = j;
	}
	// Remember to update this array every time a cut is added to the subproblem!!!
#ifdef DEBUG_FUNC
	fprintf (stdout, "Function setupSubProbMip:\n");
	for(i = 0; i < subprobPtr->nrows_mip; i++){
		fprintf(stdout, "subprobPtr->indices_rowmip[%d] = %d \n", i, subprobPtr->indices_rowmip[i]);
	}
#endif
	
	subprobPtr->ncols = CPXgetnumcols (env, lp_submip);
	
	return (0);
} // End function setupSubProbMip

//************************** // End setupSubProbMip function *****************************/

int addFeasColToSubProbLP(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr)
/**
 * This function adds a column to subprobl CPLEX lp object to enable
 * complete resourse as required by the D^2 algorithm
 * This column is penalized in the object function
 * @param env  a pointer to the CPLEX environment
 * @param lp_sub  a pointer to the CPLEX LP subproblem object
 * @param subprobPtr a pointer to the subproblem data structure
 * @return returns 0 on success, otherwise returns a nonzero value
 */
{
	int i, nzcnt;
	int status;
	
	double *obj = NULL;
	int    *cmatbeg = NULL;
	int    *cmatind = NULL;
	double *cmatval= NULL;
	double *lb = NULL;
	double *ub = NULL;
	char **colname = NULL;
	int cur_numcols;
	
	nzcnt = subprobPtr->nrows_mip;
	
	// Allocate memory
	obj = (double*)malloc(sizeof(double));
	cmatbeg = (int*)malloc(sizeof(int));
	cmatind = (int*)malloc(nzcnt*sizeof(int));
	cmatval = (double*)malloc(nzcnt*sizeof(double));
	lb = (double*)malloc(sizeof(double));
	ub = (double*)malloc(sizeof(double));
	colname = (char**)malloc(sizeof(char*));
	
	if ( obj == NULL || cmatbeg == NULL || cmatind == NULL ||
		cmatval == NULL || lb == NULL || ub == NULL ||
		colname == NULL) {
		fprintf (stderr, "\nd2algFuncs:\n");
		fprintf (stderr, "Function addFeasColToSubProbLP: \n");
		fprintf (stderr, " Failed to allocate memory to feasibility column variables\n");
		fprintf (stderr, " for subproblem lp. \n");
		return 1;
	}
	colname[0] = (char*)malloc(NAMELEN*sizeof(char*));
	if (colname == NULL) {
		fprintf (stderr, "\nd2algFuncs:\n");
		fprintf (stderr, "Function addFeasColToSubProbLP: \n");
		fprintf (stderr, " Failed to allocate memory to feasibility colname variable\n");
		fprintf (stderr, " for subproblem lp. \n");
		return 1;
	}
	
	// Add a feasibility column to subproblem lp data to enable complete resourse
	obj[0] = PENALTY; // High penalty cost: see config.h
	
	//Add feasibility col name to sub prob data
	strcpy(colname[0], "slack");
	
	cur_numcols = CPXgetnumcols (env, lp_sub);
	
	strcpy(subprobPtr->colnames[cur_numcols], "slack");
	
	cmatbeg[0] = 0;
	
	//nzcnt = subprobPtr->nrows_mip;
	
	for (i = 0; i < nzcnt; i++) {
		cmatind[i] = i;
		cmatval[i] = 1;
	}
	
	lb[0] = 0  ;
	ub[0] = CPX_INFBOUND ;
	
	status = CPXaddcols(env, lp_sub, 1, nzcnt, obj, cmatbeg, cmatind, cmatval,
						lb, ub, colname);
	if (status) {
		fprintf (stderr, "\nd2algFuncs:\n");
		fprintf (stderr, "Function addFeasColToSubProbLP: \n");
		fprintf (stderr, " Failed to add feasibility column to subproblem lp. \n");
		return status;
	}
	
	
	// Free up the memory allocated to the arrays
	if ( obj != NULL ){
		free (obj);
		obj = NULL;
	}
	if ( cmatbeg != NULL ){
		free (cmatbeg);
		cmatbeg = NULL;
	}
	if ( cmatind != NULL ){
		free (cmatind);
		cmatind = NULL;
	}
	//Added by Dinakar Gade
	if(cmatval != NULL) {
		free(cmatval);
		cmatval = NULL;
	}
	if ( lb != NULL ){
		free (lb);
		lb = NULL;
	}
	if ( ub != NULL ){
		free (ub);
		ub = NULL;
	}
	if(colname[0] != NULL) {
		free(colname[0]);
		colname[0] = NULL;
	}
	if ( colname != NULL ){
		free (colname);
		colname = NULL;
	}
	
	return status;
	
} // ************ addFeasColToSubProbLP *********************//

int loadSubProblemData(CPXENVptr env, CPXLPptr lp_sub, subproblem_t *subprobPtr, int nrows_master,
					   int ncols_master, int nrows_sub, int ncols_sub, stochfile_info *stochdataPtr,
					   FILE *fpout)
/**
 * Sets up the subproblem lp CPLEX LP obj by deleting the first stage rows and cols
 * from the core file LP. This LP object will be used for solving subproblem lps.
 * The function also extracts subproblem data from the lp and puts it in subproblem data
 * structure to be used in creating D^2 cuts. It is assumed that the subproblem data structure
 * has already been allocated memory
 * @param lp pointer to the lp model
 * @param subprobPtr  pointer to the subproblem lp data structure
 * @param nrows_master   integer indicating number of master problem rows (constraints)
 * @param ncols_mp   integer indicating number of master problem columns (vars)
 * @param nrows_sub  integer indicating number of subproblem rows (constraints)
 * @param ncols_sub   integer indicating number of subproblem columns (vars)
 * @param stochdataPtr  pointer to the subproblem lp stoch data structure
 * @param fpout output file pointer
 * @return 0 data extraction is successfully, otherwise return a nonzero integer
 */
{
	
	int status;
	int surplus; // To check for sufficiency of the allocated array size
	int i, j;
	char lpname[20];
	//int curr_row;
	//char *ub = NULL;
	//int *indices = NULL;
	int numscens = stochdataPtr->nscens;
	
	strcpy(lpname, "subproblp.lp");;
	
	
	//Create subproblem lp by deleting stage_1 cols and rows from core file lp
	
	//***************** Delete stage-1 rows (constraints) *********************************
	if (nrows_master > 0) {
		status = CPXdelrows(env, lp_sub, 0, nrows_master-1);
		if ( status ) {
			fprintf(stderr,"\nFunction loadSubProblemData(...): \n");
			fprintf (stderr,"Failure to delete suproblem lp rows, error %d.\n", status);
			return(status);
		}
	}
	// Get the sense vector for the range of contraints for the subproblem       //
	status = CPXgetsense(env, lp_sub, subprobPtr->sense, 0, subprobPtr->nrows-1);
	if ( status ) {
		fprintf (stderr, "loadSubProblemData:\n");
		fprintf (stderr, "Failure to get the senses from the subproblem lp, error %d.\n", status);
		return(status);
	}
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadSubProblemData:\n");
	fprintf(stdout, "Subproblem senses\n");
	for (j = 0; j < subprobPtr->nrows; j++) {
		fprintf(stdout, "subprobPtr->sense[%d:] %c\n", j, subprobPtr->sense[j]);
	}
#endif
	
	//////////////////// ADDED Jan 7, 2002 //////////////////
	//            Reset the <= contraints to >=            //
	/////////////////////////////////////////////////////////
	status = resetconstraints(env, lp_sub);
	if ( status ) {
		fprintf (stderr, "loadSubProblemData(...):\n");
		fprintf (stderr, "Failure to reset master MIP constraints, error %d.\n", status);
		return status;
	}
	
#ifdef DEBUG_FUNC
	fprintf (stdout, "writing subprob lp:\n");
	// Write sub prob to file in lp format
	status = CPXwriteprob(env, lp_sub, "prob3.lp", NULL);
	if ( status ) {
		fprintf (stderr, "resetconstraints:\n");
		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
		return status;
	}
#endif
	////////////////////////////////////////////////////////
	
	
	
	
	
	//**********Extract the T matrix from the core file and store it in T *****************
	status = CPXgetcols(env, lp_sub, &subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
						subprobPtr->cmatind_T, subprobPtr->cmatval_T,
						subprobPtr->cmatspace_T, &surplus, 0, ncols_master-1);
	if ( status ) {
		fprintf (stderr, "Failure to extraxt T from lp_sub from T(w) matrix, error %d.\n", status);
		fprintf (stderr, "Surplus value is: %d.\n", surplus);
		return(status);
	}
	
	//******************** Print the T matrix ******************************************
#ifdef DEBUG_FUNC
	fprintf(stderr, "\nFunction loadSubProblemData(...): \n");
	fprintf(stderr, "subprobPtr->nzcnt_T = %d\n", subprobPtr->nzcnt_T);
	fprintf(stderr, "*******The T matrix is********: \n");
	printSparseMatrix(ncols_master, subprobPtr->nzcnt_T, subprobPtr->cmatbeg_T,
					  subprobPtr->cmatbeg_T, subprobPtr->cmatind_T, subprobPtr->cmatval_T, stderr);
#endif
	
	
	//********************** Delete stage-1 columns (vars) **********************************
	status = CPXdelcols(env, lp_sub, 0, ncols_master-1);
	if ( status ) {
		fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
		fprintf (stderr, "Failure to delete lp_sub cols, error %d.\n", status);
		return(status);
	}
	
	//************* Get and store the ctype (var type) array from the subproblem lp ********
	status = CPXgetctype(env, lp_sub, subprobPtr->ctype, 0, ncols_sub-1);
	if ( status ) {
		fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
		fprintf (stderr, "*Failure to get ctype from lp_sub, error %d.\n", status);
		return (status);
	}
	//fprintf(stdout, "@@@ CTYPE SUB @@@\n");
	//for (i = 0; i < ncols_sub; i++) {
	//    fprintf(stdout, "Y[%d] = %c\n",i,subprobPtr->ctype[i]);
	//}
	//fprintf(stdout, "@@@ CTYPE SUB @@@\n");
#ifdef DEBUG_FUNC
	fprintf(fpout, "\nFunction loadSubProblemData(...): \n");
	fprintf(stdout,"Subproblem ctype: \n");
	for (j = 0; j < ncols_sub; j++)
		fprintf(stdout, "col %d: %c\n", j, subprobPtr->ctype[j]);
#endif
	
	
	
	// Store a vector of continuous ctypes for the linear relaxation
	// and ctype indices
	// Include feasibility column to be added soon
	subprobPtr->ctype[ncols_sub] = 'C';
	for (j = 0; j < ncols_sub+1; j++){
		subprobPtr->ctype_ctns[j] = 'C';
		subprobPtr->indices_ctype[j] = j;
	}
	
	
	/////////////////// Updates
	
	subprobPtr->nrows   = CPXgetnumrows (env, lp_sub);
	subprobPtr->nrows_T = CPXgetnumrows (env, lp_sub);
	stochdataPtr->nrows = CPXgetnumrows (env, lp_sub);;	// Initialize: VERY IMPORTANT!!!!!!
#ifdef DEBUG_FUNC
	printf("subprobPtr->nrows = %d \n", subprobPtr->nrows);
	printf("subprobPtr->nrows_T = %d \n", subprobPtr->nrows_T);
	printf("stochdataPtr->nrows = %d \n", stochdataPtr->nrows);
	
	//Write sub prob to file in lp format
	status = CPXwriteprob(env, lp_sub, "prob2.lp", NULL);
	if ( status ) {
		fprintf (stderr, "resetconstraints:\n");
		fprintf (stderr, "Failure to write sub problem lp to file, error %d.\n", status);
		return status;
	}
#endif
	
	
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Add explicit feasibility column to the model to enable complete recourse //
	// as require by the D^2 algorithm                                          //
	//////////////////////////////////////////////////////////////////////////////
	
	status = addFeasColToSubProbLP(env, lp_sub, subprobPtr);
	
	if ( status ) {
		fprintf (stderr, "loadSubProblemData:\n");
		fprintf (stderr, "Failure to add binary contraints to lp_sub, error %d.\n", status);
		return(status);
	}
	
	///////////////////////////////////////////////////////////////////////
	// Add explicit binary contraints to the model for the D^2 algorithm //
	// cut generation.                                                   //
	///////////////////////////////////////////////////////////////////////
	
	// status = addBinaryConstrsSubProbLP(env, lp_sub, subprobPtr);  //UnCommented by Dinakar Gade - No need to add binary constraints to subproblem
	//if ( status ) {
	//    fprintf (stderr, "loadSubProblemData:\n");
	//  	fprintf (stderr, "Failure to add binary contraints to lp_sub, error %d.\n", status);
	//  	return(status);
	// }
	
	
	// Make updates
	subprobPtr->nrows_total  = CPXgetnumrows (env, lp_sub);
	subprobPtr->nrowsW       = CPXgetnumrows (env, lp_sub);
	ncols_sub                = CPXgetnumcols (env, lp_sub);
	subprobPtr->ncols        = ncols_sub;
	//fprintf(stdout, "@@@ CTYPE SUB @@@\n");
	//for (i = 0; i < ncols_sub; i++) {
	//    fprintf(stdout, "Y[%d] = %c\n",i,subprobPtr->ctype[i]);
	//}
	
	
	//********************** Write prob to file in lp format **************************
#ifdef DEBUG_FUNC
	status = CPXwriteprob(env, lp_sub, lpname, NULL);
	if ( status ) {
		fprintf (stderr, "loadSubProblemData:\n");
		fprintf (stderr, "Failure to write lp_sub to file, error %d.\n", status);
		return(status);
	}
#endif
	
	//************* Extract the obj from the lp_subfrom the core file ****************
	status = CPXgetobj(env, lp_sub, subprobPtr->obj, 0, ncols_sub-1);
	if ( status ) {
		fprintf (stderr, "loadSubProblemData:\n");
		fprintf (stderr, "Failure to get obj from lp_sub, error %d.\n", status);
		return(status);
	}
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadSubProblemData:\n");
	fprintf(stdout, "Subproblem obj\n");
	for (j = 0; j < subprobPtr->ncols; j++) {
		fprintf(stdout, "subprobPtr->obj[%d]: %f\n", j, subprobPtr->obj[j]);
	}
#endif
	
	// Get the rhs vector for the range of contraints for the subproblem  //
	// and initialize each scenario rhs with these values                 //
	status = CPXgetrhs(env, lp_sub, subprobPtr->rhs, 0, nrows_sub-1);
	if ( status ) {
		fprintf (stderr, "Failure to get the rhs from the subproblem lp, error %d.\n", status);
		return(status);
	}
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadSubProblemData:\n");
	fprintf(stdout, "Subproblem rhs\n");
	for (j = 0; j < subprobPtr->nrows; j++) {
		fprintf(stdout, "subprobPtr->rhs[%d]: %f\n", j, subprobPtr->rhs[j]);
	}
#endif
	
	// Initialize each scenario rhs with these values:
	// Initialize the scenario subproblem rhs(w)
	for (i = 0; i < numscens; i++){
		for (j = 0; j < subprobPtr->nrows; j++)
			stochdataPtr->rhs[i][j] = subprobPtr->rhs[j];
	}
	
	
	
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadSubProblemData:\n");
	for (i = 0; i < numscens; i++){
		fprintf(stdout, "Subproblem rhs for scenario %d\n", i);
		for (j = 0; j < subprobPtr->nrows; j++) {
			fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n", i, j, stochdataPtr->rhs[i][j]);
		}
	}
#endif
	
	//Get the lower bounds on the subproblem cols (vars)
	status = CPXgetlb(env, lp_sub, subprobPtr->lb, 0, ncols_sub-1);
	if ( status ) {
		fprintf (stderr, "Failure to get lower bounds lb from lp_sub, error %d.\n", status);
		return (status);
	}
	//Get the upper bounds on the subproblem cols (vars)
	status = CPXgetub(env, lp_sub, subprobPtr->ub, 0, ncols_sub-1);
	if ( status ) {
		fprintf (stderr, "Failure to get upper bounds ub from lp_sub, error %d.\n", status);
		return (status);
	}
	
	
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadSubProblemData:\n");
	fprintf(stdout, "\n\n\n  Subproblem lp column ub and lb:\n");
	for(i = 0; i < ncols_sub; i++){
		fprintf(stdout, "lb[%d] = %f	ub[%d] = %f \n", i,subprobPtr->lb[i],i,subprobPtr->ub[i]);
	}
#endif
	
	
	//****** Initialize the indices array corresponding to the constraints *********//
	//************* for which the rhs coefs are to be changed *********************//
	for (j = 0; j < subprobPtr->nrows; j++) {
		subprobPtr->indices_row[j] = j;
	}
	// Remember to update this array every time a cut is added to the subproblem!!!
#ifdef DEBUG_FUNC
	fprintf (stdout, "loadSubProblemData:\n");
	for(i = 0; i < subprobPtr->nrows; i++){
		fprintf(stdout, "subprobPtr->indices_row[%d] = %d \n", i, subprobPtr->indices_row[i]);
	}
#endif
	
	return (0); // Successful return
	
} //************************** End loadSubProblemData function *****************************

int changeSlacksToIntegers(CPXENVptr env, CPXLPptr lp_submip) {
	int status = 0,i;
	int ncols;
	char *coltype = NULL;
	int *indices = NULL;
	
	ncols = CPXgetnumcols(env, lp_submip);
	coltype = (char *)malloc(ncols*sizeof(char));
	indices = (int *) malloc(ncols*sizeof(int));
	if (coltype == NULL || indices == NULL) {
		fprintf(stderr, "In changeSlacksToIntegers():\n\tError allocating coltype");
		return 1;
	}
	status = CPXgetctype(env, lp_submip, coltype, 0, ncols-1);
	if (status) {
		fprintf(stderr, "In changeSlacksToIntegers():\n\tCPLEX error getting ctype. %d\n",status);
		return 1;
	}
	for (i = 0; i < ncols; i++) {
		indices[i] = i;
		if (coltype[i] == 'C') {
			coltype[i] = 'I';
		}
	}
	status = CPXchgctype(env, lp_submip, ncols, indices, coltype);
	if (status) {
		fprintf(stderr, "In changeSlacksToIntegers():\n\tCPLEX error changing ctype. %d\n",status);
		return 1;
	}
	if (coltype != NULL) {
		free(coltype);
		coltype = NULL;
	}
	if (indices != NULL) {
		free(indices);
		indices = NULL;
	}
	
	return status;
}
//************************** End changeSlacksToIntegers function *****************************

int addObjConstraint2(CPXENVptr env, CPXLPptr lp) {
	/*
	 * Dinakar
	 * Add the objective function as a constrain by adding y0 column as the FIRST variable and adding y0 - cx = 0
	 * Needed for Gomory cut computation
	 * as  the FIRST constraint
	 * @param env  - CPLEX environment
	 * @param lp   - CPLEX lp
	 * @return     - 0 if everything is o.k., > 0 otherwise
	 */
	int status = 0,i, nzcnt, surplus, count,cur_colnamespace, cur_rownamespace;
	int *rmatbeg = NULL, *rmatind = NULL;
	double *rhs = NULL, *obj = NULL, *rmatval = NULL, *lb = NULL, *ub = NULL;;
	char **colnames = NULL, **rownames = NULL, *colnamestore = NULL, *rownamestore = NULL, *sense = NULL;
	int ncols, nrows, numnz, numnzobj, newnz;
	
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);
	
	//For rows and column data, first allocate one more than necessary since we'll add y0 column and y0-cx = 0 row
	
	obj          = (double *)malloc((ncols + 1)*sizeof(double));
	rhs          = (double *)malloc((nrows + 1)*sizeof(double));
	lb           = (double *)malloc((ncols + 1)*sizeof(double));
	ub           = (double *)malloc((ncols + 1)*sizeof(double));
	sense        = (char *)  malloc((nrows + 1)*sizeof(char));
	
	
	if (obj == NULL || rhs == NULL || lb == NULL || ub == NULL || sense == NULL) {
		fprintf(stderr, "In addObjConstraint():\n\tError allocating memory to obj...sense\n");
		return (1);
	}
	
	
	status = CPXgetcolname (env, lp, NULL, NULL, 0, &surplus, 0,ncols-1);
	
	if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))  {
		fprintf (stderr, "Could not determine amount of space for column names.\n");
		return (1);
	}
	cur_colnamespace = -surplus;
	if ( cur_colnamespace > 0 ) {
		colnames      = (char **)  malloc ((ncols+1)*sizeof(char *));
		colnamestore  =  (char *)  malloc ((cur_colnamespace + NAMELEN)*sizeof(char));
		if ( colnames == NULL || colnamestore == NULL   ) {
			fprintf (stderr, "Failed to get memory for column names.\n");
			status = -1;
			return status;
		}
	}
	
	//Get Column names starting from 1 since column 0 will contain "y_0"
	status = CPXgetcolname(env, lp, (colnames + 1), (colnamestore + NAMELEN), cur_colnamespace, &surplus, 0, ncols-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error getting colnames. CPLEX STATS = %d\n",status);
		return (1);
	}
	colnames[0] = colnamestore;
	strcpy(colnames[0], "y_0");
	
	status = CPXgetrowname(env, lp, NULL, NULL, 0, &surplus, 0,nrows-1);
	
	if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))  {
		fprintf (stderr, "Could not determine amount of space for row names.\n");
		return (1);
	}
	cur_rownamespace = -surplus;
	if ( cur_rownamespace > 0 ) {
		rownames      =  (char **)  malloc ((nrows+1)*sizeof(char *));
		rownamestore  =  (char *)  malloc ((cur_rownamespace + NAMELEN)*sizeof(char));
		if ( rownames == NULL || rownamestore == NULL   ) {
			fprintf (stderr, "Failed to get memory for column names.\n");
			status = -1;
			return status;
		}
	}
	
	
	//Get Row names starting from 1 since column 0 will contain "objrow"
	status = CPXgetrowname(env, lp, (rownames + 1), (rownamestore + NAMELEN), cur_rownamespace, &surplus, 0, nrows-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error getting colnames. CPLEX STATS = %d\n",status);
	}
	rownames[0] = rownamestore;
	strcpy(rownames[0], "objrow");
	
	//Get the required problem data shifted to account for the new row/column
	status = CPXgetobj(env, lp, (obj+1), 0, ncols-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error extracting obj\n. CPLEX STATUS = %d\n",status);
		return (1);
	}
	status = CPXgetrhs(env, lp, (rhs+1), 0, nrows-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error extracting rhs\n. CPLEX STATUS = %d\n",status);
		return (1);
	}
	status = CPXgetub(env, lp, ub+1, 0, ncols-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error extracting ub\n. CPLEX STATUS = %d\n",status);
		return (1);
	}
	status = CPXgetlb(env, lp, lb+1, 0, ncols-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error extracting lb\n. CPLEX STATUS = %d\n",status);
		return (1);
	}
	status = CPXgetsense(env, lp, sense+1, 0, nrows-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error getting sense. CPLEX STATS = %d\n",status);
		return (1);
	}
	/*Start counting the total number of new nonzeros that will be added to constraints
	 * These new nonzeros correspond to the constraint y0 - cx = 0;
	 */
	numnzobj = numnz + 1; //y_0 element
	for (i = 1; i < ncols + 1; i++) {
		if (myfabs(obj[i]) > EPSILON_TOLERANCE) {
			numnzobj++;
		}
	}
	newnz = numnzobj - numnz; //newnz is the new number of nonzeros added to the constraint matrix
	//fprintf(stdout, "Old Numz = %d, New Numnz = %d, new nz = %d\n",numnz, numnzobj,newnz);
	rmatbeg = (int *)calloc((nrows+1), sizeof(int));
	rmatind = (int *)calloc((numnzobj), sizeof(int));
	rmatval = (double *)calloc((numnzobj), sizeof(double));
	
	//Extract rows
	status = CPXgetrows(env, lp, &nzcnt, (rmatbeg + 1), (rmatind + newnz), (rmatval + newnz), numnzobj, &surplus, 0, nrows-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\tError extracting rows\n");
		return (1);
	}
	
	//Set y0's row/colum data
	rhs[0] = 0.0;
	rmatval[0] = 1.0;
	rmatind[0] = 0;
	lb[0] = -CPX_INFBOUND;
	ub[0] =  CPX_INFBOUND;
	sense[0] = 'E';
	obj[0] = 1.0;
	
	//Populate the sparse matrices for the new elements that are being added
	rmatbeg[0] = 0;
	for (i = 1; i < nrows + 1; i++) {
		rmatbeg[i] = rmatbeg[i] + (numnzobj - numnz);
	}
	for (i = newnz; i < numnzobj; i++) {
		rmatind[i]++;
	}
	
	count = 1;
	for (i = 0; i < ncols; i++) {
		if (myfabs(obj[i+1]) > EPSILON_TOLERANCE) {
			rmatind[count] = i+1;
			rmatval[count] = -obj[i+1];
			obj[i+1] = 0.0; //Set remaining objective coeffs to zero because objective now is \min y_0
			count++;
		}
	}
	//Delete the columns and rows
	status = CPXdelcols(env, lp, 0, ncols-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error deleting cols\n. CPLEX STATUS = %d\n",status);
		return (1);
	}
	status = CPXdelrows(env, lp, 0, nrows-1);
	if (status) {
		fprintf(stderr, "In addObjConstraint():\n\t Error deleting rows\n. CPLEX STATUS = %d\n",status);
		return (1);
	}
	//Add new columns
	status = CPXaddcols(env, lp, ncols+1, 0, obj, NULL, NULL, NULL, lb, ub, colnames);
	if (status) {
		fprintf(stderr, "Error adding cols in addObjConstraint(): CPLEX STATUS = %d\n",status);
		return (1);
	}
	status = CPXaddrows(env, lp, 0, nrows+1, numnzobj, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if (status) {
		fprintf(stderr, "Error adding rows in addObjConstraint(): CPLEX STATUS = %d\n",status);
		return (1);
	}
	//status = CPXwriteprob(env, lp, "newrows.lp", "lp");
	
	//freeup memory
	if (rmatbeg != NULL) {
		free(rmatbeg);
		rmatbeg = NULL;
	}
	if (rmatind != NULL) {
		free(rmatind);
		rmatind = NULL;
	}
	if (rmatval != NULL) {
		free(rmatval);
		rmatval = NULL;
	}
	if(rhs != NULL) {
		free(rhs);
		rhs = NULL;
	}
	if (obj != NULL) {
		free(obj);
		obj = NULL;
	}
	if (lb != NULL) {
		free(lb);
		lb = NULL;
	}
	if (ub != NULL) {
		free(ub);
	}
	if (sense != NULL) {
		free(sense);
		sense = NULL;
	}
	if (colnames != NULL) {
		free(colnames);
		colnames = NULL;
	}
	if (rownames != NULL) {
		free(rownames);
		rownames = NULL;
	}
	if (colnamestore != NULL) {
		free(colnamestore);
	}
	if (rownamestore != NULL) {
		free(rownamestore);
	}
	return status;
}

int convert_col_to_row(int *cmatbeg, int *cmatcnt, int *cmatind, double *cmatval, int *rmatbeg, int *rmatcnt, int *rmatind, double *rmatval, int nrows, int ncols, int numnz) {
	/**
	 * Convert the T-matrix in column sparse matrix into T-matrix in row sparse matrix
	 * @param cmatbeg,cmatcnt, cmatind, cmatval                 - Column sparse matrix representation of the T-matrix
	 * @param rmatbeg,rmatcnt, rmatind, rmatval                 - Row sparse matrix representation that will be populated by this method
	 * @param nrows                                             - Number of rows in the matrix
	 * @param ncols                                             - Number of columns in the matrix
	 * @param numnz                                             - Number of nonzeros in the matrix
	 * @return                                                    0 = success, nonzero = fail
	 */
	int status = 0,i,j,count;
	
	double **A = NULL;
	
	if (cmatbeg == NULL || cmatcnt == NULL || cmatind == NULL || cmatval == NULL || rmatbeg == NULL || rmatind == NULL || rmatval == NULL) {
		fprintf(stderr, "In function convert_col_to_row():\n\t NullPointer Error for input arguments\n");
		return (1);
	}
	A = (double **)malloc(nrows*sizeof(double *));
	if (A == NULL) {
		fprintf(stderr, "In function convert_col_to_row():\n\tError allocating memory to A\n");
		return (1);
	}
	for (i = 0; i < nrows; i++) {
		A[i] = (double *)malloc(ncols*sizeof(double));
		if (A[i] == NULL) {
			fprintf(stderr, "In function convert_col_to_row():\n\tError allocating memory to A[i]\n");
			return (1);
		}
		for (j = 0; j < ncols; j++) {
			A[i][j] = 0.0;
		}
	}
	//Copy col to dense
	for (j = 0; j < ncols; j++) {
		for (i = 0; i < cmatcnt[j]; i++) {
			A[cmatind[cmatbeg[j] + i]][j] = cmatval[cmatbeg[j] + i];
		}
	}
	/*
	 for (i = 0; i < nrows; i++) {
	 for (j = 0; j < ncols; j++) {
	 fprintf(stdout, "%f ",A[i][j]);
	 }
	 fprintf(stdout, "\n");
	 }*/
	//Convert back to rowsparse
	count = 0;
	for (i = 0; i < nrows; i++) {
		rmatcnt[i] = 0;
		for (j = 0; j < ncols; j++) {
			if (fabs(A[i][j]) > 0.0) {
				rmatcnt[i]++;
				rmatind[count] = j;
				rmatval[count] = A[i][j];
				count++;
			}
		}
	}
	rmatbeg[0] = 0;
	for (i = 1; i < nrows; i++) {
		rmatbeg[i] = rmatbeg[i-1] + rmatcnt[i-1];
	}
	for (i = 0; i < nrows; i++) {
		if (A != NULL && A[i] != NULL) {
			free(A[i]);
			A[i] = NULL;
		}
	}
	if (A != NULL) {
		free(A);
		A = NULL;
	}
	return status;
}

//************************** End convert_col_to_row() *****************************/
int
loadStochFile(stochfile_info *stochdataPtr, subproblem_t *subprobPtr, char *filename, int *random_T, int *random_obj, int ncols_master, FILE *fpout)
/**
 * Reads the STOCH file data and puts the data in the structure timedataPtr
 * @param timedataPtr pointer to data structure to hold STOCH file data
 * @param subprobPtr pointer to master problem data structure
 * @param filename  STOCH file name
 * @param probname  problem name as read from the TIME file
 * @param random_T a pointer to an integer value: 0 indicates that the technology
 *                 matrix T(w) = T is constant, 1 indicates that it is random T(W)
 * @param random_obj a pointer to an integer value: 0 indicates that the scenario subprob
 *                 objective is constant, 1 indicates that it is random
 * @param ncols_master number of columns in the master problem
 * @param fpout output file pointer
 * @return 0 if STOCH file is successfully read, otherwise return a nonzero integer
 */

{
	
	FILE *stoch;
	char field1[NAMELEN], field2[NAMELEN], field3[NAMELEN], field4[NAMELEN];
	char buffer[LENGTH];
	char myfilename[LENGTH];
	int i, j, col, scenario;  // counters
	double tempDouble;
	int rowcount = 0;
	//int status;
	double sum;
	int numscens = stochdataPtr->nscens;
	int scen_index;
	
	
	// Open the TIME file
	stoch = fopen(filename,"r" );
	if(stoch == NULL) {
		fprintf(stderr, "\nloadStochFile(...): \n");
		fprintf(stderr, " Could not open the TIME file %s for reading!\n", filename);
		return(1);
	}
	
	
	//************************* Read the TIME file ************************************
	// Read first line: e.g. STOCH	example
	if (fgets (buffer, LENGTH, stoch) != NULL) {
		sscanf(buffer, "%s%s", myfilename, stochdataPtr->probname);
#ifdef DEBUG_FUNC
		fprintf(fpout, "\nloadStochFile(...): \n");
		fprintf(fpout, "\n %s  %s\n", myfilename, stochdataPtr->probname);
#endif
		if ( strcmp(myfilename, "STOCH") != 0) {
			fprintf(stderr, "\nloadStochFile(...): \n");
			fprintf(stderr, " The first line of the STOCH file %s must start with ""STOCH"" \n!", filename);
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	}
	
	
	// Read the second line: e.g SCENARIOS	DISCRETE
	if (fgets (buffer, LENGTH, stoch) != NULL){
		sscanf(buffer, "%s%s", stochdataPtr->content, stochdataPtr->distn);
#ifdef DEBUG_FUNC
		fprintf(fpout, "\nloadStochFile(...): \n");
		fprintf(fpout, "%s  %s\n", stochdataPtr->content, stochdataPtr->distn);
#endif
		if ( strcmp(stochdataPtr->content, "SCENARIOS") != 0){
			fprintf(stderr, "\nloadStochFile(...): \n");
			fprintf(stderr, " The second line of file %s must start with ""SCENARIOS""!\n", filename);
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
		if ( strcmp(stochdataPtr->distn, "DISCRETE") != 0){
			fprintf(stderr, "\nloadStochFile(...): \n");
			fprintf(stderr, " The second word in the second line of file %s must be ""DISCRETE""!\n", filename);
			fprintf(stderr, "Exiting...\n");
			return(1);
		}
	}
	
	
	// Initialize the counters first
	*random_T = 0;                      // First assume constant technology matrix;
	*random_obj = 0;                    // First assume constant scenario subproblem objective;
	scen_index = -1;   	                // Scenario index
	
	
	for (j = 0; j < numscens; j++) {
		stochdataPtr->cnzcnt_T[j] = 0;
		for (i = 0; i < ncols_master-1; i++)
			stochdataPtr->cmatcnt_T[j][i] = 0;
	}
	
	// Read the rest of the file
	while (fgets (buffer, LENGTH, stoch) != NULL) {
		//fprintf(stderr, "BUFFER: %s\n", buffer);
		sscanf(buffer, "%s", field1);
		
		// printf("Number 10b\n");
		//fprintf(stdout, "%s\n",field1);
		
		if (( strcmp(field1, "ENDATA") == 0)) { // End of file
			scen_index++; // Total number of scenarios
			
			// Probabilites must add to one
			sum = 0;
			for (i = 0; i < scen_index; i++)
				sum +=  stochdataPtr->scenProb[i];
			
			if (sum < 1-NONZERO_LB || sum > 1+NONZERO_LB){
				fprintf (stderr, "loadStochFile: \n");
				fprintf(stderr, "Scenario probabities do NOT add to one!\n");
				fprintf(stderr, "Please check your STOCH file %s\n", filename);
				fprintf(stderr, "Exiting...\n");
				//exit(0);
				return (1);
			}
			
			// Make copy of probabilities for conditional prob
			for (i = 0; i < scen_index; i++) {
				stochdataPtr->scenCondProb[i] =  stochdataPtr->scenProb[i];
			}
			
#ifdef DEBUG_FUNC
			for (i = 0; i < scen_index; i++)
				printf("scenCondProb[%d] = %f \n", i, stochdataPtr->scenCondProb[i]);
#endif
			fclose(stoch);   // close the STOCH file
			return(0);      // successfull return
		}
		else if (strcmp(field1, "SC") == 0){  // Read in new scenario data T(w) amd r(w)
			rowcount = 0; 			// reset the row counter
			scen_index++; 			// count scenario
			stochdataPtr->obj_cnt = 0;	// Re-initialize scenario subproblem random objective
			// coefficients count
			
			
			sscanf(buffer, "%s%s%s%lf%s", field1, field2, field3, &tempDouble, field4);
			// Store the scenario name
			strcpy(stochdataPtr->scenName[scen_index], field2);
			
			// Read in the prob for this scenario
			if (strcmp(field3, "'ROOT'") == 0 || strcmp(field3, "ROOT") == 0){
				stochdataPtr->scenProb[scen_index] = tempDouble;
#ifdef DEBUG_FUNC
				fprintf(fpout, "loadStochFile(...):\n");
				fprintf(fpout, "ECHO: = %s  %s  %s  %f\n", field1, field2, field3,tempDouble);
#endif
			}
			else {
				fprintf(stderr, "\nloadStochFile(...): \n");
				fprintf(stderr, "The third word of line: %s\n", buffer);
				fprintf(stderr, "Must be ""'ROOT'"" and not ""%s""\n", field3);
				fprintf(stderr, "Exiting...\n");
				return (1);
			}
		} else if (strcmp(field1, "RHS") == 0){ // Read the RHS
			
			sscanf(buffer, "%s%s%lf", field1, field2, &tempDouble);
			// Get the row index for this rhs
			for (i = 0; i < subprobPtr->nrows; i++){
				if (strcmp(subprobPtr->rownames[i], field2) == 0){
					break;
				}
			} // End for loop
			
			// Multiply rhs by -1 for constraints originally with a <= sense
			if (subprobPtr->sense[i] == 'L')
				stochdataPtr->rhs[scen_index][i] = -1*tempDouble;
			else
				stochdataPtr->rhs[scen_index][i] = tempDouble;
			
#ifdef DEBUG_FUNC
			printf("buffer: %s\n", buffer);
			fprintf(stderr, "loadStochFile(...):\n");
			fprintf(stderr, "stochdataPtr->rhs[%d][%d] = %f\n", scen_index,
					rowcount, stochdataPtr->rhs[scen_index][rowcount]);
#endif
			rowcount++; // Increment row counter
			
		}
		else if (strcmp(field1, "OBJ") == 0) { // Stochastic objective function g(w) - Comment changed by Dinakar Gade
			
			printf("Random objective function for scenario: %d\n", scen_index);
			*random_obj = 1;
			
			sscanf(buffer, "%s%s%lf", field1, field2, &tempDouble);
			
			//fprintf(stdout, "field1 = %s \t field2 = %s\t  obj = %f \n", field1, field2, tempDouble);
			
			// Get the row index for this rhs
			//fprintf(stdout, "subprobPtr->ncols = %d \n", subprobPtr->ncols);
			for (i = 0; i < subprobPtr->ncols; i++) {
				if (strcmp(subprobPtr->colnames[i], field2) == 0){
					break;
				}
			} // End for loop
			
			//fprintf(stdout, "obj_index[%d] = %d \n", i, stochdataPtr->obj_index[i]);
			
			stochdataPtr->obj[scen_index][stochdataPtr->obj_cnt] = tempDouble;
			
			if (scen_index == 0) { // Store random objective indices only once
				stochdataPtr->obj_index[stochdataPtr->obj_cnt] = i;
			}
			stochdataPtr->obj_cnt++;
			
		} else { // Stochastic Technology matrix T(w)
			fprintf(stdout, "Filed 1 = %s",field1);
			printf("Stochastic Technology matrix T(w)\n");
			*random_T = 1; 		      // Technology matrix is random. More work!!;
			scenario = scen_index;     // Current scenario
			j = stochdataPtr->cnzcnt_T[scenario]; // Current nonzero count
			//printf("scenario = %d\n", scenario);
			//printf("nzcnt_T[%d] = %d\n", scenario, j);
			
			// Read in col name, row name and T(w) nonzero
			sscanf(buffer, "%s%s%lf", field1, field2, &tempDouble);
			//printf("colnames[j] = %s\n", field1);
			//printf("rownames[j] = %s\n", field2);
			
			// Store this nonzero T(w) value
			stochdataPtr->cmatval_T[scenario][j] = tempDouble;
			
			// Set start index for this colname read in
			if (j == 0) { // start of a column
				col = 0;
				stochdataPtr->cmatbeg_T[scenario][col] = 0;
				stochdataPtr->cmatcnt_T[scenario][col]++; // count nonzeros in this column
				//printf("cmatbeg_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatbeg_T[scenario][col]);
				//printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
			} else {
				// count nonzeros in this column
				if (strcmp(subprobPtr->colnames[j-1], subprobPtr->colnames[j]) == 0) {
					stochdataPtr->cmatcnt_T[scenario][col]++;
					//printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
				} else { // new column
					col++;
					stochdataPtr->cmatbeg_T[scenario][col] = j; // set beginning of this new col read
					stochdataPtr->cmatcnt_T[scenario][col]++;   // count nonzeros in this col
					//printf("cmatbeg_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatbeg_T[scenario][col]);
					//printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
				}
				
			}
			// Find row index for current rowname read and store it
			for (i = 0; i < subprobPtr->nrows; i++) {
				if (strcmp(subprobPtr->rownames[i], field2) == 0){
					stochdataPtr->cmatind_T[scenario][j] = i;
					//printf("cmatcnt_T[%d][%d] = %d\n", scenario, col, stochdataPtr->cmatcnt_T[scenario][col]);
					break;
				}
			}
			
			
			stochdataPtr->cnzcnt_T[scenario]++; // Increment nonzeros counter
			
#ifdef DEBUG_FUNC
			printf("cmatval_T[%d][%d] = %f\n", scenario, j, stochdataPtr->cmatval_T[scenario][j]);
			printf("cnzcnt_T[%d] = %d\n", scenario, stochdataPtr->cnzcnt_T[scenario]);
			fprintf(fpout, "loadStochFile(...):\n");
			fprintf(fpout, "ECHO: = %s  %s    %f\n", field1,
					field2, stochdataPtr->cmatval_T[scenario][j]);
#endif
			
		}
		
	} // end while loop
	
	return 1;
	
	fprintf(stderr, "STOCH file %s must end with ""ENDATA""!\n", filename);
	fprintf(stderr, "Exiting...\n");
	// close the TIME file
	fclose(stoch);
	return (1);
	
} //************************** End loadStochFile() function *****************************


int mmult(int *rmatbeg, int *rmatcnt, int *rmatind, double *rmatval, double *x, double *store, int nrows, int ncols, int numnz, int start) {
	/**
	 *  Calculate A*x and put in store
	 * @param rmatbeg,rmatcnt,rmatind,rmatval       - Sparse matrix representation of A matrix
	 * @param x                                     - The vector you want A multiplied with
	 * @param store                                 - The vector where you want to store A*x
	 * @param nrows                                 - The number of rows of A
	 * @param ncols                                 - The number of columns of A
	 * @param numnz                                 - Nonzeros of A
	 * @param start                                 - The index in store from which point onwards you want to store A*x
	 */
	
	int status = 0;
	int i,j,count;
	
	for (i = start; i < nrows + start; i++) {
		store[i] = 0.0;
	}
	
	count = 0;
	/*
	 if (start > 0) {
	 for (i = 0; i < nrows; i++) {
	 for (j = 0; j < rmatcnt[i]; j++) {
	 fprintf(stdout, "rmatind[%d] = %d, rmatval [%d] = %f\n",rmatbeg[i] + j,rmatind[rmatbeg[i] + j],rmatbeg[i] + j, rmatval[rmatbeg[i] + j]);
	 }
	 }
	 }*/
	
	//fprintf(stdout, "START = %d\n",start);
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < rmatcnt[i]; j++) {
			store[start + i] += rmatval[rmatbeg[i] + j]*x[rmatind[rmatbeg[i] + j]];
			//fprintf(stdout, "computing %f*%f for row = %d\n",rmatval[rmatbeg[i] + j],x[rmatind[rmatbeg[i] + j]],start + i);
		}
		//fprintf(stdout, "---\n");
	}
	return status;
}

//************************** End mmult *****************************/

int memAllocFixedRecourseStruct(CUTS_FR *frcuts, int ncols_W, int ncols_T, int nscens) {
	/**
	 * Allocate memory for the fixed recourse data structure
	 * @param
	 */
	int status = 0, scenario;
	int numnz_W = MAX_CUTS_PER_SCENARIO*(ncols_W+1);
	int numnz_T = MAX_CUTS_PER_SCENARIO*ncols_T;
	
	//Initialize the number of cuts to zero
	frcuts->ncuts = 0;
	frcuts->numnz_T = 0;
	frcuts->numnz_W = 0;
	
	//Check Input pointer is not null
	if (frcuts == NULL) {
		fprintf(stderr, "In function memAllocFixedRecourseStruct\n\t. Input pointer to frcuts is NULL\n");
		status = 1;
		return status;
	}
	//Initialize the data to Null first
	
	frcuts->rmatbeg_T = NULL; frcuts->rmatind_T = NULL; frcuts->rmatcnt_T = NULL;
	frcuts->rmatval_T = NULL; frcuts->rmatbeg_W = NULL; frcuts->rmatcnt_W = NULL;
	frcuts->rmatind_W = NULL; frcuts->rmatval_T = NULL; frcuts->cutrhs    = NULL;
	
	//Allocate memory to T sparse matrix for the cuts
	frcuts->rmatbeg_T = (int *)malloc(MAX_CUTS_PER_SCENARIO*sizeof(int));
	frcuts->rmatcnt_T = (int *)malloc(MAX_CUTS_PER_SCENARIO*sizeof(int));
	frcuts->rmatind_T = (int *)malloc(numnz_T*sizeof(int));
	frcuts->rmatval_T = (double *)malloc(numnz_T*sizeof(double));
	
	if (frcuts->rmatbeg_T == NULL || frcuts->rmatcnt_T == NULL || frcuts->rmatind_T == NULL || frcuts->rmatval_T == NULL) {
		fprintf(stderr, "In function memAllocFixedRecourseStruct\n\tError Allocating memory to cut T sparse matrices\n");
		status = 1;
		goto TERMINATE;
	}
	
	//Allocate memory to W sparse matrix for the cuts
	frcuts->rmatbeg_W = (int *)malloc(MAX_CUTS_PER_SCENARIO*sizeof(int));
	frcuts->rmatcnt_W = (int *)malloc(MAX_CUTS_PER_SCENARIO*sizeof(int));
	frcuts->rmatind_W = (int *)malloc(numnz_W*sizeof(int));
	frcuts->rmatval_W = (double *)malloc(numnz_W*sizeof(double));
	
	if (frcuts->rmatbeg_W == NULL || frcuts->rmatcnt_W == NULL || frcuts->rmatind_W == NULL || frcuts->rmatval_W == NULL) {
		fprintf(stderr, "In function memAllocFixedRecourseStruct\n\tError Allocating memory to cut W sparse matrices\n");
		status = 1;
		goto TERMINATE;
	}
	
	//Cut RHS
	frcuts->cutrhs = (double **)malloc(nscens*sizeof(double *));
	if (frcuts->cutrhs == NULL) {
		fprintf(stderr, "In function memAllocFixedRecourseStruct\n\tError allocating memory to CUT RHS\n");
		status = 1;
		goto TERMINATE;
	}
	for (scenario = 0; scenario < nscens; scenario++) {
		frcuts->cutrhs[scenario] = (double *)malloc(MAX_CUTS_PER_SCENARIO*sizeof(double));
		if (frcuts->cutrhs[scenario] == NULL) {
			fprintf(stderr, "In function memAllocFixedRecourseStruct\n\t. Error allocating memory to cut rhs for scenario %d\n", scenario);
			status = 1;
			goto TERMINATE;
		}
	}
	
TERMINATE:
	return status;
}

//************************** End memAllocFixedRecourseStruct *****************************/

void freeFixedRecourseStructs(CUTS_FR *frcuts, int nscens) {
	int s;
	if (frcuts != NULL) {
		if (frcuts->rmatbeg_T != NULL) {
			free(frcuts->rmatbeg_T);
			frcuts->rmatbeg_T = NULL;
		}
		if (frcuts->rmatcnt_T != NULL) {
			free(frcuts->rmatcnt_T);
			frcuts->rmatcnt_T = NULL;
		}
		if (frcuts->rmatind_T != NULL) {
			free(frcuts->rmatind_T);
			frcuts->rmatind_T = NULL;
		}
		if (frcuts->rmatval_T != NULL) {
			free(frcuts->rmatval_T);
			frcuts->rmatval_T = NULL;
		}
		if (frcuts->rmatbeg_W != NULL) {
			free(frcuts->rmatbeg_W);
			frcuts->rmatbeg_W = NULL;
		}
		if (frcuts->rmatcnt_W != NULL) {
			free(frcuts->rmatcnt_W);
			frcuts->rmatcnt_W = NULL;
		}
		if (frcuts->rmatind_W != NULL) {
			free(frcuts->rmatind_W);
			frcuts->rmatind_W = NULL;
		}
		for (s = 0; s < nscens; s++) {
			if (frcuts->cutrhs[s] != NULL) {
				free(frcuts->cutrhs[s]);
				frcuts->cutrhs[s] = NULL;
			}
		}
		free(frcuts->cutrhs);
		frcuts->cutrhs = NULL;
		free(frcuts);
		frcuts = NULL;
	}
	
}

//************************** End Free *****************************/


int compute_rhs_rows(subproblem_t *subprobPtr, stochfile_info *stochdataPtr, CUTS_FR *frcuts, double *rhs, int rhssize, int omega, double *solnX, int x_size) {
	
	/**
	 * Compute the right hand side r - Tx for scenario omega using row format
	 * @param subprobPtr        - Subproblem pointer containing T information
	 * @param stochdataPtr      - The stochastic data pointer
	 * @param frcuts            - The pointer to the fixed recourse cuts data structure
	 * @param rhs               - The right hand size which must be populated with r - Tx  (memory  must already alloced to rhs)
	 * @param rhssize           - The size of rhs (= #rows of T matrix)
	 * @param omega             - Scenario omega
	 * @param solnX             - The first stage solution
	 */
	int status = 0,i;
	double *Tx = NULL;
	
	//Extract r^k(\omega)
	for (i = 0; i < subprobPtr->nrows_T; i++) {
		rhs[i] = stochdataPtr->rhs[omega][i];
		//fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n",omega,i,stochdataPtr->rhs[omega][i]);
	}
	for (i = subprobPtr->nrows_T; i < frcuts->ncuts + subprobPtr->nrows_T; i++) {
		rhs[i] = frcuts->cutrhs[omega][i-subprobPtr->nrows_T];
		//fprintf(stdout, "stochdataPtr->rhs[%d][%d] = %f\n",omega,i,stochdataPtr->rhs[omega][i]);
	}
	
	Tx = (double *)calloc(subprobPtr->nrows_T + frcuts->ncuts, sizeof(double));
	if (Tx == NULL) {
		fprintf(stdout, "In compute_rhs_rows():\n\tError allocating memory to Tx\n");
		return (1);
	}
	status = mmult(subprobPtr->rmatbeg, subprobPtr->rmatcnt, subprobPtr->rmatind, subprobPtr->rmatval, solnX, Tx, subprobPtr->nrows_T, subprobPtr->ncols_T, subprobPtr->nzcnt_T,0);
	if (status) {
		fprintf(stderr, "In compute_rhs_rows():\n\tError with mmult():\n");
		status = 1;
		goto TERMINATE;
	}
	status = mmult(frcuts->rmatbeg_T, frcuts->rmatcnt_T, frcuts->rmatind_T, frcuts->rmatval_T, solnX, Tx, frcuts->ncuts, subprobPtr->ncols_T, frcuts->numnz_T, subprobPtr->nrows_T);
	if (status) {
		fprintf(stderr, "In compute_rhs_rows():\n\tError with mmult(): in the second time\n");
		if (Tx != NULL) {
			free((Tx));
			Tx = NULL;
		}
		return (1);
	}
	for (i = 0; i < rhssize; i++) {
		rhs[i] = (rhs[i] - Tx[i]);
	}
	
	
TERMINATE:
	if (Tx != NULL) {
		free(Tx);
		Tx = NULL;
	}
	return status;
}

//************************** End Compute_RHS_ROWS *****************************/


int check_sol_frac(double *x, int n) {
	/*
	 * Check whether the vector x has fractional elements
	 *@param x - Vector
	 *@param n - Size of the vector x
	 *@return  - The number of fractional variables
	 */
	int i,numfrac = 0;
	
	for(i = 0; i < n; i++) {
		if(is_element_int(x[i]) == FALSE) {
			numfrac++;
		}
	}
	return numfrac;
}

//************************** End check_sol_frac *****************************/


int is_element_int(double x) {
	/**
	 * Check whether x is integer or fractional
	 * @param x
	 * @return - 1 if x is integer, 0 if x is fractional
	 **/
	int isint = TRUE;
	if((x - floor(x) > INT_PRECISION) && (ceil(x) - x > INT_PRECISION)) {
		isint = FALSE;
	}
	
	return (isint);
	
}
//************************** End is_element_int *****************************/


int check_sorted_bhead(int *bhead, int basis_size, int ncols) {
	int sorted = TRUE,i;
	int compare1, compare2;
	for (i = 0; i < basis_size - 1; i++) {
		if (bhead[i] < 0) {
			compare1 = -(bhead[i] + 1) + ncols;
		}
		if (bhead[i] >= 0) {
			compare1 = bhead[i];
		}
		if (bhead[i+1] < 0) {
			compare2 = -(bhead[i+1] + 1) + ncols;
		}
		if (bhead[i+1] >= 0) {
			compare2 = bhead[i+1];
		}
		if (compare1 > compare2) {
			return FALSE;
		}
	}
	
	
	return sorted;
}

//************************** End check_sorted_bhead *****************************/

int sort_bhead2(int *bhead, double *barx, int *bindex, int basis_size, int ncols) {
	/**
	 *  The basis inverse returned by CPLEX could countain rows based on bhead in non-sorted order, for e.g.
	 *  the order of basic columns in bhead could be [3,-4,0,-2,1] implying 3rd column, slack corresponding to 4th row,
	 *  etc are in the basis. BINVERSE returned by CPLEX also contains rows in this order. Since, we are computing binverse T,
	 *  the rows must be sorted to get the correct values.
	 *  @param bhead        - Basis header
	 *  @param barx         - Values of the basic vars
	 *  @param binv         - Basis inverse
	 *  @param basis_size   - No. columns in the basis
	 *  @param ncols        - Number of columns
	 *
	 **/
	int status = 0,i,key1,key2,j;
	double tempx;
	int tempidx;
	int *newbhead = NULL;
	
	newbhead = (int *)malloc(basis_size*sizeof(int));
	if (newbhead == NULL) {
		fprintf(stderr, "Error in sort_binverse():\n\t Could not allocate newbhead\n");
		return 1;
	}
	//Increase the indices of slacks by ncols and making them positive
	for (i = 0; i < basis_size; i++) {
		if (bhead[i] < 0) {
			newbhead[i] = -(bhead[i] + 1) + ncols;
		}
		else {
			newbhead[i] = bhead[i];
		}
	}
	//Now sort bhead, barx, binv based on newbhead in ascending order - Uses insertion sort - could be improved!
	for (j = 1; j < basis_size; j++) {
		key1 = newbhead[j];
		key2 = bhead[j];
		tempx = barx[j];
		tempidx = bindex[j];
		i = j;
		while (i > 0 && newbhead[i-1] > key1) {
			newbhead[i] = newbhead[i-1];
			bhead[i]    = bhead[i-1];
			barx[i]     = barx[i-1];
			bindex[i]   = bindex[i-1];
			i--;
		}
		newbhead[i]     = key1;
		bhead[i]        = key2;
		barx[i]         = tempx;
		bindex[i]       = tempidx;
	}
	
	if (newbhead != NULL) {
		free(newbhead);
		newbhead = NULL;
	}
	return status;
}
//************************** End sort_binverse() *****************************/

int get_nonbasic_vars(int *cstat, int *rstat, int ncols, int nrows, int *nonbasic_vars) {
	/**
	 * Dinakar Gade
	 * Get the list of nonbasic variables
	 * @param cstat         - Basis status of the columns
	 * @param rstat         - Basis status of the rows
	 * @param ncols         - Number of columns in the problem
	 * @param nrows         - Number of rows in the problem
	 * @param nonbasiv_vars - The array where the list of nonbasic vars are stored (assumed alloced earlier)
	 **/
	
	int i,status = 0;
	int count = 0;
	
	if (cstat == NULL || rstat == NULL || nonbasic_vars == NULL) {
		fprintf(stderr, "In get_non_basic_vars():\n\t NullPointerError for rstat/cstst/nonbasic_vars\n");
		return (1);
	}
	
	for(i = 0; i < ncols; i++) {
		if(cstat[i] == CPX_AT_LOWER || cstat[i] == CPX_AT_UPPER) {
			nonbasic_vars[count] = i;
			count++;
		}
	}
	for(i = 0; i < nrows; i++) {
		if(rstat[i] == CPX_AT_LOWER || rstat[i] == CPX_AT_UPPER) {
			nonbasic_vars[count] = -(i) - 1;
			count++;
		}
	}
	return status;
}
//************************** End get_nonbasic_vars() *****************************/

int compute_binvTrow(double *binvrow, subproblem_t *subprobPtr, CUTS_FR *frcuts, int basis_size, double *binvTrow) {
	
	int status = 0;
	int i,j;
	
	//Reset binvTrow
	for (i = 0; i < subprobPtr->ncols_T; i++) {
		binvTrow[i] = 0.0;
	}
	for (i = 0; i < basis_size-1; i++) {
		if (i < subprobPtr->nrows_T) {
			for (j = 0; j < subprobPtr->rmatcnt[i]; j++) {
				binvTrow[subprobPtr->rmatind[subprobPtr->rmatbeg[i] + j]] += binvrow[i+1]*subprobPtr->rmatval[subprobPtr->rmatbeg[i] + j];
			}
		}
		if (i >= subprobPtr->nrows_T) {
			for (j = 0; j < frcuts->rmatcnt_T[i-subprobPtr->nrows_T]; j++) {
				binvTrow[frcuts->rmatind_T[frcuts->rmatbeg_T[i-subprobPtr->nrows_T] + j]] += binvrow[i+1]*frcuts->rmatval_T[frcuts->rmatbeg_T[i-subprobPtr->nrows_T] + j];
			}
		}
	}
	return status;
}

//************************** End compute_binvTrow() *****************************/


int generate_fixed_recourse_cuts(CPXENVptr env, CPXLPptr lp, subproblem_t *subprobPtr, stochfile_info *stochdataPtr, CUTS_FR *frcuts, int ncuts, double *solnX, int scen) {
	/**
	 * env              - CPLEX environment
	 * lp               - LP Pointer to subproblem
	 * subprobPtr       - Pointer to the subproblem
	 * stochdataPtr     - Stochastic data pointer
	 * frcuts           - Cuts data structure
	 * ncuts            - Number of cuts to be generated
	 * solnX            - First stage solution being passed
	 */
	int status = 0, ncols, nrows, numnz, basis_size, nzcnt, surplus,i, cutcount, headcount,j,scenario;
	//Local Arrays for storage
	int *rmatbeg_W = NULL, *rmatind_W = NULL, *rmatcnt_W = NULL;
	double *rmatval_W = NULL;
	char *sense_W = NULL;
	int *cstat = NULL, *rstat = NULL;
	
	int *bhead = NULL; double *barx = NULL;//Basis info - header and values
	
	//Calculated arrays for basis
	int *nonbasics      = NULL;
	double *binvrow     = NULL;
	double *binvn       = NULL;
	double *binva       = NULL;
	double *binvTrow    = NULL;
	
	
	double *cutcoeffs   = NULL;
	int *bindex         = NULL;
	double *T_update    = NULL;
	char *sense         = NULL;
	double *newrhs      = NULL;
	
	int old_row_idx;
	int old_nz_idx;
	
	int *temp_rmatbeg_W = NULL;
	
	double rhs, rhs_dueto_nonbasics, rhs_dueto_binvt;
	double rhsdueto_rounding;
	
	//Check input arguments
	if (stochdataPtr == NULL || subprobPtr == NULL || stochdataPtr == NULL || solnX == NULL || frcuts == NULL) {
		fprintf(stderr, "In generate_fixed_recourse_cuts():\n\tNullPointerError for input arguments\n");
		status = 1;
		goto TERMINATE;
	}
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);
	basis_size = nrows;
	
	rmatbeg_W = (int *)malloc(nrows*sizeof(int));
	rmatcnt_W = (int *)malloc(nrows*sizeof(int));
	rmatind_W = (int *)malloc(numnz*sizeof(int));
	rmatval_W = (double *)malloc(numnz*sizeof(double));
	bhead     = (int *)malloc(nrows*sizeof(int));
	barx      = (double *)malloc(nrows*sizeof(double));
	nonbasics = (int *)malloc(ncols*sizeof(int));
	sense_W   = (char *)malloc(nrows*sizeof(char));
	cstat     = (int *)malloc(ncols*sizeof(int));
	rstat     = (int *)malloc(nrows*sizeof(int));
	
	
	if (rmatbeg_W == NULL || rmatind_W == NULL || rmatval_W == NULL || bhead == NULL || barx == NULL || sense_W == NULL) {
		status = 1;
		fprintf(stderr, "In function generate_fixed_recourse_cuts():\n\tError allocating memory to row for W matrix. STATUS = %d",status);
		goto TERMINATE;
	}
	
	status = CPXgetrows(env, lp, &nzcnt, rmatbeg_W, rmatind_W, rmatval_W, numnz, &surplus, 0, nrows-1);
	if (status) {
		fprintf(stderr, "In function generate_fixed_recourse_cuts():\n\tError getting rows from CPLEX. STATUS = %d",status);
		goto TERMINATE;
	}
	
	//Convert rmatbeg to rmatcnt
	for (i = 0; i < nrows-1; i++) {
		rmatcnt_W[i] = rmatbeg_W[i+1] - rmatbeg_W[i];
	}
	rmatcnt_W[nrows-1] = numnz - rmatbeg_W[nrows-1];
	
	
	status = CPXgetsense(env, lp, sense_W, 0, nrows-1);
	if (status) {
		fprintf(stderr, "In function generate_fixed_recourse_cuts():\n\tError getting sense from CPLEX. STATUS = %d",status);
		goto TERMINATE;
		
	}
	status = CPXgetbhead(env, lp, bhead, barx);
	if (status) {
		fprintf(stderr, "In generate_fixed_recourse_cuts():\n\tError getting basis header. CPLEX STATUS = %d\n",status);
		goto TERMINATE;
	}
	status = CPXgetbase(env, lp, cstat, rstat);
	if (status) {
		fprintf(stderr, "In generate_fixed_recourse_cuts():\n\tError getting basis statuses. CPLEX STATUS = %d\n",status);
		goto TERMINATE;
	}
	
	bindex = (int *)malloc(basis_size*sizeof(int));
	for (i = 0; i < basis_size; i++) {
		bindex[i] = i;
	}
	if (check_sorted_bhead(bhead, basis_size, ncols) == FALSE) {
		//fprintf(stdout, "BASIS HEADER IS NOT SORTED. Sorting...\n");
		sort_bhead2(bhead, barx, bindex, basis_size, ncols);
	}
	
	//Extract the nonasic variables
	status = get_nonbasic_vars(cstat, rstat, ncols, nrows, nonbasics);
	if (status) {
		fprintf(stderr, "In generate_fixed_recourse_cuts():\n\tError getting nonbasic variables\n");
		goto TERMINATE;
	}
	binvrow  = (double *)malloc(nrows*sizeof(double));
	binva    = (double *)malloc(ncols*sizeof(double));
	binvn    = (double *)malloc(ncols*sizeof(double));
	T_update = (double *)malloc((subprobPtr->ncols_T)*sizeof(double));
	cutcoeffs = (double *)calloc(ncols,sizeof(double));
	binvTrow  = (double *)malloc(subprobPtr->ncols_T*sizeof(double));
	
	sense = (char *)malloc(ncuts*sizeof(char));
	newrhs = (double *)malloc(ncuts*sizeof(double));
	
	if (T_update == NULL || cutcoeffs == NULL || binvrow == NULL || binvn == NULL ||sense == NULL || newrhs == NULL ||binva == NULL
		|| binvTrow == NULL) {
		fprintf(stderr, "In generate_fixed_recourse_cuts():\n\tError allocating memory to cutcoeffs/T_update/\n");
		return (1);
	}
	for (i = 0; i < ncuts; i++) {
		sense[i] = 'G';
	}
	old_row_idx = frcuts->ncuts;
	old_nz_idx  = frcuts->numnz_W;
	
	/*****************************************************************
	 START CUT GENERATION
	 *****************************************************************/
	cutcount = 0;
	for (headcount = 0; headcount < basis_size; headcount++) {
		if (cutcount >= ncuts) {
			break;
		}
		if (bhead[headcount] >= 0 && is_element_int(barx[headcount]) == FALSE) {
			//fprintf(stdout, "BHEAD [%d] = %d, barx[%d] = %f",headcount,bhead[headcount],headcount,barx[headcount]);
			//fprintf(stdout, "=======\nCUT%d\n=======\n",cutcount);
			
			//Get BINVERSE for the source row
			status = CPXbinvrow(env, lp, bindex[headcount], binvrow);
			if (status) {
				fprintf(stderr, "In generate_fixed_recourse_cuts():\n\t Error getting binvrow,CPLEX STATUS = %d\n",status);
				goto TERMINATE;
			}
			//for(i = 0; i < basis_size; i++) {
			//    fprintf(stdout, "BINV[%d] = %f\n", i,binvrow[i]);
			//}
			//Get BINV*A for the source row
			status = CPXbinvarow(env, lp, bindex[headcount], binva);
			if (status) {
				fprintf(stderr, "In generate_fixed_recourse_cuts():\n\t Error getting binvarow. CPLEX STATUS = %d\n",status);
				goto TERMINATE;
			}
			//Calculate binvn
			for (i = 0; i < ncols; i++) {
				if (nonbasics[i] >= 0) {
					binvn[i] = binva[nonbasics[i]];
				}
				if(nonbasics[i] < 0) {
					if(sense_W[(-1)*(nonbasics[i]+1)] == 'L') {
						binvn[i] = binvrow[(-1)*(nonbasics[i]+1)];
					}
					if(sense_W[(-1)*(nonbasics[i]+1)] == 'G') {
						binvn[i] = -1.0*binvrow[(-1)*(nonbasics[i]+1)];
					}
					if(sense_W[(-1)*(nonbasics[i]+1)] == 'E') {
						binvn[i] = 0.0;//Because slacks corresponding to equality constraints are zero, they dont count to cutcoeffs
					}
				}
			}
			
			//for (i = 0; i < ncols; i++) {
			//    fprintf(stdout, "Nonbasics[%d] = %d, BINV*N[%d] = %0.24f\n",i,nonbasics[i],i,binvn[i]);
			//}
			//Calcuate binvTrow
			status = compute_binvTrow(binvrow, subprobPtr, frcuts, basis_size, binvTrow);
			if (status) {
				fprintf(stderr, "In generate_cuts2():\n\tError computing binvTrow\n");
				goto TERMINATE;
			}
			newrhs[cutcount] = myceil(barx[headcount]);
			
			rhs_dueto_nonbasics = 0.0;
			rhsdueto_rounding   = 0.0;
			for (i = 0; i < ncols; i++) {
				if (nonbasics[i] >= 0 && cstat[nonbasics[i]] == CPX_AT_UPPER) {
					//newrhs[cutcount] += binvn[i]*subprobPtr->ub[nonbasics[i]-1];
					rhs_dueto_nonbasics += (binvn[i])*subprobPtr->ub[nonbasics[i]-1];
					rhsdueto_rounding += myceil(-1.0*binvn[i])*subprobPtr->ub[nonbasics[i]-1];
				}
				if (nonbasics[i] >= 0 && cstat[nonbasics[i]] == CPX_AT_LOWER) {
					//newrhs[cutcount] += binvn[i]*subprobPtr->lb[nonbasics[i]-1];
					rhs_dueto_nonbasics += (binvn[i])*subprobPtr->lb[nonbasics[i]-1];
					rhsdueto_rounding += -1.0*myceil(binvn[i])*subprobPtr->lb[nonbasics[i]-1];
				}
			}
			for (i = 0; i < ncols; i++) {
				if (nonbasics[i] >= 0 && cstat[nonbasics[i]] == CPX_AT_UPPER) {
					
					binvn[i] = myfloor(binvn[i]);
				}
				else {
					binvn[i] = myceil(binvn[i]);
				}
			}
			
			
			//for (i = 0; i < ncols; i++) {
			//    fprintf(stdout, "Nonbasics[%d] = %d, BINV*N[%d] = %0.24f\n",i,nonbasics[i],i,binvn[i]);
			//}
			
			
			//Reset T_update
			for (i = 0; i < subprobPtr->ncols_T; i++) {
				T_update[i] = 0.0;
			}
			rhs_dueto_binvt = 0.0;
			//Compute the T_update
			for (i = 0; i < subprobPtr->ncols_T; i++) {
				if (myfabs(solnX[i] - 1.0) < EPSILON_TOLERANCE) {
					T_update[i] = myfloor(binvTrow[i]);
					rhs_dueto_binvt+= binvTrow[i];
					rhsdueto_rounding += myceil(-binvTrow[i]);
					newrhs[cutcount]+= T_update[i];
					
				}
				if (myfabs(solnX[i]) < EPSILON_TOLERANCE) {
					T_update[i] = myceil(binvTrow[i]);
				}
				//fprintf(stdout, "T_update[%d] = %f\n",i,T_update[i]);
			}
			
			//Now compute the projection by eliminating the slacks in the constraints and deriving the cut in the space of original (x,y) variables
			//First reset cutcoeffs
			for (i = 0; i < ncols; i++) {
				cutcoeffs[i] = 0.0;
			}
			cutcoeffs[bhead[headcount]] = 1.0;
			for (i = 0; i < ncols; i++) {
				if (nonbasics[i] >= 0) {
					cutcoeffs[nonbasics[i]] = binvn[i];
				}
				//Substitute out the slack variables
				if (nonbasics[i] < 0) {
					if (sense_W[-(nonbasics[i] + 1)] == 'G') {
						//Update the W-matrix
						for(j = 0; j < rmatcnt_W[-1*(nonbasics[i]+1)]; j++) {
							cutcoeffs[rmatind_W[rmatbeg_W[-1*(nonbasics[i]+1)] + j]] +=
							binvn[i]*rmatval_W[rmatbeg_W[-1*(nonbasics[i]+1)] + j];
						}
						//Update the RHS
						//rhs += binvn[i]*stochdatPtr->rhs[scenario][-1*(nonbasics[i]+1)-1];
						//Update the actual rhs for the given SolnX and Scenario - start with the current rhs and compute r - Tx
						newrhs[cutcount] += binvn[i]*stochdataPtr->rhs[scen][-1*(nonbasics[i]+1)-1];
						//Update T-coefficients
						if (-(nonbasics[i] + 1) <= subprobPtr->nrows) {
							//Slack for the original T matrix
							for (j = 0; j < subprobPtr->rmatcnt[-(nonbasics[i] + 1)-1]; j++) {
								T_update[subprobPtr->rmatind[subprobPtr->rmatbeg[-(nonbasics[i] + 1)-1]+j]]
								+= subprobPtr->rmatval[subprobPtr->rmatbeg[-(nonbasics[i] + 1)-1]+j]*binvn[i];
							}
						}
						if (-(nonbasics[i] + 1) > subprobPtr->nrows) {
							//Slack for the new T matrix
							for (j = 0; j < frcuts->rmatcnt_T[-(nonbasics[i] + 1) - subprobPtr->nrows -1]; j++) {
								T_update[frcuts->rmatind_T[frcuts->rmatbeg_T[-(nonbasics[i] + 1)-1-subprobPtr->nrows]+j]]
								+= frcuts->rmatval_T[frcuts->rmatbeg_T[-(nonbasics[i] + 1)-1-subprobPtr->nrows]+j]*binvn[i];
							}
						}
					}
					if (sense_W[-(nonbasics[i] + 1)] == 'L') {
						//Update the W-matrix
						for(j = 0; j < rmatcnt_W[-1*(nonbasics[i]+1)]; j++) {
							cutcoeffs[rmatind_W[rmatbeg_W[-1*(nonbasics[i]+1)] + j]] -=
							binvn[i]*rmatval_W[rmatbeg_W[-1*(nonbasics[i]+1)] + j];
						}
						//Update the RHS
						newrhs[cutcount] -= binvn[i]*stochdataPtr->rhs[scen][-1*(nonbasics[i]+1)-1];
						//Update T-coefficients
						if (-(nonbasics[i] + 1) <= subprobPtr->nrows) {
							//Slack for the original T matrix
							for (j = 0; j < subprobPtr->rmatcnt[-(nonbasics[i] + 1)-1]; j++) {
								T_update[subprobPtr->rmatind[subprobPtr->rmatbeg[-(nonbasics[i] + 1)-1]+j]]
								-= subprobPtr->rmatval[subprobPtr->rmatbeg[-(nonbasics[i] + 1)-1]+j]*binvn[i];
							}
						}
						if (-(nonbasics[i] + 1) > subprobPtr->nrows) {
							//Slack for the new T matrix
							for (j = 0; j < frcuts->rmatcnt_T[-(nonbasics[i] + 1) - subprobPtr->nrows -1]; j++) {
								T_update[frcuts->rmatind_T[frcuts->rmatbeg_T[-(nonbasics[i] + 1)-1-subprobPtr->nrows]+j]]
								-= frcuts->rmatval_T[frcuts->rmatbeg_T[-(nonbasics[i] + 1)-1-subprobPtr->nrows]+j]*binvn[i];
							}
						}
					}
				}
			}
			
			//If y_0 has nonzero coeff - project out the objective function variable:
			if (myfabs(cutcoeffs[0]) > INT_PRECISION) {
				cutcoeffs[0] = 0.0;
				for (i = 0; i < ncols-1; i++) {
					cutcoeffs[i+1] += subprobPtr->obj[i];
				}
			}
			
			//for (i = 0; i < subprobPtr->ncols_T; i++) {
			//    fprintf(stdout, "CUTCOEFFS_T[%d] = %f\n",i,T_update[i]);
			//}
			/*
			 for (i = 0; i < ncols; i++) {
			 fprintf(stdout, "CUTCOEFFS_W[%d] = %f\n",i,cutcoeffs[i]);
			 }
			 */
			
			//Update newrhs using the T_update
			for (i = 0; i < subprobPtr->ncols_T; i++) {
				newrhs[cutcount] -= T_update[i]*solnX[i];
			}
			
			//Update the W^k(\omega)\W matrix
			//fprintf(stdout, "ORIG_ROWS = %d\n",old_row_idx);
			//fprintf(stdout, "ORIG_NZ = %d\n",old_nz_idx);
			frcuts->rmatbeg_W[frcuts->ncuts] = frcuts->numnz_W;
			frcuts->rmatcnt_W[frcuts->ncuts] = 0;
			for (i = 0; i < ncols; i++) {
				if (myfabs(cutcoeffs[i]) > INT_PRECISION) {
					frcuts->rmatcnt_W[frcuts->ncuts]++;
					frcuts->rmatind_W[frcuts->numnz_W] = i;
					frcuts->rmatval_W[frcuts->numnz_W] = cutcoeffs[i];
					frcuts->numnz_W++;
				}
			}
			
			//Update the T^k(\omega)\W matrix
			frcuts->rmatbeg_T[frcuts->ncuts] = frcuts->numnz_T;
			frcuts->rmatcnt_T[frcuts->ncuts] = 0;
			for (i = 0; i < subprobPtr->ncols_T; i++) {
				if (myfabs(T_update[i]) > INT_PRECISION) {
					frcuts->rmatcnt_T[frcuts->ncuts ]++;
					frcuts->rmatind_T[frcuts->numnz_T] = i;
					frcuts->rmatval_T[frcuts->numnz_T] = T_update[i];
					frcuts->numnz_T++;
				}
			}
   
   
			//fprintf(stdout, "RHS_DUE_TO_NONBASICS = %f\n", rhs_dueto_nonbasics);
			for(scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
				//Compute B^-1r(\omega)
				rhs = 0.0;
				for(j = 0; j < subprobPtr->nrows + frcuts->ncuts; j++) {
					if(j < subprobPtr->nrows) {
						//fprintf(stdout, "HERE!!\n");
						rhs += binvrow[j+1]*stochdataPtr->rhs[scenario][j];
					}
					if(j>=subprobPtr->nrows) {
						rhs += binvrow[j+1]*frcuts->cutrhs[scenario][j-subprobPtr->nrows];
					}
				}
				//fprintf(stdout, "SCENARIO = %d, RHS INIT = %f\n",scenario, rhs);
				//Due to nonbasics
				rhs -= rhs_dueto_nonbasics;
				//Due to BINVT (complementation)
				rhs -= rhs_dueto_binvt;
				//fprintf(stdout, "SCENARIO = %d, RHS BEFORE ROUNDING = %f\n",scenario, rhs);
				rhs = myceil(rhs);
				//fprintf(stdout, "SCENARIO = %d, RHS AFTER ROUNDING = %f\n",scenario, rhs);
				//fprintf(stdout, "SCENARIO = %d, RHS DUE TO ROUNDING = %f\n", scenario, rhsdueto_rounding);
				//Now the part where we project out slacks
				for (i = 0; i < ncols; i++) {
					if (nonbasics[i] < 0 && -1*(nonbasics[i]+1)-1 < subprobPtr->nrows) {
						//fprintf(stdout, "HEEEERERERER\n");
						if (sense_W[-(nonbasics[i] + 1)] == 'G') {
							rhs += binvn[i]*stochdataPtr->rhs[scenario][-1*(nonbasics[i]+1)-1];
						}
						if (sense_W[-(nonbasics[i] + 1)] == 'L') {
							rhs -= binvn[i]*stochdataPtr->rhs[scenario][-1*(nonbasics[i]+1)-1];
						}
					}
					if (nonbasics[i] < 0 &&  -1*(nonbasics[i]+1)-1 >= subprobPtr->nrows) {
						//fprintf(stdout, "HEEEERERERER2\n");
						if (sense_W[-(nonbasics[i] + 1)] == 'G') {
							rhs += binvn[i]*frcuts->cutrhs[scenario][-1*(nonbasics[i]+1)-1 - subprobPtr->nrows];
						}
						if (sense_W[-(nonbasics[i] + 1)] == 'L') {
							rhs -= binvn[i]*frcuts->cutrhs[scenario][-1*(nonbasics[i]+1)-1 - subprobPtr->nrows];
						}
					}
				}
				rhs -= rhsdueto_rounding;
				frcuts->cutrhs[scenario][frcuts->ncuts] = rhs;
				//fprintf(stdout, "SCENARIO = %d, FINAL RHS AFTER PROJECTION = %f\n",scenario, rhs);
			}
			frcuts->ncuts++;
			cutcount++;
		}
	}
	
	//for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
	//    for (j = 0; j < frcuts->ncuts; j++) {
	//        fprintf(stdout, "SCENRHS[%d][%d] = %f\n", scenario, j, frcuts->cutrhs[scenario][j]);
	//
	//    }
	//}
 
	temp_rmatbeg_W = (int *)calloc(ncuts, sizeof(int));
	if (temp_rmatbeg_W == NULL) {
		fprintf(stderr, "In generate_cuts():\n\tError Allocating memory to temp_rmatbeg_w\n");
		goto TERMINATE;
	}
	for (i = 0; i < ncuts; i++) {
		temp_rmatbeg_W[i] = *(frcuts->rmatbeg_W+old_row_idx + i) - old_nz_idx;
	}
	status = CPXaddrows(env, lp, 0, ncuts, (frcuts->numnz_W - old_nz_idx), newrhs, sense, temp_rmatbeg_W,
						(frcuts->rmatind_W + old_nz_idx), (frcuts->rmatval_W + old_nz_idx), NULL, NULL);
	CPXwriteprob(env, lp, "gcuts1.lp", "lp");
	if (status) {
		fprintf(stderr, "In generate_cuts2():\n\t. Error adding cuts for scenario %d. CPLEX status = %d\n",scenario,status);
		goto TERMINATE;
	}
	
	
TERMINATE:
	
	if(rmatbeg_W != NULL) {
		free(rmatbeg_W);
		rmatbeg_W = NULL;
	}
	if(rmatcnt_W != NULL) {
		free(rmatcnt_W);
		rmatcnt_W = NULL;
	}
	if(rmatind_W != NULL) {
		free(rmatind_W);
		rmatind_W = NULL;
	}
	if(rmatval_W != NULL) {
		free(rmatval_W);
		rmatval_W = NULL;
	}
	if(bhead != NULL) {
		free(bhead);
		bhead = NULL;
	}
	if(sense_W != NULL) {
		free(sense_W);
		sense_W = NULL;
	}
	if(barx != NULL) {
		free(barx);
		barx = NULL;
	}
	if (cstat != NULL) {
		free(cstat);
		cstat = NULL;
	}
	if(rstat != NULL) {
		free(rstat);
		rstat = NULL;
	}
	if(bindex != NULL) {
		free(bindex);
		bindex = NULL;
	}
	if(binva != NULL) {
		free(binva);
		binva = NULL;
	}
	if(binvn != NULL) {
		free(binvn);
		binvn = NULL;
	}
	if(binvrow != NULL) {
		free(binvrow);
		binvrow = NULL;
	}
	if(binvTrow != NULL) {
		free(binvTrow);
		binvTrow = NULL;
	}
	if(cutcoeffs != NULL) {
		free(cutcoeffs);
		cutcoeffs = NULL;
	}
	if(sense != NULL) {
		free(sense);
		sense = NULL;
	}
	if(newrhs != NULL) {
		free(newrhs);
		newrhs = NULL;
	}
	return status;
}



void print_CPX_sparse_matrix(int *rmatbeg, int *rmatcnt, int *rmatind, double *rmatval, int nrows, int ncols, int numnz) {
	/**
	 * Prints the CPLEX sparse matrix representation as a 2D matrix to screen
	 * @param rmatbeg - CPLEX style row beginning indices
	 * @param rmatind - CPLEX style column indices
	 * @param rmatval - CPLEX style column values
	 * @param nrows   - Number of rows in the matrix (size of rmatbeg)
	 * @param ncols   - Number of columns in the matrix
	 * @param numnz   - Size of the non-zeros (size of rmatind and rmatval)
	 */
	int i,j;
	double **A = NULL;
	
	for (i = 0; i < nrows; i++) {
		//fprintf(stdout, "rmatcount[%d] = %d\n",(i),rmatcnt[i]);
	}
	
	A = (double **) calloc(nrows, sizeof(double));
	for (i = 0; i < nrows; i++) {
		A[i] = (double *)calloc(ncols, sizeof(double));
	}
	
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < rmatcnt[i]; j++) {
			A[i][rmatind[rmatbeg[i] + j]] = rmatval[rmatbeg[i] + j];
		}
	}
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			printf("%f ",A[i][j]);
		}
		printf("\n");
	}
	for (i = 0; i < nrows; i++) {
		if(A[i] != NULL) {
			free(A[i]);
			A[i] = NULL;
		}
	}
	free(A);
}

ParametricGomoryAlg::ParametricGomoryAlg()
{
	
}

bool ParametricGomoryAlg::solve(int argc, const char **argv)
{
	int status;
	int DISPLAY_LOG_INFO = 1;
	int i,j,surplus,random_T = 0, random_obj = 0, scenario, nfrac, nfractotal, ncuts;
	int NGOMORY_CUTS = 0, NMASTER_NODES = 0, NCLIQUE = 0, NCOVER = 0, NZEROHALF = 0, NMIR = 0, NGFC = 0, NFLOW = 0, mastercutcnt;
	double ubobj;
	int solnstat;
#ifdef WRITE_SP
	char scenname[30];
#endif
	FILE *fpout = NULL;      												// Pointer to DEBUG output file
	FILE *fpSolnOut = NULL;
	
	int *base_indices = NULL;
	int *miprowindices = NULL;  //Used to change rhs for subproblem mip
	int origXcols = 0;
	
	double *solnX = NULL;
	double *solnY = NULL;
	double objX;
	double objY;
	double objXlb,objmip;
	// Pointer to optimal solution output file
	char debugoutfname[NAMELEN];
	char corename_lp[NAMELEN];
	char mastername_lp[NAMELEN];
	char subprobname_lp[NAMELEN];
	char subprobname_mip[NAMELEN];
	
	char lpsubname[100];
	
	char soln_filename[NAMELEN];
	char time_filename[NAMELEN];
	char core_filename[NAMELEN];
	char stoch_filename[NAMELEN];
	int UB_FLAG = FALSE;
	int iter_no_gap_chage = 0;												// Number of iterations without gap changing
	char rowname_start[NAMELEN];  // First row name in subproblem
	char colname_start[NAMELEN];  // First column name in subproblem
	
	char    bendersense = 'G';
	
	stochfile_info *stochdataPtr;                                   // Structure to hold STOCH file info
	stochdataPtr = &struct_s;                                       // set its address
	masterproblem_t *masterprobPtr;                                 // Structure to hold master problem data
	masterprobPtr = &master_prob_t;                                 // set its address
	subproblem_t *subprobPtr;                                       // Structure to hold subproblem data
	subprobPtr = &sub_prob_t;                                       // set its address
	
	time_t wallclock_start;
	time_t wallclock_stop;
	double wall_time;
	double cycles;
	time_t masterstart, masterend;
	time_t ubounding_start, ubounding_end;
	double ubounding_time  = 0.0;
	//time_t cutgenstart, cutgenend;
	//time_t firstoptstart, firstoptend;
	double totalmaster = 0.0, totalcutgen = 0.0, totalfirstopt = 0.0;
	cycles = (double) CLOCKS_PER_SEC;
	
	//double *change_rhs = NULL;
	//char *change_sense = NULL;
	//char **cut_names = NULL;
	
	double incumbent = CPX_INFBOUND;
	double lbk = -CPX_INFBOUND;
	double ubk = CPX_INFBOUND;
	double ubmip = CPX_INFBOUND;
	
	int      bendersbeg  = 0;
	int     *bendersind = NULL;
	double  *benderscut = NULL;
	double  *store      = NULL;
	double  *dj         = NULL;//Reduced costs
	double  *pi         = NULL;//Dual variables
	int     sizedj;
	//char    bendersense = 'G';
	double  bendersrhs;
	double  scenbendersrhs;
	
	double *prevSolnX = NULL;
	int k = 0;
	int nrows_core;
	int ncols_core;
	//int nrows_master;
	int ncols_master, nrows_master, nrows_sub, ncols_sub,nrows_submip;
	int nscens;
	double *scenrhs = NULL;
	
	//int cut_max = 0;
	
	CPXENVptr     env           = NULL;								// Subproblem LPs
	CPXLPptr      lp_master     = NULL;                             // pointer to lp from core file
	CPXLPptr      lp_sub        = NULL;                             // pointer to subproblem lp to be created
	CPXLPptr      lp_submip     = NULL;
	
	double gap_old = 100000.0, gap_new = 100000.0;
	
	
	CUTS_FR *frcuts = NULL;
	
	wallclock_start = clock();
	
	strcpy(time_filename,  argv[1]); strcat(time_filename, ".tim");
	strcpy(core_filename,  argv[1]); strcat(core_filename, ".cor");
	strcpy(stoch_filename, argv[1]); strcat(stoch_filename, ".sto");
	
	strcpy(soln_filename,  "LGom");
	strcat(soln_filename, argv[1]);
	strcat(soln_filename, "_1");
	strcat(soln_filename, ".out");
	
	//************* Initialize the CPLEX environment ******************
	env = CPXopenCPLEX (&status);
	
	if ( env == NULL ) {
		char  errmsg[1024];
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		goto TERMINATE;
	}
	
	//************ Turn on output to the screen and set some environment parameters **************
	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr,
				 "Failure to turn on screen indicator, error %d.\n", status);
		goto TERMINATE;
	}
	status = CPXsetdblparam(env, CPX_PARAM_WORKMEM, 4096);
	if (status) {
		fprintf(stderr, "In function main():\n");
		fprintf(stderr, "Failure set workmem.STATUS = %d\n",status);
		goto TERMINATE;
	}
	//************* Copy file and lp names ********************
	strcpy(debugoutfname,  "debug.out");
	strcpy(corename_lp,  "core.lp");
	strcpy(mastername_lp,  "master.lp");
	strcpy(subprobname_lp, "subprob.lp");
	strcpy(subprobname_mip, "subprobmip.lp");
	
	// Open debug output file
	fpout = fopen(debugoutfname,"w" );
	if(fpout == NULL) {
		fprintf (stderr, "In function main():\n");
		fprintf(stderr, "\tCould not open default DEBUG output file %s for writing!\n", debugoutfname);
		fprintf(stderr, "\tTerminating...\n");
		return(0);
	}
	// Open solution output file
	fpSolnOut = fopen(soln_filename,"w" );
	if(fpSolnOut == NULL) {
		fprintf (stderr, "In function main():\n");
		fprintf(stderr, "Could not open solution output file %s for writing!\n", soln_filename);
		fprintf(stderr, "Terminating...\n");
		return(0);
	}
	
	//****** Load the CORE file into an MIP model for creating the master ******
	lp_master = CPXcreateprob (env, &status, "master_lp");
	if ( lp_master == NULL ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failed to create LP.\n");
		goto TERMINATE;
	}
	
	//******* Now read the file, and copy the data into the created lp ******
	status = CPXreadcopyprob (env, lp_master, core_filename, NULL);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Main: Failed to read and copy the problem data from file %s.\n", core_filename);
		goto TERMINATE;
	}
	//****** Load the CORE file into an MIP model for creating subproblem mip ******
	lp_submip = CPXcreateprob (env, &status, "subprobmip");
	if ( lp_submip == NULL ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failed to create CPX subproblem lp_submip LP.\n");
		goto TERMINATE;
	}
	//******* Now read the file, and copy the data into the created lp *******
	status = CPXreadcopyprob (env, lp_submip, core_filename, NULL);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Main: Failed to read and copy the problem data from file %s.\n", core_filename);
		fprintf (stderr, " for the lp_submip LP.\n");
		goto TERMINATE;
	}
	
	//****** Load the CORE file into an MIP model for creating subproblem lp ******
	lp_sub = CPXcreateprob (env, &status, "subprob");
	if ( lp_sub == NULL ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failed to create CPX subproblem lp_sub LP.\n");
		goto TERMINATE;
	}
	
	//* Now read the file, and copy the data into the created lp
	status = CPXreadcopyprob (env, lp_sub, core_filename, NULL);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Main: Failed to read and copy the problem data from file %s.\n", core_filename);
		fprintf (stderr, " for the lp_sub LP.\n");
		goto TERMINATE;
	}
	// *****  Load TIME file data for splitting CORE lp into master and subproblem MIP and LP *****
#ifndef DEBUG_SCR_OUT
	fprintf (stdout, "Loading TIME file info...\n");
#endif
	
	//************** Load TIME file **************
	status = loadTimeFile(rowname_start, colname_start, time_filename, fpout);
	if (status) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "Failed to read TIME file %s! \n", time_filename);
		fprintf(stderr, "Exiting...\n");
		goto TERMINATE;
	}
	
	//**************** Get the number of cols and rows in the master and subproblem *********
	nrows_core = CPXgetnumrows (env, lp_master);
	ncols_core = CPXgetnumcols (env, lp_master);
	status = CPXgetcolindex (env, lp_master, colname_start, &ncols_master);
	if ( status ) {
		fprintf (stderr, "d2algMain():\n");
		fprintf (stderr, "Failure to get column index from CORE LP, error %d.\n", status);
		goto TERMINATE;
	}
	status = CPXgetcolindex (env, lp_master, colname_start, &ncols_master);
	if ( status ) {
		fprintf (stderr, "d2algMain():\n");
		fprintf (stderr, "Failure to get column index from CORE LP, error %d.\n", status);
		goto TERMINATE;
	}
	status = CPXgetcolindex (env, lp_master, colname_start, &ncols_master);
	
	if ( status ) {
		fprintf (stderr, "d2algMain():\n");
		fprintf (stderr, "Failure to get column index from CORE LP, error %d.\n", status);
		goto TERMINATE;
	}
#ifndef DEBUG_SCR_OUT
	fprintf(stdout, " Number of Original Columns in the Master Problem = %d\n", ncols_master);
#endif
	status = CPXgetrowindex (env, lp_master, rowname_start, &nrows_master);
	if ( status ) {
		fprintf (stderr, "d2algMain():\n");
		fprintf (stderr, "Failure to get row index from CORE LP, error %d.\n", status);
		goto TERMINATE;
	}
	nrows_sub = nrows_core - nrows_master;
	ncols_sub = ncols_core - ncols_master;
#ifndef DEBUG_SCR_OUT
	fprintf(stdout, " Number of Rows    in the Master Problem = %d\n", nrows_master);
#endif
	stochdataPtr->nscens = getnumscenarios(stoch_filename);
#ifndef DEBUG_SCR_OUT
	fprintf(stdout, " Number of Scenarios in the Stochastic Program = %d\n", stochdataPtr->nscens);
#endif
	nscens = stochdataPtr->nscens;
	
	//Allocate memory for problem structures
	
	//****************** Allocate memory to STOCH file data structure ********************
	if (DISPLAY_LOG_INFO) {
		fprintf (stdout, " Allocating memory to STOCH data struct...\n");
	}
	status = memAllocStochFileStruct(&struct_s, nrows_sub+ncols_sub, ncols_master, ncols_sub);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failed to allocate memory to STOCH file data structures.\n");
		goto TERMINATE;
	}
	//**************** Allocate memory to master problem data structure *******************
	if (DISPLAY_LOG_INFO) {
		fprintf (stdout, " Allocating memory to Master data struct...\n");
	}
	status = memAllocMasterProblemStructs(masterprobPtr, nrows_master, ncols_master);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to intialize master problem data structures!\n");
		goto TERMINATE;
	}
	
	//**************** Allocate memory to sub problem data structure ***********************
	if (DISPLAY_LOG_INFO) {
		fprintf (stdout, " Allocating memory to subproblem struct...\n");
	}
	status = memAllocSubProblemStruct(subprobPtr, nrows_sub, ncols_sub,  ncols_master, stochdataPtr->nscens);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to intialize subproblem data structures!\n");
		goto TERMINATE;
	}
	
	
	//*******Extract master problem lp data and store into master lp structure*******
	status = loadMasterProblemData(env, lp_master, masterprobPtr, nrows_master,
								   ncols_master, nrows_core, ncols_core, fpout);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to load master problem data into lp data structures.\n");
		goto TERMINATE;
	}
	/**** Add the optimality column "theta" to master ****/
	status = addOptColToMasterProbLP(env, lp_master, masterprobPtr);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to add optimality column to master lp.\n");
		goto TERMINATE;
	}
	masterprobPtr->ncols = CPXgetnumcols(env, lp_master);
#ifndef DEBUG_SCR_OUT
	status = CPXwriteprob(env, lp_master, "MASTER.lp", "lp");
	if (status) {
		fprintf(stderr, "In Function main():\n");
		fprintf(stderr, "\tFailed to write master problem to file\n");
		goto TERMINATE;
	}
#endif
	//*******Set up subproblem mip data  *******
	status = setupSubProbMip(env, lp_submip, subprobPtr, nrows_master,
							 ncols_master, nrows_sub, ncols_sub, fpout);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to load subproblem data into lp data structures.\n");
		goto TERMINATE;
	}
	nrows_submip = CPXgetnumrows (env, lp_submip);
	//Setup the LP subproblem lp data and extract subproblem lp data
	fprintf(stdout, " Loading subproblemdata...\n");
	status = loadSubProblemData(env, lp_sub, subprobPtr, nrows_master, ncols_master,
								nrows_sub, ncols_sub, stochdataPtr, fpout);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to load subproblem data into lp data structures.\n");
		goto TERMINATE;
	}
	
	status = changeSlacksToIntegers(env, lp_submip);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, "\tFailed to convert slacks to integers.\n");
		goto TERMINATE;
	}
	
#ifndef DEBUG_SCR_OUT
	status = CPXwriteprob(env,lp_submip, "SUBMIP.lp", "lp");
#endif
	//****** Get master column names *****
	status = CPXgetcolname(env, lp_master, masterprobPtr->colnames, masterprobPtr->colnamestore,
						   ncols_master*FIELDLEN, &surplus, 0, ncols_master-1);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failure to get master lp colnames, error %d.\n", status);
		fprintf (stderr, " surplus: %d.\n", surplus);
		goto TERMINATE;
	}
	status = CPXgetcolname(env, lp_sub, subprobPtr->colnames, subprobPtr->colnamestore, ncols_sub*FIELDLEN,
						   &surplus, 0, ncols_sub-1);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failure to get subproblem lp colnames, error %d.\n", status);
		fprintf (stderr, " surplus: %d.\n", surplus);
		goto TERMINATE;
	}
	status = CPXgetrowname(env, lp_submip, subprobPtr->rownames, subprobPtr->rownamestore, nrows_submip*FIELDLEN,
						   &surplus, 0, nrows_submip-1);
	if ( status ) {
		fprintf (stderr, "In function main():\n");
		fprintf (stderr, " Failure to get subproblem lp rownames, error %d.\n", status);
		fprintf (stderr, " surplus: %d.\n", surplus);
		goto TERMINATE;
	}
	
	/**********************************************************************
	 Load STOCH file data into the struct stochdataPtr
	 ***********************************************************************/
#ifndef DEBUG_SCR_OUT
	fprintf(stdout, " Loading stoch file data...\n");
#endif
	status = loadStochFile(stochdataPtr, subprobPtr, stoch_filename,
						   &random_T, &random_obj, masterprobPtr->ncols, stdout);
	if (status){
		fprintf (stderr, "In function main():\n");
		fprintf(stderr, "Failed to read STOCH file %s \n", stoch_filename);
		fprintf(stderr, "Exiting...\n");
		goto TERMINATE;
	}
	if(random_T == 1) {
		fprintf(stdout, "Random Technology Matrix Present\n");
	}
	else {
		fprintf(stdout, "Fixed Technology Matrix\n");
	}
	
	//* Allocate space for master solution *
	prevSolnX = (double *)  malloc(masterprobPtr->ncols*sizeof(double));
	solnX =     (double *)  malloc (masterprobPtr->ncols*sizeof(double));
	if ( solnX == NULL || prevSolnX == NULL) {
		fprintf (stderr, "No memory for master solution values.\n");
		goto TERMINATE;
	}
	
	//Compute cmatcnt for T matrix
	origXcols = CPXgetnumcols(env, lp_master) - 1;
	subprobPtr->ncols_T = origXcols;
	fprintf(stdout, "ORIGXCOLS = %d\n",origXcols);
	for (i = 0; i < origXcols - 1; i++) {
		subprobPtr->cmatcnt_T[i] = subprobPtr->cmatbeg_T[i+1] - subprobPtr->cmatbeg_T[i];
	}
	subprobPtr->cmatcnt_T[origXcols - 1] = subprobPtr->nzcnt_T - subprobPtr->cmatbeg_T[origXcols-1];
	
	/**********************************************************************
	 Add the objective column and constraint y0 - cx = 0
	 as the first row and column of the subproblem
	 ***********************************************************************/
	status = addObjConstraint2(env, lp_sub);
	if (status) {
		fprintf(stderr, "Error adding objective constraint in main():\n");
		goto TERMINATE;
	}
	
	// Allocate space for second stage solution
	ncols_sub = CPXgetnumcols(env, lp_sub);
	solnY = (double *)malloc(ncols_sub*sizeof(double));
	if ( solnY == NULL ) {
		fprintf (stderr, "No memory for subproblem solution.\n");
		goto TERMINATE;
	}
	
	/**********************************************************************
	 Solve the master problem once to get a first stage x\in X
	 ***********************************************************************/
	masterstart = clock();
	status = CPXmipopt(env, lp_master);
	masterend   = clock();
	totalmaster   +=  ((double)(masterend - masterstart))/cycles;
	NMASTER_NODES += CPXgetnodecnt(env, lp_master);
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_COVER, &mastercutcnt);
	NCOVER += mastercutcnt;
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_CLIQUE, &mastercutcnt);
	NCLIQUE += mastercutcnt;
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_FLOWCOVER, &mastercutcnt);
	NFLOW += mastercutcnt;
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_FLOWPATH, &mastercutcnt);
	NFLOW+= mastercutcnt;
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_ZEROHALF, &mastercutcnt);
	NZEROHALF += mastercutcnt;
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_FRAC, &mastercutcnt);
	NGOMORY_CUTS += mastercutcnt;
	status = CPXgetnumcuts(env, lp_master, CPX_CUT_MIR, &mastercutcnt);
	NMIR += mastercutcnt;
	
	if (status) {
		fprintf(stderr, "In function main():\n\t Erro optimizing first master\n");
		goto TERMINATE;
	}
	solnstat = CPXgetstat(env, lp_master);
	if (solnstat != CPXMIP_OPTIMAL && solnstat != CPXMIP_OPTIMAL_TOL) {
		fprintf(stderr, "CPLEX master problem not optimal. Solution status = %d\n",solnstat);
		goto TERMINATE;
	}
	status = CPXgetmipx(env, lp_master, solnX, 0, masterprobPtr->ncols-2);
	if (status) {
		fprintf(stderr, "In main():\n\tError extracting master problem solution.CPLEX status = %d\n",status);
		goto TERMINATE;
	}
	status = CPXgetmipx(env, lp_master, prevSolnX, 0, masterprobPtr->ncols-2);
	if (status) {
		fprintf(stderr, "In main():\n\tError extracting master problem solution.CPLEX status = %d\n",status);
		goto TERMINATE;
	}
	
	//Change subproblem to lp
	status = CPXchgprobtype(env, lp_sub, CPXPROB_LP);
	if (status) {
		fprintf(stderr, "Main():\n\tCould not change subproblem to LP\n");
		goto TERMINATE;
	}
	
	//Allocate memory for converting the original T-matrix into row sparse matrix
	nrows_sub = CPXgetnumrows(env, lp_sub);
	subprobPtr->rmatbeg = (int *)malloc((nrows_sub-1)*sizeof(int));
	subprobPtr->rmatcnt = (int *)malloc((nrows_sub-1)*sizeof(int));
	subprobPtr->rmatind = (int *)malloc((subprobPtr->nzcnt_T)*sizeof(int));
	subprobPtr->rmatval = (double *)malloc((subprobPtr->nzcnt_T)*sizeof(double));
	//Convert the T-matrix into row format
	status = convert_col_to_row(subprobPtr->cmatbeg_T, subprobPtr->cmatcnt_T, subprobPtr->cmatind_T, subprobPtr->cmatval_T, subprobPtr->rmatbeg, subprobPtr->rmatcnt, subprobPtr->rmatind, subprobPtr->rmatval, nrows_sub-1, masterprobPtr->ncols-1, subprobPtr->nzcnt_T);
	objX = 0.0;
	
	for (i = 0; i < masterprobPtr->ncols-1 ; i++) {
		objX += solnX[i]*masterprobPtr->obj[i];
	}
	fprintf(stdout, "First stage objective = %f\n",objX);
	
	/******************************************************
	 Allocate memory for Benders' Cut
	 ******************************************************/
	bendersbeg  = 0;
	bendersind  = (int *)   malloc((subprobPtr->ncols_T + 1)*sizeof(int));
	benderscut  = (double *)malloc((subprobPtr->ncols_T + 1)*sizeof(double));
	store       = (double *)malloc(subprobPtr->ncols_T*sizeof(double));
	pi          = (double *)malloc((subprobPtr->nrows -1 + INIT_ROWSIZE_FACTOR*MAX_CUTS_PER_SCENARIO)*sizeof(double));
	sizedj      = CPXgetnumcols(env, lp_sub);
	sizedj--;   //No need to consider dj for y0 column
	dj          = (double *)malloc(sizedj*sizeof(double));
	fprintf(stdout, "SIZEDJ = %d, ncols_T = %d\n",sizedj,subprobPtr->ncols_T);
	
	if (bendersind == NULL || benderscut == NULL || store  == NULL || pi == NULL || dj == NULL) {
		fprintf(stderr, "In main():\n\tError allocating memeory for benders' cut calculations\n");
		goto TERMINATE;
	}
	
	for (i = 0; i < subprobPtr->ncols_T + 1; i++) {
		bendersind[i] = i;
	}
	
	//Scenario RHS Memory - Used for compute rho = r - Tx for original constraints
	// base_indices for changing the rhs values. (1..nrows)
	nrows_sub = CPXgetnumrows(env, lp_sub);
	miprowindices = (int *)malloc((nrows_sub-1)*sizeof(int));
	for (i = 0; i < nrows_sub - 1; i++) {
		miprowindices[i] = i;
	}
	base_indices = (int *)malloc((nrows_sub-1 + INIT_ROWSIZE_FACTOR*MAX_CUTS_PER_SCENARIO)*sizeof(int));
	scenrhs = (double *)calloc(nrows_sub-1 + INIT_ROWSIZE_FACTOR*MAX_CUTS_PER_SCENARIO, sizeof(double));
	if (scenrhs == NULL || base_indices == NULL) {
		fprintf(stderr, "In main():\n\t Failed to allocate memory to Scenario RHS/base_indices \n");
		goto TERMINATE;
	}
	for (i = 0; i < nrows_sub - 1 + INIT_ROWSIZE_FACTOR*MAX_CUTS_PER_SCENARIO; i++) {
		base_indices[i] = (i+1);
	}
	
	frcuts = (CUTS_FR *)malloc(sizeof(CUTS_FR));
	
	if(frcuts == NULL) {
		fprintf(stderr, "In main():\n\t Error allocating memory to frcuts\n");
		status = 1;
		goto TERMINATE;
	}
	status = memAllocFixedRecourseStruct(frcuts, subprobPtr->ncols, subprobPtr->ncols_T, stochdataPtr->nscens);
	if (status) {
		fprintf(stderr, "In main():\n\t Errror in memAllocFixedRecourseStruct()\n ");
		goto TERMINATE;
	}
	
	k = 0; //Initialize iteration index
	/**********************************************************************
	 Start the L-Shaped algorithm
	 ***********************************************************************/
	do {
#ifdef PRINT_SCEN
		fprintf(stdout, "****************************************************\n");
		fprintf(stdout, "ITERATION %d\n",(k+1));
		fprintf(stdout, "****************************************************\n");
#endif
		nfractotal = 0;
		ubobj = 0.0;
		ubmip = 0;
		for (i = 0; i < subprobPtr->ncols_T; i++) {
			benderscut[i] = 0;
		}
		//Set cut coefficient for theta column
		benderscut[subprobPtr->ncols_T] = 1.0;
		//Reset benders RHS
		bendersrhs = 0.0;
		for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
#ifdef PRINT_SCEN
			fprintf(stdout, "++++++++++++++++++++++++++++++++++++++\n");
			fprintf(stdout, "Scenario %d\n",scenario);
			fprintf(stdout, "++++++++++++++++++++++++++++++++++++++\n");
#endif
			nfrac = 0;
			ncuts = 0;
			status = compute_rhs_rows(subprobPtr, stochdataPtr, frcuts, scenrhs, subprobPtr->nrows + frcuts->ncuts, scenario, solnX, masterprobPtr->ncols-1);
			if(status) {
				fprintf(stderr, "In main():Error Computing rhs_rows for scenario %d", scenario);
				goto TERMINATE;
			}
			status = CPXchgrhs(env, lp_sub, (nrows_sub - 1 + frcuts->ncuts), base_indices, scenrhs);
			if (status) {
				fprintf(stderr, "In Main():\n\tError changing RHS for Scenario %d\n. CPLEX status = %d",scenario,status);
				goto TERMINATE;
			}
			
			sprintf(lpsubname,"scenario%d_%d.lp",scenario,k);
			//status = CPXwriteprob(env, lp_sub, lpsubname, "lp");
			if (status) {
				fprintf(stderr, "Error writing scenario subproblem\n");
				goto TERMINATE;
			}
			status = CPXdualopt(env, lp_sub);
			
			if (status) {
				fprintf(stderr, "In Main():\n\tError Solving subproblem Scenario %d. CPLEX status = %d\n",scenario,status);
				goto TERMINATE;
			}
			solnstat = CPXgetstat(env, lp_sub);
			if (solnstat != CPX_STAT_OPTIMAL) {
				fprintf(stderr, "CPLEX status not optimal for scenario %d. CPLEX SOL STAT = %d\n",scenario, solnstat);
				status = CPXwriteprob(env, lp_sub, "notoptimal.lp","lp");
				goto TERMINATE;
			}
			status = CPXgetx(env, lp_sub, solnY, 0, ncols_sub-1);
			if (status) {
				fprintf(stderr, "In main():Error getting solution for scenario %d. CPLEX Status = %d\n",scenario,status);
				goto TERMINATE;
			}
			status = CPXgetobjval(env, lp_sub, &objY);
			if (status) {
				fprintf(stderr, "In main():Error getting objective for scenario %d. CPLEX Status = %d\n",scenario,status);
				goto TERMINATE;
			}
			nfrac = check_sol_frac(solnY, ncols_sub);
			//Calculate the number of cuts to be generated for the current scenario
			(nfrac < CUTS_PER_ROUND)?(ncuts = nfrac):(ncuts = CUTS_PER_ROUND);
			
			NGOMORY_CUTS += ncuts;
			//If the Number of cuts > 0 Generate the cut coefficient and add it to the SubProblem
			if(ncuts > 0) {
				//fprintf(stdout, "Scenario with fractional solution is %d\n", scenario);
				break;
			}
		}
		//fprintf(stdout, "Number of cuts to be generated = %d\n", ncuts);
		status = generate_fixed_recourse_cuts(env, lp_sub, subprobPtr, stochdataPtr, frcuts, ncuts, solnX, scenario);
		for (scenario = 0; scenario < stochdataPtr->nscens; scenario++) {
			scenbendersrhs = 0.0;
			status = compute_rhs_rows(subprobPtr, stochdataPtr, frcuts, scenrhs, subprobPtr->nrows + frcuts->ncuts, scenario, solnX, masterprobPtr->ncols-1);
			
			if(status) {
				fprintf(stderr, "In main():Error Computing rhs_rows for scenario %d", scenario);
				goto TERMINATE;
			}
			status = CPXchgrhs(env, lp_sub, (nrows_sub - 1 + frcuts->ncuts), base_indices, scenrhs);
			if (status) {
				fprintf(stderr, "In Main():\n\tError changing RHS for Scenario %d\n. CPLEX status = %d",scenario,status);
				goto TERMINATE;
			}
			sprintf(lpsubname,"scenarioafter%d_%d.lp",scenario,k);
			//status = CPXwriteprob(env, lp_sub, lpsubname, "lp");
			
			status = CPXdualopt(env, lp_sub);
			
			if (status) {
				fprintf(stderr, "In Main():\n\tError Solving subproblem Scenario %d. CPLEX status = %d\n",scenario,status);
				goto TERMINATE;
			}
			solnstat = CPXgetstat(env, lp_sub);
			if (solnstat != CPX_STAT_OPTIMAL) {
				fprintf(stderr, "CPLEX status not optimal for scenario %d. CPLEX SOL STAT = %d\n",scenario, solnstat);
				status = CPXwriteprob(env, lp_sub, "notoptimal.lp","lp");
				goto TERMINATE;
			}
			status = CPXgetx(env, lp_sub, solnY, 0, ncols_sub-1);
			if (status) {
				fprintf(stderr, "In main():Error getting solution for scenario %d. CPLEX Status = %d\n",scenario,status);
				goto TERMINATE;
			}
			status = CPXgetobjval(env, lp_sub, &objY);
			if (status) {
				fprintf(stderr, "In main():Error getting objective for scenario %d. CPLEX Status = %d\n",scenario,status);
				goto TERMINATE;
			}
			nfrac = check_sol_frac(solnY, ncols_sub);
			nfractotal += nfrac;
			ubobj+= stochdataPtr->scenProb[scenario]*objY;
			
			/*++++++++++++++++++++++++++++++
			 Calculate Benders Cut coeffs
			 ++++++++++++++++++++++++++++++*/
			//Get the dual variables and reduced costs
			
			status = CPXgetpi(env, lp_sub, pi, 1, CPXgetnumrows(env, lp_sub)-1);
			if (status) {
				fprintf(stderr, "In main():\n\tError getting dual variables from CPLEX. SCENARIO = %d,STATUS = %d\n",scenario,status);
				goto TERMINATE;
			}
			//for (i = 0; i < CPXgetnumrows(env, lp_sub)-1; i++) {
			//    fprintf(stdout, "PI[%d] = %f\n", i,pi[i]);
			//}
			status = CPXgetdj(env, lp_sub, dj, 1, CPXgetnumcols(env, lp_sub)-1);
			if (status) {
				fprintf(stderr, "In main():\n\tError getting reduced costs from CPLEX. STATUS = %d. SCENARIO = %d",status,scenario);
				goto TERMINATE;
			}
			//Calculate the optmality cut right hand side for this scenario - first for bounds
			for (i = 0; i < CPXgetnumcols(env, lp_sub)-1; i++) {
				if (myfabs(solnY[i+1] - subprobPtr->lb[i]) < EPSILON_TOLERANCE ) {
					//printf("YELLO! %d %f\n", i+1,subprobPtr->lb[i]);
					scenbendersrhs += subprobPtr->lb[i]*dj[i];
				}
				if (myfabs(solnY[i+1] - subprobPtr->ub[i]) < EPSILON_TOLERANCE) {
					//printf("subprobPtr->ub[%d]*dj[%d] = %f*%f\n", i+1,i+1,subprobPtr->ub[i], dj[i]);
					scenbendersrhs += subprobPtr->ub[i]*dj[i];
				}
			}
			//Calculate the optimality cut rhs for this scenario - second for actual rows
			for (i = 0; i < frcuts->ncuts + subprobPtr->nrows; i++) {
				if (i < subprobPtr->nrows) {
					//fprintf(stdout, "rhs[%d][%d]*pi[%d] = %f*%f\n",scenario,i,i,stochdataPtr->rhs[scenario][i],pi[i]);
					scenbendersrhs+= stochdataPtr->rhs[scenario][i]*pi[i];
				}
				else {
					//fprintf(stdout, "rhs[%d][%d]*pi[%d] = %f*%f\n",scenario,i,i,frcuts->cutrhs[scenario][i-subprobPtr->nrows],pi[i]);
					scenbendersrhs += frcuts->cutrhs[scenario][i-subprobPtr->nrows]*pi[i];
				}
				
			}
			//Update the final benders rhs
			bendersrhs += stochdataPtr->scenProb[scenario]*scenbendersrhs;
			//Calculate the optimality cut coeffs u^T*T
			
			//Reset store to zero
			for (i = 0; i < subprobPtr->ncols_T; i++) {
				store[i] = 0.0;
			}
			
			for (i = 0; i < subprobPtr->nrows_T + frcuts->ncuts; i++) {
				if (i < subprobPtr->nrows_T) {
					for (j = 0; j < subprobPtr->rmatcnt[i]; j++) {
						store[subprobPtr->rmatind[subprobPtr->rmatbeg[i] + j]] += subprobPtr->rmatval[subprobPtr->rmatbeg[i] + j]*pi[i];
					}
				}
				if (i >= subprobPtr->nrows_T) {
					for (j = 0; j < frcuts->rmatcnt_T[i-subprobPtr->nrows_T]; j++) {
						store[frcuts->rmatind_T[frcuts->rmatbeg_T[i-subprobPtr->nrows_T] + j]]
						+= frcuts->rmatval_T[frcuts->rmatbeg_T[i-subprobPtr->nrows_T] + j]*pi[i];
					}
				}
			}
			for (i = 0; i < subprobPtr->ncols_T; i++) {
				if (myfabs(store[i]) > EPSILON_TOLERANCE) {
					benderscut[i] += stochdataPtr->scenProb[scenario]*store[i];
				}
				
			}
			if (UB_FLAG == TRUE) {
				//fprintf(stdout, "UPPERBOUNDING IN ITERATION %d\n",(k+1));
				status = CPXchgrhs(env, lp_submip, nrows_sub-1, miprowindices, scenrhs);
				if (status) {
					fprintf(stderr, "In main():\n\t Error changing RHS for lp_submip-scenario %d in iter %d. STATUS = %d\n",scenario,k,status);
					goto TERMINATE;
				}
				status = CPXmipopt(env, lp_submip);
				if (status) {
					fprintf(stderr, "IN main():\n\tError optimizing the sub problem mip for upper bounding. STAT = %d\n",status);
					goto TERMINATE;
				}
				solnstat = CPXgetstat(env, lp_submip);
				if (solnstat != CPXMIP_OPTIMAL && solnstat != CPXMIP_OPTIMAL_TOL) {
					fprintf(stderr, "Solution status Error in optimizing subproblem MIP. solstat = %d\n",solnstat);
					goto TERMINATE;
				}
				status = CPXgetmipobjval(env, lp_submip, &objmip);
				if (status) {
					fprintf(stderr, "Could not extract subproblem MIP objval. status = %d\n",status);
					goto TERMINATE;
				}
				ubmip += stochdataPtr->scenProb[scenario]*objmip;
			}
		}
		if(UB_FLAG == TRUE) {
			ubounding_end = clock();
			ubounding_time =  (double)((ubounding_end - ubounding_start)/CLOCKS_PER_SEC);
			fprintf(stdout, "UPPER BOUNDING TIME = %0.4f\n", ubounding_time);
		}
		status = CPXaddrows(env, lp_master, 0, 1, subprobPtr->ncols_T + 1, &bendersrhs, &bendersense,
							&bendersbeg, bendersind, benderscut, NULL, NULL);
		if (status) {
			fprintf(stderr, "In Main():\n\t Error adding optimality cuts to the master problem in iteration%d\n",(k+1));
			goto TERMINATE;
		}
		
		//		status = CPXsetintparam (env, CPX_PARAM_SCRIND, 1);
		//		if ( status ) {
		//			fprintf (stderr, "In function main():\n");
		//			fprintf (stderr,
		//					 "Failure to turn on screen indicator, error %d.\n", status);
		//			goto TERMINATE;
		//		}
		
		masterstart = clock();
		status = CPXmipopt(env, lp_master);
		masterend = clock();
		
		//		status = CPXsetintparam (env, CPX_PARAM_SCRIND, 0);
		//		if ( status ) {
		//			fprintf (stderr, "In function main():\n");
		//			fprintf (stderr,
		//					 "Failure to turn on screen indicator, error %d.\n", status);
		//			goto TERMINATE;
		//		}
		//
		
		if (status) {
			fprintf(stderr, "In main():\n\tError optimizing master problem in iteration %d. STATUS = %d\n",(k+1),status);
			status = CPXwriteprob(env,lp_master,"masternotopt.lp","lp");
			goto TERMINATE;
		}
		totalmaster   +=  ((double)(masterend - masterstart))/cycles;
		NMASTER_NODES += CPXgetnodecnt(env, lp_master);
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_COVER, &mastercutcnt);
		NCOVER += mastercutcnt;
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_CLIQUE, &mastercutcnt);
		NCLIQUE += mastercutcnt;
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_FLOWCOVER, &mastercutcnt);
		NFLOW += mastercutcnt;
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_FLOWPATH, &mastercutcnt);
		NFLOW+= mastercutcnt;
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_ZEROHALF, &mastercutcnt);
		NZEROHALF += mastercutcnt;
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_FRAC, &mastercutcnt);
		NGOMORY_CUTS += mastercutcnt;
		status = CPXgetnumcuts(env, lp_master, CPX_CUT_MIR, &mastercutcnt);
		NMIR += mastercutcnt;
		
		solnstat = CPXgetstat(env, lp_master);
		if (solnstat != CPXMIP_OPTIMAL && solnstat != CPXMIP_OPTIMAL_TOL) {
			fprintf(stderr, "In main():\n\tMaster problem solution status is not optimal. SOLSTAT = %d\n",solnstat);
			status = CPXwriteprob(env, lp_master, "finalmaster.lp", "lp");
			goto TERMINATE;
		}
		status = CPXgetmipx(env, lp_master, solnX, 0, subprobPtr->ncols_T-1);
		if (status) {
			fprintf(stderr, "In main():\n\tMaster problem solution could not be extracted in iteration %d\n. STAT = %d\n",(k+1),status);
			goto TERMINATE;
		}
		//for (i = 0; i < masterprobPtr->ncols-1; i++) {
		//    fprintf(stdout, "X{%d} = %f\n",i,solnX[i]);
		//}
		status = CPXgetmipobjval(env, lp_master, &objXlb);
		if (status) {
			fprintf(stderr, "In main():\n\tMaster problem obj fn value could not be extracted in iteration %d\n. STAT = %d\n",(k+1),status);
			goto TERMINATE;
		}
		
		//Update objX - Used to update UB if feasible solution is found
		objX = 0.0;
		for (i = 0; i < masterprobPtr->ncols-1 ; i++) {
			objX += prevSolnX[i]*masterprobPtr->obj[i];
		}
		if (nfractotal == 0) {
			//fprintf(stdout, "Iteration = %d, NFRACTOTAL ==== 0\n",k);
			fprintf(stdout, "+ ");
			ubk = ubobj + objX;
		}
		
		if (ubk < incumbent && nfractotal == 0) {
			incumbent = ubk;
		}
		if (UB_FLAG == TRUE && (ubmip + objX) < incumbent) {
			fprintf(stdout, "UPPERBOUNDING IN ITERATION %d\n",(k+1));
			incumbent = (ubmip + objX);
		}
		UB_FLAG = FALSE;
		lbk = objXlb;
		fprintf(stdout, "\t[LB,UB] = \t[%f,%f], \trel-gap = %0.2f%%\n", lbk,incumbent,((myfabs(incumbent - lbk)/myfabs(incumbent))*100.0));
		gap_new = (myfabs(incumbent - lbk)/(myfabs(incumbent))*100.0);
		//The the new gap is close to the old gap -
		if (myfabs(gap_new - gap_old) > ZERO_TOL) {
			iter_no_gap_chage = 0;
			gap_old = gap_new;
		}
		if (myfabs(gap_new - gap_old) <= ZERO_TOL) {
			iter_no_gap_chage++;
		}
		if (iter_no_gap_chage >= GAP_UBOUNDING || (k+1)%UB_ITER == 0) {
			fprintf(stdout, "UBOUNDING FLAG CHANGED\n");
			UB_FLAG = TRUE;
			ubounding_start = clock();
			iter_no_gap_chage = 0;
		}
		//Copy the new solution into the old solution
		for (i = 0; i < subprobPtr->ncols_T; i++) {
			prevSolnX[i] = solnX[i];
		}
		wallclock_stop = clock();
		wall_time = (double)(wallclock_stop - wallclock_start);
		wall_time = wall_time/cycles;
		if (wall_time > TIME_LIMIT) {
			fprintf(stdout,"Time Limit Exceeded\n");
			goto TERMINATE;
		}
		
		k++;
		//fprintf(stdout, "k %d\n", k);
		
	} while (k < MAX_ITER && (myfabs(incumbent - lbk))/(myfabs(incumbent)) > PERCENT_GAP && k < 10000);
	CPXwriteprob(env, lp_sub, "FINAL.lp", "lp");
TERMINATE:
	printf("---------------------------------------------------------\n");
	fprintf(stdout, "\t\tMaster Problem Iterations = %d\n",k);
	fprintf(stdout, "\t\tMaster Problem Solves = %fs\n",totalmaster);
	fprintf(stdout, "\t\tInitial Subproblem Solves = %fs\n",totalfirstopt);
	fprintf(stdout, "\t\tCut gen + reoptimization = %fs\n",totalcutgen);
	fprintf(stdout, "\t\tSTART TO END RUN TIME = %f\n",wall_time);
	fprintf(stdout, "\t\tNumber of 2nd Stage Gomory Cuts = %d\n",NGOMORY_CUTS);
	fprintf(stdout, "\t\tTotal Master problem Nodes = %d\n",NMASTER_NODES);
	fprintf(stdout, "\t\tTotal Cover Cuts in Master = %d\n",NCOVER);
	fprintf(stdout, "\t\tTotal Clique Cuts in Master = %d\n", NCLIQUE);
	fprintf(stdout, "\t\tTotal MIR cuts in Master = %d\n",NMIR);
	fprintf(stdout, "\t\tTotal Fractional cuts in Master = %d\n",NGFC);
	fprintf(stdout, "\t\tTotal Zero-half cuts in Master = %d\n",NZEROHALF);
	fprintf(stdout, "\t\tTotal Flow cuts in Master = %d\n",NFLOW);
	fprintf(stdout, "\t\tTotal Master Problem Cuts = %d\n", (NFLOW + NZEROHALF + NGFC + NMIR + NCLIQUE + NCOVER));
	printf("---------------------------------------------------------\n");
	
	if(prevSolnX != NULL) {
		free(prevSolnX);
		prevSolnX = NULL;
	}
	if (solnX != NULL) {
		free(solnX);
		solnX = NULL;
	}
	if(solnY != NULL) {
		free(solnY);
		solnY = NULL;
	}
	if(benderscut != NULL) {
		free(benderscut);
		benderscut = NULL;
	}
	if (bendersind != NULL) {
		free(bendersind);
		bendersind = NULL;
	}
	if (store != NULL) {
		free(store);
		store = NULL;
	}
	if (pi != NULL) {
		free(pi);
		pi = NULL;
	}
	if (dj != NULL) {
		free(dj);
		dj = NULL;
	}
	if (miprowindices != NULL) {
		free(miprowindices);
		miprowindices = NULL;
	}
	if (scenrhs != NULL) {
		free(scenrhs);
	}
	
	freeSubProbStruct(subprobPtr, nscens);
	free_stochDataPtr(stochdataPtr);
	
	return 0;
	
}




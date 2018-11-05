/*
**		Authors: Mingfu Zhu & David B. Goldstein
**		Organization: Center for Human Genome Variation, Duke School of Medicine
**		Memory allocation routines borrowed from the book "Numerical Recipes" 
**		by Press, Flannery, Teukolsky, and Vetterling. 																	    
*/


float *vector();
float **matrix();
float **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
float **submatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();

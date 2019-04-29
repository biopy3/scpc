double detFourByFour(double *x);
void dist_dna_raw_without_pairdel(unsigned char *x, int *n, int *s, double *d);
void dist_dna_raw_with_pairdel(unsigned char *x, int *n, int *s, double *d);
void dist_dna_JC69_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				   int *variance, double *var, int *gamma, double *alpha);
void dist_dna_JC69_with_pairdel(unsigned char *x, int *n, int *s, double *d,
				int *variance, double *var, int *gamma, double *alpha);
void dist_dna_K80_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  int *variance, double *var, int *gamma, double *alpha);
void dist_dna_K80_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       int *variance, double *var, int *gamma, double *alpha);

void dist_dna_F81_without_pairdel(unsigned char *x, int *n, int *s, double *d, double *BF,
				  int *variance, double *var, int *gamma, double *alpha);
void dist_dna_F81_with_pairdel(unsigned char *x, int *n, int *s, double *d, double *BF,
			       int *variance, double *var, int *gamma, double *alpha);

void dist_dna_K81_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  int *variance, double *var);
void dist_dna_K81_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       int *variance, double *var);
void dist_dna_F84_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  double *BF, int *variance, double *var);
void dist_dna_F84_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       double *BF, int *variance, double *var);
void dist_dna_T92_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  double *BF, int *variance, double *var);
void dist_dna_T92_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       double *BF, int *variance, double *var);
void dist_dna_TN93_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				   double *BF, int *variance, double *var,
				   int *gamma, double *alpha);
void dist_dna_TN93_with_pairdel(unsigned char *x, int *n, int *s, double *d,
				double *BF, int *variance, double *var,
				int *gamma, double *alpha);
void dist_dna_BH87(unsigned char *x, int *n, int *s, double *d,
		   int *variance, double *var);
void BaseProportion(unsigned char *x, int *n, double *BF);
void SegSites(unsigned char *x, int *n, int *s, int *seg);
void NuclearDiversity(unsigned char *x, int *n, int *s,
		      int *pairdel, double *ans);
void dist_dna(unsigned char *x, int *n, int *s, int *model, double *d,
	      double *BF, int *pairdel, int *variance, double *var,
	      int *gamma, double *alpha);
double* dist_dna_api( unsigned char *Seqs, char *model, int variance, int gamma,
                    int pairwise_deletion, double base_freq, int *num_seqs, int *num_distances,
					int len_seq, int konw_bases);

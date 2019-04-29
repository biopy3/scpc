/* Copyright 2019 */
/* This file is part of the 'SCPC' */
#include <math.h>
#include "scpc.h"
#include "dist_dna.h"
#include "parse_fasta2bin.h"
#include "opra_matrix.h"

const char models[17][11] = {"RAW", "JC69", "K80", "F81", "K81", "F84", "T92", "TN93",
                            "GG95", "LOGDET", "BH87", "PA", "N", "TS", "TV",
                            "INDEL", "INDELBLOCK"};

/* computes directly the determinant of a 4x4 matrix */
double detFourByFour(double *x)
{
    double det, a33a44, a34a43, a34a42, a32a44, a32a43, a33a42, a34a41, a31a44, a31a43, a33a41, a31a42, a32a41;
    a33a44 = x[10]*x[15]; a34a43 = x[14]*x[11];
    a34a42 = x[14]*x[7];  a32a44 = x[6]*x[15];
    a32a43 = x[6]*x[11];  a33a42 = x[10]*x[7];
    a34a41 = x[14]*x[3];  a31a44 = x[2]*x[15];
    a31a43 = x[2]*x[11];  a33a41 = x[10]*x[3];
    a31a42 = x[2]*x[7];   a32a41 = x[6]*x[3];

    det = x[0]*x[5]*(a33a44 - a34a43) + x[0]*x[9]*(a34a42 - a32a44) +
      x[0]*x[13]*(a32a43 - a33a42) + x[4]*x[9]*(a31a44 - a34a41) +
      x[4]*x[13]*(a33a41 - a31a43) + x[4]*x[1]*(a34a43 - a33a44) +       
      x[8]*x[13]*(a31a42 - a32a41) + x[8]*x[1]*(a32a44 - a34a42) +
      x[8]*x[5]*(a34a41 - a31a44) + x[12]*x[1]*(a33a42 - a32a43) +
      x[12]*x[5]*(a31a43 - a33a41) + x[12]*x[9]*(a32a41 - a31a42);

    return det;
}

#define CHECK_PAIRWISE_DELETION\
    if (KnownBase(x[s1]) && KnownBase(x[s2])) L++;\
    else continue;

void dist_dna_raw_without_pairdel(unsigned char *x, int *n, int *s, double *d)
{
    int i1, i2, s1, s2, target, Nd;

    target = 0;
    for (i1 = 1; i1 < *n; i1++)
    {
        for (i2 = i1 + 1; i2 <= *n; i2++)
        {
            Nd = 0;
            for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
                if (DifferentBase(x[s1], x[s2])) Nd++;
            d[target] = ((double) Nd / *s);
            target++;
	    }
    }
}

void dist_dna_raw_with_pairdel(unsigned char *x, int *n, int *s, double *d)
{
    int i1, i2, s1, s2, target, Nd, L;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
                CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    d[target] = ((double) Nd/L);
	    target++;
	}
    }
}

#define COMPUTE_DIST_JC69\
    p = ((double) Nd / *s);\
    if (*gamma)\
      d[target] = 0.75 * *alpha*(pow(1 - 4*p/3, -1/ *alpha) - 1);\
    else d[target] = -0.75*log(1 - 4*p/3);\
    if (*variance) {\
        if (*gamma) var[target] = p*(1 - p)/(pow(1 - 4*p/3, -2/(*alpha + 1)) * L);\
	else var[target] = p*(1 - p)/(pow(1 - 4*p/3, 2)*L);\
    }

void dist_dna_JC69_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				   int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p;

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
	      if (DifferentBase(x[s1], x[s2])) Nd++;
	    COMPUTE_DIST_JC69
	    target++;
	}
    }
}

void dist_dna_JC69_with_pairdel(unsigned char *x, int *n, int *s, double *d,
				int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    COMPUTE_DIST_JC69
	    target++;
	}
    }
}

#define COUNT_TS_TV\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if (IsPurine(x[s1]) && IsPurine(x[s2])) {\
        Ns++;\
        continue;\
    }\
    if (IsPyrimidine(x[s1]) && IsPyrimidine(x[s2])) Ns++;

#define COMPUTE_DIST_K80\
    P = ((double) Ns/L);\
    Q = ((double) (Nd - Ns)/L);\
    a1 = 1 - 2*P - Q;\
    a2 = 1 - 2*Q;\
    if (*gamma) {\
        b = -1 / *alpha;\
    	d[target] = *alpha * (pow(a1, b) + 0.5*pow(a2, b) - 1.5)/2;\
    }\
    else d[target] = -0.5 * log(a1 * sqrt(a2));\
    if (*variance) {\
        if (*gamma) {\
    	    b = -(1 / *alpha + 1);\
    	    c1 = pow(a1, b);\
    	    c2 = pow(a2, b);\
    	    c3 = (c1 + c2)/2;\
    	} else {\
    	  c1 = 1/a1;\
    	  c2 = 1/a2;\
    	  c3 = (c1 + c2)/2;\
    	}\
    	var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;\
    }

void dist_dna_K80_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, Ns, L;
    double P, Q, a1, a2, b, c1, c2, c3;

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_K80
	    target++;
	}
    }
}

void dist_dna_K80_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, Ns, L;
    double P, Q, a1, a2, b, c1, c2, c3;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		COUNT_TS_TV
	    }
	    COMPUTE_DIST_K80
	    target++;
	}
    }
}

#define COMPUTE_DIST_F81\
    p = ((double) Nd/L);\
    if (*gamma) d[target] = E * *alpha * (pow(1 - p/E, -1/ *alpha) - 1);\
    else d[target] = -E*log(1 - p/E);\
    if (*variance) {\
	if (*gamma) var[target] = p*(1 - p)/(pow(1 - p/E, -2/(*alpha + 1)) * L);\
	else var[target] = p*(1 - p)/(pow(1 - p/E, 2)*L);\
    }

void dist_dna_F81_without_pairdel(unsigned char *x, int *n, int *s, double *d, double *BF,
				  int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p, E;

    L = *s;
    E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
	      if (DifferentBase(x[s1], x[s2])) Nd++;
	    COMPUTE_DIST_F81
	    target++;
	}
    }
}

void dist_dna_F81_with_pairdel(unsigned char *x, int *n, int *s, double *d, double *BF,
			       int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p, E;

    E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    COMPUTE_DIST_F81
	    target++;
	}
    }
}

#define COUNT_TS_TV1_TV2\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if ((x[s1] | x[s2]) == 152 || (x[s1] | x[s2]) == 104) {\
        Nv1++;\
        continue;\
    }\
    if ((x[s1] | x[s2]) == 168 || (x[s1] | x[s2]) == 88) Nv2++;


#define COMPUTE_DIST_K81\
    P = ((double) (Nd - Nv1 - Nv2)/L);\
    Q = ((double) Nv1/L);\
    R = ((double) Nv2/L);\
    a1 = 1 - 2*P - 2*Q;\
    a2 = 1 - 2*P - 2*R;\
    a3 = 1 - 2*Q - 2*R;\
    d[target] = -0.25*log(a1*a2*a3);\
    if (*variance) {\
        a = (1/a1 + 1/a2)/2;\
    	b = (1/a1 + 1/a3)/2;\
    	c = (1/a2 + 1/a3)/2;\
      var[target] = (a*a*P + b*b*Q + c*c*R - pow(a*P + b*Q + c*R, 2))/2;\
    }

void dist_dna_K81_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  int *variance, double *var)
{
    int i1, i2, Nd, Nv1, Nv2, L, s1, s2, target;
    double P, Q, R, a1, a2, a3, a, b, c;

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = Nv1 = Nv2 = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV1_TV2
	    }
	    COMPUTE_DIST_K81
	    target++;
	}
    }
}

void dist_dna_K81_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       int *variance, double *var)
{
    int i1, i2, Nd, Nv1, Nv2, L, s1, s2, target;
    double P, Q, R, a1, a2, a3, a, b, c;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = Nv1 = Nv2 = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
	        COUNT_TS_TV1_TV2
	    }
	    COMPUTE_DIST_K81
	    target++;
	}
    }
}

#define PREPARE_BF_F84\
    A = (BF[0]*BF[2])/(BF[0] + BF[2]) + (BF[1]*BF[3])/(BF[1] + BF[3]);\
    B = BF[0]*BF[2] + BF[1]*BF[3];\
    C = (BF[0] + BF[2])*(BF[1] + BF[3]);

#define COMPUTE_DIST_F84\
   P = ((double) Ns/L);\
   Q = ((double) (Nd - Ns)/L);\
   d[target] = -2*A*log(1 - (P/(2*A) - (A - B)*Q/(2*A*C))) + 2*(A - B - C)*log(1 - Q/(2*C));\
   if (*variance) {\
       t1 = A*C;\
       t2 = C*P/2;\
       t3 = (A - B)*Q/2;\
       a = t1/(t1 - t2 - t3);\
       b = A*(A - B)/(t1 - t2 - t3) - (A - B - C)/(C - Q/2);\
       var[target] = (a*a*P + b*b*Q - pow(a*P + b*Q, 2))/2;\
   }

void dist_dna_F84_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, A, B, C, a, b, t1, t2, t3;

    PREPARE_BF_F84
    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_F84
	    target++;
	}
    }
}

void dist_dna_F84_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, A, B, C, a, b, t1, t2, t3;

    PREPARE_BF_F84

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		CHECK_PAIRWISE_DELETION
		COUNT_TS_TV
	    }
	    COMPUTE_DIST_F84
	    target++;
	}
    }
}

#define COMPUTE_DIST_T92\
    P = ((double) Ns/L);\
    Q = ((double) (Nd - Ns)/L);\
    a1 = 1 - P/wg - Q;\
    a2 = 1 - 2*Q;\
    d[target] = -wg*log(a1) - 0.5*(1 - wg)*log(a2);\
    if (*variance) {\
        c1 = 1/a1;\
        c2 = 1/a2;\
        c3 = wg*(c1 - c2) + c2;\
        var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;\
    }

void dist_dna_T92_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				  double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, wg, a1, a2, c1, c2, c3;

    L = *s;
    wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_T92
	    target++;
	}
    }
}

void dist_dna_T92_with_pairdel(unsigned char *x, int *n, int *s, double *d,
			       double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, wg, a1, a2, c1, c2, c3;

    wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_T92
	    target++;
	}
    }
}

/* returns 1 if one of the base is adenine and
   the other one is guanine surely, 0 otherwise */
#define AdenineAndGuanine(a, b) (a | b) == 200

/* returns 1 if one of the base is cytosine and
   the other one is thymine surely, 0 otherwise */
#define CytosineAndThymine(a, b) (a | b) == 56

#define PREPARE_BF_TN93\
    gR = BF[0] + BF[2];\
    gY = BF[1] + BF[3];\
    k1 = 2 * BF[0] * BF[2] / gR;\
    k2 = 2 * BF[1] * BF[3] / gY;\
    k3 = 2 * (gR * gY - BF[0]*BF[2]*gY/gR - BF[1]*BF[3]*gR/gY);

#define COUNT_TS1_TS2_TV\
    if (DifferentBase(x[s1], x[s2])) {\
        Nd++;\
        if (AdenineAndGuanine(x[s1], x[s2])) {\
            Ns1++;\
    	continue;\
        }\
        if (CytosineAndThymine(x[s1], x[s2])) Ns2++;\
    }

#define COMPUTE_DIST_TN93\
    P1 = ((double) Ns1/L);\
    P2 = ((double) Ns2/L);\
    Q = ((double) (Nd - Ns1 - Ns2)/L);\
    w1 = 1 - P1/k1 - Q/(2*gR);\
    w2 = 1 - P2/k2 - Q/(2*gY);\
    w3 = 1 - Q/(2*gR*gY);\
    if (*gamma) {\
        k4 = 2*(BF[0]*BF[2] + BF[1]*BF[3] + gR*gY);\
    	b = -1 / *alpha;\
    	c1 = pow(w1, b);\
    	c2 = pow(w2, b);\
    	c3 = pow(w3, b);\
    	c4 = k1*c1/(2*gR) + k2*c2/(2*gY) + k3*c3/(2*gR*gY);\
    	d[target] = *alpha * (k1*pow(w1, b) + k2*pow(w2, b) + k3*pow(w3, b) - k4);\
    } else {\
        k4 = 2*((BF[0]*BF[0] + BF[2]*BF[2])/(2*gR*gR) + (BF[2]*BF[2] + BF[3]*BF[3])/(2*gY*gY));\
    	c1 = 1/w1;\
    	c2 = 1/w2;\
    	c3 = 1/w3;\
    	c4 = k1 * c1/(2 * gR) + k2 * c2/(2 * gY) + k4 * c3;\
    	d[target] = -k1*log(w1) - k2*log(w2) - k3*log(w3);\
    }\
    if (*variance)\
      var[target] = (c1*c1*P1 + c2*c2*P2 + c4*c4*Q - pow(c1*P1 + c2*P2 + c4*Q, 2))/L;

void dist_dna_TN93_without_pairdel(unsigned char *x, int *n, int *s, double *d,
				   double *BF, int *variance, double *var,
				   int *gamma, double *alpha)
{
    int i1, i2, k, Nd, Ns1, Ns2, L, target, s1, s2;
    double P1, P2, Q, A, B, C, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, c4, b;

    L = *s;

    PREPARE_BF_TN93

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns1 = Ns2 = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		COUNT_TS1_TS2_TV
	    }
	    COMPUTE_DIST_TN93
	    target++;
	}
    }
}

void dist_dna_TN93_with_pairdel(unsigned char *x, int *n, int *s, double *d,
				double *BF, int *variance, double *var,
				int *gamma, double *alpha)
{
    int i1, i2, k, Nd, Ns1, Ns2, L, target, s1, s2;
    double P1, P2, Q, A, B, C, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, c4, b;

    PREPARE_BF_TN93

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns1 = Ns2 = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		CHECK_PAIRWISE_DELETION
		COUNT_TS1_TS2_TV
	    }
	    COMPUTE_DIST_TN93
	    target++;
	}
    }
}

#define DO_CONTINGENCY_NUCLEOTIDES\
    switch (x[s1]) {\
    case 136 : m = 0; break;\
    case 72 : m = 1; break;\
    case 40 : m = 2; break;\
    case 24 : m = 3; break;\
    }\
    switch (x[s2]) {\
    case 72 : m += 4; break;\
    case 40 : m += 8; break;\
    case 24 : m += 12; break;\
    }\
    Ntab[m]++;


void dist_dna_BH87(unsigned char *x, int *n, int *s, double *d,
		   int *variance, double *var)
/* <FIXME>
   For the moment there is no need to check for pairwise deletions
   since DO_CONTINGENCY_NUCLEOTIDES considers only the known nucleotides.
   In effect the pairwise deletion has possibly been done before.
   The sequence length(s) are used only to compute the variances, which is
   currently not available.
   </FIXME> */
{
    int i1, i2, k, kb, s1, s2, m, Ntab[16], ROWsums[4], ndim = 4, info, ipiv[16];
    double P12[16], P21[16], U[16];

    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }

            /* get the rowwise sums of Ntab */
            ROWsums[0] = Ntab[0] + Ntab[4] + Ntab[8] + Ntab[12];
            ROWsums[1] = Ntab[1] + Ntab[5] + Ntab[9] + Ntab[13];
            ROWsums[2] = Ntab[2] + Ntab[6] + Ntab[10] + Ntab[14];
            ROWsums[3] = Ntab[3] + Ntab[7] + Ntab[11] + Ntab[15];

            for (k = 0; k < 16; k++)
              P12[k] = ((double) Ntab[k]);

            /* scale each element of P12 by its rowwise sum */
            for (k = 0; k < 4; k++)
              for (kb = 0; kb < 16; kb += 4)
            	P12[k + kb] = P12[k + kb]/ROWsums[k];

            d[*n*(i2 - 1) + i1 - 1] = -log(detFourByFour(P12))/4;

            /* compute the columnwise sums of Ntab: these
               are the rowwise sums of its transpose */
            ROWsums[0] = Ntab[0] + Ntab[1] + Ntab[2] + Ntab[3];
            ROWsums[1] = Ntab[4] + Ntab[5] + Ntab[6] + Ntab[7];
            ROWsums[2] = Ntab[8] + Ntab[9] + Ntab[10] + Ntab[11];
            ROWsums[3] = Ntab[12] + Ntab[13] + Ntab[14] + Ntab[15];

            /* transpose Ntab and store the result in P21 */
            for (k = 0; k < 4; k++)
               for (kb = 0; kb < 4; kb++)
            	 P21[kb + 4*k] = Ntab[k + 4*kb];

            /* scale as above */
            for (k = 0; k < 4; k++)
              for (kb = 0; kb < 16; kb += 4)
            	P21[k + kb] = P21[k + kb]/ROWsums[k];

            d[*n*(i1 - 1) + i2 - 1] = -log(detFourByFour(P21))/4;
	}
    }
}


void BaseProportion(unsigned char *x, int *n, double *BF)
{
    int i, m;

    m = 0;
    for (i = 0; i < *n; i++) {
        if (KnownBase(x[i])) {
	    m++;
	    switch (x[i]) {
	    case 136 : BF[0]++; break;
	    case 40 : BF[1]++; break;
	    case 72 : BF[2]++; break;
	    case 24 : BF[3]++; break;
	    }
	}
    }
    for (i = 0; i < 4; i++) BF[i] /= m;
}

void SegSites(unsigned char *x, int *n, int *s, int *seg)
{
    int i, j;
    unsigned char basis;

    for (j = 0; j < *s; j++) {
        i = *n * j;
	while (!KnownBase(x[i])) i++;
	basis = x[i];
	i++;
	while (i < *n * (j + 1)) {
	    if (x[i] == basis) i++;
	    else {
	        seg[j] = 1;
		break;
	    }
	}
    }
}

void NuclearDiversity(unsigned char *x, int *n, int *s,
		      int *pairdel, double *ans)
{
    int i1, i2, s1, s2, Nd, L;

    if (!*pairdel) L = *s;

    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = 0;
	    if (*pairdel) L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
                CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    *ans += ((double) Nd/L);
	}
    }
    *ans /= (*n * (*n - 1)/2);
}

void dist_dna(unsigned char *x, int *n, int *s, int *model, double *d,
	      double *BF, int *pairdel, int *variance, double *var,
	      int *gamma, double *alpha)
{
    switch (*model) {
    case 0 : if (pairdel) dist_dna_raw_with_pairdel(x, n, s, d);
             else dist_dna_raw_without_pairdel(x, n, s, d); break;

    case 1 : if (pairdel) dist_dna_JC69_with_pairdel(x, n, s, d, variance, var, gamma, alpha);
             else dist_dna_JC69_without_pairdel(x, n, s, d, variance, var, gamma, alpha); break;

    case 2 : if (pairdel) dist_dna_K80_with_pairdel(x, n, s, d, variance, var, gamma, alpha);
             else dist_dna_K80_without_pairdel(x, n, s, d, variance, var, gamma, alpha); break;

    case 3 : if (pairdel) dist_dna_F81_with_pairdel(x, n, s, d, BF, variance, var, gamma, alpha);
             else dist_dna_F81_without_pairdel(x, n, s, d, BF, variance, var, gamma, alpha); break;

    case 4 : if (pairdel) dist_dna_K81_with_pairdel(x, n, s, d, variance, var);
             else dist_dna_K81_without_pairdel(x, n, s, d, variance, var); break;

    case 5 : if (pairdel) dist_dna_F84_with_pairdel(x, n, s, d, BF, variance, var);
             else dist_dna_F84_without_pairdel(x, n, s, d, BF, variance, var); break;

    case 6 : if (pairdel) dist_dna_T92_with_pairdel(x, n, s, d, BF, variance, var);
             else dist_dna_T92_without_pairdel(x, n, s, d, BF, variance, var); break;

    case 7 : if (pairdel) dist_dna_TN93_with_pairdel(x, n, s, d, BF, variance, var, gamma, alpha);
             else dist_dna_TN93_without_pairdel(x, n, s, d, BF, variance, var, gamma, alpha); break;

    case 10 : dist_dna_BH87(x, n, s, d, variance, var); break;
    }
}

double* dist_dna_api( unsigned char *Seqs, char *model, int variance, int gamma,
                    int pairwise_deletion, double base_freq, int *num_seqs, int *num_distances,
                    int len_seq, int konw_bases)
{
    double alpha, var, *d;
    int index_model = 0, i; 
    *num_distances = 0;

    for( i = 0; i < 17; i++)
    {
        if(strcmp(models[i], model) == 0)
        {
            index_model = i;
            break;
        }else if(i != 16){
            ;
        }else{
            printf("\nThe given model is wrong, ours supples :\
            \"RAW\", \"JC69\", \"K80\", \"F81\", \"K81\", \"F84\", \"T92\", \"TN93\", \
            \"GG95\", \"LOGDET\", \"BH87\", \"PARALIN\", \"N\", \"TS\", \"TV\", \
            \"INDEL\", \"INDELBLOCK\" models.\n");
            exit(1);
        }
    }

    if(index_model == 10 && variance)
    {
        printf("\nWarning: computing variance not available for model BH87,\nwe \
        already changed the variance parameter to false and continue.\n");
        variance = 0;
    }

    if( gamma && (index_model == 0 || ( 4 <= index_model && index_model <= 6 ) || ( 8 <= index_model && index_model <= 16 ) ) )
    {
        printf("\nWarning: gamma-correction not available for %s, we \
        already changed the gamma parameter to false and continue.\n", models[index_model]);
        gamma = 0;
    }

    if( index_model == 3 || ( 5 <= index_model && index_model <= 7 ) )
    {
        if ( (base_freq>=-0.0001) && (base_freq<=0.00001) )
        {
            if( base_freq <= 0.0 || base_freq >= 1.0 )
            {
                printf("\nWarning: base frequencies invalid, we \
                    already changed the base_freq parameter to the default, which is \
                    the base frequencies are computed from the whole set of sequences.\n");
                base_freq = konw_bases / ( *num_seqs*len_seq );
            }

        }else{
            base_freq = konw_bases / ( *num_seqs*len_seq );
        }
    }else{
        base_freq = 0;
    }

    if (15 <= index_model && index_model <= 16)
    {
        pairwise_deletion = 1;
    }

    if( pairwise_deletion == 0 )
    {
        GlobalDeletionDNA(Seqs, num_seqs, &len_seq);
    }

    if ( index_model == 10)
    {
        *num_distances = *num_seqs * *num_seqs;
        printf("num_distances:%d\n", *num_distances);
    }else{
        *num_distances = *num_seqs * (*num_seqs - 1) / 2;
        printf("num_distances:%d\n", *num_distances);
    }
    
    if( variance != 0 )
    {
        var = *num_distances;
    }else{
        var = 0;
    }

    if( gamma == 0 )
    {
        gamma = alpha = 0;
    }else{
        alpha = gamma;
        gamma = 1;
    }
 
    d = (double *)malloc((*num_distances)*sizeof(double));
    printf("The base_fre:%f\n", base_freq);
    printf("index_model:%d\n", index_model);
    dist_dna( Seqs, num_seqs, &len_seq, &index_model,
            d, &base_freq, &pairwise_deletion,
            &variance, &var, &gamma, &alpha );
    return d;
}

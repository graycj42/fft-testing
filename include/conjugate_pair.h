#ifndef _CONJUGATE_PAIR_H
#define _CONJUGATE_PAIR_H

#define N 128

#include <math.h>

#include <stdint.h>

typedef struct {
	double real;
	double imag;
} complex;

extern "C" {void cpfft_init(complex tw[N/8]);}
extern "C" {static inline void cpfft_bf4(unsigned s, complex out[N], complex w);}
extern "C" {void cpfft_dfi(complex in[N], complex out[N], complex twid[N]);}

#define CADD(Z, X, Y)  ({ \
	(Z).real = (X).real + (Y).real;   \
	(Z).imag = (X).imag + (Y).imag; \
	})

#define CSUB(Z, X, Y)  ({ \
	(Z).real = (X).real - (Y).real; \
	(Z).imag = (X).imag - (Y).imag; \
	})

#define CMUL(Z, X, Y)  ({ \
	(Z).real = (X).real * (Y).real - (X).imag * (Y).imag; \
	(Z).imag = (X).real * (Y).imag + (X).imag * (Y).real; \
	})

#define CEXP(X, Y) ({ \
    (Y).real = cos((X).imag); \
    (Y).imag = sin((X).imag); \
    })

//X + i*Y
#define CADD_X_IMUL_Y(Z, X, Y) ({ \
    (Z).real = (X).real - (Y).imag; \
    (Z).imag = (X).imag + (Y).real; \
})

//X - i*Y
#define CSUB_X_IMUL_Y(Z, X, Y) ({ \
    (Z).real = (X).real + (Y).imag; \
    (Z).imag = (X).imag - (Y).real; \
})




#endif
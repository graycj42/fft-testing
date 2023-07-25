#ifndef _CONJUGATE_PAIR_H
#define _CONJUGATE_PAIR_H

#define DATASIZE 1024

#include <math.h>
#include <stdint.h>
#include <hls_vector.h>

typedef struct {
	double real;
	double imag;
} cmplx_type;


extern "C" {void cpfft_init(cmplx_type tw[]);}
extern "C" {void cpfft_bf4(unsigned s, hls::vector<cmplx_type, DATASIZE> out, uint32_t offset, cmplx_type w, cmplx_type write_only[4]);}
extern "C" {void cpfft_dfi(hls::vector<cmplx_type, DATASIZE>& in, hls::vector<cmplx_type, DATASIZE>& out, cmplx_type twid[DATASIZE/8]);}

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
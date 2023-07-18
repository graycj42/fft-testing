// #include "header.h"
#include "constants.h"
#include "conjugate_p.h"




//Compute twiddle factors
extern "C" {
    void cpfft_init(cmplx_type tw[])
    {
        #pragma HLS inline
        cmplx_type exp;
        exp.real = 0;
        exp.imag = 0;
        for (uint32_t i = 0; i < N/8; i++) {
            exp.imag = -2 * M_PI * i/N;
            CEXP(exp, tw[i]);
        }

    }
}

extern "C" {
    void cpfft_bf4(uint32_t s, cmplx_type out[N], cmplx_type w)
    {
        #pragma HLS inline
        cmplx_type a = out[0];
        cmplx_type b = out[s];
        cmplx_type c = out[s*2];
        cmplx_type d = out[s*3];
        cmplx_type conj_w;
        conj_w.real = w.real;
        conj_w.imag = -1*w.imag;
        cmplx_type wc; 
        CMUL(wc, w, c);
        cmplx_type conj_wd;
        CMUL(conj_wd, conj_w, d);
        
        cmplx_type cpxsum;
        CADD(cpxsum, wc, conj_wd);
        //out[0] = a + (w * c + conj(w) * d);
        CADD(out[0], a, cpxsum);
        //out[s * 2] = a - (w * c + conj(w) * d);
        CSUB(out[s*2], a, cpxsum);

        CSUB(cpxsum, wc, conj_wd);
        //out[s * 3] = b + I * (w * c - conj(w) * d);
        CADD_X_IMUL_Y(out[s*3], b, cpxsum);
        //out[s ] = b - I * (w * c - conj(w) * d);
        CSUB_X_IMUL_Y(out[s], b, cpxsum);

    }
}

//depth first iterative fft algorithm
extern "C" {void cpfft_dfi(cmplx_type in[N], cmplx_type out[N], cmplx_type twid[N])
    {
        uint32_t log2_n = 31 - __builtin_clz(N);
        // uint32_t log2_n = 31 - clz(N);
        uint32_t r = 32 - log2_n;
        uint32_t p = 0; 
        uint32_t q = 0;

        //mapping input to output indices
        // uint32_t test_array[32] = {0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5};
        // uint32_t k = 0;
        outer_loop : for(uint32_t h = 0; h < N; h += 2){

            //generate the binary carry sequence
            uint32_t h2 = h + 2;
            // uint32_t c = 30 - __builtin_clz(h ^ h2);
            // k++;
            // uint32_t c = 30 - clz(h ^ h2);

            /* input indices */
            uint32_t i0 = (p - q) >> r;
            uint32_t i1 = i0 ^ (N >> 1);
            uint32_t j;
            if (carry[h/2] & 1) { // stage 1
                j = 2;
                // out[h ] = in[i0];
                out[h].real = in[i0].real;
                out[h].imag = in[i0].imag;
                // out[h + 1] = in[i1];
                out[h + 1].real = in[i1].real;
                out[h + 1].imag = in[i1].imag;
                cmplx_type one;
                one.real = 1;
                one.imag = 0;
                cpfft_bf4(1, out + h - 2, one);
            } else { 
                j = 1;
                // stage 0
                CADD(out[h], in[i0], in[i1]);
                CSUB(out[h + 1], in[i0], in[i1]);
            }

            //higher stages
            inner_loop : for (j; j < 9; j += 2) {
                if(j < carry[h/2]){
                    uint32_t s = 1 << j;
                    uint32_t z = h2 - 4 * s;
                    uint32_t t = log2_n - j - 2;

                    //butterfly blocks ------ SOURCE OF CRASH
                    butterfly_loop : for (uint32_t b = 1; b < s / 2; b++) {
                    // w = e^(-2 * M_PI * I * b / s / 4);
                        cmplx_type w = twid[b << t];
                        cmplx_type w_rev;
                        w_rev.real = -1*w.imag;
                        w_rev.imag = -1*w.real;
                        
                        cpfft_bf4(s, out + z + b, w);
                        cpfft_bf4(s, out + z + s - b, w_rev);
                    }

                    cmplx_type spec_case; //From what I understand, these cover specific twiddle factor scenarios, i.e. 2k pi and k pi/4
                    spec_case.real = 1;
                    spec_case.imag = 0;
                    //spec_case = 1+0i aka 1
                    cpfft_bf4(s, out + z, spec_case);
                    spec_case.real = M_SQRT1_2;
                    spec_case.imag = -1*M_SQRT1_2;
                    //spec_case = (1-i)*(1/sqrt(2))
                    cpfft_bf4(s, out + z + s/2, spec_case);
                }
            }

            /* next input index */
            uint32_t m2 = 0x20000000 >> carry[h/2];
            uint32_t m1 = m2 - 1;
            uint32_t m = p & m2;
            q = (q & m1) | m;
            p = (p & m1) | ((m ^ m2) << 1);
        }
    }
}

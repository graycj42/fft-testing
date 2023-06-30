#include "ap_int.h"
#include "conjugate_pair.h"

//Compute twiddle factors





extern "C" void cpfft_init(complex tw[N/8])
{
    complex exp;
    exp.real = 0;
    exp.imag = 0;

for (unsigned i = 0; i < N/8; i++) {
    exp.imag = -2 * M_PI * i/N;
    CEXP(exp, tw[i]);
}

}

extern "C" {static inline void cpfft_bf4(unsigned s, complex out[N], complex w)
    {
        complex a = out[0];
        complex b = out[s];
        complex c = out[s*2];
        complex d = out[s*3];
        complex conj_w = w;
        conj_w.imag *= -1;
        complex wc; 
        CMUL(wc, w, c);
        complex conj_wd;
        CMUL(conj_wd, conj_w, d);
        
        complex cpxsum;
        CADD(cpxsum, wc, conj_wd);
        CADD(out[0], a, cpxsum);
        
        CSUB(out[s*2], a, cpxsum);

        CSUB(cpxsum, wc, conj_wd);
        CADD_X_IMUL_Y(out[s], b, cpxsum);
        CSUB_X_IMUL_Y(out[s*3], b, cpxsum);

    }
}

//depth first iterative fft algorithm
extern "C" {void cpfft_dfi(complex in[N], complex out[N], complex twid[N])
    {
        unsigned log2_n = __builtin_clz(N);
        unsigned r = 32 - log2_n;
        uint32_t p = 0; 
        uint32_t q = 0;

        //mapping input to output indices
        for(uint32_t h, h2 = 0; h < N; h = h2){

            //generate the binary carry sequence
            h2 = h + 2;
            unsigned c = 30 - __builtin_clz(h ^ h2);

            /* input indices */
            unsigned i0 = (p - q) >> r;
            unsigned i1 = i0 ^ (N >> 1);

            if (c & 1) { // stage 1
                // out[h ] = in[i0];
                out[h].real = in[i0].real;
                out[h].imag = in[i0].imag;
                // out[h + 1] = in[i1];
                out[h + 1].real = in[i1].real;
                out[h + 1].imag = in[i1].imag;
                complex one;
                one.real = 1;
                one.imag = 0;
                cpfft_bf4(1, out + h - 2, one);
            } else { // stage 0
                CADD(out[h], in[i0], in[i1]);
                CSUB(out[h + 1], in[i0], in[i1]);
            }

            //higher stages
            for (unsigned j = 1 + (c & 1); j < c; j += 2) {
                unsigned s = 1 << j;
                unsigned r = h2 - 4 * s;
                unsigned t = log2_n - j - 2;

                // butterfly blocks
                for (unsigned b = 1; b < s / 2; b++) {
                // w = e^(-2 * M_PI * I * b / s / 4);
                    complex w = twid[b << t];
                    complex w_rev;
                    w_rev.real = -1*w.imag;
                    w_rev.imag = -1*w.real;
                    
                    cpfft_bf4(s, out + r + b, w);
                    cpfft_bf4(s, out + r + s - b, w_rev);
                }

                complex spec_case; //From what I understand, these cover specific twiddle factor scenarios, i.e. 2k pi and k pi/4
                spec_case.real = 1;
                spec_case.imag = 0;
                cpfft_bf4(s, out + r, spec_case);
                spec_case.real = M_SQRT1_2;
                spec_case.imag = -1*M_SQRT1_2;
                cpfft_bf4(s, out + r + s/2, spec_case);
            }

            /* next input index */
            uint32_t m2 = 0x20000000 >> c;
            uint32_t m1 = m2 - 1;
            uint32_t m = p & m2;
            q = (q & m1) | m;
            p = (p & m1) | ((m ^ m2) << 1);
        }

    }
}
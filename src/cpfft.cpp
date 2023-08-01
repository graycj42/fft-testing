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
    void cpfft_bf4(uint32_t s, cmplx_type out[N], cmplx_type w, cmplx_type write_only[4])
    {
        // #pragma HLS inline
        #pragma HLS pipeline II=1
        #pragma HLS INTERFACE mode=bram port=out storage_type=ram_t2p
        #pragma HLS array_reshape variable=out type=block factor=N/2

        cmplx_type compute[4];
        #pragma HLS array_partition variable=compute

        for(uint8_t i = 0; i<4; i++){
            #pragma HLS unroll
            compute[i].real = out[s*i].real;
            compute[i].imag = out[s*i].imag;
        }
        cmplx_type conj_w;
        conj_w.real = w.real;
        conj_w.imag = -1*w.imag;
        cmplx_type wc; 
        CMUL(wc, w, compute[2]);
        cmplx_type conj_wd;
        CMUL(conj_wd, conj_w, compute[3]);
        
        cmplx_type cpxsum;
        CADD(cpxsum, wc, conj_wd);
        //out[0] = a + (w * c + conj(w) * d);
        CADD(write_only[0], compute[0], cpxsum);
        //out[s * 2] = a - (w * c + conj(w) * d);
        CSUB(write_only[2], compute[0], cpxsum);

        CSUB(cpxsum, wc, conj_wd);
        //out[s * 3] = b + I * (w * c - conj(w) * d);
        CADD_X_IMUL_Y(write_only[3], compute[1], cpxsum);
        //out[s ] = b - I * (w * c - conj(w) * d);
        CSUB_X_IMUL_Y(write_only[1], compute[1], cpxsum);
    }
}

//depth first iterative fft algorithm
extern "C" {void cpfft_dfi(cmplx_type in[N], cmplx_type out[N], cmplx_type twid[N])
    {
        #pragma HLS INTERFACE mode=bram port=out storage_type=ram_t2p
        #pragma HLS array_reshape variable=out type=block factor=N/2
        #pragma HLS INTERFACE mode=bram port=in storage_type=ram_t2p
        #pragma HLS INTERFACE mode=bram port=twid storage_type=ram_t2p
        #pragma HLS array_reshape variable=in type=block factor=N/2
        #pragma HLS array_reshape variable=twid type=block factor=N/16

        uint32_t log2_n = 31 - __builtin_clz(N);

        uint32_t r = 32 - log2_n;
        uint32_t p = 0; 
        uint32_t q = 0;
        uint32_t h2 = 0;
        cmplx_type bf4_write[4];
        #pragma HLS array_partition variable=bf4_write

        //mapping input to output indices
        outer_loop : for(uint32_t h = 0; h < N; h += 2){
            #pragma HLS loop_tripcount min=512 max=512
            #pragma HLS pipeline
            //generate the binary carry sequence
            h2 = h + 2;
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
                cpfft_bf4(1, out + h - 2, one, bf4_write);
                out[h-2] = bf4_write[0];
                out[h-1] = bf4_write[1];
                out[h]   = bf4_write[2];
                out[h+1] = bf4_write[3];
            } else { 
                j = 1;
                // stage 0
                CADD(out[h], in[i0], in[i1]);
                CSUB(out[h + 1], in[i0], in[i1]);
            }

            //higher stages
            inner_loop : for (j; j < carry[h/2]; j += 2) {
                    #pragma HLS loop_tripcount min=1 max=9 avg=2
                    #pragma HLS pipeline
                    uint32_t s = 1 << j;
                    uint32_t z = h2 - 4 * s;
                    uint32_t t = log2_n - j - 2;
                    //butterfly blocks ------ SOURCE OF CRASH
                    butterfly_loop : for (uint32_t b = 1; b < s / 2; b++) {
                    // w = e^(-2 * M_PI * I * b / s / 4);

                        #pragma HLS loop_tripcount min=1 max=127 avg=3
                       #pragma HLS dataflow

                        cmplx_type w = twid[b << t];
                        cmplx_type w_rev;
                        w_rev.real = -1*w.imag;
                        w_rev.imag = -1*w.real;
                        
                        cpfft_bf4(s, out + z + b, w, bf4_write);

                        for(uint8_t i = 0; i<4; i++){
                            #pragma HLS unroll
                            out[z + b + s*i]  = bf4_write[i];
                        }
                        cpfft_bf4(s, out + z + s - b, w_rev, bf4_write);
                        for(uint8_t i = 0; i<4; i++){
                            #pragma HLS unroll
                            out[z - b + s*(i+1)]  = bf4_write[i];
                        }
                    }

                    cmplx_type spec_case; //From what I understand, these cover specific twiddle factor scenarios, i.e. 2k pi and k pi/4
                    spec_case.real = 1;
                    spec_case.imag = 0;
                    //spec_case = 1+0i aka 1
                    cpfft_bf4(s, out + z, spec_case, bf4_write);
                    out[z]       = bf4_write[0];
                    out[z + s]   = bf4_write[1];
                    out[z + s*2] = bf4_write[2];
                    out[z + s*3] = bf4_write[3];
                    spec_case.real = M_SQRT1_2;
                    spec_case.imag = -1*M_SQRT1_2;
                    //spec_case = (1-i)*(1/sqrt(2))
                    cpfft_bf4(s, out + z + s/2, spec_case, bf4_write);
                    out[z + s/2]       = bf4_write[0];
                    out[z + s/2 + s]   = bf4_write[1];
                    out[z + s/2 + s*2] = bf4_write[2];
                    out[z + s/2 + s*3] = bf4_write[3];
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

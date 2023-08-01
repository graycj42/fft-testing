#include "constants.h"
#include "conjugate_p.h"

extern "C" {
    void krnl(data_t* a, data_t* b) {
        #pragma HLS INTERFACE m_axi port = a bundle = gmem0
        #pragma HLS INTERFACE m_axi port = b bundle = gmem1

        cmplx_type fft_in[DATA_SIZE], fft_out[DATA_SIZE];
        #pragma HLS bind_storage variable=fft_out type=ram_t2p impl=bram
        #pragma HLS array_reshape variable=fft_out type=block factor=N/2
        #pragma HLS bind_storage variable=fft_in type=ram_t2p impl=bram
        #pragma HLS array_reshape variable=fft_in type=block factor=N/2


        cmplx_type twid[DATA_SIZE/8];
        // #pragma HLS bind_storage variable=twid type=ram_t2p impl=bram
        #pragma HLS array_reshape variable=twid type=block factor=N/16

        cpfft_init(twid);

        for (int i = 0; i < DATA_SIZE; i++) {
            fft_in[i].real = a[i];
            fft_in[i].imag = 0.0;
        }
        
        // pease_fft(fft_in, fft_out);
        cpfft_dfi(fft_in, fft_out, twid);
        
        for (int i = 0; i < DATA_SIZE; i++) {
            b[i] = fft_out[i].real;
            //b[i] = fft_out[i].imag;
        }

        /***********TWIDDLE FACTOR INITIALIZATION TESTING**********/
        // int j = 0;
        // for(int i = 0; i < DATA_SIZE/8; i++){
        //     b[j] = twid[i].real;
        //     b[j+1] = twid[i].imag;
        //     // b[j] = i;
        //     // b[j+1] = i+1;
        //     j+=2;
        // }
    }
}
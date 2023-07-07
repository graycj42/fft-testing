#include "constants.h"

extern "C" {
    void krnl(data_t* a, data_t* b) {
        #pragma HLS INTERFACE m_axi port = a bundle = gmem0
        #pragma HLS INTERFACE m_axi port = b bundle = gmem1

        cmplx_type fft_in[DATA_SIZE], fft_out[DATA_SIZE];
        cmplx_type twid[DATA_SIZE/8];
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
    }
}
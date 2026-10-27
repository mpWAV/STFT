/*
STFT - Short Time Fourier Transoform cpp
https://github.com/kooBH/STFT
nine4409@sogang.ac.kr
*/

#ifndef _H_STFT_
#define _H_STFT_

#include "Ooura_FFT.h"
#include "DFT.h"
#include "WindowFunction.h"
#include "OA.h"

class STFT{
  private : 
    const double MATLAB_scale = 32768;

    HannWindow *hw;

    /* WOLA synthesis window : w[n] / sum_k w^2[n - k*shift]
       Without this division the overlap-added w^2 envelope is not constant,
       which adds a periodic ripple at (sr/shift_size) Hz to the output.
       Matches ONNXSTFT.window_istft on the python side. */
    bool wola = false;
    double* syn_window = nullptr;
    inline void synthesis_window(double* buffer);
    inline void synthesis_window(double** buffer, int n_channels);
    FFTbase*fft;
    OA *ap;

    int channels;
    int frame_size;
    int shift_size;
    int ol;

    double**buf = nullptr;
    double* tmp = nullptr;

    bool is_power_of_two(int num) {
      if (num <= 0) {
        return false;
      }
      while (num > 1) {
        if (num % 2 != 0) {
          return false;
        }
        num /= 2;
      }
      return true;
    }

  public :
    inline STFT(int channels,int frame,int shift,bool wola_synthesis=false);
    inline ~STFT();

    /* in from input device or file
    
      in : raw buffer from wav or mic
      length : shift_size * channels   (for not fully occupied input)
      out : STFTed buffer [channels][frame_size + 2] (half FFT in complex)
      */
    inline void stft(short*in,int length,double**out);
    inline void istft(double**in,short*out);

    inline void stft(short*in,int length,double**out,int target_channels);

    /* 2-D raw input STFT
       in  : [channels][shift_size]   raw data in double
       out : [channels][frame_size+2]
    */
    inline void stft(double** in, double** out);
    inline void istft(double** in, double** out);
	  inline void stft(double** in, double** out,int target_channels);

    inline void stft(float** in, double** out);
    inline void istft(double** in, float** out);

    /* Single-Channel STFT   
      in : 1 x shift
      out : 1 x frame_size + 2 (half FFT in complex)
    */
    inline void stft(short* in, double* out);
    inline void stft(float* in, float* out);
    inline void stft(double* in, double* out);

    /* Single-Channel ISTFT   
      in : 1 x frame_size + 2 (half FFT in complex)
      out : 1 x shift_size     */
    inline void istft(double* in, short* out);
    inline void istft(float* in, float* out);
    inline void istft(double* in, double* out);

    //for separated 3-channels wav
    inline void stft(short* in_1, short* in_2, short* in_3, int length, double** out);
};


STFT::STFT(int channels_,int frame_,int shift_,bool wola_synthesis){
  int i;
  channels = channels_;
  frame_size = frame_;
  shift_size = shift_;
  ol = frame_size - shift_size;

  wola = wola_synthesis;
  if (wola) {
    /* analysis : plain hann, no normalization (== torch.hann_window) */
    hw = new HannWindow(frame_size, shift_size, -1, false);
    /* synthesis : w[n] / sum_k w^2[n - k*shift] */
    const double PI_ = 3.141592653589793;
    syn_window = new double[frame_size];
    double* env = new double[shift_size];
    for (int n = 0; n < shift_size; n++) {
      double e = 0.0;
      for (int j = n; j < frame_size; j += shift_size) {
        double w = 0.5 - 0.5 * cos(2.0 * PI_ * (j / (double)frame_size));
        e += w * w;
      }
      env[n] = e;
    }
    for (int n = 0; n < frame_size; n++) {
      double w = 0.5 - 0.5 * cos(2.0 * PI_ * (n / (double)frame_size));
      syn_window[n] = w / env[n % shift_size];
    }
    delete[] env;
  }
  else {
    hw = new HannWindow(frame_size, shift_size);
  }

  if (is_power_of_two(frame_size)) {
    //printf("STFT::FFT\n");
    fft = new Ooura_FFT(frame_size, channels);
  }
  else {
    //printf("STFT::DFT\n");
    fft = new DFTbrute(frame_size, channels);
  }
  ap = new OA(frame_size, shift_size, channels);

  buf =  new double*[channels];
  for(i=0;i<channels;i++){
    buf[i] = new double[frame_size];
    memset(buf[i],0,sizeof(double)*frame_size);
  }

  tmp = new double[frame_size + 2];
  memset(tmp, 0, sizeof(double) * (frame_size + 2));
}

STFT::~STFT(){
  int i;
  delete hw;
  if (syn_window) delete[] syn_window;
  delete fft;
  delete ap;
  for(i=0;i<channels;i++)
    delete[] buf[i];
  delete[] buf;
  if (tmp)delete[] tmp;

}

inline void STFT::synthesis_window(double* buffer) {
  if (!wola) { hw->Process(buffer); return; }
  for (int j = 0; j < frame_size; j++) buffer[j] *= syn_window[j];
}

inline void STFT::synthesis_window(double** buffer, int n_channels) {
  if (!wola) { hw->Process(buffer, n_channels); return; }
  for (int i = 0; i < n_channels; i++)
    for (int j = 0; j < frame_size; j++) buffer[i][j] *= syn_window[j];
}

void STFT::stft(short*in,int length,double**out){
  int i,j;
  /*** Shfit & Copy***/
  for (j = 0; j < channels; j++) {
    for (i = 0; i < ol; i++) {
      buf[j][i] = buf[j][i + shift_size];
    }
  }
  //// EOF
  if(length!=shift_size*channels){

    length = length/channels;
    for (i = 0; i < length; i++) {
      for (j = 0; j < channels; j++)
        buf[j][i + ol]
          =  (double)(in[i * channels+ j]);
    }
    for (i = length; i < shift_size; i++) {
      for (j = 0; j < channels; j++)
        buf[j][i + ol] = 0;
    }
    //// continue
  }else{
    for (i = 0; i < shift_size; i++) {
      for (j = 0; j < channels; j++){
        buf[j][i + ol] 
          = (double)(in[i * channels+ j]);
      }
    }
  }
  /*** Copy input -> hann_input buffer ***/
  for (i = 0; i < channels; i++)
    memcpy(out[i], buf[i], sizeof(double) * frame_size);

  // scaling for precision
  for (i = 0; i < channels; i++)
      for (j = 0; j < frame_size; j++)
        out[i][j] /= MATLAB_scale;

  /*** Window ***/
  hw->Process(out, channels);

  /*** FFT ***/
  fft->FFT(out);
}

void STFT::stft(short* in, int length, double** out, int target_channels) {
  int tmp = channels;
  channels = target_channels;
  stft(in, length, out);
  channels = tmp;
}


void STFT::istft(double**in,short*out){
  /*** iFFT ***/
  fft->iFFT(in);

  /*** Window ***/
  synthesis_window(in, channels);

  // scaling for precision
  for (int i = 0; i < channels; i++)
    for (int j = 0; j < frame_size; j++)
      in[i][j] *= MATLAB_scale;

  /*** Output ***/
  memcpy(out,ap->Overlap(in),sizeof(short)*shift_size*channels);
}


// Single-Channel double
void STFT::stft(short* in, double* out){
	int i;
    /*** Shfit & Copy***/
    for (i = 0; i < ol; i++) {
        buf[0][i] = buf[0][i + shift_size];
    }
    for (i = 0; i < shift_size; i++)
        buf[0][ol + i] = static_cast<double>(in[i]);


    memcpy(out, buf[0], sizeof(double) * frame_size);

    // scaling for precision
    for (int j = 0; j < frame_size; j++)
      out[j] /= MATLAB_scale;

    /*** Window ***/
    hw->Process(out);

    /*** FFT ***/
    fft->FFT(out);
}

void STFT::stft(float* in, float* out) {
  int i;
  /*** Shfit & Copy***/
  for (i = 0; i < ol; i++) {
    buf[0][i] = buf[0][i + shift_size];
  }
  for (i = 0; i < shift_size; i++)
    buf[0][ol + i] = static_cast<double>(in[i]);

  memcpy(tmp, buf[0], sizeof(double) * frame_size);

  /*** Window ***/
  hw->Process(tmp);

  /*** FFT ***/
  fft->FFT(tmp);

  for (int i = 0; i < frame_size+2; i++) {
    out[i] = static_cast<float>(tmp[i]);
  }
}

void STFT::stft(double* in, double* out) {
    int i;
	/*** Shfit & Copy***/
    for (i = 0; i < ol; i++) {
        buf[0][i] = buf[0][i + shift_size];
    }
    for (i = 0; i < shift_size; i++)
        buf[0][ol + i] = in[i];

    memcpy(out, buf[0], sizeof(double) * frame_size);

    /*** Window ***/
    hw->Process(out);

    /*** FFT ***/
    fft->FFT(out);
}

void STFT::stft(double** in, double** out) {
	/*** Shfit & Copy***/
#pragma omp parallel for
    for (int j = 0; j < channels; j++) {
      for (int i = 0; i < ol; i++) {
        buf[j][i] = buf[j][i + shift_size];
      }
      for (int i = 0; i < shift_size; i++){
        buf[j][ol + i] = in[j][i];
      }
        memcpy(out[j], buf[j], sizeof(double) * frame_size);
    }

  // scaling for precision
  for (int i = 0; i < channels; i++)

    /*** Window ***/
    hw->Process(out,channels);

    /*** FFT ***/
    fft->FFT(out);
}

void STFT::stft(double** in, double** out,int target_channels){
		/*** Shfit & Copy***/
#pragma omp parallel for
    for (int j = 0; j < target_channels; j++) {
      for (int i = 0; i < ol; i++) {
        buf[j][i] = buf[j][i + shift_size];
      }
      for (int i = 0; i < shift_size; i++){
        buf[j][ol + i] = in[j][i];
      }
        memcpy(out[j], buf[j], sizeof(double) * frame_size);
    }

  // scaling for precision
  for (int i = 0; i < target_channels; i++)
    /*** Window ***/
    hw->Process(out,target_channels);

    /*** FFT ***/
    fft->FFT(out,target_channels);
}

void STFT::istft(double* in, short* out) {
  /*** iFFT ***/
  fft->iFFT(in);

  /*** Window ***/
  synthesis_window(in);

  for (int j = 0; j < frame_size; j++)
    in[j] *= MATLAB_scale;

  /*** Output ***/
  memcpy(out,ap->Overlap(in),sizeof(short)*shift_size);
}

void STFT::istft(float* in, float* out) {
  for(int i=0;i<frame_size+2;i++)
    tmp[i] = static_cast<double>(in[i]);

  /*** iFFT ***/
  fft->iFFT(tmp);

  /*** Window ***/
  synthesis_window(tmp);

  //for (int j = 0; j < frame_size; j++)tmp[j] *= MATLAB_scale;

  /*** Output ***/
  ap->Overlap(tmp);
  //memcpy(out,ap->Overlap(in),sizeof(short)*shift_size);
  for (int i = 0; i < shift_size; i++)
    out[i] = static_cast<float>(ap->Get_buf()[0][i]);
}

void STFT::istft(double* in, double* out) {
  /*** iFFT ***/
  fft->iFFT(in);

  /*** Window ***/
  synthesis_window(in);

  //for (int j = 0; j < frame_size; j++)in[j] *= MATLAB_scale;

  /*** Output ***/
  ap->Overlap(in);
  //memcpy(out,ap->Overlap(in),sizeof(short)*shift_size);
  for(int i=0;i<shift_size;i++)
    out[i] = ap->Get_buf()[0][i];
}

void STFT::istft(double**in,double**out){
  /*** iFFT ***/
  fft->iFFT(in);

  /*** Window ***/
  synthesis_window(in, channels);

  /*** Output ***/
  ap->Overlap(in);
  for (int i = 0; i < channels; i++)
    for(int j=0;j<shift_size;j++)
      out[i][j] = ap->Get_buf()[i][j];
}

void STFT::stft(float** in, double** out) {
	/*** Shfit & Copy***/
#pragma omp parallel for
    for (int j = 0; j < channels; j++) {
      for (int i = 0; i < ol; i++) {
        buf[j][i] = buf[j][i + shift_size];
      }
      for (int i = 0; i < shift_size; i++){
        buf[j][ol + i] = static_cast<double>(in[j][i]);
      }
      
      for (int i = 0; i < frame_size; i++)
        out[j][i] = buf[j][i];
      memcpy(out[j], buf[j], sizeof(double) * frame_size);
    }

    /*** Window ***/
    hw->Process(out,channels);

    /*** FFT ***/
    fft->FFT(out);
}

void STFT::istft(double**in,float**out){
  /*** iFFT ***/
  fft->iFFT(in);

  /*** Window ***/
  synthesis_window(in, channels);

  /*** Output ***/
  for (int i = 0; i < channels; i++)
    for(int j=0;j<shift_size;j++)
      out[i][j] = static_cast<float>(ap->Get_buf()[i][j]);
}



#endif

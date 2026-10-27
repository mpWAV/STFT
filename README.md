# STFT

STFT(Short Time Fourier Transform), ISTFT(Inverse - Short Time Fourier Transform) for wav,mic input.  
And other pre-processings. 

+ Provides 25%,50% overlap STFT process.   
+ DFT, iDFT for the case that n_fft(frame_size) is not power of 2.  
+ Mel-filterbank for both HTK-style(pytorch default) and slaney-style(librosa-default)  
+ LPF for perceptual quality. 

## TODO
- ERB filterbank & normalization  
- Vorbis window  & spectrogram normalization  
- C++ wrapper for C  

## NOTE

```git clone --recursive https://github.com/kooBH/STFT.git```

To build test code you need to clone ```--recursive``` to use submodule  ```koobh/WAV.git```
If you already cloned, then,
```
git submoudle init
git submodule update
```
to use submodule.

## About FFT  
I'm currently using FFT of [Ooura](http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html).  
Since, it is fastest FFT in a single header file.   
But, sometimes (usually  not), there are errors between MATLAB FFT output and Ooura FFT output.     
If you need to perfectly same output as MATLAB, you have to use other FFT library.  

## CMake

```CMAKE
...
set(SRC)
set(INCL)

# it will append sources and paths to SRC and INCL
include(${PATH_TO_STFT}/STFT/cpp/STFT_cpp.cmake )  

...

add_executable(${PROJECT_NAME} ${SRC})
target_include_directories(${PROJECT_NAME}    PRIVATE    ${INCL}  )

```

## EXAMPLE

+ Directly from buffer  

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame_size,shift_size);
WAV input;
WAV output(channels,sample_rate);
input.OpenFile('input.wav');

short* buf_in; //[channels * shift_size]
double **data; //[channels][frame_size]
short* buf_out; //[channels * shift_size]

while(!input.IsEOF()){
    length = input.ReadUnit(buf_in,shift*channels);
    process.stft(buf_in,length,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

+ From 2D array

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame_size,shift_size);
WAV input;
WAV output(channels,sample_rate);
input.OpenFile('input.wav');

double **raw; //[channels][shift_size]
double **data; //[channels][frame_size]
short* buf_out; //[channels * shift_size]

while(!input.IsEOF()){
    input.Convert2Array(raw);
    process.stft(raw,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

# 포크 메모 (mpWAV)

이 레포는 [kooBH/STFT](https://github.com/kooBH/STFT) 의 포크다. upstream 에
푸시 권한이 없어서 2026-10-27 에 갈라놨다.

+ 갈라진 지점 : `36b3a9d`
+ 작업 브랜치 : `ko`
+ 전체 차이 보기 : `git diff 36b3a9d..ko`
+ 쓰는 쪽 : `mpWAV/mpNC`, `mpWAV/2025_irlink_mpNCw` (IR_Link)

upstream 과 다른 점은 `STFT.h` 하나뿐이다 (`e99525d`).

## `STFT.h` : WOLA 합성 윈도 옵션

생성자에 네 번째 인자가 붙었다. 기본값이 `false` 라서 기존 호출자는 영향이 없다.

```cpp
STFT(int channels, int frame, int shift, bool wola_synthesis = false);
```

`wola_synthesis = true` 일 때 :

+ 분석 윈도는 정규화 없는 순수 Hann (`HannWindow(frame, shift, -1, false)`).
  `torch.hann_window` 와 같다.
+ 합성 윈도 `w[n] / sum_k w^2[n - k*shift]` 를 생성자에서 미리 계산해
  `syn_window` 에 넣어 둔다.
+ `istft` 계열이 `hw->Process(...)` 대신 새로 추가된 `synthesis_window(...)` 를
  부른다. `wola == false` 면 그 안에서 기존 `hw->Process(...)` 를 그대로
  부르므로 동작이 같다.

왜 : 나누지 않으면 overlap-add 된 w^2 포락선이 일정하지 않아서 출력에
`sr / shift_size` Hz 주기의 리플이 생긴다. 파이썬 쪽
`ONNXSTFT.window_istft` 와 출력을 맞추려면 필요하다.

`mpNC` 는 이 옵션을 켜고 쓴다 (`mpNC_impl.cpp` 의
`new STFT(n_channel, n_fft, n_hop, true)`). 즉 **`ko` 브랜치 없이는 mpNC 가
컴파일되지 않는다** - 3인자 생성자에 기본값이 없기 때문이다.

## upstream 과 다시 합칠 때

`origin/main` 이 `36b3a9d` 보다 두 커밋 앞서 있다 (2026-10-27 기준) :

```
32c1147 string.h
c3a4332 [Add]Clear()
```

아직 가져오지 않았다. 둘 다 `STFT.h` 를 건드리므로 위 변경과 충돌할 수 있다
(`FFT.h`, `OA.h` 도 바뀐다). 머지한 뒤에는 `wola = true` 로 돌려서 출력이
그대로인지 확인할 것.

upstream 이 받아야 할 수정이므로, 권한이 생기면 PR 로 올려서 포크를 없애는
쪽이 낫다.

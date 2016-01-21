# ECGSYN version for Arduino

>ECGSYN generates a synthesized ECG signal with user-settable mean heart rate, number of beats, sampling frequency, waveform morphology (P, Q, R, S, and T timing, amplitude,and duration), standard deviation of the RR interval, and LF/HF ratio (a measure of the relative contributions of the low and high frequency components of the RR time series to total heart rate variability). Using a model based on three coupled ordinary differential equations, ECGSYN reproduces many of the features of the human ECG, including beat-to-beat variation in morphology and timing, respiratory sinus arrhythmia, QT dependence on heart rate, and R-peak amplitude modulation. The output of ECGSYN may be employed to assess biomedical signal processing techniques which are used to compute clinical statistics from the ECG.

**ecgsyn-arduino** is a software implementation of [ECGSYN](http://www.physionet.org/physiotools/ecgsyn/) for Arduino. It has been tested on Arduino Uno R3, but it will work on any compatible variant.

## Limitations

As Arduino's memory is limited, this implementation has the following limitations:

- Maximum sampling frequency of 64.
- It does not use the peak detection function.

## License
Distributed under the MIT license. See *LICENSE*.

## Version History
- 1.0 Initial version (October 2013).
	Complete software implementation for Arduino.

## Acknowledgements
This implementation was developed during an internship at SAPYC Group, Universidad de las Fuerzas Armadas ESPE.

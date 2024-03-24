# mel-spectrum-calculator

## [Javadocs](https://jasgarcha.github.io/mel-spectrum-calculator/)

## Logic
`MelSpectrum` contains the Mel Filter bank frequency table as a constant two-dimensional array, indexed from lower bound to upper bound of the frequency range (`lower bound =  center frequency-(band width/2)` & `upper bound = center frequency+(bandwidth/2)`).

The filter function `M` takes the parameter `l`, the `l`th filter of the Mel Filter bank, and a continuous frequency in Hertz. `M` returns 1 if the frequency is in the range of the `l`th filter, 0 if not, due to the assumption of a square versus triangular filter (throws the frequency into the right “bucket”).

`computeMelSpectrum` computes the Mel spectrum of a set of sampled data. Its parameters are `N` (the number of sampled data points), sampling frequency, and arrays in which the real and imaginary terms of the result `S(l)` are stored. An overridden method returns `S(l)` as a two-dimensional array of size 20x1, analogous to a 20x1 column vector. 

`printMelSpectrumWithSteps` computes the Mel Spectrum and prints the intermediate calculations. 

## Usage
`MelSpectrumExample` demonstrates the `computeMelSpectrum` and `printMelSpectrumWithSteps` methods using an input of `n = 10` data points and `f = frequency of 2000 Hz`.

Manually tweak `n` and `f` in the code as desired.

## Build
The `discreteFourierTransform` package is required for operations in `MelSpectrum`. 
It is included here, but can also be found in my [discrete-fourier-transform-calculator-verbose](https://github.com/jasgarcha/discrete-fourier-transform-calculator-verbose) repo.
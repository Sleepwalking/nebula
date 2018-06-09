Nebula
===

Zero-data (yet trainable) probabilistic fundamental frequency estimator.

[Paper](https://arxiv.org/abs/1710.11317) available on arXiv, to be presented in Interspeech 2018.

Proudly implemented in GNU Octave. Probably won't work in MATLAB.

License: GPL v3

What is Nebula
---

Nebula is a low-gross-error pitch estimator designed for speech synthesis systems. The main idea is to define a very rough probabilistic distribution of speech signals, generate synthetic signals from the distribution by means of Monte-Carlo simulation, and train an estimator on the synthetic data.

The difficulty with this type of data-free approaches is the poor generalization from synthetic speech to real speech data. This is tackled in Nebula by factorizing the problem into a lot of time-frequency local pieces and training a frequency-dependent model for each band. The prediction from all frequency-dependent models are fused together by taking the average log posterior.

The said factorization is made possible by re-using Hideki Kawahara's SNR and instantaneous frequency feature extractors, which have very good time-frequency resolution. The models are simply Gaussian Mixture Models (GMM). These low-dimensional GMMs can be efficiently converted into conditional forms when all except one of the variables are known.

> H. Kawahara, Y. Agiomyrgiannakis, and H. Zen, "Using instantaneous frequency and aperiodicity detection to estimate F0 for high-quality speech synthesis," in 9th ISCA Workshop on Speech Synthesis, Sunnyvale, 2016.

### Why the name Nebula

The Monte-Carlo simulation generates lots of points in a 6-dimensional space. After some dimensionality reduction and visualization by 3D scatter plot, the resulting image looks like some glowing celestial objects and hence the algorithm was given the name Nebula.

![fig2](https://user-images.githubusercontent.com/4531595/41190370-17881fdc-6c18-11e8-91a6-a013db1a0543.jpg)

How to use
---

Like many of my other projects on Github, some C-implemented functions depend on [`ciglet`](https://github.com/Sleepwalking/ciglet) (a library of DSP snippets) and [`libgvps`](https://github.com/Sleepwalking/libgvps) (a library implementing many variants of Viterbi algorithm).

* `cd` into `ciglet`, run `make single-file`. This creates `ciglet.h` and `ciglet.c` under `ciglet/single-file/`. Copy and rename this directory to `nebula/external/ciglet`.
* Put `libgvps` under `nebula/external/` and run `make` from `nebula/external/libgvps/`.
* Finally launch Octave from `nebula/`, run `startup` in Octave.

### Inference

To estimate F0 from a given audio signal,

```matlab
[x fs] = audioread('./test.wav');
M = load_model('./model/'); % load from directory ./model
f0 = nebula_est(M, x, fs, 0.005); % estimate F0 at a 0.005s interval
```

A pretrained model is included in `nebula/model`. It contains 36 GMMs for all frequency bands and a calibration file `Lcal`.

You can also let Nebula output the F0 posterior map,

```matlab
[f0 v pv lmap] = nebula_est(M, x, fs, 0.005);
imagesc(log(lmap));
```

![fig-lmap-1](https://user-images.githubusercontent.com/4531595/41190368-12fdfeaa-6c18-11e8-87cb-78b610b000e6.png)

### Training

The training part depends on [SPTK](http://sp-tk.sourceforge.net/).

```matlab
make_random_dataset; % this is going to take a while
train_gmm; % the actual GMM training and calibration
```

How to cite this work
---

K. Hua, "Nebula: F0 estimation and voicing detection by modeling the statistical properties of feature extractors," in *Interspeech*, Hyderabad, 2018 [to be presented].

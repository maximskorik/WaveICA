# WaveICA
WaveICA is an R package designed to remove batch effects from metabolomics untargeted MS data based on wavelet analysis.
The WaveICA method uses the time trend of samples over the injection order, decomposes the original data into 
multi-scale data with different features, extracts and removes the batch effect information in multi-scale data,
and obtains clean data.

This WaveICA is has been slightly refactored and modified by RECETOX: is has fewer dependencies, and includes both original
WaveICA and its successor WaveICA 2.0.

To use WaveICA, first install the package:
```
devtools::install_github("RECETOX/WaveICA",host="https://api.github.com")
```
Then you can use either of two batch effects removal functions provided by the package. If your data has been measured
in multiple batches and contains batch information, we recommend you to use `WaveICA` function. You can use it as follows:
```
library(WaveICA)
features <- WaveICA(features_table, wf, batch, factorization, group, K, t, t2, alpha) 
```
+ `features_table` is a sample-by-feature matrix in the following format:

| feature#1                   | feature#2                    | ... | feature#M                    |
|-----------------------------|------------------------------|-----|------------------------------|
| sample#1_feature#1 intensity| ...                          | ... | ...                          |
| ...                         | sample#2_feature#2 intensity | ... | ...                          |
| ...                         | ...                          | ... | ...                          |
| ...                         | ...                          | ... | sample#N_feature#M intensity |

+ `wf` is the wavelet function to be used. It is a character string, consisting of a filter prefix and its length.
The supported filters (and its prefixes) are: Daubechies (`d`), Least Asymetric (`la`), Best Localized (`bl`), and
Coiflet (`c`). The supported length are as follows:
   + **Daubechies** 2,4,6,8,10,12,14,16,18,20;
   + **Least Asymetric** 8,10,12,14,16,18,20;
   + **Best Localized** 14,18,20;
   + **Coiflet** 6,12,18,24,30;
  
  Thus, to obtain **Coiflet** filter with length **12**, you should use `wf="c12"`. The only exception is the **Daubechies**
with length **2** which is, by convention, named `haar`.
+ `batch` is a numeric vector with samples batch numbers. It's size has to be *number of samples* **x** 1.
+ `factorization` is used to specify data-matrix factorization method. It can be either `"svd"` for **singular value
decomposition** or `"stICA"` for **unbiased ICA**.
+ `group` is a numeric vector representing the type of sample. The accepted values are **0** for **blank**, **1** for 
**sample** and **2** for **QC**. The vector's size has to be *number of samples* **x** 1.
+ `K` is the maximal number of independent components if `factorization="stICA"` or the number of singular vectors 
computed by **singular value decomposition** if `factorization="svd"`.
+ `t` and `t2` are floating point numbers between 0 and 1 specifying the association of components with groups.
+ `alpha` is a floating point number specifying the trade-off between spacial (features) and temporal (samples) ICA.

Alternatively, if your data has been measured in a single batch or the batch information is not available, we recommend
using `WaveICA_nonbatchwise` function. You can use it as follows:
```
library(WaveICA)
features <- WaveICA_nonbatchwise(features_table, wf, injection_order, alpha, cutoff, K) 
```
+ `features_table`, `wf`, `alpha`, and `K` are the same as above.
+ `injection_order` is a numeric vector representing the injection order of samples. It's size has to be *number of samples*
**x** 1.
+ `cutoff` is a floating point number between 0 and 1, which specifies the thershold of the variation explained by the
injection order for independent components.

Both functions return a sample-by-feature matrix with the same size as the input `features` matrix.

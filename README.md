#### 2024-01-04 

## Overview

This README file provides information about the data and the computer code used to generate the results presented in Ardia et al. (2024), **Is it alpha or beta? Decomposing hedge fund returns when models are misspecified**, *Journal of Financial Economics*, 154, 103805. [https://doi.org/10.1016/j.jfineco.2024.103805](https://doi.org/10.1016/j.jfineco.2024.103805)

By using the code, you agree to the following rules:

- You must cite the paper in working papers and published papers that use the code.
- You must place the [DOI](https://doi.org/10.5281/zenodo.10459612) of this data/code repository in a footnote to help others find it.
- You assume all risk for the use of the code.

All datasets are proprietary. We do not have the rights to share any of the data. We provide pseudo data to illustrate the code usage.

The computer code is written in R. We provide the code to reproduce the tables and figures in the paper.

## Data Availability and Provenance Statements

All datasets are proprietary. Factor data are available from other researchers' websites, while hedge fund data and mutual fund data come from commercial databases. See below for details.

### Data Sources for the Factors

Please refer to Section 4.2 in the paper and Section III.B in the Supplementary Materials.

### Data Sources for the Hedge Funds

Please refer to Section 4.1 in the paper and Section III.A in the Supplementary Materials.

### Data Sources for the Mutual Funds

Please refer to Section III.A in Barras et al. (2022) and Section V of its Supplementary Material.

### Summary of Availability

**No data can be made** publicly available. Data used in this paper and not provided as part of the public replication package will be preserved for one year after publication. 

## Datasets

We provide two pseudo datasets for running the code: `Db_Factors_Pseudo.rda` and `Db_HF_Pseudo.rda` located in the folder `data`.

| Data File | Description | 
|--------|----------------|
| `Db_Factors_Pseudo.rda` | Matrix `Factors` of size 324 x 28 that contains monthly returns from 1994-01 to 2020-12 of 28 fake pseudo factors |
| | `FF_EMKT` market factor |
| | `FF_HML` value factor |
| | `FF_SMB` size factor |
| | `FF_MOM` momentum factor |
| | `FF_CMA` investment factor |
| | `FF_RMW` profitability factor |
| | `FH_TERM_FH` term factor |
| | `FH_DFLT_FH` default factor |
| | `FH_PTFSBD` straddle on bonds factor |
| | `FH_PTFSCOM` straddle on commodities factor |
| | `FH_PTFSFX` straddle on currencies factor |
| | `AQR_VME_VAL` global value factor |
| | `AQR_VME_MOM` global momentum factor |
| | `AQR_TSMOM` time-series momentum factor |
| | `AQR_BAB_USA` BAB factor |
| | `BH_VARVIX` variance factor |
| | `KMPV_GCF` carry factor |
| | `PS_LIQ` illiquidity factor |
| | `KNS_dindrrevlv`, `KNS_dindmomrev`, `KNS_dindrrev`, `KNS_dseason`, and `KNS_dsue` machine learning portfolios factors |
| | `KNS_pc_d5`, `KNS_pc_d1`, `KNS_pc_d6`, `KNS_pc_d12`, and `KNS_pc_d2` machine learning principal components factors |
| `Db_HF_Pseudo.rda` | List `Db_HF` containing information for 2000 fake pseudo hedge funds |
|  | `id` vector of size 2000 of fund id |
|  | `dates` vector of size 324 of dates |
|  | `ret` matrix of size 324 x 2000 of monthly returns [%] |
|  | `aum` matrix of size 324 x 2000 of monthly aums [usd] |
|  | `strat` vector of size 2000 of fund main strategy |
|  | `substrat` vector of size 2000 of fund sub-strategy |
|  | `mf` vector of size 2000 of fund management fees [%] |
|  | `pf` vector of size 2000 of fund performance fees [%] |
|  | `hwm` vector of size 2000 of fund highwater mark [true/false] |
|  | `hr` vector of size 2000 of fund hurdle rate [true/false] |
|  | `np` vector of size 2000 of fund notice period [month] |
|  | `lp` vector of size 2000 of fund lookup period [year] |

## Computational Requirements

You must have R installed and a C++ compiler (such as GCC) to run Rcpp/RcppArmadillo.

### Software Requirements

The file `run_install_packages.R` will install all dependencies (latest version) and should be run once before running other programs. See the file `session_info.txt` in the folder `outputs` to see the exact setup that generated the results in the paper.

## Description of Programs/Code

- The file `run_install_packages.R` will install all dependencies (latest version) and should be run once before running other programs. 
- The file `run_results.R` will generate all tables and figures in the paper using the pseudo data sets. Be careful, as it will overwrite the content of the `outputs` folder.

### License for Code

The code is under the GPL-3 license; see the file `LICENSE.txt`. Moreover, if you use part of the computer code in your research, you must cite Ardia et al. (202x) and add a footnote pointing to the code/data repository. 

## Instructions to Replicators

- Run `run_install_packages.R` to install the missing packages. 
- Run `run_results.R` to generate all tables (except Tables 8 and 9) and all figures in the papers. They will be saved in the folder `outputs`.

## List of Tables/Figures and Programs

The provided code reproduces all tables and figures in the paper.

| Figure/Table #    | Line Number | Output file   | Note          |
|-------------------|-------------|---------------|---------------|
| Table 1           |    27       | `Table1.txt`    |               |      
| Table 2           |    86       | `Table2.txt`    |               |    
| Table 4           |   153       | `Table4r2.txt`  | R2 in Table 4 |  
| Table 6           |   198       | `Table6a.txt`, `Table6b.txt`  | |  
| Table 5           |   267       | `Table5.txt`  | |  
| Table 4           |   309       | `Table4.txt`  | |  
| Table 7a          |   336       | `Table7a.txt`  | |  
| Table 7b          |   371       | `Table7b.txt`  | | 
| Figure 3a-4        |   385       | `Figure3a.pdf`, `Figure4.pdf`  | | 
| Table 11        |   505       | `Table11a.txt`, `Table11b.txt`, `Table11c.txt`  | | 
| Figure 2        |   628       | `Figure2a.pdf`,`Figure2b.pdf`  | | 
| Table 10        |   683       | `Table10.txt`  | Generate results for management fees. To generate each panel, uncomment other parts of the code | 
| Figure 1        |   739       | `Figure1.pdf` | | 

Table 3 does not require any code to be run. Tables 8 and 9 are category-specific cases of Table 5 and Table 7. Please proceed as follows to generate them.

| Table #    | Note        |
|-------------------|-------------|
| Table 8           |  Uncomment line 248 and run the code for Table 5 above |       
| Table 9           |  Uncomment line 248 and run the code for Table 7 above | 

## References

Ardia D., Barras L., Gagilardini P., and O. Scaillet. 2024. Is it alpha or beta? Decomposing hedge fund returns when models are misspecified. Journal of Financial Economics 154:103805. https://doi.org/10.1016/j.jfineco.2024.103805

Barras L., Gagilardini P., and O. Scaillet. 2022. Skill, scale, and value creation in the mutual fund industry. Journal of Finance 77:601-638. https://doi.org/10.1111/jofi.13096

Carhart, M. 1997. On persistence in mutual fund performance. Journal of Finance 52:57-82. https://doi.org/10.1111/j.1540-6261.1997.tb03808.x

Fama, E. F., and K. R. French. 2015. A five-factor asset pricing model. Journal of Financial Economics 116:1-22. https://doi.org/10.1016/j.jfineco.2014.10.010

Frazzini, A., and L. H. Pedersen. 2014. Betting against beta. Journal of Financial Economics 111:1-25. https://doi.org/10.1016/j.jfineco.2013.10.005

Fung, W., and D. A. Hsieh. 2004. Hedge fund benchmarks: A risk-based approach. Financial Analysts Journal 60:65-80. https://doi.org/10.2469/faj.v60.n5.2657

Koijen, R. S. J., T. J. Moskowitz, L. H. Pedersen, and E. B. Vrugt. 2018. Carry. Journal of Financial Economics 127:197-225. https://doi.org/10.1016/j.jfineco.2017.11.002

Kozak, S., S. Nagel, and S. Santosh. 2020. Shrinking the cross-section. Journal of Financial Economics 135:271-92. https://doi.org/10.1016/j.jfineco.2019.06.008

Moskowitz, T. J., Y. H. Ooi, and L. H. Pedersen. 2012. Time series momentum. Journal of Financial Economics 104:228-50. https://doi.org/10.1016/j.jfineco.2011.11.003

Pastor, L., and R. F. Stambaugh. 2003. Liquidity risk and expected stock returns. Journal of Political Economy 111:642-85. https://doi.org/10.1086/374184

## Acknowledgements

Some content on this page was copied from [Hindawi](https://www.hindawi.com/research.data/#statement.templates). Other content was adapted from [Fort (2016)](https://doi.org/10.1093/restud/rdw057), Supplementary data, with the author's permission.

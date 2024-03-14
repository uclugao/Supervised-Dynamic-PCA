# Supervised-Dynamic-PCA

This repository consists of the code to reproduce all the results in the Supervised Dynamic PCA method of Gao and Tsay (2024). See https://arxiv.org/abs/2307.07689.

We consider the macroeconomic variables studied by Stock and Watson (2002b), McCracken and Ng (2016), and Huang, Jiang, Li, Tong, and Zhou (2022), among many others. The data are obtained from the FRED-MD data base which are maintained by St. Louis Fed (https://research.stlouisfed.org/econ/mccracken/fred-databases/). The detailed variables and transformation codes to ensure the stationarity of each macro variable are provided in the online data appendix. There are 127 variables in the online data set, but 4 of them are removed due to missing values therein. The S&P 500 index volatility and return data are obtained from the online data appendix of Welch and Goyal (2008) (https://sites.google.com/view/agoyal145).

![image](https://github.com/uclugao/Supervised-Dynamic-PCA/assets/77175209/30373c32-1817-40c9-8822-a017ccc1ab96)

Table A.1 can be produced by sdPCA.R and sdPCAR1+PLS.R; Table 1 can be produced by sdPCA-out.R and sdPCAR1+PLS.R; Table A.2 can be produced by sdPCA-out-1-singular.R and PLS-singular.R;Table A.3 can be produced by lasso.R; Table A.4 can be produced by sdPCAR1PLASSO.R;Table A.5 can be produced by selectionQ.R; Figures A.1-A.7 can be produced by the macro-updt.R, R2.R, and R2aic.R; Tables 2-4 and Tables A.6-A7 can be produced by macro-updt.R and macro-q3-9ft.R; Table A.8 can be produced by macroR1PLS.R; Table A.9 can be produced by Lasso-Rescaled-data.R.![image](https://github.com/uclugao/Supervised-Dynamic-PCA/assets/77175209/55def2de-1775-4c37-b6a6-513726141319)

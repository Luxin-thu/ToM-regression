# ToM-regression
Our program:

1. Violin_plot.R produces Figure 1, i,e,  the plot of calibrated weights, leverage scores, and sample influence curve.
2. crt_gen_nolasso.R produces Figure 2-4, Figure S.11-S.17, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in completely randomized experiments, relative bias
3. 3_2stage_gen.R produces Figure 5-7, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in completely randomized survey experiments.
4. real_data_point_error_bar.R produces Figure S.1, i.e., the bar plots of the second real data analysis.
5. 1_crt_smallsample.R produces Figure S.2-S.7, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in completely randomized experiments with a small sample.
6. 1_crt_gen_lasso.R produces Figure S.8-S.10, i.e., compares the performance of Lasso and other methods in completely randomized experiments.
7. 2_stra_gen.R produces Figure S.18-S.26, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in stratified  randomized experiments.
8. 4_cluster_gen.R produces Figure S.27-S.29, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in cluster randomized experiments.
9. 3_Stra_real_data.R and 3_crs_real_data.R are two real data analyses corresponding to Section 6 and Section A (in the Supplementary Material).
10. The same .sh file with the same name as the .R file for parallel computing.
11. The real datasets are included in the folder ``data".
12. The simulation summary data are in the folder ``output"
13. The output figures are in the folder ``figures"


Steps to run Violin_plot.R, real_data_point_error_bar.R, 3_Stra_real_data.R, 3_crs_real_data.R：

Violin_plot.R: Run the Rscript directly.
real_data_point_error_bar.R: Run the Rscript directly.
3_Stra_real_data.R: Run the Rscript directly.
3_crs_real_data.R: Run the Rscript directly.

Steps to run crt_gen_nolasso.R, 3_2stage_gen.R, 1_crt_smallsample.R, 1_crt_gen_lasso.R, 2_stra_gen.R, 4_cluster_gen.R, 3_Stra_real_data.R:

1. Those Rscripts include modules for parallel computing and figure drawing.
2. You need parallel computing  to run those Rscripts.  You can run the .sh file with the same name as the Rscript to replicate our simulation results.
3. The summary data are in the folder of ``output". You can load them on your computer and use the code of the figure drawing module in the Rscript to draw those figures.
4. stra_MS+FL.Rdata and stra_FL.Rdata are summary data of 2_stra_gen.R.
5. crt_n_20.Rdata and crt_n_40.Rdata are summary data of 1_crt_smallsample.R.
6. crs.Rdata is summary data of 3_2stage_gen.R.
7. cluster.Rdata is summary data of 4_cluster_gen.R.
8. nolasso_crt_r_0.5.Rdata, nolasso_crt_r_0.4.Rdata, nolasso_crt_r_0.3.Rdata are summary data of crt_gen_nolasso.R.
9. lasso_crt_1se_r_0.3 is summary data of 1_crt_gen_lasso.R

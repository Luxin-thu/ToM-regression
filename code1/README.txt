Our program:

1. Violin_plot.R produces Figure 1, i,e,  the plot of calibrated weights, leverage scores and sample influence curve.
2. crt_gen_nolasso.R produces Figure 2-4, Figure S.11-S.17, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in completely randomized experiments, relative bias
3. 3_2stage_gen.R produces Figure 5-7, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in completely randomized survey experiments.
4. real_data_point_error_bar.R produces Figure S.1, i.e., the barplots of the second real data analysis.
5. 1_crt_smallsample.R produces Figure S.2-S.7, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in completely randomized experiments with small sample.
6. 1_crt_gen_lasso.R produces Figure S.8-S.10, i.e., compare the performance of Lasso and other methods in completely randomized experiments .
7. 2_stra_gen.R produces Figure S.18-S.26, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in stratified  randomized experiments.
8. 4_cluster_gen.R produces Figure S.27-S.29, i.e., the plots of relative RMSE, relative confidence interval length, coverage probabilities
in cluster randomized experiments.
9. 3_Stra_real_data.R and 3_crs_real_data.R are two real data analysis corresponding to Section 6 and Section A (in the Supplementary Material).
10. The same .sh file with the same name as the .R file for for parallel computing.
11. The real datasets are included in the folder ``data".


Steps to run Violin_plot.R, real_data_point_error_bar.R, 3_Stra_real_data.R, 3_crs_real_data.R：
Violin_plot.R: Run the Rscript directly.
real_data_point_error_bar.R: Run the Rscript directly.
3_Stra_real_data.R: Run the Rscript directly.
3_crs_real_data.R: Run the Rscript directly.

Steps to run crt_gen_nolasso.R, 3_2stage_gen.R, 1_crt_smallsample.R, 1_crt_gen_lasso.R, 2_stra_gen.R, 4_cluster_gen.R, 3_Stra_real_data.R:

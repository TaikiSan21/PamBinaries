#ifndef _RS_INTERFACE_RSMutils_
#define _RS_INTERFACE_RSMutils_

#include "ut_RS.h"

// RS_fra_dim.c
EXTERN_R SEXP RS_fractal_embed(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_dimension_information(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_dimension_correlation_summation(
 SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_dimension_false_nearest_neighbors(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_dimension_false_nearest_strands(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_poincare_map(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_space_time_separation_plot(
 SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra_fdp.c
EXTERN_R SEXP RS_wavelets_fdp_simulate_weights(SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_fdp_simulate(SEXP, SEXP);

// RS_fra_filt.c
EXTERN_R SEXP RS_fractal_filter_median(SEXP, SEXP);
EXTERN_R SEXP RS_fractal_filter_nonlinear_local_projection(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra_kde.c
EXTERN_R SEXP RS_fractal_kernel_density_estimate(SEXP, SEXP);

// RS_fra_lyap.c
EXTERN_R SEXP RS_fractal_local_lyapunov_spectrum(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra_modl.c
EXTERN_R SEXP RS_fractal_determinism_delta_epsilon(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra.neig.c
EXTERN_R SEXP RS_fractal_neighbor_find(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra_scale.c
EXTERN_R SEXP RS_fractal_piecwise_linear_segmentation(
 SEXP, SEXP, SEXP, SEXP);

// RS_fra_sdf.c
EXTERN_R SEXP RS_fractal_spectral_density_function_direct(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_spectral_density_function_lag_window(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_spectral_density_function_wosa(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_spectral_density_function_multitaper(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra_stat.c
EXTERN_R SEXP RS_fractal_stationarity_priestley_subba_rao(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_fra_surr.c
EXTERN_R SEXP RS_fractal_bootstrap_theiler(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_fractal_bootstrap_davison_hinkley(SEXP, SEXP);
EXTERN_R SEXP RS_fractal_bootstrap_circulant_embedding(SEXP, SEXP);

// RS_fra_tdmi.c
EXTERN_R SEXP RS_fractal_time_delayed_mutual_information(SEXP, SEXP);

// wmtsa\src\RS_mth_var.c or sapa\src\RS_math_acvs.c
EXTERN_R SEXP RS_math_acvs(SEXP, SEXP, SEXP);

// RS_sig_win.c
EXTERN_R SEXP RS_signal_taper(SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_wav_boot.c
EXTERN_R SEXP RS_wavelets_transform_packet_whitest(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_bootstrap(SEXP, SEXP, SEXP, SEXP);

// RS_wav_fdp.c
EXTERN_R SEXP RS_wavelets_fdp_estimator_instantaneous(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_fdp_estimator_block(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_fdp_bandpass_variance(SEXP, SEXP, SEXP);

// RS_wav_filt.c
EXTERN_R SEXP RS_wavelets_filters_daubechies(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_filters_daubechies_gain(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_filter_zero_phase(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_filters_continuous(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_coefficient_boundaries(
 SEXP, SEXP, SEXP, SEXP);

// RS_wav_mrd.c
EXTERN_R SEXP RS_wavelets_transform_packet_detail(
 SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_wav_shrk.c
EXTERN_R SEXP RS_wavelets_shrink(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_shrink_threshold(
 SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_wav_var.c
EXTERN_R SEXP RS_wavelets_variance(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_variance_confidence(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_variance_edof(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_wav_wtmm.c
EXTERN_R SEXP RS_wavelets_transform_continuous_wavelet_modulus_maxima(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_continuous_wavelet_modulus_maxima_tree(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// RS_wav_xform.c
EXTERN_R SEXP RS_wavelets_transform_continuous_wavelet(
 SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_discrete_wavelet_convolution(
 SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_packet(
 SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_discrete_wavelet_convolution_inverse(
 SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_packet_convert_indices(SEXP);
EXTERN_R SEXP RS_wavelets_transform_packet_basis(SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_packet_inverse(
 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_maximum_overlap(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_maximum_overlap_packet(SEXP, SEXP, SEXP);
EXTERN_R SEXP RS_wavelets_transform_maximum_overlap_inverse(SEXP, SEXP);


#endif

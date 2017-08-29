#include <R.h>

#ifdef USING_R

#include <R_ext/Rdynload.h>

#include "RSInterface.h"

void R_init_rsmutils(DllInfo *dll)
{
	static const R_CallMethodDef CallEntries[] = 
	{
		/**********************
			For package SAPA
		***********************/
		// RS_fra_sdf.c
		{"RS_fractal_spectral_density_function_direct", 
			(DL_FUNC) &RS_fractal_spectral_density_function_direct, 6},
		{"RS_fractal_spectral_density_function_lag_window", 
			(DL_FUNC) &RS_fractal_spectral_density_function_lag_window, 7},
		{"RS_fractal_spectral_density_function_wosa", 
			(DL_FUNC) &RS_fractal_spectral_density_function_wosa, 7},
		{"RS_fractal_spectral_density_function_multitaper", 
			(DL_FUNC) &RS_fractal_spectral_density_function_multitaper, 6},
		// RS_math_acvs.c	
		{"RS_math_acvs", (DL_FUNC) &RS_math_acvs, 3},
		// RS_sig_win.c    
		{"RS_signal_taper", (DL_FUNC) &RS_signal_taper, 5},
		
		/*************************
			For package wmtsa
		**************************/
		
		// RS_wav_boot.c
		{"RS_wavelets_transform_packet_whitest", 
			(DL_FUNC) &RS_wavelets_transform_packet_whitest, 3},
		{"RS_wavelets_bootstrap", (DL_FUNC) &RS_wavelets_bootstrap, 4},
		// RS_wav_fdp.c
		{"RS_wavelets_fdp_estimator_instantaneous", 
			(DL_FUNC) &RS_wavelets_fdp_estimator_instantaneous, 8},
		{"RS_wavelets_fdp_estimator_block", 
			(DL_FUNC) &RS_wavelets_fdp_estimator_block, 9},
		{"RS_wavelets_fdp_bandpass_variance", 
			(DL_FUNC) &RS_wavelets_fdp_bandpass_variance, 3},
		// RS_wav_filt.c
		{"RS_wavelets_filters_daubechies", 
			(DL_FUNC) &RS_wavelets_filters_daubechies, 3},
		{"RS_wavelets_filters_daubechies_gain", 
			(DL_FUNC) &RS_wavelets_filters_daubechies_gain, 6},
		{"RS_wavelets_filter_zero_phase", 
			(DL_FUNC) &RS_wavelets_filter_zero_phase, 3},
		{"RS_wavelets_filters_continuous", 
			(DL_FUNC) &RS_wavelets_filters_continuous, 3},
		{"RS_wavelets_transform_coefficient_boundaries", 
			(DL_FUNC) &RS_wavelets_transform_coefficient_boundaries, 4},
		// RS_wav_mrd.c
		{"RS_wavelets_transform_packet_detail", 
			(DL_FUNC) &RS_wavelets_transform_packet_detail, 5},
		// RS_wav_shrk.c
		{"RS_wavelets_shrink", 
			(DL_FUNC) &RS_wavelets_shrink, 9},
		{"RS_wavelets_shrink_threshold", 
			(DL_FUNC) &RS_wavelets_shrink_threshold, 5},
		// RS_wav_var.c
		{"RS_wavelets_variance", (DL_FUNC) &RS_wavelets_variance, 6},
		{"RS_wavelets_variance_confidence",
			(DL_FUNC) &RS_wavelets_variance_confidence, 3},
		{"RS_wavelets_variance_edof", (DL_FUNC) &RS_wavelets_variance_edof, 7},
		// RS_wav_xform.c
		{"RS_wavelets_transform_continuous_wavelet", 
			(DL_FUNC) &RS_wavelets_transform_continuous_wavelet, 5},
		{"RS_wavelets_transform_discrete_wavelet_convolution", 
			(DL_FUNC) &RS_wavelets_transform_discrete_wavelet_convolution, 3},
		{"RS_wavelets_transform_packet", 
			(DL_FUNC) &RS_wavelets_transform_packet, 3},
		{"RS_wavelets_transform_discrete_wavelet_convolution_inverse", 
			(DL_FUNC) &RS_wavelets_transform_discrete_wavelet_convolution_inverse, 2},
		{"RS_wavelets_transform_packet_convert_indices", 
			(DL_FUNC) &RS_wavelets_transform_packet_convert_indices, 1},
		{"RS_wavelets_transform_packet_basis", 
			(DL_FUNC) &RS_wavelets_transform_packet_basis, 2},
		{"RS_wavelets_transform_packet_inverse", 
			(DL_FUNC) &RS_wavelets_transform_packet_inverse, 6},
		{"RS_wavelets_transform_maximum_overlap", 
			(DL_FUNC) &RS_wavelets_transform_maximum_overlap, 3},
		{"RS_wavelets_transform_maximum_overlap_packet", 
			(DL_FUNC) &RS_wavelets_transform_maximum_overlap_packet, 3},
		{"RS_wavelets_transform_maximum_overlap_inverse", 
			(DL_FUNC) &RS_wavelets_transform_maximum_overlap_inverse, 2},
			
		/***********************
			For package fractal
		***********************/
		
		// RS_fra_dim.c
		{"RS_fractal_embed", 
			(DL_FUNC) &RS_fractal_embed, 3},
		{"RS_fractal_dimension_information", 
			(DL_FUNC) &RS_fractal_dimension_information, 8},
		{"RS_fractal_dimension_correlation_summation", 
			(DL_FUNC) &RS_fractal_dimension_correlation_summation, 5},
		{"RS_fractal_dimension_false_nearest_neighbors", 
			(DL_FUNC) &RS_fractal_dimension_false_nearest_neighbors, 6},
		{"RS_fractal_dimension_false_nearest_strands", 
			(DL_FUNC) &RS_fractal_dimension_false_nearest_strands, 6},
		{"RS_fractal_poincare_map", 
			(DL_FUNC) &RS_fractal_poincare_map, 3},
		{"RS_fractal_space_time_separation_plot", 
			(DL_FUNC) &RS_fractal_space_time_separation_plot, 5},
		// RS_fra_fdp.c
		{"RS_wavelets_fdp_simulate_weights", 
			(DL_FUNC) &RS_wavelets_fdp_simulate_weights, 2},
		{"RS_wavelets_fdp_simulate", 
			(DL_FUNC) &RS_wavelets_fdp_simulate, 2},
		// RS_fra_filt.c
		{"RS_fractal_filter_median", 
			(DL_FUNC) &RS_fractal_filter_median, 2},
		{"RS_fractal_filter_nonlinear_local_projection", 
			(DL_FUNC) &RS_fractal_filter_nonlinear_local_projection, 8},
		// RS_fra_kde.c
		{"RS_fractal_kernel_density_estimate", 
			(DL_FUNC) &RS_fractal_kernel_density_estimate, 2},
		// RS_fra_lyap.c
		{"RS_fractal_local_lyapunov_spectrum", 
			(DL_FUNC) &RS_fractal_local_lyapunov_spectrum, 10},
		// RS_fra_modl.c
		{"RS_fractal_determinism_delta_epsilon", 
			(DL_FUNC) &RS_fractal_determinism_delta_epsilon, 8},
		// RS_fra.neig.c
		{"RS_fractal_neighbor_find", 
			(DL_FUNC) &RS_fractal_neighbor_find, 6},
		// RS_fra_scale.c
		{"RS_fractal_piecwise_linear_segmentation", 
			(DL_FUNC) &RS_fractal_piecwise_linear_segmentation, 4},
		// RS_fra_stat.c
		{"RS_fractal_stationarity_priestley_subba_rao", 
			(DL_FUNC) &RS_fractal_stationarity_priestley_subba_rao, 7},
		// RS_fra_surr.c
		{"RS_fractal_bootstrap_theiler", 
			(DL_FUNC) &RS_fractal_bootstrap_theiler, 3},
		{"RS_fractal_bootstrap_davison_hinkley", 
			(DL_FUNC) &RS_fractal_bootstrap_davison_hinkley, 2},
		{"RS_fractal_bootstrap_circulant_embedding", 
			(DL_FUNC) &RS_fractal_bootstrap_circulant_embedding, 2},
		// RS_fra_tdmi.c
		{"RS_fractal_time_delayed_mutual_information", 
			(DL_FUNC) &RS_fractal_time_delayed_mutual_information, 2},
		// RS_wav_wtmm.c
		{"RS_wavelets_transform_continuous_wavelet_modulus_maxima", 
			(DL_FUNC) &RS_wavelets_transform_continuous_wavelet_modulus_maxima, 5},
		{"RS_wavelets_transform_continuous_wavelet_modulus_maxima_tree", 
			(DL_FUNC) &RS_wavelets_transform_continuous_wavelet_modulus_maxima_tree, 5}
	};
    	
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}

#endif


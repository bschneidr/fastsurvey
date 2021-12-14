//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat arma_onestage(arma::mat Y,
                        arma::colvec samp_unit_ids,
                        arma::colvec strata_ids,
                        arma::colvec strata_samp_sizes,
                        arma::colvec strata_pop_sizes,
                        Rcpp::CharacterVector singleton_method,
                        Rcpp::LogicalVector use_singleton_method_for_domains,
                        int stage) {
  
  // Determine dimensions of result
  size_t n_col_y = Y.n_cols;
  arma::mat result(n_col_y, n_col_y, arma::fill::zeros);
  
  // Get distinct strata ids and their length, H
  arma::colvec distinct_strata_ids = unique(strata_ids);
  arma::uword H = distinct_strata_ids.n_elem;
  
  // Check for singleton strata
  bool any_singleton_strata = min(strata_samp_sizes) < 2;
  arma::uword n_singleton_strata = 0;
  
  // If `singleton_method = "adjust", get mean of all sampling units
  arma::rowvec Y_means;
  if (any_singleton_strata | use_singleton_method_for_domains[0]) {
    if (singleton_method[0] == "adjust") {
      
      // Get number of distinct sampling units
      int n = 0;
      for (arma::uword h = 0; h < H; ++h) {
        arma::uvec h_indices = arma::find(strata_ids==distinct_strata_ids[h]);
        arma::uvec first_h_index = h_indices.head(1);
        n += min(strata_samp_sizes.elem(first_h_index));
      }
      // Calculate average across all sampling units of all strata
      Y_means = sum(Y, 0) / n;
    }
  }
  
  // Get information from each stratum:
  // - Number of sampling units, in sample and population
  // - Sampling fraction, if applicable
  // - Contribution to sampling variance
  for (arma::uword h = 0; h < H; ++h) {
    
    // Determine which rows of data correspond to the current stratum
    arma::uvec h_indices = arma::find(strata_ids==distinct_strata_ids[h]);
    arma::uvec first_h_index = h_indices.head(1);
    int h_num_observations = h_indices.n_elem;
    
    // Get counts of sampling units in stratum, and corresponding sampling rate
    int n_h = min(strata_samp_sizes.elem(first_h_index));
    double N_h = static_cast<double>(min(strata_pop_sizes.elem(first_h_index)));
    double f_h;
    if (arma::is_finite(N_h)) {
      f_h = static_cast<double>(n_h) /  N_h;
    } else {
      f_h = 0.0;
    }
    
    // Determine whether there's clustering
    bool h_has_clusters = false;
    arma::colvec h_distinct_samp_unit_ids = unique(samp_unit_ids.elem(h_indices));
    int h_num_distinct_samp_unit_ids = h_distinct_samp_unit_ids.n_elem;
    h_has_clusters = h_num_observations != h_num_distinct_samp_unit_ids;
    
    // Increment count of singleton strata
    // and determine denominator to use for
    // estimating variance of PSU totals
    arma::uword df;
    if (n_h < 2) {
      n_singleton_strata += 1;
      df = 1;
      if (singleton_method[0] == "fail") {
        Rcpp::String error_msg("At least one stratum contains only one PSU at stage ");
        error_msg += stage;
        Rcpp::stop(error_msg);
      }
    } else {
      df = n_h - 1;
      if (use_singleton_method_for_domains[0] & h_num_observations == 1) {
        n_singleton_strata += 1;
        any_singleton_strata = TRUE;
      }
    }
    
    if ((n_h > 1) | (singleton_method[0] == "adjust")) {
      // Subset variables of interest to stratum
      // and calculate means for stratum
      arma::mat Y_h = Y.rows(h_indices);
      arma::rowvec mean_Yhi = arma::sum(Y_h, 0) / n_h;
      
      // Initialize variance-covariance of PSU totals
      arma::mat cov_mat(n_col_y, n_col_y, arma::fill::zeros);
      
      // If there's no clustering, simply get each row's contribution
      // to stratum's variance-covariance matrix
      if (!h_has_clusters) {
        
        for (int i=0; i < h_num_observations; ++i ) {
          arma::rowvec Yhi = Y_h.row(i);
          if (n_h > 1) {
            Yhi.each_row() -= mean_Yhi;
          } else {
            Yhi.each_row() -= Y_means;
          }
          
          cov_mat += (arma::trans(Yhi)*Yhi);
        }
      }
      
      // If there is clustering, sum up values in each cluster
      // and add this sum's contribution to the between sampling-units
      // variance-covariance matrix
      if (h_has_clusters) {
        
        for (arma::uword i=0; i < h_distinct_samp_unit_ids.n_elem; ++i ) {
          arma::uvec unit_indices = arma::find(samp_unit_ids.elem(h_indices) == h_distinct_samp_unit_ids[i]);
          arma::rowvec Yhi = sum(Y_h.rows(unit_indices), 0);
          if (n_h > 1) {
            Yhi.each_row() -= mean_Yhi;
          } else {
            Yhi.each_row() -= Y_means;
          }
          
          cov_mat += (arma::trans(Yhi)*Yhi);
        }
      }
      
      // If the data were subsetted, some sampling units
      // may not have rows of data that appear in inputs.
      // Make sure these units contribute to the variance.
      int n_h_missing = n_h - h_num_distinct_samp_unit_ids;
      if (n_h_missing > 0) {
        cov_mat += n_h_missing*(arma::trans(mean_Yhi)*mean_Yhi);
      }
      
      cov_mat = cov_mat / df;
      
      // Add variance contribution
      result += ((1.0 - f_h) * n_h) * cov_mat;
    }
  }
  
  if (any_singleton_strata & (singleton_method[0] == "average")) {
    int n_nonsingleton_strata = H - n_singleton_strata;
    double scaling_factor;
    if (n_nonsingleton_strata > 0) {
      scaling_factor = static_cast<double>(H)/static_cast<double>(n_nonsingleton_strata);
    } else {
      scaling_factor = 1;
    }
    result *= scaling_factor;
  }
  
  return result;
}

// [[Rcpp::export]]
arma::mat arma_multistage(arma::mat Y,
                          arma::mat samp_unit_ids,
                          arma::mat strata_ids,
                          arma::mat strata_samp_sizes,
                          arma::mat strata_pop_sizes,
                          Rcpp::CharacterVector singleton_method,
                          Rcpp::LogicalVector use_singleton_method_for_domains,
                          Rcpp::LogicalVector use_only_first_stage,
                          int stage) {
  
  size_t n_stages = samp_unit_ids.n_cols;
  
  // If there are later stages of sampling,
  // obtain the necessary columns from inputs,
  // which will be used recursively
  
  arma::mat later_stage_ids;
  arma::mat later_stage_strata;
  arma::mat later_stage_strata_samp_sizes;
  arma::mat later_stage_strata_pop_sizes;
  
  if ((n_stages > 1) & !use_only_first_stage[0]) {
    later_stage_ids = samp_unit_ids.tail_cols(n_stages - 1);
    later_stage_strata = strata_ids.tail_cols(n_stages - 1);
    later_stage_strata_samp_sizes = strata_samp_sizes.tail_cols(n_stages-1);
    later_stage_strata_pop_sizes = strata_pop_sizes.tail_cols(n_stages-1);
  }
  
  // Obtain first stage information
  arma::colvec first_stage_ids = samp_unit_ids.col(0);
  arma::colvec first_stage_strata = strata_ids.col(0);
  arma::colvec first_stage_strata_samp_sizes = strata_samp_sizes.col(0);
  arma::colvec first_stage_strata_pop_sizes = strata_pop_sizes.col(0);
  
  // Calculate first-stage variance
  arma::mat V = arma_onestage(Y = Y,
                              samp_unit_ids = first_stage_ids,
                              strata_ids = first_stage_strata,
                              strata_samp_sizes = first_stage_strata_samp_sizes,
                              strata_pop_sizes = first_stage_strata_pop_sizes,
                              singleton_method = singleton_method,
                              use_singleton_method_for_domains = use_singleton_method_for_domains,
                              stage = stage);
  
  // For each first-stage unit, get variance contribution from next stage
  if ((n_stages > 1) & !use_only_first_stage[0]) {
    
    // Get distinct first-stage strata ids and their length, H
    arma::colvec distinct_strata_ids = unique(first_stage_strata);
    arma::uword H = distinct_strata_ids.n_elem;
    
    for (arma::uword h = 0; h < H; ++h) {
      
      // Determine which rows of data correspond to the current first-stage stratum
      arma::uvec h_indices = arma::find(first_stage_strata==distinct_strata_ids(h));
      
      // Get submatrices of inputs corresponding to the current first-stage stratum
      arma::mat Y_h = Y.rows(h_indices);
      
      arma::mat h_samp_unit_ids = later_stage_ids.rows(h_indices);
      
      arma::mat h_strata = later_stage_strata.rows(h_indices);
      
      arma::mat h_strata_samp_sizes = later_stage_strata_samp_sizes.rows(h_indices);
      arma::mat h_strata_pop_sizes = later_stage_strata_pop_sizes.rows(h_indices);
      
      // Get count of first-stage sampling units in first-stage stratum
      // and finite population correction, based on all first-stage units in sample design
      arma::uword n_h = min(first_stage_strata_samp_sizes.elem(h_indices));
      double N_h = static_cast<double>(min(strata_pop_sizes.elem(h_indices)));
      
      double f_h;
      if (arma::is_finite(N_h)) {
        f_h = static_cast<double>(n_h) /  N_h;
      } else {
        f_h = 0.0;
      }
      
      // Get list of first-stage untis in the current subset of data, and count them
      arma::colvec h_first_stage_units = first_stage_ids.elem(h_indices);
      arma::colvec h_unique_first_stage_units = unique(h_first_stage_units);
      arma::uword n_h_subset = h_unique_first_stage_units.n_elem;
      
      for (arma::uword i=0; i < n_h_subset; ++i ) {
        // Create subsets of inputs specific to current first-stage sampling unit
        arma::uvec unit_indices = arma::find(first_stage_ids.elem(h_indices) == h_unique_first_stage_units(i));
        arma::mat Y_hi = Y_h.rows(unit_indices);
        arma::mat hi_samp_unit_ids = h_samp_unit_ids.rows(unit_indices);
        arma::mat hi_strata = h_strata.rows(unit_indices);
        arma::mat hi_strata_samp_sizes = h_strata_samp_sizes.rows(unit_indices);
        arma::mat hi_strata_pop_sizes = h_strata_pop_sizes.rows(unit_indices);
        
        // Estimate later-stage variance contribution
        arma::mat V_hi = f_h * arma_multistage(Y_hi,
                                               hi_samp_unit_ids,
                                               hi_strata,
                                               hi_strata_samp_sizes,
                                               hi_strata_pop_sizes,
                                               singleton_method,
                                               use_singleton_method_for_domains,
                                               use_only_first_stage,
                                               stage = stage + 1);
        V += V_hi;
        
      }
    }
  }
  return V;
}

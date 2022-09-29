//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(rng=FALSE)]]
arma::mat arma_onestage(arma::mat Y,
                        arma::colvec samp_unit_ids,
                        arma::colvec strata_ids,
                        arma::colvec strata_samp_sizes,
                        arma::colvec strata_pop_sizes,
                        Rcpp::CharacterVector singleton_method,
                        Rcpp::LogicalVector use_singleton_method_for_domains,
                        int stage) { 
  
  arma::uword number_of_data_rows = samp_unit_ids.n_elem;
  arma::uvec row_indices = arma::linspace<arma::uvec>(0L, number_of_data_rows - 1L, number_of_data_rows);
  
  // Determine dimensions of result
  size_t n_col_y = Y.n_cols;
  arma::mat result(n_col_y, n_col_y, arma::fill::zeros);
  
  // Check for singleton strata
  bool any_singleton_strata = min(strata_samp_sizes) < 2;
  
  // If `singleton_method = "adjust", get mean of all sampling units
  arma::rowvec Y_means;
  Y_means = Y_means.zeros(n_col_y);
  if (any_singleton_strata | use_singleton_method_for_domains[0]) {
    if (singleton_method[0] == "adjust") {
      
      int n = 0;
      bool at_end_of_stratum;
      arma::uword next_row_index = 0;
      for (arma::uvec::iterator row_index = row_indices.begin(); row_index != row_indices.end(); ++row_index) {
        Y_means += Y.row(*row_index);
        if (next_row_index < number_of_data_rows) {
          at_end_of_stratum = strata_ids(*row_index) != strata_ids(next_row_index);
        } else {
          at_end_of_stratum = true;
        }
        if (at_end_of_stratum) {
          n += strata_samp_sizes(*row_index);
        }
      }
      // Calculate average across all sampling units of all strata
      Y_means = Y_means / n;
    }
  }
  
  // Initialize count of singleton strata
  int n_singleton_strata = 0;
  
  // Initialize stratum-level summaries
  int H = 0;
  arma::rowvec Ybar_h;
  Ybar_h = Ybar_h.zeros(n_col_y);
  double n_h_in_data = 0;
  double scale;
  arma::mat h_sum_of_squares(n_col_y, n_col_y, arma::fill::zeros);
  arma::mat cov_h(n_col_y, n_col_y, arma::fill::zeros);
  
  // Initialize sample unit total
  arma::rowvec Yhi;
  Yhi = Yhi.zeros(n_col_y);
  
  // Initialize checks for end of stratum or end of samp unit
  
  bool at_end_of_stratum = true;
  bool at_end_of_samp_unit = true;
  
  // Iterate over each row and its following row,
  // in the process iterating over strata and sampling units within strata
  arma::uword next_row_index = 0;
  
  for (arma::uvec::iterator row_index = row_indices.begin(); row_index != row_indices.end(); ++row_index) {
    
    // Determine whether the current row is the last observation in a stratum or sampling unit
    next_row_index = (*(row_index+1));
    
    if ((*(row_index)) == (number_of_data_rows - 1)) {
      at_end_of_stratum = true;
      at_end_of_samp_unit = true;
    } else {
      at_end_of_stratum = strata_ids(*row_index) != strata_ids(next_row_index);
      if (!at_end_of_stratum) {
        at_end_of_samp_unit = samp_unit_ids(*row_index) != samp_unit_ids(next_row_index);
      } else {
        at_end_of_samp_unit = true;
      }
    }
    
    // Get contribution to sampling unit's total
    Yhi += Y.row(*row_index);
    
    // Add contribution of sampling unit
    // to the stratum's sum of squares
    if (at_end_of_samp_unit) {
      Ybar_h += Yhi;
      n_h_in_data += 1;
      h_sum_of_squares += (arma::trans(Yhi)*Yhi);
      Yhi = Yhi.zeros();
    }
    
    if (at_end_of_stratum) {
      
      H += 1;
      
      // Determine sampling fraction
      double n_h = static_cast<double>(strata_samp_sizes(*row_index));
      double N_h = static_cast<double>(strata_pop_sizes(*row_index));
      double f_h;
      if (arma::is_finite(N_h)) {
        f_h = static_cast<double>(n_h) /  N_h;
      } else {
        f_h = 0.0;
      }
      
      // Increment count of singleton strata
      bool h_is_singleton = false;
      if (n_h == 1) {
        h_is_singleton = true;
        n_singleton_strata += 1;
        if (singleton_method[0] == "fail") {
          Rcpp::String error_msg("At least one stratum contains only one PSU at stage ");
          error_msg += stage;
          Rcpp::stop(error_msg);
        }
      } else if (n_h_in_data == 1) {
        if (use_singleton_method_for_domains[0] & ((singleton_method[0] == "adjust") | (singleton_method[0] == "average"))) {
          h_is_singleton = true;
          n_singleton_strata += 1;
        }
      }
      
      // Determine scaling factor to use for normalizing sum of squares
      // and handling finite population correction
      scale = ((1.0 - f_h) * n_h);
      
      // Determine what to do if the stratum is a singleton (has only one sampling unit)
      if (h_is_singleton) {
        if (singleton_method[0] == "adjust") {
          Ybar_h = Y_means;
        } else {
          // Reset stratum count and sum variables so they can be used for the next stratum,
          // then move on to next stratum and row of data
          Ybar_h = Ybar_h.zeros();
          h_sum_of_squares = h_sum_of_squares.zeros();
          cov_h = cov_h.zeros();
          n_h = 0;
          n_h_in_data = 0;
          continue;
        }
      } else {
        Ybar_h = Ybar_h / n_h;
        scale /= (n_h - 1);
      }
      
      // Get variance-covariance contribution of stratum
      cov_h = scale * (h_sum_of_squares - n_h*(arma::trans(Ybar_h)*Ybar_h));
      result += cov_h;
      
      // Reset stratum count and sum variables so they can be used for the next stratum
      Ybar_h = Ybar_h.zeros();
      h_sum_of_squares = h_sum_of_squares.zeros();
      cov_h = cov_h.zeros();
      n_h = 0;
      n_h_in_data = 0;
    }
    
  }
  
  // If the user specified 'average' method for handling singleton strata,
  // scale the total variance estimate
  any_singleton_strata = n_singleton_strata > 0;
  if ((singleton_method[0] == "average") & any_singleton_strata) {
    int n_nonsingleton_strata = H - n_singleton_strata;
    double scaling_factor;
    if (n_nonsingleton_strata > 0) {
      scaling_factor = static_cast<double>(H)/static_cast<double>(n_nonsingleton_strata);
    } else {
      scaling_factor = R_NaN;
    }
    result *= scaling_factor;
  }
  
  return result;
}

// [[Rcpp::export(rng=FALSE)]]
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
  
  // First reorder inputs by first-stage sample unit IDs
  arma::uvec samp_unit_id_order = arma::stable_sort_index(samp_unit_ids.col(0), "ascend");
  Y = Y.rows(samp_unit_id_order);
  samp_unit_ids = samp_unit_ids.rows(samp_unit_id_order);
  strata_ids = strata_ids.rows(samp_unit_id_order);
  strata_samp_sizes = strata_samp_sizes.rows(samp_unit_id_order);
  strata_pop_sizes = strata_pop_sizes.rows(samp_unit_id_order);
  
  // Next reorder inputs by first-stage strata IDs
  arma::uvec strata_id_order = arma::stable_sort_index(strata_ids.col(0), "ascend");
  Y = Y.rows(strata_id_order);
  samp_unit_ids = samp_unit_ids.rows(strata_id_order);
  strata_ids = strata_ids.rows(strata_id_order);
  strata_samp_sizes = strata_samp_sizes.rows(strata_id_order);
  strata_pop_sizes = strata_pop_sizes.rows(strata_id_order);
  
  // Obtain first stage information
  arma::colvec first_stage_ids = samp_unit_ids.col(0);
  arma::colvec first_stage_strata = strata_ids.col(0);
  arma::colvec first_stage_strata_samp_sizes = strata_samp_sizes.col(0);
  arma::colvec first_stage_strata_pop_sizes = strata_pop_sizes.col(0);
  
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

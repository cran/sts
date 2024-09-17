#include <RcppArmadillo.h>
#include <armadillo>
// [[Rcpp::depends(RcppArmadillo)]]

// Log-Sum-Exp function
//double logSumExp(const arma::vec &x)
//{
//  double max_elem = x.max();
//  return std::log(arma::accu(arma::exp(x - max_elem))) + max_elem;
//}

// Softmax function
//arma::vec softmax(const arma::vec &x)
//{
//  double lse = logSumExp(x);
//  return arma::exp(x - lse);
//}
// Optimized Softmax function
//arma::vec softmax(const arma::vec &x)
//{
//    double max_elem = x.max();
//    arma::vec exp_x = arma::exp(x - max_elem);
//    double sum_exp = arma::accu(exp_x);
//
//    return exp_x / sum_exp;
//}

// [[Rcpp::export]]
double lpbdcpp(const arma::vec& alpha_d,
               const arma::mat& Sigma_Inv,
               const arma::mat& kappa_t,
               const arma::mat& kappa_s,
               const arma::vec& mu_d,
               const arma::umat& docs,
               int V,
               const arma::vec& mv)
{
    int K = (alpha_d.n_elem + 1) / 2;

    // Precompute residual and term1
    arma::vec resid = alpha_d - mu_d;
    double term1 = -0.5 * arma::as_scalar(resid.t() * Sigma_Inv * resid);

    // Convert document matrix to zero-based indices for C++
    const arma::uvec cdv_ind = docs.row(0).t() - 1;
    const arma::vec cdv = arma::conv_to<arma::vec>::from(docs.row(1).t());

    // Precompute etaexp
    arma::vec etaexp = arma::exp(arma::join_vert(alpha_d.subvec(0, K - 2), arma::zeros<arma::vec>(1)));

    // Precompute linpred_base and the softmax calculations
    arma::mat linpred_base = mv * arma::ones<arma::rowvec>(K) + kappa_t;
    arma::mat beta(V, K, arma::fill::zeros);

    for (int k = 0; k < K; ++k) {
        arma::vec linpred = linpred_base.col(k) + kappa_s.col(k) * alpha_d(K + k - 1);
        beta.col(k) = arma::exp(linpred - arma::max(linpred));  // Stabilized softmax computation
        beta.col(k) /= arma::accu(beta.col(k));  // Normalize to get probabilities
    }

    // Compute log_probs and term2 efficiently
    const arma::vec log_probs = arma::log(beta.rows(cdv_ind) * etaexp);
    double term2 = arma::dot(cdv, log_probs) - arma::accu(cdv) * std::log(arma::accu(etaexp));

    return term1 + term2;
}


// [[Rcpp::export]]
arma::vec lgaecpp(const arma::vec& alpha_d,
               const arma::mat& Sigma_Inv,
               const arma::mat& kappa_t,
               const arma::mat& kappa_s,
               const arma::vec& mu_d,
               const arma::umat& docs,
               int V,
               const arma::vec& mv)
{
    int K = (alpha_d.n_elem + 1) / 2;

    // Calculate beta and beta_bar
    arma::vec etaexp = arma::exp(arma::join_vert(alpha_d.subvec(0, K-2), arma::zeros<arma::vec>(1)));
    // Precompute linpred_base and the softmax calculations
    arma::mat linpred_base = mv * arma::ones<arma::rowvec>(K) + kappa_t;
    arma::mat beta(V, K, arma::fill::zeros);
    arma::mat beta_bar(V, K, arma::fill::zeros);

    for (int k = 0; k < K; ++k) {
//        arma::vec linpred = mv + kappa_t.col(k) + kappa_s.col(k) * alpha_d(K + k - 1);
        arma::vec linpred = linpred_base.col(k) + kappa_s.col(k) * alpha_d(K + k - 1);
        beta.col(k) = arma::exp(linpred - arma::max(linpred));  // Stabilized softmax computation
        beta.col(k) /= arma::accu(beta.col(k));  // Normalize to get probabilities
        beta_bar.col(k) = etaexp[k] * beta.col(k);
    }


    arma::vec rs = arma::sum(beta_bar, 1);
    beta_bar.each_col() /= rs;

    // Convert document matrix to zero-based indices for C++
    const arma::uvec cdv_ind = docs.row(0).t() - 1;
    const arma::vec cdv = arma::conv_to<arma::vec>::from(docs.row(1).t());

    // Expand cdv to match dimensions of beta_bar.rows(cdv_ind)
    arma::mat cdv_expanded = arma::repmat(cdv, 1, K);

    // Compute g1 efficiently
    arma::vec term1 = arma::sum(beta_bar.rows(cdv_ind) % cdv_expanded, 0).t();
    double sum_cdv = arma::accu(cdv);
   
    arma::vec term2 = sum_cdv * etaexp / arma::accu(etaexp);
    arma::vec g1 = term1.subvec(0, K-2) - term2.subvec(0, K-2);

    // Compute g2 efficiently
    arma::rowvec kap_beta_sum_vec = arma::sum(kappa_s % beta, 0);
    arma::mat kap_beta_sum = arma::repmat(kap_beta_sum_vec, V, 1);
    arma::mat beta_bar_cdv_mult = cdv_expanded % beta_bar.rows(cdv_ind);
    arma::mat kappa_diff = kappa_s.rows(cdv_ind) - kap_beta_sum.rows(cdv_ind);
    arma::vec g2 = arma::sum(beta_bar_cdv_mult % kappa_diff, 0).t();

    // Compute term3
    arma::vec resid = alpha_d - mu_d;
    arma::vec term3 = Sigma_Inv * resid;

    return arma::join_vert(g1, g2) - term3;
}


// [[Rcpp::export]]
arma::mat esthcpp(const arma::vec& alpha_d, 
                  const arma::mat& kappa_t,
                  const arma::mat& kappa_s,
                  const arma::mat& Sigma_Inv, 
                  const arma::mat& doc, 
                  int V, 
                  const arma::vec& mv) 
{
    int K = (alpha_d.n_elem + 1) / 2;

    arma::mat H(2 * K - 1, 2 * K - 1, arma::fill::zeros);
    
    arma::uvec cdv_ind = arma::conv_to<arma::uvec>::from(doc.row(0) - 1);
    arma::vec cdv = doc.row(1).t();

    arma::vec etaexp = arma::exp(arma::join_cols(alpha_d.subvec(0, K - 2), arma::zeros<arma::vec>(1)));
        
    arma::mat beta(V, K, arma::fill::zeros);
    arma::mat beta_bar(V, K, arma::fill::zeros);
    
    for (int k = 0; k < K; ++k) {
        beta.col(k) = arma::exp(mv + kappa_t.col(k) + kappa_s.col(k) * alpha_d(K + k - 1));
        beta.col(k) /= arma::accu(beta.col(k));  // In-place softmax
        beta_bar.col(k) = etaexp(k) * beta.col(k);
    }
    
    arma::vec rs = arma::sum(beta_bar, 1);
    beta_bar.each_col() /= rs;
    
    arma::vec theta = etaexp / arma::sum(etaexp);
    
    arma::mat kap_beta_sum = arma::sum(kappa_s % beta, 0);
    kap_beta_sum = arma::repmat(kap_beta_sum, V, 1);
    arma::mat kappa_bar = kappa_s - kap_beta_sum;
    
    arma::mat cdv_rep = arma::repmat(cdv, 1, K);
    arma::vec term1 = arma::sum(cdv_rep % beta_bar.rows(cdv_ind), 0).t();
    
    double term2_sum = arma::sum(etaexp);
    arma::vec term2 = arma::sum(cdv) * etaexp / term2_sum;
    
    arma::vec g1 = term1 - term2;
    
    arma::vec g2 = arma::sum(cdv_rep % beta_bar.rows(cdv_ind) % kappa_bar.rows(cdv_ind), 0).t();
    
    arma::mat h_pp = arma::diagmat(g1) - beta_bar.rows(cdv_ind).t() * arma::diagmat(cdv) * beta_bar.rows(cdv_ind) + arma::sum(cdv) * theta * theta.t();
    
    arma::mat h_ps = arma::diagmat(g2) - (beta_bar.rows(cdv_ind) % kappa_bar.rows(cdv_ind)).t() * arma::diagmat(cdv) * beta_bar.rows(cdv_ind);
    
    arma::mat kap_beta_sum_rows = arma::sum(kappa_bar % kappa_s % beta, 0);
    kap_beta_sum_rows = arma::repmat(kap_beta_sum_rows, V, 1);
    
    arma::vec hssterm1 = arma::sum(cdv_rep % beta_bar.rows(cdv_ind) % (arma::pow(kappa_bar.rows(cdv_ind), 2) - kap_beta_sum_rows.rows(cdv_ind)), 0).t();
    arma::mat hssterm2 = (beta_bar.rows(cdv_ind) % kappa_bar.rows(cdv_ind)).t() * arma::diagmat(cdv) * (beta_bar.rows(cdv_ind) % kappa_bar.rows(cdv_ind));
    
    arma::mat h_ss = arma::diagmat(hssterm1) - hssterm2;
    
    H.submat(0, 0, K - 2, K - 2) = h_pp.submat(0, 0, K - 2, K - 2);
    H.submat(K - 1, 0, 2 * K - 2, K - 2) = h_ps.cols(0, K - 2);
    H.submat(0, K - 1, K - 2, 2 * K - 2) = h_ps.cols(0, K - 2).t();
    H.submat(K - 1, K - 1, 2 * K - 2, 2 * K - 2) = h_ss;
    
    return -H + Sigma_Inv;
}



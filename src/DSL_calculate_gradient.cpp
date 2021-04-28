#include "math.h"
#include <Rcpp.h>
#include <RcppEigen.h>

#include <omp.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]]

//[[Rcpp::export]]
Eigen::ArrayXd DSL_calculate_gradient(Eigen::Map<Eigen::MatrixXd> alpha, Eigen::Map<Eigen::MatrixXd> Xdata, Eigen::Map<Eigen::MatrixXd> Kmatrix, double sigma2)
{
  // omp_set_num_threads(5);
  Eigen::MatrixXd Alpha, X;
  X = Xdata;
  int N = X.rows();
  int P = X.cols();
  //   Alpha.resize(1,N);
  Alpha = alpha;

  Eigen::ArrayXXd K = Kmatrix;

  Eigen::MatrixXd I_1byn(1, N);
  for (int i = 0; i < N; i++)
  {
    I_1byn(0, i) = 1;
  }
  Eigen::ArrayXd derivative = ArrayXd::Zero(P);
#pragma omp parallel for
  for (int k = 0; k < P; k++)
  {
    Eigen::MatrixXd Xk = X.col(k);
    Eigen::ArrayXXd temp = kroneckerProduct(Xk, I_1byn).eval();
    Eigen::ArrayXXd temp2 = (temp - temp.transpose()) / sigma2;
    derivative[k] = (Alpha.transpose() * (-1.0) * ((K * temp2).matrix())).array().square().mean();
    derivative[k] = sqrt(derivative[k]);
  }
  // VectorXd rresult = derivative.sqrt();
  return derivative;
}

//[[Rcpp::export]]
Eigen::ArrayXXd Calculate_Gradient_order2_Cpp(Eigen::Map<Eigen::MatrixXd> alpha, Eigen::Map<Eigen::MatrixXd> Xdata, Eigen::Map<Eigen::MatrixXd> Kmatrix, double sigma2, Eigen::Map<Eigen::ArrayXi> Ind1, Eigen::Map<Eigen::ArrayXi> Ind2)
{

  Eigen::MatrixXd Alpha, X;
  X = Xdata;
  int N = X.rows();
  // int P = X.cols();
  int d1 = Ind1.size();
  int d2 = Ind2.size();
  // for(int i1=0;i1<d1;i1++){
  //   Ind1[i1] = Ind1[i1]-1;
  //  }
  // for(int i1=0;i1<d1;i1++){
  //   Ind1[i1] = Ind1[i1]-1;
  //  }
  Alpha = alpha;
  Eigen::ArrayXXd K = Kmatrix;

  Eigen::MatrixXd I_1byn(1, N);
  for (int i = 0; i < N; i++)
  {
    I_1byn(0, i) = 1;
  }


  Eigen::ArrayXXd derivative = ArrayXXd::Zero(d1, d2);

  for (int i1 = 0; i1 < d1; i1++)
  {
    int k = Ind1[i1]-1;
    Eigen::MatrixXd Xk = X.col(k);
    Eigen::ArrayXXd tempk = kroneckerProduct(Xk, I_1byn).eval();
    Eigen::ArrayXXd Xk_cross = (tempk - tempk.transpose()) / sigma2;

    #pragma omp parallel for
    for (int i2 = 0; i2 < d2; i2++)
    {
      int l = Ind2[i2]-1;
      Eigen::MatrixXd Xl = X.col(l);
      Eigen::ArrayXXd templ = kroneckerProduct(Xl, I_1byn).eval();
      Eigen::ArrayXXd Xl_cross = (templ - templ.transpose()) / sigma2;
      if (k == l)
      {
        derivative(i1, i2) = (Alpha.transpose() * ((-1.0) * (K * (Xk_cross * Xl_cross - (1.0 / sigma2))).matrix())).array().square().mean();
        derivative(i1, i2) = sqrt(derivative(i1, i2));
      }
      else
      {
        derivative(i1, i2) = (Alpha.transpose() *  ((-1.0) * (K * Xk_cross * Xl_cross).matrix())).array().square().mean();
        derivative(i1, i2) = sqrt(derivative(i1, i2));
      }
    }
  }
  // VectorXd rresult = derivative.sqrt();
  return derivative;
}
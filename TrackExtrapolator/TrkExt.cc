#include <iostream>

#include "TrackExtrapolator.hh"

using namespace std;

void printMatrix(const TMatrixD& m, std::string str)
{  
  int nRow = m.GetNrows();
  int nCol = m.GetNcols();
  
  std::cout << "Printing the content of matrix: " << str << std::endl;  
  std::cout << "The matrix has " << nRow << " rows and " << nCol << " columns." << std::endl;
  for(int i = 0; i < nRow; i++)
    {
      std::cout << "Line " << i << ":  "; 
      for(int j = 0; j < nCol; j++)
	{
	  std::cout << m[i][j] << "  ";
	}
      std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
  TrackExtrapolator j;

  j.init("geometry_R997", true);

  TMatrixD state_i(5, 1), state_f(5, 1);
  //state_i[0][0] = 0.0197863;
  //state_i[1][0] = 3.14749/sqrt(1./0.0197863/0.0197863-3.14749*3.14749-0.563128*0.563128);
  //state_i[2][0] = 0.563128/sqrt(1./0.0197863/0.0197863-3.14749*3.14749-0.563128*0.563128);
  //state_i[3][0] = 55.3971;
  //state_i[4][0] = 15.865;

  state_i[0][0] = 0.0203687;
  state_i[1][0] = 0.067167;
  state_i[2][0] = 0.0111264;
  state_i[3][0] = 55.2596;
  state_i[4][0] = 15.3771;

  TMatrixD cov_i(5, 5), cov_f(5, 5);
  for(Int_t i = 0; i < 5; i++)
    {
      state_f[i][0] = 0.;
      for(Int_t j = 0; j < 5; j++)
	{
	  cov_i[i][j] = 0.;
	  cov_f[i][j] = 0.;
	}
    }

  cov_i[0][0] = 1.5E-5;
  cov_i[1][1] = 0.008;
  cov_i[2][2] = 0.009;
  cov_i[3][3] = 0.195;
  cov_i[4][4] = 2.738;

  //cov_i.Zero();

  for(int i = 0; i < 1; i++) {
  j.setInitialStateWithCov(atof(argv[1]), state_i, cov_i);
  j.extrapolateTo(atof(argv[2]));
  j.getFinalStateWithCov(state_f, cov_f);
  }
  
  TMatrixD prop(5, 5);
  j.getPropagator(prop);
    
  TMatrixD state_calc = prop*state_i;
  
  printMatrix(prop, "propagator");
  printMatrix(state_f, "final state vector");
  printMatrix(prop*state_i, "calculated state vector");
  printMatrix(cov_f, "error matrix");
  printMatrix(cov_i, "original error matrix");

  TMatrixD prop_T = prop; prop_T.T();
  printMatrix(prop*cov_i*prop_T, "propagated error matrix");

  printMatrix(cov_f - prop*cov_i*prop_T, "extra added on err matrix");

  return 1;
}

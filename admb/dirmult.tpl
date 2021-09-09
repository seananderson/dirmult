DATA_SECTION
  init_int nyrs
  init_int nages
  init_matrix paa_obs(1,nyrs,1,nages)
  init_vector Neff(1,nyrs)
PARAMETER_SECTION
  init_number log_phi(1)
  init_bounded_vector p(1,nages,0,1)
  objective_function_value objfun
PROCEDURE_SECTION
  eval_objective_function();

FUNCTION double ddirmultinom(const dvar_vector& obs, const dvar_vector& p, const dvariable& phi)
  RETURN_ARRAYS_INCREMENT();
  dvariable N = sum(obs);
  dvariable ll = gammln(N + 1.0) +
                 gammln(phi); -
                 gammln(N + phi);
  for(int a = obs.indexmin(); a <= obs.indexmax(); a++){
    ll += -gammln(obs(a) + 1.0) +
           gammln(obs(a) + phi * (p(a) + 1.0e-15)) -
           gammln(phi * (p(a) + 1.0e-15));
  }
  RETURN_ARRAYS_DECREMENT();
  return(value(ll));

FUNCTION eval_objective_function
  dvar_vector paa_pred = mfexp(p) / sum(mfexp(p));
  for(int i = 1; i < nyrs; i++){
    objfun -= ddirmultinom(Neff(i) * paa_obs(i), paa_pred, exp(log_phi));
  }

REPORT_SECTION
  report << "p " << endl;
  report << p << endl;

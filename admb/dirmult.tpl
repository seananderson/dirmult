DATA_SECTION
  init_int nyrs
  init_int nages
  init_matrix paa_obs(1,nyrs,1,nages)
  init_vector Neff(1,nyrs)

PARAMETER_SECTION
  init_number log_phi
  init_bounded_vector p(1,nages,0,1)
  objective_function_value objfun

PRELIMINARY_CALCS_SECTION
  random_number_generator rng(10293);
  log_phi = log(200) + randn(rng) * 0.5;
  //p.fill_seqadd(1.0, 0.0);

PROCEDURE_SECTION
  eval_objective_function();

FUNCTION double ddirmultinom(const dvar_vector& obs, const dvar_vector& p, const dvariable& phi, int do_log)
  //cout<<"In ddirmultinom\n";
  dvariable N = sum(obs);
  dvariable ll = gammln(N + 1.0) +
                 gammln(phi); -
                 gammln(N + phi);
  for(int a = obs.indexmin(); a <= obs.indexmax(); a++){
    ll += -gammln(obs(a) + 1.0) +
           gammln(obs(a) + phi * (p(a) + 1.0e-15)) -
           gammln(phi * (p(a) + 1.0e-15));
  }
  //return(do_log?value(ll):exp(value(ll)));
  return(value(ll));

FUNCTION eval_objective_function
  dvar_vector paa_pred = mfexp(p) / sum(mfexp(p));
  for(int i = 1; i <= paa_obs.indexmax(); i++){
    objfun -= ddirmultinom(Neff(i) * paa_obs(i), paa_pred, exp(log_phi), 1);
  }

REPORT_SECTION
  report<<"paa_obs.indexmin() = "<<paa_obs.indexmin()<<"\n";
  report<<"paa_obs.indexmax() = "<<paa_obs.indexmax()<<"\n";
  report<<"log_phi\n";
  report<<log_phi<<"\n\n";
  report<<"p\n";
  report<<p<<"\n";
  report<<"objfun = "<<objfun<<"\n";

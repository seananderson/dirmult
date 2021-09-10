DATA_SECTION
  init_int nyrs
  init_int nages
  init_matrix paa_obs(1,nyrs,1,nages)
  init_vector Neff(1,nyrs)

PARAMETER_SECTION
  init_bounded_number log_phi(0.00001,10)
  init_bounded_vector p(1,nages,-10,10)
  vector paa_pred(1,nages)
  vector temp_n(1,nages);
  objective_function_value objfun

PRELIMINARY_CALCS_SECTION
  //random_number_generator rng(10293);
  //log_phi = log(200) + randn(rng) * 0.5;
  //p.fill_seqadd(1.0, 0.0);

PROCEDURE_SECTION
  paa_pred = mfexp(p) / sum(mfexp(p));
  for(int i = 1; i <= paa_obs.indexmax(); i++){
    temp_n = Neff(i) * paa_obs(i);
    objfun -= ddirmultinom(temp_n, paa_pred, log_phi);
  }

FUNCTION dvariable ddirmultinom(dvar_vector& obs, dvar_vector& p, dvariable& log_phi)
  dvariable phi = exp(log_phi);
  dvariable N = sum(obs);
  dvariable ll = gammln(N + 1.0) +
                 gammln(phi) -
                 gammln(N + phi);
  for(int a = obs.indexmin(); a <= obs.indexmax(); a++){
    ll += -gammln(obs(a) + 1.0) +
           gammln(obs(a) + phi * (p(a) + 1.0e-15)) -
           gammln(phi * (p(a) + 1.0e-15));
  }
  return(ll);

REPORT_SECTION
  report<<"log_phi\n";
  report<<log_phi<<"\n";
  report<<"p\n";
  report<<p<<"\n";
  report<<"paa_pred\n";
  report<<paa_pred<<"\n";

#include "SyphUtils.hpp"

#define R_NO_REMAP

using namespace boost;
using namespace std;

extern "C" {
  
  SEXP checkBoostAsserts()
  {
	#ifndef BOOST_DISABLE_ASSERTS
    Rprintf("BOOST ASSERTS ENABLED\n");
	#endif
    return R_NilValue;
  } //end checkBoostAsserts()
  
  SEXP RunSyphProj(SEXP mendist, SEXP mentotal, SEXP menCondomUse, SEXP menScreeningRate, SEXP menNumPart, SEXP menNumSex,
	  SEXP womendist, SEXP womentotal, SEXP womenCondomUse, SEXP womenScreeningRate, SEXP womenNumPart, SEXP womenNumSex, SEXP simtimes, 
	  SEXP mixingMatrix, SEXP ReferralProbMatrix, SEXP turnoverFSW, SEXP fNatHistParameters, SEXP finitInci)
  {
	  Rcpp::List result = Rcpp::List::create();
	  Rcpp::NumericVector in_simtimes(simtimes);
	  Rcpp::NumericMatrix in_mendist(mendist);
	  Rcpp::NumericVector in_mentotal(mentotal);
	  Rcpp::NumericMatrix in_menCondomUse(menCondomUse);
	  Rcpp::NumericMatrix in_menScreeningRate(menScreeningRate);
	  Rcpp::NumericMatrix in_menNumPart(menNumPart);
	  Rcpp::NumericMatrix in_menNumSex(menNumSex);

	  Rcpp::NumericMatrix in_womendist(womendist);
	  Rcpp::NumericVector in_womentotal(womentotal);
	  Rcpp::NumericMatrix in_womenCondomUse(womenCondomUse);
	  Rcpp::NumericMatrix in_womenScreeningRate(womenScreeningRate);
	  Rcpp::NumericMatrix in_womenNumPart(womenNumPart);
	  Rcpp::NumericMatrix in_womenNumSex(womenNumSex);

	  Rcpp::List in_NatHistParameters(fNatHistParameters);
	  Rcpp::NumericVector in_initInci(finitInci);

	  Rcpp::NumericVector in_turnoverFSW(turnoverFSW);

	  int len_simtimes = in_simtimes.size();
	  int len_mentotal = in_mentotal.size();
	  int len_womentotal = in_womentotal.size();
	  int nr_mendist = in_mendist.nrow();
	  int nc_mendist = in_mendist.ncol();
	  int nr_womendist = in_womendist.nrow();
	  int nc_womendist = in_womendist.ncol();

	  int perturbM = in_menScreeningRate.nrow()/5;
	  int perturbW = in_womenScreeningRate.nrow()/5;
	  if ((perturbM != SyphNumRisksMen) | (perturbW != SyphNumRisksWom))
	  {
		  Rcpp::Rcout << ">Error: Not run. Check Screening Rates" << std::endl;
		  len_womentotal++;
		  len_mentotal++;
	  }

	  if ((len_simtimes == nc_mendist) &(len_simtimes == len_mentotal) & (len_simtimes == nc_womendist) &
		  (len_simtimes == len_womentotal) & (len_simtimes >= 2) & (nc_mendist == nc_womendist) & (nc_mendist >= 1) & (nr_mendist == SyphNumRisksMen) &(nr_womendist== SyphNumRisksWom))
	  {
		  //Time
		  Array1d in_time;
		  in_time.resize(boost::extents[len_simtimes]);
		  //Men
		  Array2d arrdist_in_men;
		  arrdist_in_men.resize(boost::extents[nc_mendist][len_simtimes]);
		  Array1d arrtotal_in_men;
		  arrtotal_in_men.resize(boost::extents[len_simtimes]);

		  //Women
		  Array2d arrdist_in_women;
		  arrdist_in_women.resize(boost::extents[nc_womendist][len_simtimes]);
		  Array1d arrtotal_in_women;
		  arrtotal_in_women.resize(boost::extents[len_simtimes]);
		  		 
		  in_time[0] = std::round(in_simtimes[0]);
		  arrtotal_in_men[0] = in_mentotal[0];
		  arrtotal_in_women[0] = in_womentotal[0];
		  		  
		  for (int rg = 0; rg < nr_mendist; rg++)
		  {
			  arrdist_in_men[rg][0] = in_mendist(rg, 0);
		  }

		  for (int rg = 0; rg < nr_womendist; rg++)
		  {
			  arrdist_in_women[rg][0] = in_womendist(rg, 0);
		  }

		  int t0 = 0; 
		  int count = 1;
		  while (TRUE)
		  {
			  if (count >= len_simtimes) break;
			  t0++;
			  in_time[count] = in_time[count - 1] + 1;
			  arrtotal_in_men[count] = in_mentotal[count];
			  arrtotal_in_women[count] = in_womentotal[count];

			  for (int rg = 0; rg < nr_mendist; rg++)
			  {
				  arrdist_in_men[rg][count] = in_mendist(rg, count);
			  }

			  for (int rg = 0; rg < nr_womendist; rg++)
			  {
				  arrdist_in_men[rg][count] = in_mendist(rg, count);
				  arrdist_in_women[rg][count] = in_womendist(rg, count);
			  }
			  count++;
		  }
		  	
		  //Creating Syphilis Population
		  SyphilisPopulation SyphPop= SyphilisPopulation();
		  SyphPop.init(in_time[t0], in_time[0], 1); //Data are given for each year but progression will be performed on a weekly basis
		  SyphPop.SetSubPopMen(arrdist_in_men, arrtotal_in_men, in_time);//Create Men objects
		  SyphPop.SetSubPopWoMen(arrdist_in_women, arrtotal_in_women, in_time);//Create WoMen objects
		  
		  int len_trunoverfsw = in_turnoverFSW.length();
		  if (len_trunoverfsw >= 1)
		  {
			  SyphPop.WoMen.allwomen[SyphHig].turnover.resize(SyphPop.WoMen.allwomen[SyphHig].numsteps);
			  if (len_trunoverfsw == len_simtimes)
			  {
				  for (int li = 0; li < len_simtimes; li++)
				  {
					  SyphPop.WoMen.allwomen[SyphHig].turnover[li] = (in_turnoverFSW[li]>=0) ? in_turnoverFSW[li] :100;
				  }
			  }
			  else
			  {
				  for (int li = 0; li < len_simtimes; li++)
				  {
					  SyphPop.WoMen.allwomen[SyphHig].turnover[li] = (in_turnoverFSW[0]>=0) ? in_turnoverFSW[0]:100;
				  }
				  Rcpp::Rcout << ">Warning! turn over does not have the same length as simtimes, only the first element will be used." << std::endl;
			  }
		  }
		  else
		  {
			  SyphPop.WoMen.allwomen[SyphHig].turnover.resize(SyphPop.WoMen.allwomen[SyphHig].numsteps);
			  for (int li = 0; li < len_simtimes; li++)
			  {
				  SyphPop.WoMen.allwomen[SyphHig].turnover[li] = 100;
			  }
			  Rcpp::Rcout << ">Warning! Turn over rates was non entered for FSW. It will be assumed that they never exit the group." << std::endl;
		  }
		  SyphPop.WoMen.allwomen[SyphHig].is_turnover_on = true;

		  Rcpp::NumericMatrix in_mixingMatrix(mixingMatrix);
		  int numelmmr = in_mixingMatrix.nrow();
		  int numelmmc = in_mixingMatrix.ncol();
		  if ((numelmmr == numelmmc) &(numelmmr == (SyphNumRisksMen+ SyphNumRisksWom)))
		  {
			  std::vector<std::vector<double>> part_mixmat;
			  part_mixmat.resize(numelmmr);
			  for (int ii = 0; ii < numelmmr; ii++)
			  {
				  part_mixmat[ii].resize(numelmmr);
				  for (int jj = 0; jj < numelmmr; jj++)
				  {
					  part_mixmat[ii][jj] = in_mixingMatrix(ii, jj);
				  }
			  }
			  SyphPop.mixingmat.SetIt(part_mixmat);
		  }
		  else
		  {
			  Rcpp::Rcout << "Warning! Wrong dimension for the mixing matrix; default values will be used instead" << std::endl;
		  }

		  Rcpp::NumericMatrix in_ReferralProbMatrix(ReferralProbMatrix);
		  int numelrmr = in_ReferralProbMatrix.nrow();
		  int numelrmc = in_ReferralProbMatrix.ncol();
		  if ((numelrmr == numelrmc) &(numelrmr == (SyphNumRisksMen + SyphNumRisksWom)))
		  {
			  std::vector<std::vector<double>> part_refmat;
			  part_refmat.resize(numelrmr);
			  for (int ii = 0; ii < numelrmr; ii++)
			  {
				  part_refmat[ii].resize(numelrmr);
				  for (int jj = 0; jj < numelrmr; jj++)
				  {
					  part_refmat[ii][jj] = in_ReferralProbMatrix(ii, jj);
				  }
			  }
			  SyphPop.referralmat.SetIt(part_refmat);
		  }
		  else
		  {
			  Rcpp::Rcout << "Warning! Wrong dimension for the mixing matrix; default values will be used instead" << std::endl;
		  }
		  
		  //Set the number of partners and the number of sexacts per partner
		  for (int t = 0; t < len_simtimes; t++)
		  {
			  //Sexual Behavior by risk group. Low
			  SyphPop.Men.allmen[SyphLow].MeanNumPart[t] = in_menNumPart(SyphLow, t);
			  SyphPop.Men.allmen[SyphLow].MeanNumSexByPart[t] = in_menNumSex(SyphLow,t); 
			  SyphPop.Men.allmen[SyphLow].CondomUse[t] = in_menCondomUse(SyphLow, t);
			  SyphPop.Men.allmen[SyphLow].nh_param.screen_rate[t] = in_menScreeningRate(SyphLow, t);
			  SyphPop.Men.allmen[SyphLow].nh_param.screen_rateTPHA[t] = in_menScreeningRate(SyphNumRisksMen+SyphLow, t);
			  SyphPop.Men.allmen[SyphLow].nh_param.vect_tau[t] = in_menScreeningRate(2*SyphNumRisksMen+ SyphLow, t);
			  SyphPop.Men.allmen[SyphLow].nh_param.prop_part_referred[t] = in_menScreeningRate(3 * SyphNumRisksMen + SyphLow, t);
			  SyphPop.Men.allmen[SyphLow].pinfectdetect[t] = in_menScreeningRate(4 * SyphNumRisksMen + SyphLow, t);
			  //Medium
			  SyphPop.Men.allmen[SyphMed].MeanNumPart[t] = in_menNumPart(SyphMed, t);
			  SyphPop.Men.allmen[SyphMed].MeanNumSexByPart[t] = in_menNumSex(SyphMed, t); 
			  SyphPop.Men.allmen[SyphMed].CondomUse[t] = in_menCondomUse(SyphMed, t);
			  SyphPop.Men.allmen[SyphMed].nh_param.screen_rate[t] = in_menScreeningRate(SyphMed, t);
			  SyphPop.Men.allmen[SyphMed].nh_param.screen_rateTPHA[t] = in_menScreeningRate(SyphNumRisksMen+SyphMed, t);
			  SyphPop.Men.allmen[SyphMed].nh_param.vect_tau[t] = in_menScreeningRate(2*SyphNumRisksMen + SyphMed, t);
			  SyphPop.Men.allmen[SyphMed].nh_param.prop_part_referred[t] = in_menScreeningRate(3 * SyphNumRisksMen + SyphMed, t);
			  SyphPop.Men.allmen[SyphMed].pinfectdetect[t] = in_menScreeningRate(4 * SyphNumRisksMen + SyphMed, t);
			  //High
			  SyphPop.Men.allmen[SyphHig].MeanNumPart[t] = in_menNumPart(SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].MeanNumSexByPart[t] = in_menNumSex(SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].CondomUse[t] = in_menCondomUse(SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].nh_param.screen_rate[t] = in_menScreeningRate(SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].nh_param.screen_rateTPHA[t] = in_menScreeningRate(SyphNumRisksMen+SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].nh_param.vect_tau[t] = in_menScreeningRate(2*SyphNumRisksMen + SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].nh_param.prop_part_referred[t] = in_menScreeningRate(3 * SyphNumRisksMen + SyphHig, t);
			  SyphPop.Men.allmen[SyphHig].pinfectdetect[t] = in_menScreeningRate(4 * SyphNumRisksMen + SyphHig, t);
			  
			  //MSM
			  SyphPop.Men.allmen[SyphMSM].MeanNumPart[t] = in_menNumPart(SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].MeanNumSexByPart[t] = in_menNumSex(SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].CondomUse[t] = in_menCondomUse(SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].nh_param.screen_rate[t] = in_menScreeningRate(SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].nh_param.screen_rateTPHA[t] = in_menScreeningRate(SyphNumRisksMen+SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].nh_param.vect_tau[t] = in_menScreeningRate(2*SyphNumRisksMen + SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].nh_param.prop_part_referred[t] = in_menScreeningRate(3 * SyphNumRisksMen + SyphMSM, t);
			  SyphPop.Men.allmen[SyphMSM].pinfectdetect[t] = in_menScreeningRate(4 * SyphNumRisksMen + SyphMSM, t);
			  //No Sex
			  SyphPop.Men.allmen[SyphNoSex].MeanNumPart[t] = in_menNumPart(SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].MeanNumSexByPart[t] = in_menNumSex(SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].CondomUse[t] = in_menCondomUse(SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].nh_param.screen_rate[t] = in_menScreeningRate(SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].nh_param.screen_rateTPHA[t] = in_menScreeningRate(SyphNumRisksMen+SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].nh_param.vect_tau[t] = in_menScreeningRate(2*SyphNumRisksMen + SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].nh_param.prop_part_referred[t] = in_menScreeningRate(3 * SyphNumRisksMen + SyphNoSex, t);
			  SyphPop.Men.allmen[SyphNoSex].pinfectdetect[t] = in_menScreeningRate(4 * SyphNumRisksMen + SyphNoSex, t);
			  //Women. Low
			  SyphPop.WoMen.allwomen[SyphLow].MeanNumPart[t] = in_womenNumPart(SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].MeanNumSexByPart[t] = in_womenNumSex(SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].CondomUse[t] = in_womenCondomUse(SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].nh_param.screen_rate[t] = in_womenScreeningRate(SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].nh_param.screen_rateTPHA[t] = in_womenScreeningRate(SyphNumRisksWom+ SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].nh_param.vect_tau[t] = in_womenScreeningRate(2*SyphNumRisksWom + SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].nh_param.prop_part_referred[t] = in_womenScreeningRate(3 * SyphNumRisksWom + SyphLow, t);
			  SyphPop.WoMen.allwomen[SyphLow].pinfectdetect[t] = in_womenScreeningRate(4 * SyphNumRisksWom + SyphLow, t);
			  //Medium
			  SyphPop.WoMen.allwomen[SyphMed].MeanNumPart[t] = in_womenNumPart(SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].MeanNumSexByPart[t] = in_womenNumSex(SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].CondomUse[t] = in_womenCondomUse(SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].nh_param.screen_rate[t] = in_womenScreeningRate(SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].nh_param.screen_rateTPHA[t] = in_womenScreeningRate(SyphNumRisksWom+ SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].nh_param.vect_tau[t] = in_womenScreeningRate(2*SyphNumRisksWom + SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].nh_param.prop_part_referred[t] = in_womenScreeningRate(3 * SyphNumRisksWom + SyphMed, t);
			  SyphPop.WoMen.allwomen[SyphMed].pinfectdetect[t] = in_womenScreeningRate(4 * SyphNumRisksWom + SyphMed, t);
			  //High
			  SyphPop.WoMen.allwomen[SyphHig].MeanNumPart[t] = in_womenNumPart(SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].MeanNumSexByPart[t] = in_womenNumSex(SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].CondomUse[t] = in_womenCondomUse(SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].nh_param.screen_rate[t] = in_womenScreeningRate(SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].nh_param.screen_rateTPHA[t] = in_womenScreeningRate(SyphNumRisksWom+SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].nh_param.vect_tau[t] = in_womenScreeningRate(2*SyphNumRisksWom + SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].nh_param.prop_part_referred[t] = in_womenScreeningRate(3 * SyphNumRisksWom + SyphHig, t);
			  SyphPop.WoMen.allwomen[SyphHig].pinfectdetect[t] = in_womenScreeningRate(4 * SyphNumRisksWom + SyphHig, t);
			  //No sex.
			  SyphPop.WoMen.allwomen[SyphNoSexWom].MeanNumPart[t] = in_womenNumPart(SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].MeanNumSexByPart[t] = in_womenNumSex(SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].CondomUse[t] = in_womenCondomUse(SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].nh_param.screen_rate[t] = in_womenScreeningRate(SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].nh_param.screen_rateTPHA[t] = in_womenScreeningRate(SyphNumRisksWom+ SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].nh_param.vect_tau[t] = in_womenScreeningRate(2*SyphNumRisksWom + SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].nh_param.prop_part_referred[t] = in_womenScreeningRate(3 * SyphNumRisksWom + SyphNoSexWom, t);
			  SyphPop.WoMen.allwomen[SyphNoSexWom].pinfectdetect[t] = in_womenScreeningRate(4 * SyphNumRisksWom + SyphNoSexWom, t);
		  }

		  int len_init_inci = in_initInci.size();
		  Rcpp::StringVector nomdeNat= in_NatHistParameters.names();
		  int len_nomdeNat = in_NatHistParameters.length();
		  std::vector<std::string> cpnomdeNat;
		  cpnomdeNat.resize(len_nomdeNat);
		  for (int li = 0; li < len_nomdeNat; li++)
		  {
			  cpnomdeNat[li] = nomdeNat[li];
		  }

		  //MEN
		  for (int pop = 0; pop < SyphNumRisksMen; pop++)
		  {
			  //Make sure that all the natural history parameters are available
			  if (len_nomdeNat >= 8)
			  {
				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration1") != cpnomdeNat.end())
				  {
					  double Duration1 = in_NatHistParameters["Duration1"];
					  SyphPop.Men.allmen[pop].nh_param.sigma1 = (Duration1 > 0) ? 1.0 / Duration1  : SyphPop.Men.allmen[pop].nh_param.sigma1;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration2") != cpnomdeNat.end())
				  {
					  double Duration2 = in_NatHistParameters["Duration2"];
					  SyphPop.Men.allmen[pop].nh_param.sigma2 = (Duration2 > 0) ? 1.0 / Duration2 : SyphPop.Men.allmen[pop].nh_param.sigma2;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration3") != cpnomdeNat.end())
				  {
					  double Duration3 = in_NatHistParameters["Duration3"];
					  SyphPop.Men.allmen[pop].nh_param.sigma3 = (Duration3 > 0) ? 1.0 / Duration3  : SyphPop.Men.allmen[pop].nh_param.sigma3;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration4") != cpnomdeNat.end())
				  {
					  double Duration4 = in_NatHistParameters["Duration4"];
					  SyphPop.Men.allmen[pop].nh_param.sigma4 = (Duration4 > 0) ? 1.0 / Duration4  : SyphPop.Men.allmen[pop].nh_param.sigma4;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration5") != cpnomdeNat.end())
				  {
					  double Duration5 = in_NatHistParameters["Duration5"];
					  SyphPop.Men.allmen[pop].nh_param.sigma5 = (Duration5 > 0) ? 1.0 / Duration5  : SyphPop.Men.allmen[pop].nh_param.sigma5;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "nu") != cpnomdeNat.end())
				  {
					  double nu = in_NatHistParameters["nu"];
					  SyphPop.Men.allmen[pop].nh_param.nu = ((nu >= 0) & (nu <= 1)) ? nu : SyphPop.Men.allmen[pop].nh_param.nu;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "phi") != cpnomdeNat.end())
				  {
					  double phi = in_NatHistParameters["phi"];
					  SyphPop.Men.allmen[pop].nh_param.phi = ((phi >= 0) & (phi <= 1)) ? phi : SyphPop.Men.allmen[pop].nh_param.phi;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "psi") != cpnomdeNat.end())
				  {
					  double psi = in_NatHistParameters["psi"];
					  SyphPop.Men.allmen[pop].nh_param.psi = ((psi >= 0) & (psi <= 1)) ? psi : SyphPop.Men.allmen[pop].nh_param.psi;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "MtoMtp") != cpnomdeNat.end())
				  {
					  double mtomtp = in_NatHistParameters["MtoMtp"];
					  SyphPop.Men.allmen[pop].nh_param.MtoM_tp = ((mtomtp >= 0) & (mtomtp <= 1)) ? mtomtp : SyphPop.Men.allmen[pop].nh_param.MtoM_tp;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "MtoFtp") != cpnomdeNat.end())
				  {
					  double mtoftp = in_NatHistParameters["MtoFtp"];
					  SyphPop.Men.allmen[pop].nh_param.MtoF_tp = ((mtoftp >= 0) & (mtoftp <= 1)) ? mtoftp : SyphPop.Men.allmen[pop].nh_param.MtoF_tp;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "FtoMtp") != cpnomdeNat.end())
				  {
					  double ftomtp = in_NatHistParameters["FtoMtp"];
					  SyphPop.Men.allmen[pop].nh_param.FtoM_tp = ((ftomtp >= 0) & (ftomtp <= 1)) ? ftomtp : SyphPop.Men.allmen[pop].nh_param.FtoM_tp;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "CondomEfficacy") != cpnomdeNat.end())
				  {
					  double CondomEfficacy = in_NatHistParameters["CondomEfficacy"];
					  SyphPop.Men.allmen[pop].nh_param.condomeff = ((CondomEfficacy >= 0) & (CondomEfficacy <=1)) ? CondomEfficacy : SyphPop.Men.allmen[pop].nh_param.condomeff;
				  }
			  }

			  SyphPop.Men.allmen[pop].nh_param.updateTreatment();

			  if (len_init_inci > (pop))
			  {
				  SyphPop.Men.allmen[pop].nh_param.ValinitInci = in_initInci[pop];
			  }
			  SyphPop.Men.allmen[pop].SeedSteadyState();
		  } //end for (int pop = 0; pop < SyphNumRisksMen; pop++)

		  //WOMEN
		  for (int pop = 0; pop < SyphNumRisksWom; pop++)
		  {
			  if (len_nomdeNat >= 8)
			  {
				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration1") != cpnomdeNat.end())
				  {
					  double Duration1 = in_NatHistParameters["Duration1"];
					  SyphPop.WoMen.allwomen[pop].nh_param.sigma1 = (Duration1 > 0) ? 1.0 / Duration1  : SyphPop.WoMen.allwomen[pop].nh_param.sigma1;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration2") != cpnomdeNat.end())
				  {
					  double Duration2 = in_NatHistParameters["Duration2"];
					  SyphPop.WoMen.allwomen[pop].nh_param.sigma2 = (Duration2 > 0) ? 1.0 / Duration2  : SyphPop.WoMen.allwomen[pop].nh_param.sigma2;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration3") != cpnomdeNat.end())
				  {
					  double Duration3 = in_NatHistParameters["Duration3"];
					  SyphPop.WoMen.allwomen[pop].nh_param.sigma3 = (Duration3 > 0) ? 1.0 / Duration3 : SyphPop.WoMen.allwomen[pop].nh_param.sigma3;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration4") != cpnomdeNat.end())
				  {
					  double Duration4 = in_NatHistParameters["Duration4"];
					  SyphPop.WoMen.allwomen[pop].nh_param.sigma4 = (Duration4 > 0) ? 1.0 / Duration4  : SyphPop.WoMen.allwomen[pop].nh_param.sigma4;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "Duration5") != cpnomdeNat.end())
				  {
					  double Duration5 = in_NatHistParameters["Duration5"];
					  SyphPop.WoMen.allwomen[pop].nh_param.sigma5 = (Duration5 > 0) ? 1.0 / Duration5  : SyphPop.WoMen.allwomen[pop].nh_param.sigma5;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "nu") != cpnomdeNat.end())
				  {
					  double nu = in_NatHistParameters["nu"];
					  SyphPop.WoMen.allwomen[pop].nh_param.nu = ((nu >= 0) & (nu <= 1)) ? nu : SyphPop.WoMen.allwomen[pop].nh_param.nu;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "phi") != cpnomdeNat.end())
				  {
					  double phi = in_NatHistParameters["phi"];
					  SyphPop.WoMen.allwomen[pop].nh_param.phi = ((phi >= 0) & (phi <= 1)) ? phi : SyphPop.WoMen.allwomen[pop].nh_param.phi;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "psi") != cpnomdeNat.end())
				  {
					  double psi = in_NatHistParameters["psi"];
					  SyphPop.WoMen.allwomen[pop].nh_param.psi = ((psi >= 0) & (psi <= 1)) ? psi : SyphPop.WoMen.allwomen[pop].nh_param.psi;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "MtoMtp") != cpnomdeNat.end())
				  {
					  double mtomtp = in_NatHistParameters["MtoMtp"];
					  SyphPop.WoMen.allwomen[pop].nh_param.MtoM_tp = ((mtomtp >= 0) & (mtomtp <= 1)) ? mtomtp : SyphPop.WoMen.allwomen[pop].nh_param.MtoM_tp;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "MtoFtp") != cpnomdeNat.end())
				  {
					  double mtoftp = in_NatHistParameters["MtoFtp"];
					  SyphPop.WoMen.allwomen[pop].nh_param.MtoF_tp = ((mtoftp >= 0) & (mtoftp <= 1)) ? mtoftp : SyphPop.WoMen.allwomen[pop].nh_param.MtoF_tp;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "FtoMtp") != cpnomdeNat.end())
				  {
					  double ftomtp = in_NatHistParameters["FtoMtp"];
					  SyphPop.WoMen.allwomen[pop].nh_param.FtoM_tp = ((ftomtp >= 0) & (ftomtp <= 1)) ? ftomtp : SyphPop.WoMen.allwomen[pop].nh_param.FtoM_tp;
				  }

				  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "CondomEfficacy") != cpnomdeNat.end())
				  {
					  double CondomEfficacy = in_NatHistParameters["CondomEfficacy"];
					  SyphPop.WoMen.allwomen[pop].nh_param.condomeff = ((CondomEfficacy >= 0) & (CondomEfficacy <= 1)) ? CondomEfficacy : SyphPop.WoMen.allwomen[pop].nh_param.condomeff;
				  }
			  }

			  SyphPop.WoMen.allwomen[pop].nh_param.updateTreatment();

			  if (len_init_inci > (SyphNumRisksMen+pop))
			  {
				  SyphPop.WoMen.allwomen[pop].nh_param.ValinitInci = in_initInci[SyphNumRisksMen+pop];
			  }
			  
			  SyphPop.WoMen.allwomen[pop].SeedSteadyState();
		  } //end for (int pop = 0; pop < SyphNumRisksWom; pop++)
		 
		  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "w_assortativeness") != cpnomdeNat.end())
		  {
			  double w_assortativeness = in_NatHistParameters["w_assortativeness"];
			  SyphPop.w_assortativeness = ((w_assortativeness >= 0) & (w_assortativeness <= 1)) ? w_assortativeness : SyphPop.w_assortativeness;
		  }

		  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "assortativity") != cpnomdeNat.end())
		  {
			  double assortativity = in_NatHistParameters["assortativity"];
			  SyphPop.assortativity = ((assortativity >= 0) & (assortativity <= 1)) ? assortativity : SyphPop.assortativity;
		  }

		  if (std::find(cpnomdeNat.begin(), cpnomdeNat.end(), "MRtoHR_assortativ") != cpnomdeNat.end())
		  {
			  double MRtoHR_assortativ = in_NatHistParameters["MRtoHR_assortativ"];
			  SyphPop.MRtoHR_assortativ = ((MRtoHR_assortativ >= 0) & (MRtoHR_assortativ <= 1)) ? MRtoHR_assortativ : SyphPop.MRtoHR_assortativ;
		  }
		  else
		  {
			  Rcpp::Rcout << "Medium risk to high risk assortativity parameter missing: the default value 1 will be used!" << std::endl;
		  }

		  int in_DT = 52; // Projections scheduled every week
		  doSyphilisCalc(&SyphPop, in_DT);
		  SyphilisPopulationforR resultcpp = SyphilisPopulationforR(SyphPop);
		  result = resultcpp.outforR(in_time[0]);
	  }

	  return Rcpp::wrap(result);
  }

  //Calculates prevalence (TPHA+) at equilibrium
  SEXP calcRPRPrevCpp(SEXP in_ValinitInci, SEXP in_nu, SEXP in_tau, SEXP in_phi, SEXP in_psi, SEXP in_sigma1, SEXP in_sigma2, SEXP in_sigma3, SEXP in_sigma4,
	  SEXP in_sigma5, SEXP in_screen_rate, SEXP in_screen_rateTPHA)
  {
	  double lambda = REAL(in_ValinitInci)[0];
	  double nu = REAL(in_nu)[0];
	  double tau = REAL(in_tau)[0];
	  double phi = REAL(in_phi)[0];
	  double psi = REAL(in_psi)[0];
	  double sigma1 = REAL(in_sigma1)[0];
	  double sigma2 = REAL(in_sigma2)[0];
	  double sigma3 = REAL(in_sigma3)[0];
	  double sigma4 = REAL(in_sigma4)[0];
	  double sigma5 = REAL(in_sigma5)[0];
	  double screen_rate = REAL(in_screen_rate)[0];
	  double screen_rateTPHA = REAL(in_screen_rateTPHA)[0];

	  nu = nu * tau * sigma2;
	  double r_ActiveToSus = (nu + screen_rate + screen_rateTPHA)*psi*phi;
	  double r_ActiveRecov = (nu + screen_rate + screen_rateTPHA)*psi*(1 - phi);
	  double x = r_ActiveToSus;
	  double y = r_ActiveRecov;
	  double z = (screen_rate + screen_rateTPHA)*psi;
	  double mu = 1.0 / 35.0;
	  double numSus = 1.0;

	  double be = (lambda + mu)*numSus;
	  double cteR3e = 1.0 - lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*sigma3 / (sigma3 + z + mu)*sigma5 / (sigma5 + mu)
		  - lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*sigma4 / (sigma4 + mu)
		  - lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*sigma4 / (sigma4 + mu)
		  - lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*x / (x + y + sigma2 + mu);

	  double R3e = 1.0 / cteR3e;
	  R3e *= be * (1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*sigma3 / (sigma3 + z + mu)*sigma5 / (sigma5 + mu)
		  + 1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*sigma4 / (sigma4 + mu)
		  + 1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*sigma4 / (sigma4 + mu)
		  + 1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*x / (x + y + sigma2 + mu));

	  double R2e = R3e * (lambda*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*1.0 / (sigma4 + mu)
		  + lambda * sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*1.0 / (sigma4 + mu));
	  R2e += be * (lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*1.0 / (sigma4 + mu)
		  + lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*1.0 / (sigma4 + mu));

	  double R1e = lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*sigma3 / (sigma3 + z + mu)*1.0 / (sigma5 + mu)*be
		  + lambda * sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*sigma3 / (sigma3 + z + mu)*1.0 / (sigma5 + mu)*R3e;

	  double Le = lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*1.0 / (sigma3 + z + mu)*be;
	  Le += lambda * sigma1 / (sigma1 + mu) *sigma2 / (x + y + sigma2 + mu)*1.0 / (z + sigma3 + mu)*R3e;

	  double Ae = lambda / (lambda + mu)* sigma1 / (sigma1 + mu) *1.0 / (x + y + sigma2 + mu)*be;
	  Ae += lambda * sigma1 / (sigma1 + mu)*1.0 / (x + y + sigma2 + mu)*R3e;

	  double Ie = lambda / (lambda + mu)*1.0 / (sigma1 + mu) *be + lambda / (sigma1 + mu)*R3e;
	  double total = numSus + Ie + Ae + Le + R1e + R2e + R3e;

	  double prevRPRPos = (Ae + Le + R1e + R2e ) / total;
	  Rcpp::NumericVector result(1);
	  result[0] = prevRPRPos;
	  return(Rcpp::wrap(result));
  }

  SEXP expitCpp(SEXP x)
  {
	  Rcpp::NumericVector nx(x);
	  int len_x = nx.size();
	  Rcpp::NumericVector Result = clone(nx);
	  for(int ii=0; ii< len_x; ii++)
	  {
		  Result[ii] = expit(Result[ii]);
	  }
	  return(Rcpp::wrap(Result));
  }

  SEXP logitCpp(SEXP x)
  {
	  Rcpp::NumericVector nx(x);
	  int len_x = nx.size();
	  Rcpp::NumericVector Result=Rcpp::clone(nx);
	  for (int ii = 0; ii < len_x; ii++)
	  {
		  Result[ii] = logit(Result[ii]);
	  }
	  return(Rcpp::wrap(Result));
  }

}

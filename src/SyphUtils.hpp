#include "SyphConsts.hpp"

/**************************DECLARATIONS*******************************/
class sub_population //This is the atom of the population to be simulated
{
public:
	int tinit, tfinal, numsteps;//tinit: intial, tfinal final time; numsteps:number of time steps
	double dt;//Time Step
	Array1d SSneg, IIncub, IActive;//Negative, Incubation and Prim+Sec
	Array1d EarlyRecov, ILatent, LateRecov, RecovSSneg, Total;//Box 4, Latent, Box 5, box 6 and Total population
	Array1d MeanNumPart, MeanNumSexByPart;//Number of partners and number of se by partner
	Array1d CondomUse;//condom Use
	Array1d Incidence, Prevalence;//Incidence and Prevalence
	Array2d NumScreened, NumTreatedaftScreened;//Number screened, number treated and screened
	Array1d NumIncidentCases, NumSymptomaticClinTreat;//Number of incident cases; Number of symptomatic clinically treated
	Array1d NumReferred;//Number referred to clinic (not to number of partner traced!)
	Array2d NumReferredTestedPositive;//
	Array2d NumContactTraced;
	std::vector<double> referral_rate; // Referral rate, will be calculated as a function of the number of partners with symptomatic Syphilis
	std::vector<double> pinfectdetect; //propensity to trace infectious partner
	Array2d Prevalence2d;//Prevalence (TPHA and RPR)
	std::vector<double> turnover;//Only turned on for FSW
	bool is_turnover_on;//Only tru for FSW

	NaturalHistParam nh_param;

	sub_population& operator = (const sub_population &t)
	{
		this->tinit=t.tinit;
		this->tfinal=t.tfinal;
		this->numsteps=t.numsteps;
		this->dt=t.dt;
		
		this->SSneg.resize(boost::extents[t.numsteps]);
		this->SSneg = boost::multi_array_ref<double,1>(t.SSneg);
		this->IIncub.resize(boost::extents[t.numsteps]);
		this->IIncub = boost::multi_array_ref<double, 1>(t.IIncub);

		this->IActive.resize(boost::extents[t.numsteps]);
		this->IActive = boost::multi_array_ref<double, 1>(t.IActive);

		this->RecovSSneg.resize(boost::extents[t.numsteps]);
		this->RecovSSneg = boost::multi_array_ref<double, 1>(t.RecovSSneg);

		this->EarlyRecov.resize(boost::extents[t.numsteps]);
		this->EarlyRecov = boost::multi_array_ref<double, 1>(t.EarlyRecov);

		this->ILatent.resize(boost::extents[t.numsteps]);
		this->ILatent = boost::multi_array_ref<double, 1>(t.ILatent);
		
		this->LateRecov.resize(boost::extents[t.numsteps]);
		this->LateRecov = boost::multi_array_ref<double, 1>(t.LateRecov);
		
		this->Total.resize(boost::extents[t.numsteps]);
		this->Total = boost::multi_array_ref<double, 1>(t.Total);
		
		this->MeanNumSexByPart.resize(boost::extents[t.numsteps]);
		this->MeanNumSexByPart = boost::multi_array_ref<double, 1>(t.MeanNumSexByPart);
		
		this->MeanNumPart.resize(boost::extents[t.numsteps]);
		this->MeanNumPart = boost::multi_array_ref<double, 1>(t.MeanNumPart);

		this->Prevalence.resize(boost::extents[t.numsteps]);
		this->Prevalence = boost::multi_array_ref<double, 1>(t.Prevalence);

		this->Prevalence2d.resize(boost::extents[numPrevTypes][t.numsteps]);
		this->Prevalence2d = boost::multi_array_ref<double, 2>(t.Prevalence2d);

		this->Incidence.resize(boost::extents[t.numsteps]);
		this->Incidence = boost::multi_array_ref<double, 1>(t.Incidence);
		this->poptime = t.poptime;

		this->CondomUse.resize(boost::extents[t.numsteps]);
		this->CondomUse = boost::multi_array_ref<double, 1>(t.CondomUse);

		this->NumScreened.resize(boost::extents[numTests][t.numsteps]);
		this->NumScreened = boost::multi_array_ref<double, 2>(t.NumScreened);
		
		this->NumTreatedaftScreened.resize(boost::extents[numTests][t.numsteps]);
		this->NumTreatedaftScreened = boost::multi_array_ref<double, 2>(t.NumTreatedaftScreened);
		
		this->NumIncidentCases.resize(boost::extents[t.numsteps]);
		this->NumIncidentCases = boost::multi_array_ref<double, 1>(t.NumIncidentCases);

		this->NumSymptomaticClinTreat.resize(boost::extents[t.numsteps]);
		this->NumSymptomaticClinTreat = boost::multi_array_ref<double, 1>(t.NumSymptomaticClinTreat);

		this->NumReferred.resize(boost::extents[t.numsteps]);
		this->NumReferred = boost::multi_array_ref<double, 1>(t.NumReferred);

		this->NumReferredTestedPositive.resize(boost::extents[numTests][t.numsteps]);
		this->NumReferredTestedPositive = boost::multi_array_ref<double, 2>(t.NumReferredTestedPositive);

		this->NumContactTraced.resize(boost::extents[sizeContactTraced][t.numsteps]);
		this->NumContactTraced = boost::multi_array_ref<double, 2>(NumContactTraced);
		
		this->referral_rate = t.referral_rate;

		this->pinfectdetect = t.pinfectdetect;

		this->is_turnover_on = t.is_turnover_on;
		this->turnover = t.turnover;
		this->nh_param = t.nh_param;

		return *this;
	}

	void init(int t_final, int t_init, int in_dt);
	void calcPrevalence();
	void SeedSteadyState();
	void ProgressPopulation(int n, const double foi);
	void copyfrom_sub_population(int n1, sub_population subpop, int n2, bool cpnath);
	std::vector<double> poptime;
	sub_population();
	~sub_population();
};

class Men_Pop //Men populations
{
public:
	std::vector<sub_population> allmen; //vector of sub-populations
	std::vector<double> total;
	Array2d distributions;

	Men_Pop& operator = (const Men_Pop &t)
	{
		int len_pop = t.allmen.size();
		this->allmen.resize(len_pop);
		for (int pop = 0; pop < len_pop; pop++)
		{
			this->allmen[pop] = t.allmen[pop];
		}
		this->total = t.total;
		if (len_pop > 0)
		{
			this->distributions.resize(boost::extents[SyphNumRisksMen][t.allmen[0].numsteps]);
			this->distributions = boost::multi_array_ref<double, 2>(t.distributions);
		}
		return *this;
	}

	void init(int t_final, int t_init, int in_dt);
	void SetSubPop(Array2d in_dist, Array1d in_total, Array1d in_time);
	void ProgressPopulation(int n, const std::vector<double> vect_foi);
	void copyfromMen(int n1, Men_Pop popmen2, int n2, bool cpnath);
	Men_Pop();
	~Men_Pop();
}; //end definition of the class Men_Pop

class WoMen_Pop //Women populations
{
public:
	std::vector<sub_population> allwomen; //vector of sub-populations
	std::vector<double> total;
	Array2d distributions;

	WoMen_Pop& operator = (const WoMen_Pop &t)
	{
		int len_pop = t.allwomen.size();
		this->allwomen.resize(len_pop);
		for (int pop = 0; pop < len_pop; pop++)
		{
			this->allwomen[pop] = t.allwomen[pop];
		}
		this->total = t.total;
		if (len_pop > 0)
		{
			this->distributions.resize(boost::extents[SyphNumRisksWom][t.allwomen[0].numsteps]);
			this->distributions = boost::multi_array_ref<double, 2>(t.distributions);
		}

		return *this;
	}

	void init(int t_final, int t_init, int in_dt);
	void SetSubPop(Array2d in_dist, Array1d in_total, Array1d in_time);
	void ProgressPopulation(int n, const std::vector<double> vect_foi);
	void copyfromWoMen(int n1, WoMen_Pop popwomen2, int n2, bool cpnath);
	void turnFSW_Over(int n);
	WoMen_Pop();
	~WoMen_Pop();
}; //end class WoMen_Pop; data for women

class SyphilisPopulation
{
public:
	Men_Pop Men;//Men population, by risk group
	WoMen_Pop WoMen; //Women population, by risk group
	MixingMatrixParam mixingmat;//Mixing matrix
	void calcMixingBalancing(const int&n);
	//Array2d BalancedMixing_ExtraMarrital;
	Array2d BalancedMixing_Marrital;
	Array2d Balancing;//Partner choice
	Array2d mixingGoals;//Mixing Goals
	//Array2d mixingProportional;//Mixing Proportionnal

	std::vector< Array2d> all_mixing_t1;
	std::vector< Array2d> all_mixing_tf;

	void calcMixingProportional(const int&n);
	Array2d Mixing_Proportional;
	double assortativity;
	double w_assortativeness;
	double MRtoHR_assortativ;
	ReferralParam referralmat; //Referral Probabilities
	std::vector<double> FOIMen;//Force of infection for men
	std::vector<double> FOIWoMen; // Force of infection for women
	SyphilisPopulation& operator = (const SyphilisPopulation &t)
	{
		this->Men = t.Men;
		this->WoMen = t.WoMen;
		this->mixingmat = t.mixingmat;
		this->FOIMen = t.FOIMen;
		this->FOIWoMen = t.FOIWoMen;

		return *this;
	}

	void init(int t_final, int t_init, int in_dt);
	void SetSubPopMen(Array2d in_dist, Array1d in_total, Array1d in_time);
	void SetSubPopWoMen(Array2d in_dist, Array1d in_total, Array1d in_time);
	void calcFOI(int n);
	void calcinitreferral();
	void ProgressPopulation(int n);
	void copyMWPops(int n1, SyphilisPopulation mwpop2, int n2, bool cpnath);
	void copyMWPops(int n1, SyphilisPopulation *mwpop2, int n2, bool cpnath);
	SyphilisPopulation();
	~SyphilisPopulation();
};//end of declaration of SyphilisPopulation

//a class for interaction with R
class SyphilisPopulationforR
{
public:
	Rcpp::List Men;
	Rcpp::List WoMen;
	Rcpp::NumericMatrix mixingmat;
	Rcpp::List InitMixingMat;
	Rcpp::List FinalMixingMat;
	Rcpp::List InitBalancedPartners;
	Rcpp::List FinalBalancedPartners;
	Rcpp::List nathistory;
	Rcpp::List outforR(double firstyear);
	Rcpp::List outforR();
	SyphilisPopulationforR(const SyphilisPopulation sampop);
	~SyphilisPopulationforR();
}; //end declaration of SyphilisPopulationforR; 

//////////////////////////////////////////////////////////////
/*********************IMPLEMENTATION*************************/
//////////////////////////////////////////////////////////////

void copytosub_population(sub_population *pop1, int n1, sub_population pop2, int n2) // Copy a slice of pop2 into a slice of pop1
{
	if ((n1 < pop1->numsteps) & (n1 > 0) &(n2 < pop2.numsteps) & (n2 > 0))
	{
		pop1->SSneg[n1] = pop2.SSneg[n2];
		pop1->IIncub[n1] = pop2.IIncub[n2];
		pop1->IActive[n1] = pop2.IActive[n2];
		pop1->RecovSSneg[n1] = pop2.RecovSSneg[n2];
		pop1->ILatent[n1] = pop2.ILatent[n2];
		pop1->EarlyRecov[n1] = pop2.EarlyRecov[n2];
		pop1->LateRecov[n1] = pop2.LateRecov[n2];
		pop1->MeanNumPart[n1] = pop2.MeanNumPart[n2];
		pop1->MeanNumSexByPart[n1] = pop2.MeanNumSexByPart[n2];
		pop1->Prevalence[n1] = pop2.Prevalence[n2];
		pop1->Incidence[n1] = pop2.Incidence[n2];
		pop1->CondomUse[n1] = pop2.CondomUse[n2];
		pop1->Total[n1] = pop2.Total[n2];
		pop1->NumScreened[RPRTest][n1] = pop2.NumScreened[RPRTest][n2];
		pop1->NumScreened[TPHATest][n1] = pop2.NumScreened[TPHATest][n2];
		pop1->NumTreatedaftScreened[RPRTest][n1] = pop2.NumTreatedaftScreened[RPRTest][n2];
		pop1->NumTreatedaftScreened[TPHATest][n1] = pop2.NumTreatedaftScreened[TPHATest][n2];

		pop1->NumIncidentCases[n1] = pop2.NumIncidentCases[n2];
		pop1->NumSymptomaticClinTreat[n1] = pop2.NumSymptomaticClinTreat[n2];

		pop1->NumReferred[n1] = pop2.NumReferred[n2];
		pop1->NumReferredTestedPositive[RPRTest][n1] = pop2.NumReferredTestedPositive[RPRTest][n2];
		pop1->NumReferredTestedPositive[TPHATest][n1] = pop2.NumReferredTestedPositive[TPHATest][n2];

		pop1->NumContactTraced[all_contactTraced][n1] = pop2.NumContactTraced[all_contactTraced][n2];
		pop1->NumContactTraced[prim_and_sec_contactTraced][n1] = pop2.NumContactTraced[prim_and_sec_contactTraced][n2];
		pop1->NumContactTraced[prim_and_sec_contactTraced][n1] = pop2.NumContactTraced[prim_and_sec_contactTraced][n2];

		//referral_rate is not copied, it should be calculated!
		pop1->is_turnover_on = pop2.is_turnover_on;

		if (pop2.is_turnover_on)
		{
			pop1->turnover[n1] = pop2.turnover[n2];
		}
	}
} //copytosub_population(sub_population *pop1, int n1, sub_population pop2, int n2)

void doSyphilisCalc(SyphilisPopulation *Mypop, int DT) //Projection of Syphilis Population
{
	SyphilisPopulation interm_pop_step = SyphilisPopulation(); //temporary object for simulation by time steps. 
	Mypop->calcinitreferral();//Calculate initial force of infection. In fact, this will help grabing referral rates

	double unittime = Mypop->Men.allmen[0].dt;
	int numsteps = Mypop->Men.allmen[0].numsteps;
	int in_DT = (DT > 0) ? DT : 1;
	interm_pop_step.init(unittime, 0, in_DT);
	interm_pop_step.assortativity = Mypop->assortativity;
	interm_pop_step.w_assortativeness = Mypop->w_assortativeness;
	interm_pop_step.MRtoHR_assortativ = Mypop->MRtoHR_assortativ;

	int len_men = Mypop->Men.allmen.size();
	int len_women = Mypop->WoMen.allwomen.size();

	interm_pop_step.WoMen.allwomen[SyphHig].turnover.resize(interm_pop_step.WoMen.allwomen[SyphHig].numsteps);
	interm_pop_step.WoMen.allwomen[SyphHig].is_turnover_on = true;

	interm_pop_step.copyMWPops(0, Mypop, 0, true);
	interm_pop_step.mixingmat = Mypop->mixingmat;
	interm_pop_step.referralmat = Mypop->referralmat;

	double new_dt = 1.0 / in_DT;// converting rates for usage by time step
	for (int pop = 0; pop < len_men; pop++)//Men
	{
		interm_pop_step.Men.allmen[pop].dt = unittime / in_DT;
		interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[DT - 1] = 1.0 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[0] * new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[DT - 1] = 1.0 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[0] * new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[DT - 1] = 1.0 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[0] * new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.screen_rate[DT - 1] = 1.0 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.screen_rate[0] * new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[DT - 1] = 1.0 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[0] * new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.prop_part_referred[DT - 1] = Mypop->Men.allmen[pop].nh_param.prop_part_referred[0];

		interm_pop_step.Men.allmen[pop].nh_param.vect_nu[DT - 1] = 1.0 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.vect_nu[0] * new_dt);

		interm_pop_step.Men.allmen[pop].pinfectdetect[DT - 1] = Mypop->Men.allmen[pop].pinfectdetect[0];

		interm_pop_step.Men.allmen[pop].nh_param.sigma1 = 1 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.sigma1*new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.sigma2 = 1 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.sigma2*new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.sigma3 = 1 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.sigma3*new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.sigma4 = 1 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.sigma4*new_dt);
		interm_pop_step.Men.allmen[pop].nh_param.sigma5 = 1 - std::exp(-interm_pop_step.Men.allmen[pop].nh_param.sigma5*new_dt);
	}

	for (int pop = 0; pop < len_women; pop++)//Women
	{
		interm_pop_step.WoMen.allwomen[pop].dt = unittime / in_DT;
		interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[DT - 1] = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[0] * new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[DT - 1] = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[0] * new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[DT - 1] = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[0] * new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[DT - 1] = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[0] * new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[DT - 1] = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[0] * new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.prop_part_referred[DT - 1] = Mypop->WoMen.allwomen[pop].nh_param.prop_part_referred[0];

		interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[DT - 1] = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[0] * new_dt);
		interm_pop_step.WoMen.allwomen[pop].pinfectdetect[DT - 1] = Mypop->WoMen.allwomen[pop].pinfectdetect[0];

		interm_pop_step.WoMen.allwomen[pop].nh_param.sigma1 = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.sigma1*new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.sigma2 = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.sigma2*new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.sigma3 = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.sigma3*new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.sigma4 = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.sigma4*new_dt);
		interm_pop_step.WoMen.allwomen[pop].nh_param.sigma5 = 1.0 - std::exp(-interm_pop_step.WoMen.allwomen[pop].nh_param.sigma5*new_dt);
	}

	interm_pop_step.WoMen.allwomen[SyphHig].turnover[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[SyphHig].turnover[0] * new_dt);

	int ti = 1;
	while (TRUE)//Simulation of the whole epidemic
	{
		if (ti >= numsteps) break;
		interm_pop_step.copyMWPops(0, Mypop, ti - 1, false);//copying the parameters in Mypop for ti-1 into the temporary object interm_pop_step 

		for (int pop = 0; pop < len_men; pop++) //Converting rates to the appropriate time step, for men
		{
			interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[0] = interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[DT - 1];
			interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[0] = interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[DT - 1];
			interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[0] = interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[DT - 1];
			interm_pop_step.Men.allmen[pop].nh_param.screen_rate[0] = interm_pop_step.Men.allmen[pop].nh_param.screen_rate[DT - 1];
			interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[0] = interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[DT - 1];
			interm_pop_step.Men.allmen[pop].nh_param.vect_nu[0] = interm_pop_step.Men.allmen[pop].nh_param.vect_nu[DT - 1];
			interm_pop_step.Men.allmen[pop].nh_param.prop_part_referred[0] = interm_pop_step.Men.allmen[pop].nh_param.prop_part_referred[DT - 1];
			
			interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[DT - 1] = 1.0 - std::exp(-Mypop->Men.allmen[pop].nh_param.r_ActiveToSus[ti] * new_dt);
			interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[DT - 1] = 1.0 - std::exp(-Mypop->Men.allmen[pop].nh_param.r_ActiveRecov[ti] * new_dt);
			interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[DT - 1] = 1.0 - std::exp(-Mypop->Men.allmen[pop].nh_param.r_LatentRecov[ti] * new_dt);
			interm_pop_step.Men.allmen[pop].nh_param.screen_rate[DT - 1] = 1.0 - std::exp(-Mypop->Men.allmen[pop].nh_param.screen_rate[ti] * new_dt);
			interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[DT - 1] = 1.0 - std::exp(-Mypop->Men.allmen[pop].nh_param.screen_rateTPHA[ti] * new_dt);
			interm_pop_step.Men.allmen[pop].nh_param.vect_nu[DT - 1] = 1.0 - std::exp(-Mypop->Men.allmen[pop].nh_param.vect_nu[ti] * new_dt);
			interm_pop_step.Men.allmen[pop].nh_param.prop_part_referred[DT - 1] = Mypop->Men.allmen[pop].nh_param.prop_part_referred[ti];
			interm_pop_step.Men.allmen[pop].pinfectdetect[DT - 1] = Mypop->Men.allmen[pop].pinfectdetect[ti];
		}

		for (int pop = 0; pop < len_women; pop++) //Converting rates to the appropriate time step, for women
		{
			interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[DT - 1];;
			interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[DT - 1];
			interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[DT - 1];
			interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[DT - 1];
			interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[DT - 1];
			interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[DT - 1];
			interm_pop_step.WoMen.allwomen[pop].nh_param.prop_part_referred[0] = interm_pop_step.WoMen.allwomen[pop].nh_param.prop_part_referred[DT - 1];

			interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[pop].nh_param.r_ActiveToSus[ti] * new_dt);
			interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[pop].nh_param.r_ActiveRecov[ti] * new_dt);
			interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[pop].nh_param.r_LatentRecov[ti] * new_dt);
			interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[pop].nh_param.screen_rate[ti] * new_dt);
			interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[pop].nh_param.screen_rateTPHA[ti] * new_dt);
			interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[pop].nh_param.vect_nu[ti] * new_dt);
			interm_pop_step.WoMen.allwomen[pop].nh_param.prop_part_referred[DT - 1] = Mypop->WoMen.allwomen[pop].nh_param.prop_part_referred[ti];

			interm_pop_step.WoMen.allwomen[pop].pinfectdetect[DT - 1] = Mypop->WoMen.allwomen[pop].pinfectdetect[ti];
		}

		interm_pop_step.WoMen.allwomen[SyphHig].turnover[0] = interm_pop_step.WoMen.allwomen[SyphHig].turnover[DT - 1];
		interm_pop_step.WoMen.allwomen[SyphHig].turnover[DT - 1] = 1.0 - std::exp(-Mypop->WoMen.allwomen[SyphHig].turnover[ti] * new_dt);//converting turnover rates

		for (int tj = 0; tj <= in_DT; tj++)
		{
			//Men//
			interm_pop_step.Men.total[tj] = Mypop->Men.total[ti];
			for (int pop = 0; pop < len_men; pop++)
			{
				//Smoothing is being disabled
				//Copying the data for usage by time step
				interm_pop_step.Men.allmen[pop].Total[tj] = Mypop->Men.allmen[pop].Total[ti];
				interm_pop_step.Men.allmen[pop].MeanNumPart[tj] = Mypop->Men.allmen[pop].MeanNumPart[ti];
				interm_pop_step.Men.allmen[pop].MeanNumSexByPart[tj] = Mypop->Men.allmen[pop].MeanNumSexByPart[ti];
				interm_pop_step.Men.allmen[pop].CondomUse[tj] = Mypop->Men.allmen[pop].CondomUse[ti];

				interm_pop_step.Men.allmen[pop].nh_param.screen_rate[tj] = interm_pop_step.Men.allmen[pop].nh_param.screen_rate[DT - 1];
				interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[tj] = interm_pop_step.Men.allmen[pop].nh_param.screen_rateTPHA[DT - 1];
				interm_pop_step.Men.allmen[pop].nh_param.prop_part_referred[tj] = interm_pop_step.Men.allmen[pop].nh_param.prop_part_referred[DT - 1];
				interm_pop_step.Men.allmen[pop].nh_param.vect_nu[tj] = interm_pop_step.Men.allmen[pop].nh_param.vect_nu[DT - 1];

				interm_pop_step.Men.allmen[pop].pinfectdetect[tj] = interm_pop_step.Men.allmen[pop].pinfectdetect[DT - 1];

				interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[tj] = interm_pop_step.Men.allmen[pop].nh_param.r_ActiveToSus[DT - 1];
				interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[tj] = interm_pop_step.Men.allmen[pop].nh_param.r_ActiveRecov[DT - 1];
				interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[tj] = interm_pop_step.Men.allmen[pop].nh_param.r_LatentRecov[DT - 1];
			}
			//Women
			interm_pop_step.WoMen.total[tj] = Mypop->WoMen.total[ti];
			for (int pop = 0; pop < len_women; pop++)
			{
				//Smoothing is being disabled.
				//Copying the data for usage by time step
				interm_pop_step.WoMen.allwomen[pop].Total[tj] = Mypop->WoMen.allwomen[pop].Total[ti];
				interm_pop_step.WoMen.allwomen[pop].MeanNumPart[tj] = Mypop->WoMen.allwomen[pop].MeanNumPart[ti];
				interm_pop_step.WoMen.allwomen[pop].MeanNumSexByPart[tj] = Mypop->WoMen.allwomen[pop].MeanNumSexByPart[ti];
				interm_pop_step.WoMen.allwomen[pop].CondomUse[tj] = Mypop->WoMen.allwomen[pop].CondomUse[ti];

				interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rate[DT - 1];
				interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.screen_rateTPHA[DT - 1];
				interm_pop_step.WoMen.allwomen[pop].nh_param.prop_part_referred[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.prop_part_referred[DT - 1];
				interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.vect_nu[DT - 1];
				
				interm_pop_step.WoMen.allwomen[pop].pinfectdetect[tj] = interm_pop_step.WoMen.allwomen[pop].pinfectdetect[DT - 1];

				interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveToSus[DT - 1];
				interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.r_ActiveRecov[DT - 1];
				interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[tj] = interm_pop_step.WoMen.allwomen[pop].nh_param.r_LatentRecov[DT - 1];
			}
		}

		for (int tj = 0; tj <= in_DT; tj++)
		{
			interm_pop_step.calcFOI(tj);//Update the force of infection for each time step
			interm_pop_step.ProgressPopulation(tj);//Progress the population/update the state of the population, given the force of infection 
		}

		Mypop->copyMWPops(ti, interm_pop_step, in_DT - 1, false);
		//Get the total numbers of screened, treated, contacts traced and treated individuals together with the number of incident cases
		for (int pop = 0; pop < len_men; pop++) //Men
		{
			for (int ntes = 0; ntes < numTests; ntes++)
			{
				Mypop->Men.allmen[pop].NumScreened[ntes][ti] = 0;
				Mypop->Men.allmen[pop].NumTreatedaftScreened[ntes][ti] = 0;
				Mypop->Men.allmen[pop].NumReferredTestedPositive[ntes][ti] = 0;
				for (int tj = 0; tj < in_DT; tj++)
				{
					Mypop->Men.allmen[pop].NumScreened[ntes][ti] += interm_pop_step.Men.allmen[pop].NumScreened[ntes][tj];
					Mypop->Men.allmen[pop].NumTreatedaftScreened[ntes][ti] += interm_pop_step.Men.allmen[pop].NumTreatedaftScreened[ntes][tj];
					Mypop->Men.allmen[pop].NumReferredTestedPositive[ntes][ti] += interm_pop_step.Men.allmen[pop].NumReferredTestedPositive[ntes][tj];
				}
			}

			for (int ntrac = 0; ntrac < sizeContactTraced; ntrac++)
			{
				Mypop->Men.allmen[pop].NumContactTraced[ntrac][ti] = 0;
				for (int tj = 0; tj < in_DT; tj++)
				{
					Mypop->Men.allmen[pop].NumContactTraced[ntrac][ti] += interm_pop_step.Men.allmen[pop].NumContactTraced[ntrac][tj];
				}
			}

			Mypop->Men.allmen[pop].NumIncidentCases[ti] = 0;
			Mypop->Men.allmen[pop].NumSymptomaticClinTreat[ti] = 0;
			Mypop->Men.allmen[pop].NumReferred[ti] = 0;

			for (int tj = 0; tj < in_DT; tj++)
			{
				Mypop->Men.allmen[pop].NumIncidentCases[ti] += interm_pop_step.Men.allmen[pop].NumIncidentCases[tj];
				Mypop->Men.allmen[pop].NumSymptomaticClinTreat[ti] += interm_pop_step.Men.allmen[pop].NumSymptomaticClinTreat[tj];
				Mypop->Men.allmen[pop].NumReferred[ti] += interm_pop_step.Men.allmen[pop].NumReferred[tj];
			}
		}

		for (int pop = 0; pop < len_women; pop++)//Women
		{
			for (int ntes = 0; ntes < numTests; ntes++)
			{
				Mypop->WoMen.allwomen[pop].NumScreened[ntes][ti] = 0;
				Mypop->WoMen.allwomen[pop].NumTreatedaftScreened[ntes][ti] = 0;
				Mypop->WoMen.allwomen[pop].NumReferredTestedPositive[ntes][ti] = 0;
				for (int tj = 0; tj < in_DT; tj++)
				{
					Mypop->WoMen.allwomen[pop].NumScreened[ntes][ti] += interm_pop_step.WoMen.allwomen[pop].NumScreened[ntes][tj];
					Mypop->WoMen.allwomen[pop].NumTreatedaftScreened[ntes][ti] += interm_pop_step.WoMen.allwomen[pop].NumTreatedaftScreened[ntes][tj];
					Mypop->WoMen.allwomen[pop].NumReferredTestedPositive[ntes][ti] += interm_pop_step.WoMen.allwomen[pop].NumReferredTestedPositive[ntes][tj];
				}
			}

			for (int ntrac = 0; ntrac < sizeContactTraced; ntrac++)
			{
				Mypop->WoMen.allwomen[pop].NumContactTraced[ntrac][ti] = 0;
				for (int tj = 0; tj < in_DT; tj++)
				{
					Mypop->WoMen.allwomen[pop].NumContactTraced[ntrac][ti] += interm_pop_step.WoMen.allwomen[pop].NumContactTraced[ntrac][tj];
				}
			}
			//Rcpp::Rcout << std::endl << "END Num Contact Traced" << std::endl;

			Mypop->WoMen.allwomen[pop].NumIncidentCases[ti] = 0;
			Mypop->WoMen.allwomen[pop].NumSymptomaticClinTreat[ti] = 0;
			Mypop->WoMen.allwomen[pop].NumReferred[ti] = 0;
			for (int tj = 0; tj < in_DT; tj++)
			{
				Mypop->WoMen.allwomen[pop].NumIncidentCases[ti] += interm_pop_step.WoMen.allwomen[pop].NumIncidentCases[tj];
				Mypop->WoMen.allwomen[pop].NumSymptomaticClinTreat[ti] += interm_pop_step.WoMen.allwomen[pop].NumSymptomaticClinTreat[tj];
				Mypop->WoMen.allwomen[pop].NumReferred[ti] += interm_pop_step.WoMen.allwomen[pop].NumReferred[tj];
			}
		}

		if (ti == 1)
		{
			Mypop->all_mixing_t1 = interm_pop_step.all_mixing_tf;
		}

		//break;

		//Rcpp::Rcout << "ti=" << ti << std::endl;
		ti++;
	}

	Mypop->all_mixing_tf = interm_pop_step.all_mixing_tf;

	for (int pop = 0; pop < len_men; pop++)
	{
		Mypop->Men.allmen[pop].calcPrevalence();
	}

	for (int pop = 0; pop < len_women; pop++)
	{
		Mypop->WoMen.allwomen[pop].calcPrevalence();
	}
}//doSyphilisCalc(SyphilisPopulation *Mypop, const NaturalHistParam param, int DT); run Syphilis Model on a population given the natural history and the time step

/*Sub Populations*/
void sub_population::init(int t_final, int t_init, int in_dt)//Initialize the population
{
	dt = (in_dt <= 0) ? 1.0 / 12.0 : 1.0/ in_dt;
	tinit = (t_init < 0) ? 0 : t_init;
	tfinal = (t_final < (t_init + dt)) ? (t_init + dt) : t_final;
	double currt = tinit;
	numsteps = 0;
	poptime.resize(0);
	while (TRUE)
	{
		numsteps++;
		poptime.push_back(currt);
		currt += dt;
		if (currt > tfinal) break;
	}

	SSneg.resize(boost::extents[numsteps]);
	IIncub.resize(boost::extents[numsteps]);
	ILatent.resize(boost::extents[numsteps]);
	EarlyRecov.resize(boost::extents[numsteps]);
	LateRecov.resize(boost::extents[numsteps]);
	Total.resize(boost::extents[numsteps]);
	for (int ii = 0; ii < numsteps; ii++)
	{
		Total[ii] = 1;
	}

	MeanNumPart.resize(boost::extents[numsteps]);
	MeanNumSexByPart.resize(boost::extents[numsteps]);
	IActive.resize(boost::extents[numsteps]);
	RecovSSneg.resize(boost::extents[numsteps]);

	Incidence.resize(boost::extents[numsteps]);
	Prevalence.resize(boost::extents[numsteps]);
	Prevalence2d.resize(boost::extents[numPrevTypes][numsteps]);
	CondomUse.resize(boost::extents[numsteps]);
	NumScreened.resize(boost::extents[numTests][numsteps]); 
	NumTreatedaftScreened.resize(boost::extents[numTests][numsteps]);

	NumIncidentCases.resize(boost::extents[numsteps]);
	NumSymptomaticClinTreat.resize(boost::extents[numsteps]);

	NumReferred.resize(boost::extents[numsteps]);
	NumReferredTestedPositive.resize(boost::extents[numTests][numsteps]);

	NumContactTraced.resize(boost::extents[sizeContactTraced][numsteps]);

	referral_rate.resize(numsteps);
	std::fill(referral_rate.begin(), referral_rate.end(), 0);

	pinfectdetect.resize(numsteps);
	std::fill(pinfectdetect.begin(), pinfectdetect.end(), 0);

	for (int ii = 0; ii < numsteps; ii++)
	{
		CondomUse[ii] = 0;
	}

	is_turnover_on = false;

	nh_param = NaturalHistParam();
	nh_param.resize(numsteps);
}

sub_population::sub_population()
{
	init(1, 0, 12);
}

sub_population::~sub_population()
{
}

void sub_population::calcPrevalence() //Calculate the prevalence
{
	Prevalence.resize(boost::extents[numsteps]);
	for (int ii = 0; ii < numsteps; ii++)
	{
		Prevalence[ii] = (Total[ii] > 0) ? (Total[ii] - SSneg[ii]) / Total[ii] : 0;
	}
	Prevalence2d.resize(boost::extents[numPrevTypes][numsteps]);
	for (int ii = 0; ii < numsteps; ii++)
	{
		Prevalence2d[RPRnegTPHApos][ii] = (Total[ii] > 0) ? (RecovSSneg[ii])/ Total[ii] : 0;
		Prevalence2d[RPRposTPHApos][ii] = (Total[ii] > 0) ? (IActive[ii]+ EarlyRecov[ii]+ LateRecov[ii]+ILatent[ii]) / Total[ii] : 0;

		Prevalence2d[TPHApos][ii] = (Total[ii] > 0) ? (IActive[ii] + EarlyRecov[ii] + LateRecov[ii] + ILatent[ii]+ RecovSSneg[ii]) / Total[ii] : 0;
		Prevalence2d[RPRpos][ii] = (Total[ii] > 0) ? (IActive[ii] + EarlyRecov[ii] + LateRecov[ii] + ILatent[ii]) / Total[ii] : 0;
	}
}//end sub_population::calcPrevalence

void sub_population::SeedSteadyState()//Get initial state of the population, assuming it is as equilibrium
{
	nh_param.ValinitInci = (nh_param.ValinitInci < 0) ? 0 : nh_param.ValinitInci;
	nh_param.CalcSteadyState();

	if (numsteps >= 1)
	{
		SSneg[0] = nh_param.SteadyState[SSnegative];
		IIncub[0]= nh_param.SteadyState[IIncubates];
		IActive[0]= nh_param.SteadyState[IActiveSSy];
		EarlyRecov[0] = nh_param.SteadyState[RecovEarly];
		ILatent[0] = nh_param.SteadyState[ILatentSSy];
		LateRecov[0] = nh_param.SteadyState[RecovLater];
		RecovSSneg[0] = nh_param.SteadyState[RecovSusce];
		if (Total[0] > 0)
		{
			SSneg[0] *= Total[0];
			IIncub[0] *= Total[0];
			IActive[0] *= Total[0];
			EarlyRecov[0] *= Total[0];
			ILatent[0] *= Total[0];
			LateRecov[0] *= Total[0];
			RecovSSneg[0] *= Total[0];
			Prevalence[0] = IActive[0]/ Total[0];
		}
		Incidence[0] = nh_param.ValinitInci;

		NumScreened[RPRTest][0] = Total[0] * nh_param.screen_rate[0];//We should exclude those who have been treated recently
		NumScreened[TPHATest][0] = Total[0] * nh_param.screen_rateTPHA[0];//We should exclude those who have been treated recently  
		NumTreatedaftScreened[RPRTest][0] = (IActive[0] + EarlyRecov[0] + ILatent[0] + LateRecov[0]) * nh_param.screen_rate[0];//We should exclude those who have been treated recently
		NumTreatedaftScreened[TPHATest][0] = (IActive[0] + EarlyRecov[0] + ILatent[0] + LateRecov[0] + RecovSSneg[0]) * nh_param.screen_rateTPHA[0];//We should exclude those who have been treated recently  
	
		NumIncidentCases[0]= (SSneg[0]+ RecovSSneg[0])*nh_param.ValinitInci;
		NumSymptomaticClinTreat[0] = IActive[0] * nh_param.vect_nu[0];// *nh_param.psi*nh_param.phi;
	}
}//end sub_population::SeedSteadyState()

void sub_population::ProgressPopulation(int n, const double foi)
{
	//Nota: IN THE FOLLOWING, IT IS ASSUMED THAT THE RATES HAVE BEEN CONVERTED TO THE CORRECT TIME UNIT
	int nextime = n + 1;
	if ((nextime < numsteps) &(n >= 0))
	{
		double inc = foi * dt;//To do: try foi*dt
		double newinfection = SSneg[n] * inc;//Newly infected
		double re_infections = RecovSSneg[n] * inc; //Get re-infected after previous exposure and recovery
		double enter_active = IIncub[n] * nh_param.sigma1; //Enter the active phase after the incupation period
		//New, to accommodate contact tracing
		double enter_earlyrecov1 = IActive[n] * (nh_param.r_ActiveRecov[n] + nh_param.psi*(1-nh_param.phi)*referral_rate[n]); //Recoved after treatment ("failure"?)
		double enter_recovsus1 = IActive[n] * (nh_param.r_ActiveToSus[n] + nh_param.psi*nh_param.phi*referral_rate[n]);//Become susceptible after successful treatment
		double enter_latent = IActive[n] * nh_param.sigma2;//Progress to latent phase after being active
		
		//New, to accomodate those in latent phase who are treated
		double enter_earlyrecov2 = ILatent[n] * (nh_param.r_LatentRecov[n] + nh_param.psi*referral_rate[n]);//Recovered from latent compartment, after treatment
		
		double enter_laterecov = ILatent[n] * nh_param.sigma3; // Enter recovery compartment after latent the phase
		double enter_recovsus2 = EarlyRecov[n] * nh_param.sigma4;//Re-enter susceptible population after recovery
		double enter_recovsus3 = LateRecov[n] * nh_param.sigma5;// Recovered after latent phase and no treatment
		double simtotal = 0;

		NumReferred[n] = Total[n]*referral_rate[n];
		NumReferredTestedPositive[RPRTest][n] = (IActive[n] + ILatent[n])*referral_rate[n];// Number referred;
		NumReferredTestedPositive[TPHATest][n] = (IActive[n] + EarlyRecov[n] + ILatent[n] + LateRecov[n] + RecovSSneg[n])*referral_rate[n];// Number referred;

		NumScreened[RPRTest][n] = Total[n]*nh_param.screen_rate[n];//Number sceened with RPR
		NumScreened[TPHATest][n] = Total[n]*nh_param.screen_rateTPHA[n];//Number sceened with TPHA 

		NumTreatedaftScreened[RPRTest][n] = (IActive[n] + EarlyRecov[n] + ILatent[n] + LateRecov[n]) * nh_param.screen_rate[n];//Number Treated after screened with RPR diagnosis 
		NumTreatedaftScreened[TPHATest][n] = (IActive[n] + EarlyRecov[n]+ILatent[n] + LateRecov[n] + RecovSSneg[n]) * nh_param.screen_rateTPHA[n];// Number Treated after screened with TPHA diagnosis
		
		NumIncidentCases[n] = newinfection + re_infections;//Store the number of new cases
		NumSymptomaticClinTreat[n] = IActive[n] * nh_param.vect_nu[n];// Store the number that are treated

		/*Save incidence rate*/
		Incidence[nextime] = foi;

		/*Susceptibles*/
		double interm = SSneg[n] - newinfection;
		SSneg[nextime] = interm;
		simtotal += interm;

		double ageing = 1.0-1.0 / 35 * dt;

		/*Incubation Period*/
		interm = IIncub[n] * ageing + newinfection+ re_infections - enter_active;
		IIncub[nextime] = interm;
		simtotal += interm;

		/*Active Syphilis*/
		interm = IActive[n]* ageing + enter_active - enter_latent - enter_recovsus1- enter_earlyrecov1;
		IActive[nextime] = interm;
		simtotal += interm;

		/*Latent*/
		interm = ILatent[n] * ageing + enter_latent - enter_laterecov- enter_earlyrecov2;
		ILatent[nextime] = interm;
		simtotal += interm;

		/*Early recovery:Recovery after treatment*/
		interm = EarlyRecov[n] * ageing + enter_earlyrecov1+ enter_earlyrecov2- enter_recovsus2;
		EarlyRecov[nextime] = interm;
		simtotal += interm;

		/*Recovery after latent phase*/
		interm = LateRecov[n] * ageing + enter_laterecov - enter_recovsus3;
		LateRecov[nextime] = interm;
		simtotal += interm;

		/*Susceptible after recovery*/
		interm = RecovSSneg[n]* ageing + enter_recovsus1 + enter_recovsus2 + enter_recovsus3-re_infections;
		RecovSSneg[nextime] = interm;
		simtotal += interm;

		if (simtotal > 0)
		{
			double scaleit = (Total[nextime] > simtotal) ? 1: Total[nextime] / simtotal;// 
			double new_sus = (Total[nextime] > simtotal) ? Total[nextime] - simtotal:0;

			SSneg[nextime] = SSneg[nextime] * scaleit+ new_sus; //
			IIncub[nextime] = IIncub[nextime] * scaleit;
			IActive[nextime] = IActive[nextime] * scaleit;
			ILatent[nextime] = ILatent[nextime] * scaleit;
			EarlyRecov[nextime] = EarlyRecov[nextime] * scaleit;
			LateRecov[nextime] = LateRecov[nextime] * scaleit;
			RecovSSneg[nextime] = RecovSSneg[nextime] * scaleit;
		}
		Prevalence[nextime] = (Total[nextime] > 0) ? IActive[nextime] / Total[nextime] : 0;
	}
} // end sub_population::ProgressPopulation(int n, const double foi)

void sub_population::copyfrom_sub_population(int n1, sub_population subpop, int n2, bool cpnath) // copy a slice of subpop into a slice of the subpopulation simulated
{
	if ((n1 < this->numsteps) & (n1 >= 0) &(n2 < subpop.numsteps) & (n2 >= 0))
	{
		this->SSneg[n1] = subpop.SSneg[n2];
		this->IIncub[n1] = subpop.IIncub[n2];
		this->IActive[n1] = subpop.IActive[n2];
		this->RecovSSneg[n1] = subpop.RecovSSneg[n2];
		this->ILatent[n1] = subpop.ILatent[n2];
		this->EarlyRecov[n1] = subpop.EarlyRecov[n2];
		this->LateRecov[n1] = subpop.LateRecov[n2];
		this->MeanNumPart[n1] = subpop.MeanNumPart[n2];
		this->MeanNumSexByPart[n1] = subpop.MeanNumSexByPart[n2];
		this->Prevalence[n1] = subpop.Prevalence[n2];
		this->Incidence[n1] = subpop.Incidence[n2];
		this->CondomUse[n1] = subpop.CondomUse[n2];
		this->Total[n1] = subpop.Total[n2];
		this->NumScreened[RPRTest][n1] = subpop.NumScreened[RPRTest][n2];
		this->NumScreened[TPHATest][n1] = subpop.NumScreened[TPHATest][n2];

		this->NumTreatedaftScreened[RPRTest][n1] = subpop.NumTreatedaftScreened[RPRTest][n2];
		this->NumTreatedaftScreened[TPHATest][n1] = subpop.NumTreatedaftScreened[TPHATest][n2];

		this->NumIncidentCases[n1] = subpop.NumIncidentCases[n2];
		this->NumSymptomaticClinTreat[n1] = subpop.NumSymptomaticClinTreat[n2];

		this->NumContactTraced[all_contactTraced][n1] = subpop.NumContactTraced[all_contactTraced][n2];
		this->NumContactTraced[prim_and_sec_contactTraced][n1] = subpop.NumContactTraced[prim_and_sec_contactTraced][n2];
		this->NumContactTraced[latent_contactTraced][n1] = subpop.NumContactTraced[latent_contactTraced][n2];

		this->is_turnover_on = subpop.is_turnover_on;

		if (subpop.is_turnover_on)
		{
			this->turnover[n1] = subpop.turnover[n2];
		}

		//Natural History and Screening rates
		if (cpnath)
		{
			this->nh_param.sigma1 = subpop.nh_param.sigma1;
			this->nh_param.sigma2 = subpop.nh_param.sigma2;
			this->nh_param.sigma3 = subpop.nh_param.sigma3;
			this->nh_param.sigma4 = subpop.nh_param.sigma4;
			this->nh_param.sigma5 = subpop.nh_param.sigma5;
			this->nh_param.nu = subpop.nh_param.nu;
			this->nh_param.phi = subpop.nh_param.phi;
			this->nh_param.psi = subpop.nh_param.psi;
			this->nh_param.MtoF_tp = subpop.nh_param.MtoF_tp;
			this->nh_param.FtoM_tp = subpop.nh_param.FtoM_tp;
			this->nh_param.MtoM_tp = subpop.nh_param.MtoM_tp;

			//I should check this, I do not necessarily need to do this because the time steps will generally be different
			this->nh_param.screen_rate[n1] = subpop.nh_param.screen_rate[n2];
			this->nh_param.r_ActiveToSus[n1] = subpop.nh_param.r_ActiveToSus[n2];
			this->nh_param.r_ActiveRecov[n1] = subpop.nh_param.r_ActiveRecov[n2];
			this->nh_param.r_LatentRecov[n1] = subpop.nh_param.r_LatentRecov[n2];
			this->nh_param.vect_nu[n1] = subpop.nh_param.vect_nu[n2];
		}
	}
} //end sub_population::copyfrom_sub_population(int n1, sub_population subpop, int n2)

/*Men*/
void Men_Pop::init(int t_final, int t_init, int in_dt)
{
	allmen.resize(SyphNumRisksMen);
	for (int ii = SyphLow; ii < SyphNumRisksMen; ii++)
	{
		allmen[ii].init(t_final, t_init, in_dt);
	}
	
	distributions.resize(boost::extents[SyphNumRisksMen][allmen[0].numsteps]);
	total.resize(allmen[0].numsteps);

	for (int ii = 0; ii < allmen[0].numsteps; ii++)
	{
		total[ii] = 100;
		distributions[SyphLow][ii] = 90;
		distributions[SyphMed][ii] = 5;
		distributions[SyphHig][ii] = 2;
		distributions[SyphMSM][ii] = 3;
		distributions[SyphNoSex][ii] = 0;
		double scale = 0;
		for (int pop = SyphLow; pop < SyphNumRisksMen; pop++)
		{
			scale += distributions[pop][ii];
		}
		for (int pop = SyphLow; pop < SyphNumRisksMen; pop++)
		{
			allmen[pop].Total[ii] = distributions[pop][ii] / scale * total[ii];
		}
	}
} //end Men_Pop::init(int t_final, int t_init, int in_dt)

void Men_Pop::SetSubPop(Array2d in_dist, Array1d in_total, Array1d in_time)
{
	//TO DO: Sort the time
	for (int ii = 0; ii < allmen[0].numsteps; ii++)
	{
		double curr_time = allmen[0].poptime[ii];
		int ind_min=0;
		int ind_max = in_time.size();

		for (int jj = 0; jj < ind_max; jj++)
		{
			ind_min = jj;
			if (curr_time < in_time[jj])
			{
				break;
			}
		}

		ind_min = (ind_min >= 1) ? (ind_min - 1) : ind_min;
		ind_max--;
		ind_max = ((ind_min + 1) < ind_max) ? (ind_min + 1) : ind_max;
		if ((ind_max == ind_min) | (in_time[ind_max] == in_total[ind_min]))
		{
			total[ii] = in_total[ind_min];
		}
		else
		{
			total[ii] = in_total[ind_max] * (curr_time - in_time[ind_min]) + in_total[ind_min] * (in_time[ind_max] - curr_time);
			total[ii] /= (in_time[ind_max] - in_time[ind_min]);
		}

		//TO DO: CHeck if this should be smoothed out
		distributions[SyphLow][ii] = in_dist[SyphLow][ind_min];
		distributions[SyphMed][ii] = in_dist[SyphMed][ind_min];
		distributions[SyphHig][ii] = in_dist[SyphHig][ind_min];
		distributions[SyphMSM][ii] = in_dist[SyphMSM][ind_min];
		distributions[SyphNoSex][ii] = in_dist[SyphNoSex][ind_min];
		double scale = 0;
		for (int pop = SyphLow; pop < SyphNumRisksMen; pop++)
		{
			scale += distributions[pop][ii];
		}
		for (int pop = SyphLow; pop < SyphNumRisksMen; pop++)
		{
			allmen[pop].Total[ii] = distributions[pop][ii] / scale * total[ii];
		}
	}
} //end Men_Pop::SetSubPop(Array2d in_dist, Array1d in_total, Array1d in_time)?

void Men_Pop::ProgressPopulation(int n, const std::vector<double> vect_foi)
{
	int len_men = allmen.size();
	for (int pop = 0; pop< len_men; pop++)
	{
		allmen[pop].ProgressPopulation(n, vect_foi[pop]);
	}
} //end Men_Pop::ProgressPopulation(int n, const std::vector<double> vect_foi)

void Men_Pop::copyfromMen(int n1, Men_Pop popmen2, int n2, bool cpnath)
{
	int len_c = allmen.size();
	int len_out = popmen2.allmen.size();
	if (len_c == len_out)
	{
		for (int ispop = 0; ispop < len_c; ispop++)
		{
			allmen[ispop].copyfrom_sub_population(n1, popmen2.allmen[ispop], n2, cpnath);
		}
	}
} //end Men_Pop::copyfromMen(int n1, Men_Pop popmen2, int n2)

Men_Pop::Men_Pop()
{
	init(1, 0, 12);
} //end Men_Pop::Men_Pop() constructor

Men_Pop::~Men_Pop()
{
	;
} //end Men_Pop::~Men_Pop(); destructor

/*WoMen*/
void WoMen_Pop::init(int t_final, int t_init, int in_dt)
{
	allwomen.resize(SyphNumRisksWom);
	for (int ii = SyphLow; ii < SyphNumRisksWom; ii++)
	{
		allwomen[ii].init(t_final, t_init, in_dt);
	}

	distributions.resize(boost::extents[SyphNumRisksWom][allwomen[0].numsteps]);
	total.resize(allwomen[0].numsteps);
	for (int ii = 0; ii < allwomen[0].numsteps; ii++)
	{
		total[ii] = 100;
		distributions[SyphLow][ii] = 90;
		distributions[SyphMed][ii] = 8;
		distributions[SyphHig][ii] = 2;
		distributions[SyphNoSexWom][ii] = 0;
		double scale = 0;
		for (int pop = SyphLow; pop < SyphNumRisksWom; pop++)
		{
			scale += distributions[pop][ii];
		}
		for (int pop = SyphLow; pop < SyphNumRisksWom; pop++)
		{
			allwomen[pop].Total[ii] = distributions[pop][ii] / scale * total[ii];
		}
	}
} //end WoMen_Pop::init(int t_final, int t_init, int in_dt)

void WoMen_Pop::SetSubPop(Array2d in_dist, Array1d in_total, Array1d in_time)
{
	//TO DO: Sort the time ?? It can be done in R
	for (int ii = 0; ii < allwomen[0].numsteps; ii++)
	{
		double curr_time = allwomen[0].poptime[ii];
		int ind_min=0;
		int ind_max = in_time.size();

		for (int jj = 0; jj < ind_max; jj++)
		{
			ind_min = jj;
			if (curr_time < in_time[jj])
			{
				break;
			}
		}

		ind_min = (ind_min >= 1) ? (ind_min - 1) : ind_min;
		ind_max--;
		ind_max = ((ind_min + 1) < ind_max) ? (ind_min + 1) : ind_max;
		if ((ind_max == ind_min) | (in_time[ind_max] == in_total[ind_min]))
		{
			total[ii] = in_total[ind_min];
		}
		else
		{
			total[ii] = in_total[ind_max] * (curr_time - in_time[ind_min]) + in_total[ind_min] * (in_time[ind_max] - curr_time);
			total[ii] /= (in_time[ind_max] - in_time[ind_min]);
		}

		//Copying population distribution
		distributions[SyphLow][ii] = in_dist[SyphLow][ind_min];
		distributions[SyphMed][ii] = in_dist[SyphMed][ind_min];
		distributions[SyphHig][ii] = in_dist[SyphHig][ind_min];
		distributions[SyphNoSexWom][ii] = in_dist[SyphNoSexWom][ind_min];
		double scale = 0;
		for (int pop = SyphLow; pop < SyphNumRisksWom; pop++)
		{
			scale += distributions[pop][ii];
		}
		for (int pop = SyphLow; pop < SyphNumRisksWom; pop++)
		{
			allwomen[pop].Total[ii] = distributions[pop][ii] / scale * total[ii];
		}
	}
} // end WoMen_Pop::SetSubPop(Array2d in_dist, Array1d in_total, Array1d in_time)

void WoMen_Pop::turnFSW_Over(int n) //Turnover for FSW
{
	if (allwomen[SyphHig].is_turnover_on)
	{
		double interm = allwomen[SyphHig].SSneg[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].SSneg[n] += interm;
		allwomen[SyphHig].SSneg[n] -= interm;

		interm = allwomen[SyphHig].IIncub[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].IIncub[n] += interm;
		allwomen[SyphHig].IIncub[n] -= interm;

		interm = allwomen[SyphHig].IActive[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].IActive[n] += interm;
		allwomen[SyphHig].IActive[n] -= interm;

		interm = allwomen[SyphHig].EarlyRecov[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].EarlyRecov[n] += interm;
		allwomen[SyphHig].EarlyRecov[n] -= interm;

		interm = allwomen[SyphHig].ILatent[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].ILatent[n] += interm;
		allwomen[SyphHig].ILatent[n] -= interm;

		interm = allwomen[SyphHig].LateRecov[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].LateRecov[n] += interm;
		allwomen[SyphHig].LateRecov[n] -= interm;

		interm = allwomen[SyphHig].RecovSSneg[n] * allwomen[SyphHig].turnover[n];
		allwomen[SyphMed].RecovSSneg[n] += interm;
		allwomen[SyphHig].RecovSSneg[n] -= interm;
	}
}

void WoMen_Pop::ProgressPopulation(int n, const std::vector<double> vect_foi)
{
	int len_women = allwomen.size();
	turnFSW_Over(n);//
	for (int pop = 0; pop< len_women; pop++) //Loop over sub-populations
	{
		allwomen[pop].ProgressPopulation(n, vect_foi[pop]);
	}
} //end WoMen_Pop::ProgressPopulation(int n, const std::vector<double> vect_foi)

void WoMen_Pop::copyfromWoMen(int n1, WoMen_Pop popwomen2, int n2, bool cpnath)
{
	int len_c = allwomen.size();
	int len_out = popwomen2.allwomen.size();
	if (len_c == len_out)
	{
		for (int ispop = 0; ispop < len_c; ispop++)
		{
			allwomen[ispop].copyfrom_sub_population(n1, popwomen2.allwomen[ispop], n2, cpnath);
		}
	}
} // end WoMen_Pop::copyfromWoMen(int n1, WoMen_Pop popwomen2, int n2)

WoMen_Pop::WoMen_Pop()
{
	init(1, 0, 12);
} //end WoMen_Pop::WoMen_Pop() constructor

WoMen_Pop::~WoMen_Pop()
{
	;
} //end WoMen_Pop::~WoMen_Pop() destructor

///////////////////////
/*Syphilis Population*/
void SyphilisPopulation::init(int t_final, int t_init, int in_dt)
{
	Men.init(t_final, t_init,  in_dt);//Initialize men risk groups
	WoMen.init(t_final, t_init, in_dt);//Initialize women risk groups
	FOIMen.resize(Men.allmen.size()); //Initialize FOI vector for men
	FOIWoMen.resize(WoMen.allwomen.size());//Initialize FOI vector for women
	mixingmat = MixingMatrixParam();//Default parameter for the mixing matrix
	assortativity = 0.7;
	w_assortativeness = 0;
	MRtoHR_assortativ = 1;//Default marital
	Mixing_Proportional.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(Mixing_Proportional.data(), Mixing_Proportional.data() + Mixing_Proportional.num_elements(), 0);
	BalancedMixing_Marrital.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(BalancedMixing_Marrital.data(), BalancedMixing_Marrital.data() + BalancedMixing_Marrital.num_elements(), 0);

	mixingGoals.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(mixingGoals.data(), mixingGoals.data() + mixingGoals.num_elements(), 0);

	all_mixing_t1.resize(4);
	all_mixing_tf.resize(4);
	for (int jj = 0; jj < 4; jj++)
	{
		all_mixing_t1[jj].resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
		all_mixing_tf[jj].resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	}
} // end SyphilisPopulation::init(int t_final, int t_init, int in_dt)

void SyphilisPopulation::SetSubPopMen(Array2d in_dist, Array1d in_total, Array1d in_time)
{
	Men.SetSubPop(in_dist,in_total, in_time); //Initialize population of men
} //end SyphilisPopulation::SetSubPopMen(Array2d in_dist, Array1d in_total, Array1d in_time)

void SyphilisPopulation::SetSubPopWoMen(Array2d in_dist, Array1d in_total, Array1d in_time)
{
	WoMen.SetSubPop(in_dist, in_total, in_time); //Initialize population of women
}//end SyphilisPopulation::SetSubPopWoMen(Array2d in_dist, Array1d in_total, Array1d in_time)

void SyphilisPopulation::ProgressPopulation(int n)
{
	Men.ProgressPopulation(n, FOIMen); //Initialize Progress all men risk groups
	WoMen.ProgressPopulation(n, FOIWoMen); //Initialize Progress all women risk groups
}//end SyphilisPopulation::ProgressPopulation(int n)

SyphilisPopulation::SyphilisPopulation()
{
	init(1, 0, 12);
} // end SyphilisPopulation::SyphilisPopulation(); constructor

SyphilisPopulation::~SyphilisPopulation()
{

}//end SyphilisPopulation::~SyphilisPopulation(); destructor

void SyphilisPopulation::calcMixingBalancing(const int&n)
{
	//Men.allmen[SyphLow].MeanNumPart[n] = 1;
	//Men.allmen[SyphMed].MeanNumPart[n] = 2.0- 52.70/100;
	//Men.allmen[SyphHig].MeanNumPart[n] = 11.0-59/100;
	//Men.allmen[SyphMSM].MeanNumPart[n] = 5.0-52.7/100;

	//Men.allmen[SyphLow].Total[n] = 1028356;
	//Men.allmen[SyphMed].Total[n] = 914094;
	//Men.allmen[SyphHig].Total[n] = 285655;
	//Men.allmen[SyphMSM].Total[n] = 57137;

	//WoMen.allwomen[SyphLow].MeanNumPart[n] = 1.0;
	//WoMen.allwomen[SyphMed].MeanNumPart[n] = 2.10-43.4/100;
	//WoMen.allwomen[SyphHig].MeanNumPart[n] = 70.0-0.1;

	//WoMen.allwomen[SyphLow].Total[n] = 1374196;
	//WoMen.allwomen[SyphMed].Total[n] = 893227;
	//WoMen.allwomen[SyphHig].Total[n] = 22903;

	std::vector<double> populationsizes(SyphNumRisksMen + SyphNumRisksWom);
	std::vector<double> nom_numbpart(SyphNumRisksMen + SyphNumRisksWom);
	nom_numbpart[SyphLow] = Men.allmen[SyphLow].MeanNumPart[n];
	for (int mrg = SyphLow+1; mrg < SyphNumRisksMen; mrg++)
	{
		nom_numbpart[mrg] = Men.allmen[mrg].MeanNumPart[n] + mixingmat.Mix_MenWomen[mrg][SyphLow];
	}

	nom_numbpart[SyphNumRisksMen + SyphLow] = WoMen.allwomen[SyphLow].MeanNumPart[n];
	for (int wrg = SyphLow + 1; wrg < SyphNumRisksWom; wrg++)
	{
		nom_numbpart[SyphNumRisksMen + wrg] = WoMen.allwomen[wrg].MeanNumPart[n] + mixingmat.Mix_WomenMen[wrg][SyphLow];
	}

	//Rcpp::Rcout << "Num part" << std::endl;
	//for (int wrg = 0; wrg < (SyphNumRisksWom+ SyphNumRisksMen); wrg++)
	//{
	//	Rcpp::Rcout << nom_numbpart[wrg] << "; ";
	//}
	//Rcpp::Rcout << "**" << std::endl;



	nom_numbpart[SyphNoSex] = 0;
	nom_numbpart[SyphNumRisksMen + SyphNoSexWom] = 0;

	std::vector<double> percent_married(SyphNumRisksMen + SyphNumRisksWom);
	std::vector<double> nom_overall_partnerships(SyphNumRisksMen + SyphNumRisksWom);

	double sumMarriagesMen = 0;
	for (int mrg = SyphLow; mrg < SyphNumRisksMen; mrg++)
	{
		populationsizes[mrg] = Men.allmen[mrg].Total[n];
		percent_married[mrg] = (mrg == SyphLow) ? 1.0 : mixingmat.Mix_MenWomen[mrg][SyphLow];
		nom_overall_partnerships[mrg] = populationsizes[mrg]* nom_numbpart[mrg];
		sumMarriagesMen+= populationsizes[mrg]*percent_married[mrg];
	}
	      
	double sumMarriagesWomen = 0.0;
	for (int wrg = SyphLow; wrg < SyphNumRisksWom; wrg++)
	{
		populationsizes[SyphNumRisksMen + wrg] = WoMen.allwomen[wrg].Total[n];
		percent_married[SyphNumRisksMen + wrg] = (wrg == SyphLow) ? 1.0 : mixingmat.Mix_WomenMen[wrg][SyphLow];
		nom_overall_partnerships[SyphNumRisksMen + wrg] = populationsizes[SyphNumRisksMen + wrg] * nom_numbpart[SyphNumRisksMen + wrg];
		sumMarriagesWomen += populationsizes[SyphNumRisksMen + wrg]* percent_married[SyphNumRisksMen + wrg];
	}

	//Array2d mixing;// (SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);
	mixingGoals.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(mixingGoals.data(), mixingGoals.data() + mixingGoals.num_elements(), 0);
	Array2d supply;// (SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);
	supply.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(supply.data(), supply.data() + supply.num_elements(), 0);
	Array2d demand;// (SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);
	demand.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(demand.data(), demand.data() + demand.num_elements(), 0);
	//Array2d Balancing;// (SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);
	Balancing.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom][SyphNumRisksMen + SyphNumRisksWom]);
	std::fill(Balancing.data(), Balancing.data() + Balancing.num_elements(), 0);
	if (sumMarriagesWomen > 0)
	{
		for (int wrg = SyphLow; wrg < SyphNumRisksWom; wrg++)
		{
			mixingGoals[SyphLow][SyphNumRisksMen + wrg] = populationsizes[SyphNumRisksMen + wrg] * percent_married[SyphNumRisksMen + wrg] / sumMarriagesWomen;
			demand[SyphLow][SyphNumRisksMen + wrg] = (1-w_assortativeness)*mixingGoals[SyphLow][SyphNumRisksMen + wrg] * nom_overall_partnerships[SyphLow];
			demand[SyphLow][SyphNumRisksMen + wrg] += w_assortativeness * Mixing_Proportional[SyphLow][SyphNumRisksMen + wrg] * nom_overall_partnerships[SyphLow];
		}
	}

	for (int mrg = SyphMed; mrg < SyphNumRisksMen; mrg++)
	{
		double numpart = nom_numbpart[mrg];
		mixingGoals[mrg][ SyphNumRisksMen + SyphLow] = (numpart > 0) ? percent_married[mrg] / numpart:0.0;
		double ct2 = (numpart > 0) ? ((1.0 - percent_married[mrg])*numpart + percent_married[mrg] * (numpart - 1)) / numpart : 0;

		//double MRHRConst = (numpart > 0) ? ((1.0 - percent_married[SyphNumRisksMen + wrg]) * numpart + percent_married[SyphNumRisksMen + wrg] * (numpart - 1)) / numpart : 0;
		if (mrg != SyphMSM)
		{
			double temp_Supp1 = (nom_numbpart[SyphNumRisksMen + SyphMed] - percent_married[SyphNumRisksMen + SyphMed]) * populationsizes[SyphNumRisksMen + SyphMed];
			double temp_Supp2 = (nom_numbpart[SyphNumRisksMen + SyphHig] - percent_married[SyphNumRisksMen + SyphHig]) * populationsizes[SyphNumRisksMen + SyphHig];
			double temp_Supp_all = temp_Supp1 + temp_Supp2; temp_Supp_all = (temp_Supp_all <= 0) ? 1 : temp_Supp_all;
			double num_supp = (mrg == SyphMed) ? temp_Supp1 / temp_Supp_all : ((mrg == SyphHig) ? temp_Supp2 / temp_Supp_all : 0);
			//mixingGoals[SyphNumRisksMen + wrg][wrg] = MRHRConst * ((1 - MRtoHR_assortativ) * num_supp + MRtoHR_assortativ);

			mixingGoals[mrg][SyphNumRisksMen + mrg] = ct2*((1- MRtoHR_assortativ)* num_supp + MRtoHR_assortativ);
			demand[mrg][SyphNumRisksMen + mrg] = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + mrg] * nom_overall_partnerships[mrg];
			demand[mrg][SyphNumRisksMen + mrg] += w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + mrg] * nom_overall_partnerships[mrg];

			if (mrg == SyphMed)
			{
				mixingGoals[mrg][SyphNumRisksMen + SyphHig] = ct2 * (1 - MRtoHR_assortativ) *(1-num_supp);
			}
			else if (mrg == SyphHig)
			{
				mixingGoals[mrg][SyphNumRisksMen + SyphMed] = ct2 * (1 - MRtoHR_assortativ) * (1 - num_supp);
			}

		}
		else
		{
			mixingGoals[mrg][mrg] = ct2;
			demand[mrg][mrg] = (1 - w_assortativeness) * mixingGoals[mrg][mrg] * nom_overall_partnerships[mrg];
			demand[mrg][mrg] += w_assortativeness* Mixing_Proportional[mrg][mrg] * nom_overall_partnerships[mrg];

			//New
			demand[mrg][SyphNumRisksMen + SyphHig] = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphHig] * nom_overall_partnerships[mrg];
			demand[mrg][SyphNumRisksMen + SyphHig] += w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphHig] * nom_overall_partnerships[mrg];
			demand[mrg][SyphNumRisksMen + SyphMed] = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphMed] * nom_overall_partnerships[mrg];
			demand[mrg][SyphNumRisksMen + SyphMed] += w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphMed] * nom_overall_partnerships[mrg];
		}
		
		demand[mrg][SyphNumRisksMen + SyphLow] = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphLow] * nom_overall_partnerships[mrg];
		demand[mrg][SyphNumRisksMen + SyphLow] += w_assortativeness* Mixing_Proportional[mrg][SyphNumRisksMen + SyphLow] * nom_overall_partnerships[mrg];

		//New
		if (mrg == SyphMed)
		{
			demand[mrg][SyphNumRisksMen + SyphHig] = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphHig] * nom_overall_partnerships[mrg];
			demand[mrg][SyphNumRisksMen + SyphHig] += w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphHig] * nom_overall_partnerships[mrg];
		}

		if (mrg == SyphHig)
		{
			demand[mrg][SyphNumRisksMen + SyphMed] = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphMed] * nom_overall_partnerships[mrg];
			demand[mrg][SyphNumRisksMen + SyphMed] += w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphMed] * nom_overall_partnerships[mrg];
		}
	}
	
	if (sumMarriagesMen > 0)
	{
		for (int mrg = SyphLow; mrg < SyphNumRisksMen; mrg++)
		{
			mixingGoals[SyphNumRisksMen + SyphLow][mrg] = populationsizes[mrg] * percent_married[mrg] / sumMarriagesMen;
			demand[SyphNumRisksMen + SyphLow][mrg] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + SyphLow][mrg] * nom_overall_partnerships[SyphNumRisksMen + SyphLow];// *nom_overall_partnerships[SyphNumRisksMen + SyphLow];
			demand[SyphNumRisksMen + SyphLow][mrg] += w_assortativeness* Mixing_Proportional[SyphNumRisksMen + SyphLow][mrg] * nom_overall_partnerships[SyphNumRisksMen + SyphLow];
		}
	}
	//Proportionnal mixing
	for (int mrg = SyphLow; mrg < SyphNumRisksMen; mrg++)
	{
		demand[SyphNumRisksMen + SyphLow][mrg] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + SyphLow][mrg] * nom_overall_partnerships[SyphNumRisksMen + SyphLow];
	}

	for (int wrg = SyphMed; wrg < SyphNumRisksWom; wrg++)
	{
		double numpart = nom_numbpart[SyphNumRisksMen + wrg];
		mixingGoals[SyphNumRisksMen + wrg][SyphLow] = (numpart > 0) ? percent_married[SyphNumRisksMen + wrg] / numpart : 0.0;
		
		double MRHRConst= (numpart > 0) ? ((1.0 - percent_married[SyphNumRisksMen + wrg]) * numpart + percent_married[SyphNumRisksMen + wrg] * (numpart - 1)) / numpart : 0;
		double temp_Supp1 = (nom_numbpart[SyphMed] - percent_married[SyphMed]) * populationsizes[SyphMed];
		double temp_Supp2 = (nom_numbpart[SyphHig] - percent_married[SyphHig]) * populationsizes[SyphHig];
		double temp_Supp_all = temp_Supp1 + temp_Supp2; temp_Supp_all = (temp_Supp_all <= 0) ? 1 : temp_Supp_all;
		double num_supp = (wrg == SyphMed) ? temp_Supp1 / temp_Supp_all : ((wrg == SyphHig) ? temp_Supp2 / temp_Supp_all : 0);
		mixingGoals[SyphNumRisksMen + wrg][wrg] = MRHRConst * ((1- MRtoHR_assortativ)* num_supp +MRtoHR_assortativ);
		
		demand[SyphNumRisksMen + wrg][SyphLow] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphLow]*nom_overall_partnerships[SyphNumRisksMen + wrg];
		demand[SyphNumRisksMen + wrg][SyphLow] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphLow] * nom_overall_partnerships[SyphNumRisksMen + wrg];

		demand[SyphNumRisksMen + wrg][wrg] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][wrg] * nom_overall_partnerships[SyphNumRisksMen + wrg];
		demand[SyphNumRisksMen + wrg][wrg] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][wrg] * nom_overall_partnerships[SyphNumRisksMen + wrg];
	
		if (wrg == SyphMed)
		{
			mixingGoals[SyphNumRisksMen + wrg][SyphHig] = MRHRConst * ((1 - MRtoHR_assortativ) * (1-num_supp));

			demand[SyphNumRisksMen + wrg][SyphHig] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphHig] * nom_overall_partnerships[SyphNumRisksMen + wrg];
			demand[SyphNumRisksMen + wrg][SyphHig] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphHig] * nom_overall_partnerships[SyphNumRisksMen + wrg];
		
			demand[SyphNumRisksMen + wrg][SyphMSM] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphMSM] * nom_overall_partnerships[SyphNumRisksMen + wrg];
			demand[SyphNumRisksMen + wrg][SyphMSM] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphMSM] * nom_overall_partnerships[SyphNumRisksMen + wrg];

		}

		if (wrg == SyphHig)
		{
			mixingGoals[SyphNumRisksMen + wrg][SyphMed] = MRHRConst * ((1 - MRtoHR_assortativ) * (1 - num_supp));

			demand[SyphNumRisksMen + wrg][SyphMed] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphMed] * nom_overall_partnerships[SyphNumRisksMen + wrg];
			demand[SyphNumRisksMen + wrg][SyphMed] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphMed] * nom_overall_partnerships[SyphNumRisksMen + wrg];
		
			demand[SyphNumRisksMen + wrg][SyphMSM] = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphMSM] * nom_overall_partnerships[SyphNumRisksMen + wrg];
			demand[SyphNumRisksMen + wrg][SyphMSM] += w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphMSM] * nom_overall_partnerships[SyphNumRisksMen + wrg];
		}
	}

    ///
	if (sumMarriagesWomen > 0)
	{
		for (int wrg = SyphLow; wrg < SyphNumRisksWom; wrg++)
		{
			double temp = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphLow] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphLow];
			supply[SyphLow][SyphNumRisksMen + wrg] = nom_overall_partnerships[SyphNumRisksMen + wrg]* temp;// *
			Balancing[SyphLow][SyphNumRisksMen + wrg] = (demand[SyphLow][SyphNumRisksMen + wrg] > 0) ? std::sqrt(supply[SyphLow][SyphNumRisksMen + wrg] / demand[SyphLow][SyphNumRisksMen + wrg]) : 0;
			BalancedMixing_Marrital[SyphLow][SyphNumRisksMen + wrg] = Balancing[SyphLow][SyphNumRisksMen + wrg]*nom_numbpart[SyphLow] * temp;
		}
	}

	for (int mrg = SyphMed; mrg < SyphNumRisksMen; mrg++)
	{
		double temp1 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + SyphLow][mrg] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + SyphLow][mrg];
		double temp2 = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphLow] + w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphLow];
		
		supply[mrg][SyphNumRisksMen + SyphLow] = nom_overall_partnerships[SyphNumRisksMen + SyphLow]* temp1;
		Balancing[mrg][SyphNumRisksMen + SyphLow] = (demand[mrg][SyphNumRisksMen + SyphLow] > 0) ? std::sqrt(supply[mrg][SyphNumRisksMen + SyphLow] / demand[mrg][SyphNumRisksMen + SyphLow]) : 0.0;
		BalancedMixing_Marrital[mrg][SyphNumRisksMen + SyphLow] = Balancing[mrg][SyphNumRisksMen + SyphLow]*nom_numbpart[mrg] * temp2;
	
		if (mrg != SyphMSM)
		{
			temp1 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + mrg][mrg] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + mrg][mrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + mrg] + w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + mrg];
			supply[mrg][SyphNumRisksMen + mrg] = nom_overall_partnerships[SyphNumRisksMen + mrg] * temp1;
			Balancing[mrg][SyphNumRisksMen + mrg] = (demand[mrg][SyphNumRisksMen + mrg] > 0) ? std::sqrt(supply[mrg][SyphNumRisksMen + mrg] / demand[mrg][SyphNumRisksMen + mrg]) : 0.0;
			BalancedMixing_Marrital[mrg][SyphNumRisksMen + mrg] = Balancing[mrg][SyphNumRisksMen + mrg] * nom_numbpart[mrg] * temp2;
		}
		else
		{
			temp1 = (1 - w_assortativeness) * mixingGoals[mrg][mrg] + w_assortativeness * Mixing_Proportional[mrg][mrg];
			supply[mrg][mrg] = nom_overall_partnerships[mrg] * temp1;
			Balancing[mrg][mrg] = (demand[mrg][mrg] > 0) ? std::sqrt(supply[mrg][mrg] / demand[mrg][mrg]) : 0.0;
			BalancedMixing_Marrital[mrg][mrg] = Balancing[mrg][mrg] * nom_numbpart[mrg] * temp1;

		}

		if ((mrg == SyphMed) | (mrg == SyphMSM))
		{
			temp1 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + SyphHig][mrg] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + SyphHig][mrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphHig] + w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphHig];
			supply[mrg][SyphNumRisksMen + SyphHig] = nom_overall_partnerships[mrg] * temp1;
			Balancing[mrg][SyphNumRisksMen + SyphHig] = (demand[mrg][SyphNumRisksMen + SyphHig] > 0) ? std::sqrt(supply[mrg][SyphNumRisksMen + SyphHig] / demand[mrg][SyphNumRisksMen + SyphHig]) : 0.0;
			BalancedMixing_Marrital[mrg][SyphNumRisksMen + SyphHig] = Balancing[mrg][SyphNumRisksMen + SyphHig] * nom_numbpart[mrg] * temp2;
		}

		if ((mrg == SyphHig) | (mrg == SyphMSM))
		{
			temp1 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + SyphMed][mrg] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + SyphMed][mrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphMed] + w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphMed];
			supply[mrg][SyphNumRisksMen + SyphMed] = nom_overall_partnerships[mrg] * temp1;
			Balancing[mrg][SyphNumRisksMen + SyphMed] = (demand[mrg][SyphNumRisksMen + SyphMed] > 0) ? std::sqrt(supply[mrg][SyphNumRisksMen + SyphMed] / demand[mrg][SyphNumRisksMen + SyphMed]) : 0.0;
			BalancedMixing_Marrital[mrg][SyphNumRisksMen + SyphMed] = Balancing[mrg][SyphNumRisksMen + SyphMed] * nom_numbpart[mrg] * temp2;
		}
	}

	if (sumMarriagesMen > 0)
	{
		for (int mrg = SyphLow; mrg < SyphNumRisksMen; mrg++)
		{
			double temp1 = (1 - w_assortativeness) * mixingGoals[mrg][SyphNumRisksMen + SyphLow] + w_assortativeness * Mixing_Proportional[mrg][SyphNumRisksMen + SyphLow];
			double temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + SyphLow][mrg] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + SyphLow][mrg];
			supply[SyphNumRisksMen + SyphLow][mrg] = nom_overall_partnerships[mrg] * temp1;// *
			Balancing[SyphNumRisksMen + SyphLow][mrg] = (demand[SyphNumRisksMen + SyphLow][mrg]>0) ? std::sqrt(supply[SyphNumRisksMen + SyphLow][mrg]/ demand[SyphNumRisksMen + SyphLow][mrg] ):0.0;
			BalancedMixing_Marrital[SyphNumRisksMen + SyphLow][mrg] = Balancing[SyphNumRisksMen + SyphLow][mrg] * nom_numbpart[SyphNumRisksMen + SyphLow] * temp2;
		}
	}

	for (int wrg = SyphMed; wrg < SyphNumRisksWom; wrg++)
	{
		double temp1 = (1 - w_assortativeness) * mixingGoals[SyphLow][SyphNumRisksMen + wrg] + w_assortativeness * Mixing_Proportional[SyphLow][SyphNumRisksMen + wrg];
		double temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphLow] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphLow];
		supply[SyphNumRisksMen + wrg][SyphLow] = nom_overall_partnerships[SyphLow]* temp1;
		Balancing[SyphNumRisksMen + wrg][SyphLow] = (supply[SyphNumRisksMen + wrg][SyphLow]>0) ? std::sqrt(supply[SyphNumRisksMen + wrg][SyphLow]/ demand[SyphNumRisksMen + wrg][SyphLow]):0.0;
		BalancedMixing_Marrital[SyphNumRisksMen + wrg][SyphLow] = Balancing[SyphNumRisksMen + wrg][SyphLow]*nom_numbpart[SyphNumRisksMen + SyphLow] * temp2;
	
		temp1 = (1 - w_assortativeness) * mixingGoals[wrg][SyphNumRisksMen + wrg] + w_assortativeness * Mixing_Proportional[wrg][SyphNumRisksMen + wrg];
		temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][wrg] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][wrg];
		supply[SyphNumRisksMen + wrg][wrg] = nom_overall_partnerships[wrg] * temp1;
		Balancing[SyphNumRisksMen + wrg][wrg] = (supply[SyphNumRisksMen + wrg][wrg] > 0) ? std::sqrt(supply[SyphNumRisksMen + wrg][wrg] / demand[SyphNumRisksMen + wrg][wrg]) : 0.0;
		BalancedMixing_Marrital[SyphNumRisksMen + wrg][wrg] = Balancing[SyphNumRisksMen + wrg][wrg] * nom_numbpart[SyphNumRisksMen + wrg] * temp2;

		if (wrg == SyphMed)
		{
			temp1 = (1 - w_assortativeness) * mixingGoals[SyphHig][SyphNumRisksMen + wrg] + w_assortativeness * Mixing_Proportional[SyphHig][SyphNumRisksMen + wrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphHig] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphHig];
			supply[SyphNumRisksMen + wrg][SyphHig] = nom_overall_partnerships[SyphHig] * temp1;
			Balancing[SyphNumRisksMen + wrg][SyphHig] = (supply[SyphNumRisksMen + wrg][SyphHig] > 0) ? std::sqrt(supply[SyphNumRisksMen + wrg][SyphHig] / demand[SyphNumRisksMen + wrg][SyphHig]) : 0.0;
			BalancedMixing_Marrital[SyphNumRisksMen + wrg][SyphHig] = Balancing[SyphNumRisksMen + wrg][SyphHig] * nom_numbpart[SyphNumRisksMen + SyphHig] * temp2;
		
			temp1 = (1 - w_assortativeness) * mixingGoals[SyphMSM][SyphNumRisksMen + wrg] + w_assortativeness * Mixing_Proportional[SyphMSM][SyphNumRisksMen + wrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphMSM] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphHig];
			supply[SyphNumRisksMen + wrg][SyphMSM] = nom_overall_partnerships[SyphMSM] * temp1;
			Balancing[SyphNumRisksMen + wrg][SyphMSM] = (supply[SyphNumRisksMen + wrg][SyphMSM] > 0) ? std::sqrt(supply[SyphNumRisksMen + wrg][SyphMSM] / demand[SyphNumRisksMen + wrg][SyphMSM]) : 0.0;
			BalancedMixing_Marrital[SyphNumRisksMen + wrg][SyphMSM] = Balancing[SyphNumRisksMen + wrg][SyphMSM] * nom_numbpart[SyphNumRisksMen + wrg] * temp2;
		}

		if (wrg == SyphHig)
		{
			temp1 = (1 - w_assortativeness) * mixingGoals[SyphMed][SyphNumRisksMen + wrg] + w_assortativeness * Mixing_Proportional[SyphMed][SyphNumRisksMen + wrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphMed] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphMed];
			supply[SyphNumRisksMen + wrg][SyphMed] = nom_overall_partnerships[SyphMed] * temp1;
			Balancing[SyphNumRisksMen + wrg][SyphMed] = (supply[SyphNumRisksMen + wrg][SyphMed] > 0) ? std::sqrt(supply[SyphNumRisksMen + wrg][SyphMed] / demand[SyphNumRisksMen + wrg][SyphMed]) : 0.0;
			BalancedMixing_Marrital[SyphNumRisksMen + wrg][SyphMed] = Balancing[SyphNumRisksMen + wrg][SyphMed] * nom_numbpart[SyphNumRisksMen + SyphMed] * temp2;

			temp1 = (1 - w_assortativeness) * mixingGoals[SyphMSM][SyphNumRisksMen + wrg] + w_assortativeness * Mixing_Proportional[SyphMSM][SyphNumRisksMen + wrg];
			temp2 = (1 - w_assortativeness) * mixingGoals[SyphNumRisksMen + wrg][SyphMSM] + w_assortativeness * Mixing_Proportional[SyphNumRisksMen + wrg][SyphMSM];
			supply[SyphNumRisksMen + wrg][SyphMSM] = nom_overall_partnerships[SyphMSM] * temp1;
			Balancing[SyphNumRisksMen + wrg][SyphMSM] = (supply[SyphNumRisksMen + wrg][SyphMSM] > 0) ? std::sqrt(supply[SyphNumRisksMen + wrg][SyphMSM] / demand[SyphNumRisksMen + wrg][SyphMSM]) : 0.0;
			BalancedMixing_Marrital[SyphNumRisksMen + wrg][SyphMSM] = Balancing[SyphNumRisksMen + wrg][SyphMSM] * nom_numbpart[SyphNumRisksMen + wrg] * temp2;
		}
	}//end for (int wrg = SyphMed; wrg < SyphNumRisksWom; wrg++)

	all_mixing_tf[0] = mixingGoals;
	all_mixing_tf[1] = Mixing_Proportional;
	for (int wrg = SyphLow; wrg < (SyphNumRisksMen+SyphNumRisksWom); wrg++)
	{
		for (int mrg = SyphLow; mrg < (SyphNumRisksMen + SyphNumRisksWom); mrg++)
		{
			all_mixing_tf[2][wrg][mrg]= (1 - w_assortativeness) * mixingGoals[wrg][mrg] + w_assortativeness * Mixing_Proportional[wrg][mrg];
		}
	}
	all_mixing_tf[3] = BalancedMixing_Marrital;

	/*if (n == 1)
	{
		Rcpp::Rcout << "*************************" << std::endl;

		Rcpp::Rcout << "MRtoHR_assortativ=" << MRtoHR_assortativ << "; w_assortativeness=" << w_assortativeness <<
			"; assortativity=" << assortativity << std::endl;

		Rcpp::Rcout << "Goals Mixing" << std::endl;
		for (int mrg = 0; mrg < (SyphNumRisksMen + SyphNumRisksWom); mrg++)
		{
			for (int wrg = 0; wrg < (SyphNumRisksMen + SyphNumRisksWom); wrg++)
			{
				Rcpp::Rcout << mixingGoals[mrg][wrg] << "; ";
			}
			Rcpp::Rcout << std::endl;
		}
		Rcpp::Rcout << std::endl;

		Rcpp::Rcout << "Assortative Mixing" << std::endl;
		for (int mrg = 0; mrg < (SyphNumRisksMen + SyphNumRisksWom); mrg++)
		{
			for (int wrg = 0; wrg < (SyphNumRisksMen + SyphNumRisksWom); wrg++)
			{
				Rcpp::Rcout << Mixing_Proportional[mrg][wrg] << "; ";
			}
			Rcpp::Rcout << std::endl;
		}
		Rcpp::Rcout << std::endl;

		Rcpp::Rcout << "Combined Mixing" << std::endl;
		for (int mrg = 0; mrg < (SyphNumRisksMen + SyphNumRisksWom); mrg++)
		{
			for (int wrg = 0; wrg < (SyphNumRisksMen + SyphNumRisksWom); wrg++)
			{
				Rcpp::Rcout << all_mixing_tf[2][mrg][wrg] << "; ";
			}
			Rcpp::Rcout << std::endl;
		}
		Rcpp::Rcout << std::endl;

		Rcpp::Rcout << "Balanced Partners Change" << std::endl;
		for (int mrg = 0; mrg < (SyphNumRisksMen + SyphNumRisksWom); mrg++)
		{
			for (int wrg = 0; wrg < (SyphNumRisksMen + SyphNumRisksWom); wrg++)
			{
				Rcpp::Rcout << BalancedMixing_Marrital[mrg][wrg] << "; ";
			}
			Rcpp::Rcout << std::endl;
		}
		Rcpp::Rcout << std::endl;
	}*/

}//end SyphilisPopulation::calcMixingBalancing(const int&n)

void SyphilisPopulation::calcMixingProportional(const int& n)
{
	//Men.allmen[SyphLow].MeanNumPart[n] = 1;
	//Men.allmen[SyphMed].MeanNumPart[n] = 2.0;// - 52.70 / 100;
	//Men.allmen[SyphHig].MeanNumPart[n] = 11.0;// - 59 / 100;
	//Men.allmen[SyphMSM].MeanNumPart[n] = 5.0;// - 52.7 / 100;

	//Men.allmen[SyphLow].Total[n] = 1028356;
	//Men.allmen[SyphMed].Total[n] = 914094;
	//Men.allmen[SyphHig].Total[n] = 285655;
	//Men.allmen[SyphMSM].Total[n] = 57137;

	//WoMen.allwomen[SyphLow].MeanNumPart[n] = 1.0;
	//WoMen.allwomen[SyphMed].MeanNumPart[n] = 2.10;// -43.4 / 100;
	//WoMen.allwomen[SyphHig].MeanNumPart[n] = 70.0;// - 0.1;

	//WoMen.allwomen[SyphLow].Total[n] = 1374196;
	//WoMen.allwomen[SyphMed].Total[n] = 893227;
	//WoMen.allwomen[SyphHig].Total[n] = 22903;

	Array1d supplydemand;
	supplydemand.resize(boost::extents[SyphNumRisksMen + SyphNumRisksWom]);
	double supplymen = 0;
	for (int mrg = SyphLow; mrg <= SyphMSM; mrg++)
	{
		double temp = Men.allmen[mrg].MeanNumPart[n] * Men.allmen[mrg].Total[n];
		supplydemand[mrg] = temp;
		supplymen += temp;
	}//end for (int mrg = SyphLow; mrg <= SyphMSM; mrg++)

	double supplywomen = 0;
	for (int wrg = SyphLow; wrg <= SyphFSW; wrg++)
	{
		double temp = WoMen.allwomen[wrg].MeanNumPart[n] * WoMen.allwomen[wrg].Total[n];
		supplydemand[SyphNumRisksMen+wrg] = temp;
		supplywomen += temp;
	}//end for (int wrg = SyphLow; wrg <= SyphFSW; wrg++)

	for (int mrg = SyphLow; mrg < SyphMSM; mrg++)
	{
		double fact_m = (supplymen > 0) ? supplydemand[mrg] / supplymen : 0;
		for (int wrg = SyphLow; wrg <= SyphFSW; wrg++)
		{
			double fact_w = (supplywomen > 0) ? supplydemand[SyphNumRisksMen + wrg] / supplywomen : 0;
			if (mrg == wrg)
			{
				Mixing_Proportional[mrg][SyphNumRisksMen + wrg] = (1 - assortativity) * fact_w + assortativity;
				Mixing_Proportional[SyphNumRisksMen + wrg][mrg] = (1 - assortativity) * fact_m + assortativity;
			}
			else
			{
				Mixing_Proportional[mrg][SyphNumRisksMen + wrg] = (1 - assortativity) * fact_w;
				Mixing_Proportional[SyphNumRisksMen + wrg][mrg] = (1 - assortativity) * fact_m;
			}
		}//end for (int wrg = SyphLow; wrg <= SyphFSW; wrg++)
	}//end for (int mrg = SyphLow; mrg < SyphMSM; mrg++)
	
	//Case of MSM
	double fact_m = (supplymen > 0) ? supplydemand[SyphMSM] / supplymen : 0;
	supplywomen += supplydemand[SyphMSM];
	for (int wrg = SyphLow; wrg <= SyphFSW; wrg++)
	{
		double fact_w = (supplywomen > 0) ? supplydemand[SyphNumRisksMen + wrg] / supplywomen : 0;
		Mixing_Proportional[SyphMSM][SyphNumRisksMen + wrg] = (1 - assortativity) * fact_w;
		Mixing_Proportional[SyphNumRisksMen + wrg][SyphMSM] = (1 - assortativity) * fact_m;
	}//end for (int wrg = SyphLow; wrg <= SyphFSW; wrg++)

	fact_m = supplydemand[SyphMSM] / supplywomen;
	Mixing_Proportional[SyphMSM][SyphMSM] = (1 - assortativity) * fact_m + assortativity;
}//end calcMixingProportional(const int& n)

void SyphilisPopulation::calcFOI(int n) //Calculation of the force of infection
{
	if (n < WoMen.allwomen[0].numsteps)
	{
		calcMixingProportional(n);//Calculate mixing under the assumption of proportionnal mixing
		calcMixingBalancing(n);//Performing mixing/balancing 
		for (int mrg = SyphLow; mrg < SyphNumRisksMen; mrg++)
		{
			double numsex = Men.allmen[mrg].MeanNumSexByPart[n];

			double cu_men = Men.allmen[mrg].CondomUse[n];
			double condomeff = Men.allmen[mrg].nh_param.condomeff;

			FOIMen[mrg] = 0;
			Men.allmen[mrg].referral_rate[n] = 0;
			Men.allmen[mrg].NumContactTraced[all_contactTraced][n] = 0;
			Men.allmen[mrg].NumContactTraced[prim_and_sec_contactTraced][n] = 0;
			Men.allmen[mrg].NumContactTraced[latent_contactTraced][n] = 0;

			double epsiM = Men.allmen[mrg].pinfectdetect[n];
			double infectiousM = Men.allmen[mrg].IActive[n]+Men.allmen[mrg].EarlyRecov[n];//Active individuals and individuals who have been diagnosed and treated
			double allprop2_0 = Men.allmen[mrg].IActive[n] / (epsiM*infectiousM + (1 - epsiM)*Men.allmen[mrg].Total[n]);

			for (int pop = SyphLow; pop < SyphNumRisksWom; pop++)
			{

				//if (n == 0)
				//{
				//	Rcpp::Rcout << "*** mrg="<< mrg <<"; pop="<< pop << "; " << mixingmat.Mix_MenWomen[mrg][pop] << "; " << mixingmat.Mix_WomenMen[pop][mrg] << "; " << Men.allmen[mrg].Total[n] << "; " << WoMen.allwomen[pop].Total[n] <<
				//		WoMen.allwomen[pop].IActive[n] << std::endl;
				//}



				//if ((mixingmat.Mix_MenWomen[mrg][pop] > 0) & (mixingmat.Mix_WomenMen[pop][mrg] > 0) & (Men.allmen[mrg].Total[n] > 0) & (WoMen.allwomen[pop].Total[n] > 0) &(WoMen.allwomen[pop].IActive[n]>0))
				if ((Men.allmen[mrg].Total[n] > 0) & (WoMen.allwomen[pop].Total[n] > 0) & (WoMen.allwomen[pop].IActive[n] > 0))
				{
					double new_cu_men = (mrg== SyphLow) ? cu_men : ((pop== SyphLow) ? WoMen.allwomen[pop].CondomUse[n] : 0.5*(cu_men+WoMen.allwomen[pop].CondomUse[n]));
					new_cu_men *= condomeff;// nathistory.condomeff;
					new_cu_men = 1.0 - new_cu_men; //reduction in transmission by condom use

					double pi_j =  WoMen.allwomen[pop].IActive[n] / WoMen.allwomen[pop].Total[n]; //Only active syphilis can be transmitted
					double temp_numsex = (pop == SyphLow) ? WoMen.allwomen[pop].MeanNumSexByPart[n] : numsex;

					double Ctrans_prob = std::pow(1.0 - Men.allmen[mrg].nh_param.FtoM_tp*new_cu_men, temp_numsex);
					double temp_denom = BalancedMixing_Marrital[mrg][SyphNumRisksMen + pop];

					double Prob_ij = 1.0 - std::pow(1 - pi_j + pi_j * Ctrans_prob, temp_denom);
					FOIMen[mrg] += Prob_ij;


					//if (n == 0)
					//{
					//	Rcpp::Rcout << Prob_ij << "; " << pi_j << "; " << Ctrans_prob << "; " << temp_denom <<
					//		"; " << temp_numsex << "; FOIMen[mrg]="<< FOIMen[mrg] <<std::endl;
					//}

					//Referral rates
					double epsiW = WoMen.allwomen[pop].pinfectdetect[n];//propensity to identify an infected partner
					double pi_jRef = WoMen.allwomen[pop].IActive[n] / (epsiW*(WoMen.allwomen[pop].IActive[n]+ WoMen.allwomen[pop].EarlyRecov[n]) + (1 - epsiW)* WoMen.allwomen[pop].Total[n]); //Have necessary mixed with infectious people

					Men.allmen[mrg].referral_rate[n] += temp_denom* pi_jRef*referralmat.Refer_WomenMen[pop][mrg] * WoMen.allwomen[pop].nh_param.vect_nu[n] * WoMen.allwomen[pop].nh_param.prop_part_referred[n];
					
					if (Men.allmen[mrg].IActive[n] > 0)
					{
						double allprop = Men.allmen[mrg].IActive[n]* Men.allmen[mrg].nh_param.vect_nu[n];//proportion of men in this risk group with active syphilis
						allprop *= temp_denom * referralmat.Refer_MenWomen[mrg][pop] * Men.allmen[mrg].nh_param.prop_part_referred[n];

						double  allprop2 = allprop2_0* Men.allmen[mrg].nh_param.vect_nu[n] * temp_denom* referralmat.Refer_MenWomen[mrg][pop] * Men.allmen[mrg].nh_param.prop_part_referred[n];//proportion of men in this risk group with active syphilis

						Men.allmen[mrg].NumContactTraced[all_contactTraced][n] += allprop;// *WoMen.allwomen[pop].Total[n];

						double numact2 = allprop2 * WoMen.allwomen[pop].IActive[n];
						numact2 = (numact2 < allprop) ? numact2 : allprop;
						Men.allmen[mrg].NumContactTraced[prim_and_sec_contactTraced][n] += numact2;
						Men.allmen[mrg].NumContactTraced[latent_contactTraced][n] += allprop2 * WoMen.allwomen[pop].ILatent[n]; //I should change this
					}
				}
			}

			if (mrg == SyphMSM)
			{
				if (Men.allmen[mrg].Total[n] > 0)
				{
					double pi_j = Men.allmen[mrg].IActive[n] / Men.allmen[mrg].Total[n]; //Only active syphilis can be transmitted
					double new_cu_men = cu_men * condomeff;
					new_cu_men = 1.0 - new_cu_men;//reduction of transmission
					double Ctrans_prob = std::pow(1.0 - Men.allmen[mrg].nh_param.MtoM_tp*new_cu_men, numsex);
					
					double temp_denom = BalancedMixing_Marrital[mrg][mrg];
					double Prob_ij = 1.0 - std::pow(1.0 - pi_j + pi_j * Ctrans_prob, temp_denom);
					FOIMen[mrg] += Prob_ij;

					if (Men.allmen[mrg].IActive[n] > 0)
					{
						Men.allmen[mrg].referral_rate[n] += temp_denom * allprop2_0* referralmat.Refer_MenMen[mrg][mrg] * Men.allmen[mrg].nh_param.vect_nu[n] * Men.allmen[mrg].nh_param.prop_part_referred[n];
						
						double allprop = Men.allmen[mrg].IActive[n] / Men.allmen[mrg].Total[n] * Men.allmen[mrg].nh_param.vect_nu[n];//proportion of men in this risk group with active syphilis
						allprop *= temp_denom * referralmat.Refer_MenMen[mrg][mrg] * Men.allmen[mrg].nh_param.prop_part_referred[n];
					
						Men.allmen[mrg].NumContactTraced[all_contactTraced][n] += allprop * Men.allmen[mrg].Total[n];

                    	double allprop2 = allprop2_0 * Men.allmen[mrg].nh_param.vect_nu[n] * temp_denom* referralmat.Refer_MenMen[mrg][mrg] * Men.allmen[mrg].nh_param.prop_part_referred[n];//proportion of men in this risk group with active syphilis
						double numact2 = allprop2 * Men.allmen[mrg].IActive[n];

						numact2 = (numact2 < allprop *  Men.allmen[mrg].Total[n]) ? numact2 : allprop * Men.allmen[mrg].Total[n];
						Men.allmen[mrg].NumContactTraced[prim_and_sec_contactTraced][n] += numact2;
						Men.allmen[mrg].NumContactTraced[latent_contactTraced][n] += allprop2 * Men.allmen[mrg].ILatent[n];
					}
				}
			}
		}/*end Men*/

		/*Women force of infection*/
		for (int wrg = SyphLow; wrg < SyphNumRisksWom; wrg++)
		{
			double numsex = WoMen.allwomen[wrg].MeanNumSexByPart[n];

			FOIWoMen[wrg] = 0;
			WoMen.allwomen[wrg].referral_rate[n] = 0;
			WoMen.allwomen[wrg].NumContactTraced[all_contactTraced][n] = 0;
			WoMen.allwomen[wrg].NumContactTraced[prim_and_sec_contactTraced][n] = 0;
			WoMen.allwomen[wrg].NumContactTraced[latent_contactTraced][n] = 0;

			double cu_women = WoMen.allwomen[wrg].CondomUse[n];
			double condomeff = WoMen.allwomen[wrg].nh_param.condomeff;
			double epsiW = WoMen.allwomen[wrg].pinfectdetect[n];
			double infectiousW = WoMen.allwomen[wrg].IActive[n]+WoMen.allwomen[wrg].EarlyRecov[n];
			double allprop2_0 = WoMen.allwomen[wrg].IActive[n] / (epsiW*infectiousW + (1 - epsiW)*WoMen.allwomen[wrg].Total[n]);

			for (int pop = SyphLow; pop < SyphNumRisksMen; pop++)
			{
				//Will only get here if infections are possible
				//if ((mixingmat.Mix_WomenMen[wrg][pop] > 0) & (mixingmat.Mix_MenWomen[pop][wrg] > 0) & (WoMen.allwomen[wrg].Total[n] > 0) & (Men.allmen[pop].Total[n] > 0) &(Men.allmen[pop].IActive[n]>0))
				if ((WoMen.allwomen[wrg].Total[n] > 0) & (Men.allmen[pop].Total[n] > 0) & (Men.allmen[pop].IActive[n] > 0))
				{
					double pi_j = Men.allmen[pop].IActive[n] / Men.allmen[pop].Total[n]; //Only active syphilis can be transmitted

					double new_cu_women = (wrg== SyphLow) ? cu_women : ((pop== SyphLow) ? Men.allmen[pop].CondomUse[n] : 0.5*(cu_women + Men.allmen[pop].CondomUse[n]));
					new_cu_women *= condomeff;
					new_cu_women = 1.0 - new_cu_women;

					double temp_numsex = (pop == SyphLow) ? Men.allmen[pop].MeanNumSexByPart[n] : numsex;					
					double Ctrans_prob = std::pow(1.0 - WoMen.allwomen[wrg].nh_param.MtoF_tp*new_cu_women, temp_numsex);
					
					double temp_denom = BalancedMixing_Marrital[SyphNumRisksMen + wrg][pop];
					double Prob_ij = 1 - std::pow(1 - pi_j + pi_j * Ctrans_prob, temp_denom);
					FOIWoMen[wrg] += Prob_ij; 

					double epsiM = Men.allmen[pop].pinfectdetect[n];
					double pi_jRef = Men.allmen[pop].IActive[n] / (epsiM*(Men.allmen[pop].IActive[n]+Men.allmen[pop].EarlyRecov[n]) + (1 - epsiM)*Men.allmen[pop].Total[n]); //For Referral, we assume infected individuals had at least one infected partner

					WoMen.allwomen[wrg].referral_rate[n] += temp_denom * pi_jRef * referralmat.Refer_MenWomen[pop][wrg] * Men.allmen[pop].nh_param.vect_nu[n] * Men.allmen[pop].nh_param.prop_part_referred[n];
					if (WoMen.allwomen[wrg].IActive[n] > 0)
					{
						double allprop = WoMen.allwomen[wrg].IActive[n] * WoMen.allwomen[wrg].nh_param.vect_nu[n];//proportion of women in this risk group with active syphilis
						allprop *= temp_denom * referralmat.Refer_WomenMen[wrg][pop] * WoMen.allwomen[wrg].nh_param.prop_part_referred[n];

						WoMen.allwomen[wrg].NumContactTraced[all_contactTraced][n] += allprop;

						double allprop2 = allprop2_0 * WoMen.allwomen[wrg].nh_param.vect_nu[n] * temp_denom* referralmat.Refer_WomenMen[wrg][pop] * WoMen.allwomen[wrg].nh_param.prop_part_referred[n];//proportion of men in this risk group with active syphilis
						double numact2 = allprop2 * Men.allmen[pop].IActive[n];
						numact2 = (numact2 < allprop) ? numact2 : allprop;
						WoMen.allwomen[wrg].NumContactTraced[prim_and_sec_contactTraced][n] += numact2;
						WoMen.allwomen[wrg].NumContactTraced[latent_contactTraced][n] += allprop2 * Men.allmen[pop].ILatent[n];
					}
				}//end if 
			}//end for men
		}//end for women
	}//end if the number of steps is larger than 0
} //end SyphilisPopulation::calcFOI(int n)

void SyphilisPopulation::calcinitreferral()//???
{
	calcFOI(0);
	int len_men = Men.allmen.size();
	for (int pop = 0; pop < len_men; pop++)
	{
		Men.allmen[pop].NumReferred[0] = Men.allmen[pop].Total[0] * Men.allmen[pop].referral_rate[0];
		Men.allmen[pop].NumReferredTestedPositive[RPRTest][0] = (Men.allmen[pop].IActive[0] + Men.allmen[pop].EarlyRecov[0] + Men.allmen[pop].ILatent[0] + Men.allmen[pop].LateRecov[0])*Men.allmen[pop].referral_rate[0];
		Men.allmen[pop].NumReferredTestedPositive[TPHATest][0] = (Men.allmen[pop].IActive[0] + Men.allmen[pop].EarlyRecov[0] + Men.allmen[pop].ILatent[0] + Men.allmen[pop].LateRecov[0] + Men.allmen[pop].RecovSSneg[0])*Men.allmen[pop].referral_rate[0];
	}

	int len_women = WoMen.allwomen.size();
	for (int pop = 0; pop < len_women; pop++)
	{
		WoMen.allwomen[pop].NumReferred[0] = WoMen.allwomen[pop].Total[0] * WoMen.allwomen[pop].referral_rate[0];
		WoMen.allwomen[pop].NumReferredTestedPositive[RPRTest][0] = (WoMen.allwomen[pop].IActive[0] + WoMen.allwomen[pop].EarlyRecov[0] + WoMen.allwomen[pop].ILatent[0] + WoMen.allwomen[pop].LateRecov[0])*WoMen.allwomen[pop].referral_rate[0];
		WoMen.allwomen[pop].NumReferredTestedPositive[TPHATest][0] = (WoMen.allwomen[pop].IActive[0] + WoMen.allwomen[pop].EarlyRecov[0] + WoMen.allwomen[pop].ILatent[0] + WoMen.allwomen[pop].LateRecov[0] + WoMen.allwomen[pop].RecovSSneg[0])*WoMen.allwomen[pop].referral_rate[0];
	}
}

void SyphilisPopulation::copyMWPops(int n1, SyphilisPopulation mwpop2, int n2, bool cpnath)
{
	Men.copyfromMen(n1, mwpop2.Men,n2, cpnath);
	WoMen.copyfromWoMen(n1, mwpop2.WoMen, n2, cpnath);
}

void SyphilisPopulation::copyMWPops(int n1, SyphilisPopulation *mwpop2, int n2, bool cpnath)
{
	Men.copyfromMen(n1, mwpop2->Men, n2, cpnath);
	WoMen.copyfromWoMen(n1, mwpop2->WoMen, n2, cpnath);
} //end SyphilisPopulation::copyMWPops(int n1, SyphilisPopulation mwpop2, int n2)

/////////////////////////////
/*Syphilis Population for R*/
SyphilisPopulationforR::SyphilisPopulationforR(const SyphilisPopulation sampop)
{
	int len_men = sampop.Men.allmen.size();
	Men = Rcpp::List::create(Rcpp::Named("LowRisk") = Rcpp::List::create(), Rcpp::Named("MediumRisk") = Rcpp::List::create(),
		Rcpp::Named("HighRisk") = Rcpp::List::create(), Rcpp::Named("MSM") = Rcpp::List::create(), Rcpp::Named("NoSex") = Rcpp::List::create());//Rcpp::List(len_men)
	int len_time = sampop.Men.allmen[0].numsteps;

	for (int pop = 0; pop < len_men; pop++)
	{
		Rcpp::NumericVector SSneg(len_time);
		Rcpp::NumericVector IIncub(len_time);
		Rcpp::NumericVector IActive(len_time);
		Rcpp::NumericVector RecovSSneg(len_time);
		Rcpp::NumericVector ILatent(len_time);
		Rcpp::NumericVector EarlyRecov(len_time);
		Rcpp::NumericVector LateRecov(len_time);
		Rcpp::NumericVector Total(len_time);
		Rcpp::NumericVector MeanNumPart(len_time);
		Rcpp::NumericVector MeanNumSexByPart(len_time);
		Rcpp::NumericVector Incidence(len_time);
		Rcpp::NumericVector Prevalence(len_time);
		Rcpp::NumericMatrix allPrevalence(numPrevTypes, len_time);
		Rcpp::NumericMatrix NumScreened(numTests,len_time);
		Rcpp::NumericMatrix NumTreatedaftScreened(numTests, len_time);

		Rcpp::NumericVector NumIncidentCases(len_time);
		Rcpp::NumericVector NumSymptomaticClinTreat(len_time);
		
		Rcpp::NumericMatrix NumReferredTestedPositive(numTests, len_time);
		Rcpp::NumericVector NumReferred(len_time);

		Rcpp::NumericMatrix NumContactTraced(sizeContactTraced, len_time);

		for (int ti = 0; ti < len_time; ti++)
		{
			SSneg[ti] = sampop.Men.allmen[pop].SSneg[ti];
			IIncub[ti] = sampop.Men.allmen[pop].IIncub[ti];
			IActive[ti] = sampop.Men.allmen[pop].IActive[ti];
			RecovSSneg[ti] = sampop.Men.allmen[pop].RecovSSneg[ti];
			ILatent[ti] = sampop.Men.allmen[pop].ILatent[ti];
			EarlyRecov[ti] = sampop.Men.allmen[pop].EarlyRecov[ti];
			LateRecov[ti] = sampop.Men.allmen[pop].LateRecov[ti];
			Total[ti] = sampop.Men.allmen[pop].Total[ti];
			MeanNumPart[ti] = sampop.Men.allmen[pop].MeanNumPart[ti];
			MeanNumSexByPart[ti] = sampop.Men.allmen[pop].MeanNumSexByPart[ti];
			NumScreened(TPHATest,ti) = sampop.Men.allmen[pop].NumScreened[TPHATest][ti];
			NumScreened(RPRTest, ti) = sampop.Men.allmen[pop].NumScreened[RPRTest][ti];
			NumTreatedaftScreened(TPHATest, ti) = sampop.Men.allmen[pop].NumTreatedaftScreened[TPHATest][ti];;
			NumTreatedaftScreened(RPRTest, ti) = sampop.Men.allmen[pop].NumTreatedaftScreened[RPRTest][ti];
			NumIncidentCases[ti] = sampop.Men.allmen[pop].NumIncidentCases[ti];
			NumSymptomaticClinTreat[ti] = sampop.Men.allmen[pop].NumSymptomaticClinTreat[ti];

			Incidence[ti] = sampop.Men.allmen[pop].Incidence[ti];
			Prevalence[ti] = sampop.Men.allmen[pop].Prevalence[ti];

			NumReferredTestedPositive(TPHATest, ti) = sampop.Men.allmen[pop].NumReferredTestedPositive[TPHATest][ti];
			NumReferredTestedPositive(RPRTest, ti) = sampop.Men.allmen[pop].NumReferredTestedPositive[RPRTest][ti];
			NumReferred[ti] = sampop.Men.allmen[pop].NumReferred[ti];

			NumContactTraced(all_contactTraced, ti) = sampop.Men.allmen[pop].NumContactTraced[all_contactTraced][ti];
			NumContactTraced(prim_and_sec_contactTraced, ti) = sampop.Men.allmen[pop].NumContactTraced[prim_and_sec_contactTraced][ti];
			NumContactTraced(latent_contactTraced, ti) = sampop.Men.allmen[pop].NumContactTraced[latent_contactTraced][ti];

			for (int jj = RPRnegTPHApos; jj < numPrevTypes; jj++)
			{
				allPrevalence(jj, ti) = sampop.Men.allmen[pop].Prevalence2d[jj][ti];
			}
		}

		Rcpp::List riskgroup = Rcpp::List::create(Rcpp::Named("SSneg") = SSneg, Rcpp::Named("IIncub") = IIncub,
			Rcpp::Named("IActive") = IActive, Rcpp::Named("RecovSSneg") = RecovSSneg, Rcpp::Named("ILatent") = ILatent, Rcpp::Named("EarlyRecov") = EarlyRecov,
			Rcpp::Named("LateRecov") = LateRecov, Rcpp::Named("Total") = Total, Rcpp::Named("MeanNumPart") = MeanNumPart, Rcpp::Named("MeanNumSexByPart") = MeanNumSexByPart, 
			Rcpp::Named("Incidence") = Incidence, Rcpp::Named("Prevalence") = Prevalence, Rcpp::Named("allPrevalenceTypes")= allPrevalence,
			Rcpp::Named("NumberScreened")= NumScreened, Rcpp::Named("NumTreatedaftScreened") = NumTreatedaftScreened, Rcpp::Named("NumIncidentCases")= NumIncidentCases,
			Rcpp::Named("NumSymptomaticClinTreat")= NumSymptomaticClinTreat, Rcpp::Named("NumReferred")= NumReferred,
			Rcpp::Named("NumReferredTestedPositive") = NumReferredTestedPositive, Rcpp::Named("NumContactTraced") = NumContactTraced);
		Men[pop] = riskgroup;
	}
	
	//Same thing for women
	int len_women = sampop.WoMen.allwomen.size();

	WoMen = Rcpp::List::create(Rcpp::Named("LowRisk") = Rcpp::List::create(), Rcpp::Named("MediumRisk") = Rcpp::List::create(),
		Rcpp::Named("HighRisk") = Rcpp::List::create(), Rcpp::Named("NoSex") = Rcpp::List::create());//Rcpp::List(len_men);
	len_time = sampop.WoMen.allwomen[0].numsteps;
	for (int pop = 0; pop < len_women; pop++)
	{
		Rcpp::NumericVector SSneg(len_time);
		Rcpp::NumericVector IIncub(len_time);
		Rcpp::NumericVector IActive(len_time);
		Rcpp::NumericVector RecovSSneg(len_time);
		Rcpp::NumericVector ILatent(len_time);
		Rcpp::NumericVector EarlyRecov(len_time);
		Rcpp::NumericVector LateRecov(len_time);
		Rcpp::NumericVector Total(len_time);
		Rcpp::NumericVector MeanNumPart(len_time);
		Rcpp::NumericVector MeanNumSexByPart(len_time);
		Rcpp::NumericVector Incidence(len_time);
		Rcpp::NumericVector Prevalence(len_time);
		Rcpp::NumericMatrix allPrevalence(numPrevTypes, len_time);
		Rcpp::NumericMatrix NumScreened(numTests, len_time);
		Rcpp::NumericMatrix NumTreatedaftScreened(numTests, len_time);

		Rcpp::NumericVector NumIncidentCases(len_time);
		Rcpp::NumericVector NumSymptomaticClinTreat(len_time);

		Rcpp::NumericMatrix NumReferredTestedPositive(numTests, len_time);
		Rcpp::NumericVector NumReferred(len_time);

		Rcpp::NumericMatrix NumContactTraced(sizeContactTraced, len_time);

		for (int ti = 0; ti < len_time; ti++)
		{
			SSneg[ti] = sampop.WoMen.allwomen[pop].SSneg[ti];
			IIncub[ti] = sampop.WoMen.allwomen[pop].IIncub[ti];
			IActive[ti] = sampop.WoMen.allwomen[pop].IActive[ti];
			RecovSSneg[ti] = sampop.WoMen.allwomen[pop].RecovSSneg[ti];
			ILatent[ti] = sampop.WoMen.allwomen[pop].ILatent[ti];
			EarlyRecov[ti] = sampop.WoMen.allwomen[pop].EarlyRecov[ti];
			LateRecov[ti] = sampop.WoMen.allwomen[pop].LateRecov[ti];
			Total[ti] = sampop.WoMen.allwomen[pop].Total[ti];
			MeanNumPart[ti] = sampop.WoMen.allwomen[pop].MeanNumPart[ti];
			MeanNumSexByPart[ti] = sampop.WoMen.allwomen[pop].MeanNumSexByPart[ti];
			NumScreened(TPHATest,ti) = sampop.WoMen.allwomen[pop].NumScreened[TPHATest][ti];
			NumScreened(RPRTest, ti) = sampop.WoMen.allwomen[pop].NumScreened[RPRTest][ti];
			NumTreatedaftScreened(TPHATest, ti) = sampop.WoMen.allwomen[pop].NumTreatedaftScreened[TPHATest][ti];
			NumTreatedaftScreened(RPRTest, ti) = sampop.WoMen.allwomen[pop].NumTreatedaftScreened[RPRTest][ti];

			NumIncidentCases[ti] = sampop.WoMen.allwomen[pop].NumIncidentCases[ti];
			NumSymptomaticClinTreat[ti] = sampop.WoMen.allwomen[pop].NumSymptomaticClinTreat[ti];
			Incidence[ti] = sampop.WoMen.allwomen[pop].Incidence[ti];
			Prevalence[ti] = sampop.WoMen.allwomen[pop].Prevalence[ti];

			NumReferredTestedPositive(TPHATest, ti) = sampop.WoMen.allwomen[pop].NumReferredTestedPositive[TPHATest][ti];
			NumReferredTestedPositive(RPRTest, ti) = sampop.WoMen.allwomen[pop].NumReferredTestedPositive[RPRTest][ti];
			NumReferred[ti] = sampop.WoMen.allwomen[pop].NumReferred[ti];

			NumContactTraced(all_contactTraced, ti) = sampop.WoMen.allwomen[pop].NumContactTraced[all_contactTraced][ti];
			NumContactTraced(prim_and_sec_contactTraced, ti) = sampop.WoMen.allwomen[pop].NumContactTraced[prim_and_sec_contactTraced][ti];
			NumContactTraced(latent_contactTraced, ti) = sampop.WoMen.allwomen[pop].NumContactTraced[latent_contactTraced][ti];

			for (int jj = RPRnegTPHApos; jj < numPrevTypes; jj++)
			{
				allPrevalence(jj, ti) = sampop.WoMen.allwomen[pop].Prevalence2d[jj][ti];
			}
		}

		Rcpp::List riskgroup = Rcpp::List::create(Rcpp::Named("SSneg") = SSneg, Rcpp::Named("IIncub") = IIncub,
			Rcpp::Named("IActive") = IActive, Rcpp::Named("RecovSSneg") = RecovSSneg, Rcpp::Named("ILatent") = ILatent, Rcpp::Named("EarlyRecov") = EarlyRecov,
			Rcpp::Named("LateRecov") = LateRecov, Rcpp::Named("Total") = Total, Rcpp::Named("MeanNumPart") = MeanNumPart, Rcpp::Named("MeanNumSexByPart") = MeanNumSexByPart,
			Rcpp::Named("Incidence") = Incidence, Rcpp::Named("Prevalence") = Prevalence, Rcpp::Named("allPrevalenceTypes") = allPrevalence,
			Rcpp::Named("NumberScreened") = NumScreened, Rcpp::Named("NumTreatedaftScreened") = NumTreatedaftScreened, Rcpp::Named("NumIncidentCases") = NumIncidentCases,
			Rcpp::Named("NumSymptomaticClinTreat") = NumSymptomaticClinTreat, Rcpp::Named("NumReferred") = NumReferred,
			Rcpp::Named("NumReferredTestedPositive") = NumReferredTestedPositive, Rcpp::Named("NumContactTraced") = NumContactTraced);

		WoMen[pop] = riskgroup;
	}

	mixingmat = Rcpp::NumericMatrix(SyphNumRisksMen+ SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom) ;
	nathistory = Rcpp::List::create();

	Rcpp::NumericMatrix CalcMixingt1(SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);
	Rcpp::NumericMatrix CalcMixingtf(SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);

	Rcpp::NumericMatrix BalancedPartt1(SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);
	Rcpp::NumericMatrix BalancedParttf(SyphNumRisksMen + SyphNumRisksWom, SyphNumRisksMen + SyphNumRisksWom);


	for (int mrg = 0; mrg < (SyphNumRisksMen + SyphNumRisksWom); mrg++)
	{
		for (int wrg = 0; wrg < (SyphNumRisksMen + SyphNumRisksWom); wrg++)
		{
			CalcMixingt1(mrg, wrg) = sampop.all_mixing_t1[2][mrg][wrg];
			CalcMixingtf(mrg, wrg) = sampop.all_mixing_tf[2][mrg][wrg];
			BalancedPartt1(mrg, wrg) = sampop.all_mixing_t1[3][mrg][wrg];
			BalancedParttf(mrg, wrg) = sampop.all_mixing_tf[3][mrg][wrg];
		}
	}

	InitMixingMat = Rcpp::List::create(Rcpp::Named("CombinedMat") = CalcMixingt1);
	FinalMixingMat = Rcpp::List::create(Rcpp::Named("CombinedMat") = CalcMixingtf);

	InitBalancedPartners = Rcpp::List::create(Rcpp::Named("BalancedPart") = BalancedPartt1);
	FinalBalancedPartners = Rcpp::List::create(Rcpp::Named("BalancedPart") = BalancedParttf);
} //end SyphilisPopulationforR::SyphilisPopulationforR(const SyphilisPopulation sampop)

SyphilisPopulationforR::~SyphilisPopulationforR()
{

}// end SyphilisPopulationforR::~SyphilisPopulationforR(); constructor

Rcpp::List SyphilisPopulationforR::outforR(double firstyear)
{
	Rcpp::List Result;
	Result = Rcpp::List::create(Rcpp::Named("Men") = Men, Rcpp::Named("WoMen") = WoMen, Rcpp::Named("MixingMat") = mixingmat, Rcpp::Named("nathistory")= nathistory, Rcpp::Named("firstyear") = firstyear,
		Rcpp::Named("InitMixingMat")= InitMixingMat, Rcpp::Named("FinalMixingMat")= FinalMixingMat, 
		Rcpp::Named("InitBalancedPartners")= InitBalancedPartners, Rcpp::Named("FinalBalancedPartners")= FinalBalancedPartners);
	return Result;
}//end Rcpp::List SyphilisPopulationforR::outforR();

Rcpp::List SyphilisPopulationforR::outforR()
{
	Rcpp::List Result=outforR(0);
	return Result;
}

/*Auxiliary functions*/
double expit(const double x)
{
	double result;
	if (x <= 0)
	{
		result = std::exp(x);
		result = result / (1 + result);
	}
	else
	{
		result = std::exp(-x);
		result = 1.0 / (1 + result);
	}
	return result;
}

std::vector<double> expit(const std::vector<double> x)
{
	std::vector<double> result(x);
	int len_x = x.size();
	for (int ii = 0; ii < len_x; ii++)
	{
		result[ii] = expit(x[ii]);
	}
	return result;
}

double logit(const double x)
{
	double result = NAN;
	if ((x < 1) & (x > 0))
	{
		result = log(x / (1 - x));
	}
	return result;
}

std::vector<double> logit(const std::vector<double> x)
{
	std::vector<double> result(x);
	int len_x = x.size();
	for (int ii = 0; ii < len_x; ii++)
	{
		result[ii] = logit(x[ii]);
	}
	return result;
}
////////////////////////////////END IMPLEMENTATION///////////////////////////////////////////////////////

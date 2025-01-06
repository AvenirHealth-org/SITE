#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define BOOST_DISABLE_ASSERTS
#define _ITERATOR_DEBUG_LEVEL 0

#include <boost/multi_array.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <fstream>  
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <stdint.h>

#define ARMA_DONT_PRINT_ERRORS

#define SSnegative 0 //Susceptible sero-negative
#define IIncubates 1 //Infected, incubation stage
#define IActiveSSy 2 //Infected, primary syphilis
#define ILatentSSy 3 //Infected, latent
#define RecovEarly 4 //Recovery after treatment
#define RecovLater 5 //Recovery without treatment
#define RecovSusce 6 //Recovery without treatment
#define numbStages 7 //Number of stages
#define Syph_Active 0//diagnosed/treated individuals 
#define Syph_NonAct 1//Undiagnosed/untreated individuals
//Sexual beharviour Constants
#define SyphLow 0//Low risk/1 partner
#define SyphMed 1//Medium risk
#define SyphHig 2//High risk
#define SyphNoSex 4//Men  Not Sexually active. Replacing IDU?
#define SyphNoSexWom 3//Women  Not Sexually active. Replacing IDU?
/*Men specific constants*/
#define SyphMSM 3//MSM
/*WoMen specific constants*/
#define SyphFSW 2//Female sex workers: Same as high risk women
/*Other constants*/
#define SyphNumRisksMen 5//Number of risk groups for men
#define SyphNumRisksWom 4//Number of risk groups for women
#define SyphWomen 0
#define SyphMen 1
/*Syphilis Stages Constants*/
#define RPRnegTPHApos 0 // RPR- TPHA+
#define RPRposTPHApos 1 // RPR+ TPHA+
#define TPHApos 2 //TPHA+
#define RPRpos 3  //RPR+
#define numPrevTypes 4//Prevalence depending on the test?

#define TPHATest 0
#define RPRTest  1
#define numTests 2
/*Contact Tracing Constants*/
#define all_contactTraced 0
#define prim_and_sec_contactTraced 1
#define latent_contactTraced 2
#define sizeContactTraced 3
/*Some useful objects*/
typedef boost::multi_array<int, 1> Array1d_int;
typedef boost::multi_array<int, 2> Array2d_int;
typedef boost::multi_array<double, 1> Array1d;
typedef boost::multi_array<double, 2> Array2d;
typedef boost::multi_array<double, 3> Array3d;
typedef boost::multi_array<double, 4> Array4d;
typedef boost::multi_array<double, 5> Array5d;

class NaturalHistParam
{
public:
	double sigma1;
	double sigma2;
	double sigma3;
	double sigma4;
	double sigma5;
	double phi;//proportion of prim_sec (Active) cases treated that become RPR-negative immediately
	double psi;// Probability of cure when treated
	double nu;// annual rate of treatment based on symptoms
	std::vector<double> vect_tau;
	std::vector<double> vect_nu;
	std::vector<double> screen_rate;//alpha in John's code
	std::vector<double> screen_rateTPHA;//some people will be screened with TPHA only?
	std::vector<double> r_ActiveToSus;
	std::vector<double> r_ActiveRecov;
	std::vector<double> r_LatentRecov;
	std::vector<double> prop_part_referred;
	double MtoF_tp;
	double FtoM_tp;
	double MtoM_tp;

	double condomeff;
	double ValinitInci;
	std::vector<double> SteadyState;

	NaturalHistParam& operator = (const NaturalHistParam &t)
	{
		this->sigma1 = t.sigma1;
		this->sigma2 = t.sigma2;
		this->sigma3 = t.sigma3;
		this->sigma4 = t.sigma4;
		this->sigma5 = t.sigma5;
		this->nu = t.nu;
		this->phi = t.phi;
		this->psi = t.psi;
		this->MtoF_tp = t.MtoF_tp;
		this->FtoM_tp = t.FtoM_tp;
		this->MtoM_tp = t.MtoM_tp;
		//this->OR_nu = t.OR_nu;
		this->vect_tau = t.vect_tau;
		this->vect_nu = t.vect_nu;
		this->SteadyState = t.SteadyState;
		this->screen_rate = t.screen_rate;
		this->screen_rateTPHA = t.screen_rateTPHA;
		this->r_ActiveToSus=t.r_ActiveToSus;
		this->r_ActiveRecov=t.r_ActiveRecov;
		this->r_LatentRecov= t.r_LatentRecov;
		this->condomeff = t.condomeff;
		this->ValinitInci = t.ValinitInci;
		this->prop_part_referred = t.prop_part_referred;

		return *this;
	}

	//void SetIt(double fsigma1, double fsigma2, double fsigma3, double fsigma4, double fsigma5, double in_nu, double fscreen_rate, double in_phi, double in_psi, double in_MtoF_tp, double in_FtoM_tp, double in_MtoM_tp);
	void CalcSteadyState();//Estimates initial proportions fo the population in each compartment
	void resize(int numSteps);
	void updateTreatment();
	NaturalHistParam();
	~NaturalHistParam();
};

NaturalHistParam::NaturalHistParam()
{
	sigma1 = 1.0/(4.4 / 52.0); //1/sigma1 is the duration of the incubation period
	sigma2 = 1.0 / ((6.6 + 15.6)/52); // 52; //1/sigma2 is the duration of the active syphilis phase
	sigma3 = 1.0 / (520.0 / 52.0); //1/sigma3 is the duration of the latent  phase
	sigma4 = 1.0 / (26.0 / 52.0); //1/sigma4 is the duration of the recovery phase, after treatment
	sigma5 = 1.0 / (52.0 / 52.0); //1/sigma5 is the duration of the recovery phase, after the latent phase
	psi = 0.99; //effectiveness of treatment
	phi = 0.6; // proportion sero-negative after treatment
	nu = 0.63; //treatment seeking

	ValinitInci = 1e-6;

	screen_rate.resize(1); 
	screen_rate[0] = 0;

	screen_rateTPHA.resize(1);
	screen_rateTPHA[0] = 0;

	vect_tau.resize(1);
	vect_tau[0] = 0.35;

	vect_nu.resize(1);
	vect_nu[0] = nu * vect_tau[0];

	r_ActiveToSus.resize(1);
	r_ActiveToSus[0] = (vect_nu[0] + screen_rate[0] + screen_rateTPHA[0])*psi*phi;

	r_ActiveRecov.resize(1);
	r_ActiveRecov[0] = (vect_nu[0] + screen_rate[0] + screen_rateTPHA[0])*psi*(1-phi);

	r_LatentRecov.resize(1);
	r_LatentRecov[0] = (screen_rate[0] + screen_rateTPHA[0]) * psi;

	prop_part_referred.resize(1);
	prop_part_referred[0] = 0;

	condomeff = 0.8;

	MtoF_tp = 0.25; //Male to female transmission probability
	FtoM_tp = 0.15; //Female to male transmission probability
	MtoM_tp = 0.25;
	SteadyState.resize(numbStages);
}

NaturalHistParam::~NaturalHistParam()
{

}

void NaturalHistParam::CalcSteadyState()
{
	double lambda = ValinitInci;
	double x = r_ActiveToSus[0];
	double y = r_ActiveRecov[0];
	double z = (screen_rate[0] + screen_rateTPHA[0] )* psi;
	double mu = 1.0 / 35.0;
	double numSus = 1.0;

	double be = (lambda + mu)*numSus;
	double cteR3e = 1.0 - lambda/(lambda+mu)*sigma1/(sigma1+mu)*sigma2/(x+y+sigma2+mu)*sigma3/(sigma3+z+mu)*sigma5/(sigma5+mu)
		- lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*sigma4 / (sigma4 + mu)
		- lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*sigma4 / (sigma4 + mu)
		- lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*x / (x + y + sigma2 + mu);

	double R3e = 1.0 / cteR3e;
	R3e *= be * (1.0/(lambda+mu)*lambda/(lambda+mu)*sigma1/(sigma1+mu)*sigma2/(x+y+sigma2+mu)*sigma3/(sigma3+z+mu)*sigma5/(sigma5+mu)
		+ 1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*sigma4 / (sigma4 + mu)
		+ 1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*sigma4 / (sigma4 + mu)
		+ 1.0 / (lambda + mu)*lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*x / (x + y + sigma2 + mu));

	double R2e = R3e * (lambda*sigma1/(sigma1+mu)*y/(x+y+sigma2+mu)*1.0/(sigma4+mu)
		+lambda*sigma1/(sigma1+mu)*sigma2/(x+y+sigma2+mu)*z/(sigma3+z+mu)*1.0/(sigma4+mu));
	R2e += be * (lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*y / (x + y + sigma2 + mu)*1.0 / (sigma4 + mu)
		+ lambda / (lambda + mu)*sigma1 / (sigma1 + mu)*sigma2 / (x + y + sigma2 + mu)*z / (sigma3 + z + mu)*1.0 / (sigma4 + mu));

	double R1e = lambda/(lambda+mu)*sigma1/(sigma1+mu)*sigma2/(x+y+sigma2+mu)*sigma3/(sigma3+z+mu)*1.0/(sigma5+mu)*be
		+lambda*sigma1/(sigma1+mu)*sigma2/(x+y+sigma2+mu)*sigma3/(sigma3+z+mu)*1.0/(sigma5+mu)*R3e;

	double Le = lambda/(lambda+mu)*sigma1/(sigma1+mu)*sigma2/(x+y+sigma2+mu)*1.0/(sigma3+z+mu)*be;
	Le += lambda * sigma1/ (sigma1 + mu) *sigma2 / (x + y + sigma2 + mu)*1.0 / (z + sigma3 + mu)*R3e;

	double Ae =  lambda / (lambda + mu)* sigma1 / (sigma1 + mu) *1.0/ (x + y + sigma2 + mu)*be;
	Ae += lambda * sigma1 / (sigma1 + mu)*1.0 / (x + y + sigma2 + mu)*R3e;

	double Ie = lambda / (lambda + mu)*1.0/ (sigma1 + mu) *be + lambda / (sigma1 + mu)*R3e;
	double total = numSus + Ie + Ae + Le + R1e + R2e + R3e;

	SteadyState[SSnegative] = numSus / total;
	SteadyState[IIncubates] = Ie / total;
	SteadyState[IActiveSSy] = Ae / total;
	SteadyState[ILatentSSy] = Le / total;
	SteadyState[RecovEarly] = R1e / total;
	SteadyState[RecovLater] = R2e / total;
	SteadyState[RecovSusce] = R3e / total;
}

void NaturalHistParam::resize(int numSteps)
{
	double scr = screen_rate[0];
	double scrTPHA = screen_rateTPHA[0];
	double AToSus = r_ActiveToSus[0];
	double ARecov = r_ActiveRecov[0];
	double LRecov = r_LatentRecov[0];
	double vnu = vect_nu[0];
	double vtau = vect_tau[0];
	double vprop_part_ref = prop_part_referred[0];

	vect_nu.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		vect_nu[ii] = vnu;
	}

	vect_tau.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		vect_tau[ii] = vtau;
	}

	screen_rate.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		screen_rate[ii] = scr;
	}

	screen_rateTPHA.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		screen_rateTPHA[ii] = scrTPHA;
	}

	r_ActiveToSus.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		r_ActiveToSus[ii] = AToSus;
	}

	r_ActiveRecov.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		r_ActiveRecov[ii] = ARecov;
	}

	r_LatentRecov.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		r_LatentRecov[ii] = LRecov;
	}

	prop_part_referred.resize(numSteps);
	for (int ii = 0; ii < numSteps; ii++)
	{
		prop_part_referred[ii] = vprop_part_ref;
	}
}

void NaturalHistParam::updateTreatment()
{
	int len_sr = screen_rate.size();

	for (int ii = 0; ii < len_sr; ii++)
	{
		vect_nu[ii] = nu * vect_tau[ii]*sigma2;
		r_ActiveToSus[ii] = (vect_nu[ii] + screen_rate[ii]+ screen_rateTPHA[ii])*psi*phi;
		r_ActiveRecov[ii] = (vect_nu[ii] + screen_rate[ii] + screen_rateTPHA[ii])*psi*(1 - phi);
		r_LatentRecov[ii] = (screen_rate[ii] + screen_rateTPHA[ii]) * psi;
	}
}

//*Mixing Matrix*//
class MixingMatrixParam
{
public:
	Array2d Mix_WomenMen;
	Array2d Mix_MenMen;
	Array2d Mix_MenWomen;
	Array2d Mix_WoMenWoMen;

	MixingMatrixParam& operator = (const MixingMatrixParam &t)
	{
		this->Mix_WomenMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
		this->Mix_WomenMen = boost::multi_array_ref<double, 2>(t.Mix_WomenMen);

		this->Mix_MenWomen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksWom]);
		this->Mix_MenWomen = boost::multi_array_ref<double, 2>(t.Mix_MenWomen);

		this->Mix_MenMen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksMen]);
		this->Mix_MenMen = boost::multi_array_ref<double, 2>(t.Mix_MenMen);

		this->Mix_WoMenWoMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksWom]);
		this->Mix_WoMenWoMen = boost::multi_array_ref<double, 2>(t.Mix_WoMenWoMen);

		return *this;
	}

	void SetIt(std::vector<std::vector<double>> in_matrix);
	MixingMatrixParam();
	~MixingMatrixParam();
};

MixingMatrixParam::MixingMatrixParam()
{
	/*Choice of partners by women*/
	Mix_WomenMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
	Mix_WomenMen[SyphLow][SyphLow] = 0.90;
	Mix_WomenMen[SyphLow][SyphMed] = 0.04;
	Mix_WomenMen[SyphLow][SyphHig] = 0.03;
	Mix_WomenMen[SyphLow][SyphMSM] = 0.03;
	Mix_WomenMen[SyphLow][SyphNoSex] = 0.0;

	Mix_WomenMen[SyphMed][SyphLow] = 0.25;
	Mix_WomenMen[SyphMed][SyphMed] = 0.4;
	Mix_WomenMen[SyphMed][SyphHig] = 0.3;
	Mix_WomenMen[SyphMed][SyphMSM] = 0.05;
	Mix_WomenMen[SyphMed][SyphNoSex] = 0.0;

	Mix_WomenMen[SyphHig][SyphLow] = 0.05;
	Mix_WomenMen[SyphHig][SyphMed] = 0.30;
	Mix_WomenMen[SyphHig][SyphHig] = 0.6;
	Mix_WomenMen[SyphHig][SyphMSM] = 0.05;
	Mix_WomenMen[SyphHig][SyphNoSex] = 0.0;

	/*Choice of partners by men*/
	Mix_MenWomen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksWom]);
	Mix_MenWomen[SyphLow][SyphLow] = 0.9;
	Mix_MenWomen[SyphLow][SyphMed] = 0.04;
	Mix_MenWomen[SyphLow][SyphHig] = 0.06;
	Mix_MenWomen[SyphLow][SyphNoSexWom] = 0.0;

	Mix_MenWomen[SyphMed][SyphLow] = 0.25;
	Mix_MenWomen[SyphMed][SyphMed] = 0.4;
	Mix_MenWomen[SyphMed][SyphHig] = 0.35;
	Mix_MenWomen[SyphMed][SyphNoSexWom] = 0.0;

	Mix_MenWomen[SyphHig][SyphLow] = 0.1;
	Mix_MenWomen[SyphHig][SyphMed] = 0.5;
	Mix_MenWomen[SyphHig][SyphHig] = 0.4;
	Mix_MenWomen[SyphHig][SyphNoSex] = 0.0;

	Mix_MenWomen[SyphMSM][SyphLow] = 0.01; //These men don't choose, they take everything
	Mix_MenWomen[SyphMSM][SyphMed] = 0.01;
	Mix_MenWomen[SyphMSM][SyphHig] = 0.08;
	Mix_MenWomen[SyphMSM][SyphNoSexWom] = 0;

	/*Choice of Partners of MSM: they can only choose other MSM!*/
	Mix_MenMen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksMen]);
	Mix_MenMen[SyphLow][SyphLow] = 0;
	Mix_MenMen[SyphLow][SyphMed] = 0;
	Mix_MenMen[SyphLow][SyphHig] = 0;
	Mix_MenMen[SyphLow][SyphMSM] = 0;
	Mix_MenMen[SyphLow][SyphNoSex] = 0;

	Mix_MenMen[SyphMed][SyphLow] = 0;
	Mix_MenMen[SyphMed][SyphMed] = 0;
	Mix_MenMen[SyphMed][SyphHig] = 0;
	Mix_MenMen[SyphMed][SyphMSM] = 0;
	Mix_MenMen[SyphMed][SyphNoSex] = 0;

	Mix_MenMen[SyphHig][SyphLow] = 0;
	Mix_MenMen[SyphHig][SyphMed] = 0;
	Mix_MenMen[SyphHig][SyphHig] = 0;
	Mix_MenMen[SyphHig][SyphMSM] = 0;
	Mix_MenMen[SyphHig][SyphNoSex] = 0;

	Mix_MenMen[SyphMSM][SyphLow] = 0;
	Mix_MenMen[SyphMSM][SyphMed] = 0;
	Mix_MenMen[SyphMSM][SyphHig] = 0;
	Mix_MenMen[SyphMSM][SyphMSM] = 1 - Mix_MenWomen[SyphMSM][SyphLow] - Mix_MenWomen[SyphMSM][SyphMed] - Mix_MenWomen[SyphMSM][SyphHig];// -Mix_MenWomen[SyphMSM][SyphFSW];
	Mix_MenMen[SyphMSM][SyphMSM] = (Mix_MenMen[SyphMSM][SyphMSM] <= 0) ? 0 : Mix_MenMen[SyphMSM][SyphMSM];
	Mix_MenMen[SyphMSM][SyphNoSex] = 0;

	Mix_MenMen[SyphNoSex][SyphLow] = 0;
	Mix_MenMen[SyphNoSex][SyphMed] = 0;
	Mix_MenMen[SyphNoSex][SyphHig] = 0;
	Mix_MenMen[SyphNoSex][SyphMSM] = 0;
	Mix_MenMen[SyphNoSex][SyphNoSex] = 0;

	/*Mixing of women and women*/
	Mix_WoMenWoMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
	Mix_WoMenWoMen[SyphLow][SyphLow] = 0;
	Mix_WoMenWoMen[SyphLow][SyphMed] = 0;
	Mix_WoMenWoMen[SyphLow][SyphHig] = 0;
	Mix_WoMenWoMen[SyphLow][SyphNoSexWom] = 0;

	Mix_WoMenWoMen[SyphMed][SyphLow] = 0;
	Mix_WoMenWoMen[SyphMed][SyphMed] = 0;
	Mix_WoMenWoMen[SyphMed][SyphHig] = 0;
	Mix_WoMenWoMen[SyphMed][SyphFSW] = 0;
	Mix_WoMenWoMen[SyphMed][SyphNoSexWom] = 0;

	Mix_WoMenWoMen[SyphHig][SyphLow] = 0;
	Mix_WoMenWoMen[SyphHig][SyphMed] = 0;
	Mix_WoMenWoMen[SyphHig][SyphHig] = 0;
	Mix_WoMenWoMen[SyphHig][SyphNoSexWom] = 0;
}

MixingMatrixParam::~MixingMatrixParam()
{
	;
}

void MixingMatrixParam::SetIt(std::vector<std::vector<double>> in_matrix)
{
	bool checker = true;
	int len_men = in_matrix.size();
	int len_women = 0;
	if (len_men == (SyphNumRisksMen + SyphNumRisksWom))
	{
		len_women = in_matrix[0].size();
		if (len_women != (SyphNumRisksMen + SyphNumRisksWom))
		{
			checker = false;
		}
	}
	else
	{
		checker = false;
	}
	if (checker)
	{
		
		std::vector<std::vector<double>> new_in_matrix= in_matrix;
		for (int ii = 0; ii < SyphNumRisksMen; ii++)
		{
			for (int jj = 0; jj < SyphNumRisksMen; jj++)
			{
				if ((jj != SyphMSM) |(ii!=jj))
				{
					new_in_matrix[ii][jj] = 0; //Only MSM can mix with MSM. Other men mix with women exclusively
				}
			}

			for (int jj = SyphNumRisksMen; jj < (SyphNumRisksMen + SyphNumRisksWom); jj++)
			{
				if (new_in_matrix[ii][jj]<0)
				{
					new_in_matrix[ii][jj] = 0;
					Rcpp::Rcout << "Warning! The coefficient i=" << ii << ", j=" << jj << " of the mixing matrix was negative, it will be set to zero." << std::endl;
				}
			}
		}

		for (int ii = SyphNumRisksMen; ii < (SyphNumRisksMen + SyphNumRisksWom); ii++)
		{
			for (int jj = 0; jj < SyphNumRisksMen; jj++)
			{
				if (new_in_matrix[ii][jj] < 0)
				{
					new_in_matrix[ii][jj] = 0;
					Rcpp::Rcout << "Warning! The coefficient i=" << ii << ", j=" << jj << " of the mixing matrix was negative, it will be set to zero." << std::endl;
				}
			}

			for (int jj = SyphNumRisksMen; jj < (SyphNumRisksMen + SyphNumRisksWom); jj++)
			{
				new_in_matrix[ii][jj] = 0;//Women don't mix with women because we don't have IDUs
			}
		}
		
		/*Choice of partners by women*/
		Mix_WomenMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
		Mix_WomenMen[SyphLow][SyphLow] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphLow];
		Mix_WomenMen[SyphLow][SyphLow] = 1;//
		Mix_WomenMen[SyphLow][SyphMed] = 1;//new_in_matrix[SyphNumRisksMen + SyphLow][SyphMed];
		Mix_WomenMen[SyphLow][SyphHig] = 1;//new_in_matrix[SyphNumRisksMen + SyphLow][SyphHig];
		Mix_WomenMen[SyphLow][SyphMSM] = 1;//new_in_matrix[SyphNumRisksMen + SyphLow][SyphMSM];
		Mix_WomenMen[SyphLow][SyphNoSex] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphNoSex];

		Mix_WomenMen[SyphMed][SyphLow] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphLow];
		//Mix_WomenMen[SyphMed][SyphMed] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphMed];
		Mix_WomenMen[SyphMed][SyphMed] = 1;
		Mix_WomenMen[SyphMed][SyphHig] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphHig];
		Mix_WomenMen[SyphMed][SyphMSM] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphMSM];
		Mix_WomenMen[SyphMed][SyphNoSex] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphNoSex];

		Mix_WomenMen[SyphHig][SyphLow] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphLow];
		Mix_WomenMen[SyphHig][SyphMed] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphMed];
		//Mix_WomenMen[SyphHig][SyphHig] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphHig];
		Mix_WomenMen[SyphHig][SyphHig] = 1.0;
		Mix_WomenMen[SyphHig][SyphMSM] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphMSM];
		Mix_WomenMen[SyphHig][SyphNoSex] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphNoSex];

		/*Choice of partners by men*/
		Mix_MenWomen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksWom]);
		//Mix_MenWomen[SyphLow][SyphLow] = new_in_matrix[SyphLow][SyphNumRisksMen + SyphLow];
		Mix_MenWomen[SyphLow][SyphLow] = 1.0;//
		Mix_MenWomen[SyphLow][SyphMed] = 1.0;// new_in_matrix[SyphLow][SyphNumRisksMen + SyphMed];
		Mix_MenWomen[SyphLow][SyphHig] = 1.0;//new_in_matrix[SyphLow][SyphNumRisksMen + SyphHig];
		Mix_MenWomen[SyphLow][SyphNoSexWom] = new_in_matrix[SyphLow][SyphNumRisksMen + SyphNoSexWom];

		Mix_MenWomen[SyphMed][SyphLow] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphLow];
		Mix_MenWomen[SyphMed][SyphMed] = 1.0;// new_in_matrix[SyphMed][SyphNumRisksMen + SyphMed];
		Mix_MenWomen[SyphMed][SyphHig] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphHig];
		Mix_MenWomen[SyphMed][SyphNoSexWom] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphNoSexWom];

		Mix_MenWomen[SyphHig][SyphLow] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphLow];
		Mix_MenWomen[SyphHig][SyphMed] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphMed];
		Mix_MenWomen[SyphHig][SyphHig] = 1.0;// new_in_matrix[SyphHig][SyphNumRisksMen + SyphHig];
		Mix_MenWomen[SyphHig][SyphNoSexWom] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphNoSexWom];

		Mix_MenWomen[SyphMSM][SyphLow] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphLow]; //These men don't choose, they take everything
		Mix_MenWomen[SyphMSM][SyphMed] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphMed];
		Mix_MenWomen[SyphMSM][SyphHig] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphHig];
		Mix_MenWomen[SyphMSM][SyphNoSexWom] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphNoSexWom];

		/*Choice of Partners of MSM: they can only choose other MSM!*/
		Mix_MenMen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksMen]);
		Mix_MenMen[SyphLow][SyphLow] = 0;
		Mix_MenMen[SyphLow][SyphMed] = 0;
		Mix_MenMen[SyphLow][SyphHig] = 0;
		Mix_MenMen[SyphLow][SyphMSM] = 0;
		Mix_MenMen[SyphLow][SyphNoSex] = 0;

		Mix_MenMen[SyphMed][SyphLow] = 0;
		Mix_MenMen[SyphMed][SyphMed] = 0;
		Mix_MenMen[SyphMed][SyphHig] = 0;
		Mix_MenMen[SyphMed][SyphMSM] = 0;
		Mix_MenMen[SyphMed][SyphNoSex] = 0;

		Mix_MenMen[SyphHig][SyphLow] = 0;
		Mix_MenMen[SyphHig][SyphMed] = 0;
		Mix_MenMen[SyphHig][SyphHig] = 0;
		Mix_MenMen[SyphHig][SyphMSM] = 0;
		Mix_MenMen[SyphHig][SyphNoSex] = 0;

		Mix_MenMen[SyphMSM][SyphLow] = 0;
		Mix_MenMen[SyphMSM][SyphMed] = 0;
		Mix_MenMen[SyphMSM][SyphHig] = 0;
		Mix_MenMen[SyphMSM][SyphMSM] = 1;// new_in_matrix[SyphMSM][SyphMSM];
		Mix_MenMen[SyphMSM][SyphNoSex] = 0;

		Mix_MenMen[SyphNoSex][SyphLow] = 0;
		Mix_MenMen[SyphNoSex][SyphMed] = 0;
		Mix_MenMen[SyphNoSex][SyphHig] = 0;
		Mix_MenMen[SyphNoSex][SyphMSM] = 0;
		Mix_MenMen[SyphNoSex][SyphNoSex] = 0;

		/*Mixing of women and women*/
		Mix_WoMenWoMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
		Mix_WoMenWoMen[SyphLow][SyphLow] = 0;
		Mix_WoMenWoMen[SyphLow][SyphMed] = 0;
		Mix_WoMenWoMen[SyphLow][SyphHig] = 0;
		Mix_WoMenWoMen[SyphLow][SyphNoSexWom] = 0;

		Mix_WoMenWoMen[SyphMed][SyphLow] = 0;
		Mix_WoMenWoMen[SyphMed][SyphMed] = 0;
		Mix_WoMenWoMen[SyphMed][SyphHig] = 0;
		Mix_WoMenWoMen[SyphMed][SyphFSW] = 0;
		Mix_WoMenWoMen[SyphMed][SyphNoSexWom] = 0;

		Mix_WoMenWoMen[SyphHig][SyphLow] = 0;
		Mix_WoMenWoMen[SyphHig][SyphMed] = 0;
		Mix_WoMenWoMen[SyphHig][SyphHig] = 0;
		Mix_WoMenWoMen[SyphHig][SyphNoSexWom] = 0;
	}
}

///*Referral Probabilities*///
class ReferralParam
{
public:
	Array2d Refer_WomenMen;
	Array2d Refer_MenMen;
	Array2d Refer_MenWomen;
	Array2d Refer_WoMenWoMen;

	ReferralParam& operator = (const ReferralParam &t)
	{
		this->Refer_WomenMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
		this->Refer_WomenMen = boost::multi_array_ref<double, 2>(t.Refer_WomenMen);

		this->Refer_MenWomen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksWom]);
		this->Refer_MenWomen = boost::multi_array_ref<double, 2>(t.Refer_MenWomen);

		this->Refer_MenMen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksMen]);
		this->Refer_MenMen = boost::multi_array_ref<double, 2>(t.Refer_MenMen);

		this->Refer_WoMenWoMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksWom]);
		this->Refer_WoMenWoMen = boost::multi_array_ref<double, 2>(t.Refer_WoMenWoMen);

		return *this;
	}

	void SetIt(std::vector<std::vector<double>> in_matrix);
	ReferralParam();
	~ReferralParam();
};

ReferralParam::ReferralParam()
{
	/*Referral of partners following symptoms by women*/
	Refer_WomenMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
	Refer_WomenMen[SyphLow][SyphLow] = 1;
	Refer_WomenMen[SyphLow][SyphMed] = 0.0;
	Refer_WomenMen[SyphLow][SyphHig] = 0.0;
	Refer_WomenMen[SyphLow][SyphMSM] = 0.0;
	Refer_WomenMen[SyphLow][SyphNoSex] = 0;

	Refer_WomenMen[SyphMed][SyphLow] = 0;
	Refer_WomenMen[SyphMed][SyphMed] = 0;
	Refer_WomenMen[SyphMed][SyphHig] = 0;
	Refer_WomenMen[SyphMed][SyphMSM] = 0;
	Refer_WomenMen[SyphMed][SyphNoSex] = 0;

	Refer_WomenMen[SyphHig][SyphLow] = 0.0;
	Refer_WomenMen[SyphHig][SyphMed] = 0.0;
	Refer_WomenMen[SyphHig][SyphHig] = 0.0;
	Refer_WomenMen[SyphHig][SyphMSM] = 0.0;
	Refer_WomenMen[SyphHig][SyphNoSex] = 0.0;

	/*Referral of partners following symptoms by men*/
	Refer_MenWomen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksWom]);
	Refer_MenWomen[SyphLow][SyphLow] = 1;
	Refer_MenWomen[SyphLow][SyphMed] = 0.0;
	Refer_MenWomen[SyphLow][SyphHig] = 0.0;
	Refer_MenWomen[SyphLow][SyphNoSexWom] = 0.0;

	Refer_MenWomen[SyphMed][SyphLow] = 0.0;
	Refer_MenWomen[SyphMed][SyphMed] = 0.0;
	Refer_MenWomen[SyphMed][SyphHig] = 0.0;
	Refer_MenWomen[SyphMed][SyphNoSexWom] = 0.0;

	Refer_MenWomen[SyphHig][SyphLow] = 0.0;
	Refer_MenWomen[SyphHig][SyphMed] = 0.0;
	Refer_MenWomen[SyphHig][SyphHig] = 0.0;
	Refer_MenWomen[SyphHig][SyphNoSex] = 0.0;

	Refer_MenWomen[SyphMSM][SyphLow] = 1; //
	Refer_MenWomen[SyphMSM][SyphMed] = 0.0;
	Refer_MenWomen[SyphMSM][SyphHig] = 0.0;
	Refer_MenWomen[SyphMSM][SyphNoSexWom] = 0;

	/*Referral of Partners of MSM follwoing symptomatic Syphilis: they can only choose other MSM!*/
	Refer_MenMen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksMen]);
	Refer_MenMen[SyphLow][SyphLow] = 0;
	Refer_MenMen[SyphLow][SyphMed] = 0;
	Refer_MenMen[SyphLow][SyphHig] = 0;
	Refer_MenMen[SyphLow][SyphMSM] = 0;
	Refer_MenMen[SyphLow][SyphNoSex] = 0;

	Refer_MenMen[SyphMed][SyphLow] = 0;
	Refer_MenMen[SyphMed][SyphMed] = 0;
	Refer_MenMen[SyphMed][SyphHig] = 0;
	Refer_MenMen[SyphMed][SyphMSM] = 0;
	Refer_MenMen[SyphMed][SyphNoSex] = 0;

	Refer_MenMen[SyphHig][SyphLow] = 0;
	Refer_MenMen[SyphHig][SyphMed] = 0;
	Refer_MenMen[SyphHig][SyphHig] = 0;
	Refer_MenMen[SyphHig][SyphMSM] = 0;
	Refer_MenMen[SyphHig][SyphNoSex] = 0;

	Refer_MenMen[SyphMSM][SyphLow] = 0;
	Refer_MenMen[SyphMSM][SyphMed] = 0;
	Refer_MenMen[SyphMSM][SyphHig] = 0;
	Refer_MenMen[SyphMSM][SyphMSM] = 0.5;
	Refer_MenMen[SyphMSM][SyphNoSex] = 0;

	Refer_MenMen[SyphNoSex][SyphLow] = 0;
	Refer_MenMen[SyphNoSex][SyphMed] = 0;
	Refer_MenMen[SyphNoSex][SyphHig] = 0;
	Refer_MenMen[SyphNoSex][SyphMSM] = 0;
	Refer_MenMen[SyphNoSex][SyphNoSex] = 0;

	/*Referral of women by, following asymptomatic Syphilis*/
	Refer_WoMenWoMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
	Refer_WoMenWoMen[SyphLow][SyphLow] = 0;
	Refer_WoMenWoMen[SyphLow][SyphMed] = 0;
	Refer_WoMenWoMen[SyphLow][SyphHig] = 0;
	Refer_WoMenWoMen[SyphLow][SyphNoSexWom] = 0;

	Refer_WoMenWoMen[SyphMed][SyphLow] = 0;
	Refer_WoMenWoMen[SyphMed][SyphMed] = 0;
	Refer_WoMenWoMen[SyphMed][SyphHig] = 0;
	Refer_WoMenWoMen[SyphMed][SyphFSW] = 0;
	Refer_WoMenWoMen[SyphMed][SyphNoSexWom] = 0;

	Refer_WoMenWoMen[SyphHig][SyphLow] = 0;
	Refer_WoMenWoMen[SyphHig][SyphMed] = 0;
	Refer_WoMenWoMen[SyphHig][SyphHig] = 0;
	Refer_WoMenWoMen[SyphHig][SyphNoSexWom] = 0;
}

ReferralParam::~ReferralParam()
{
	;
}

void ReferralParam::SetIt(std::vector<std::vector<double>> in_matrix)
{
	bool checker = true;
	int len_men = in_matrix.size();
	int len_women = 0;
	if (len_men == (SyphNumRisksMen + SyphNumRisksWom))
	{
		len_women = in_matrix[0].size();
		if (len_women != (SyphNumRisksMen + SyphNumRisksWom))
		{
			checker = false;
		}
	}
	else
	{
		checker = false;
	}
	if (checker)
	{
		std::vector<std::vector<double>> new_in_matrix = in_matrix;
		for (int ii = 0; ii < SyphNumRisksMen; ii++)
		{
			for (int jj = 0; jj < SyphNumRisksMen; jj++)
			{
				if ((jj != SyphMSM) | (ii != jj))
				{
					new_in_matrix[ii][jj] = 0; //Only MSM can mix with MSM. Other men mix with women exclusively
				}
			}

			for (int jj = SyphNumRisksMen; jj < (SyphNumRisksMen + SyphNumRisksWom); jj++)
			{
				if (new_in_matrix[ii][jj] < 0)
				{
					new_in_matrix[ii][jj] = 0;
					Rcpp::Rcout << ">Warning! The coefficient i=" << ii << ", j=" << jj << " of the referral matrix was out of range, it will be set to zero." << std::endl;
				}
			}
		}

		for (int ii = SyphNumRisksMen; ii < (SyphNumRisksMen + SyphNumRisksWom); ii++)
		{
			for (int jj = 0; jj < SyphNumRisksMen; jj++)
			{
				if ((new_in_matrix[ii][jj] < 0) |(new_in_matrix[ii][jj] >1000))
				{
					new_in_matrix[ii][jj] = 0;
					Rcpp::Rcout << ">Warning! The coefficient i=" << ii << ", j=" << jj << " of the referral matrix was out of range, it will be set to zero." << std::endl;
				}
			}

			for (int jj = SyphNumRisksMen; jj < (SyphNumRisksMen + SyphNumRisksWom); jj++)
			{
				new_in_matrix[ii][jj] = 0;//Women don't mix with women because we don't have IDUs
			}
		}

		/*Rescaling the mixing matrix*/
		for (int ii = 0; ii < (SyphNumRisksMen + SyphNumRisksWom); ii++)
		{
			//Now, rescaling
			for (int jj = 0; jj < (SyphNumRisksMen + SyphNumRisksWom); jj++)
			{
				new_in_matrix[ii][jj] = (new_in_matrix[ii][jj] > 1) ? new_in_matrix[ii][jj] = 1.0 / 100.0 : new_in_matrix[ii][jj];
			}
		}

		/*Choice of partners by women*/
		Refer_WomenMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
		Refer_WomenMen[SyphLow][SyphLow] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphLow];
		Refer_WomenMen[SyphLow][SyphMed] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphMed];
		Refer_WomenMen[SyphLow][SyphHig] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphHig];
		Refer_WomenMen[SyphLow][SyphMSM] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphMSM];
		Refer_WomenMen[SyphLow][SyphNoSex] = new_in_matrix[SyphNumRisksMen + SyphLow][SyphNoSex];

		Refer_WomenMen[SyphMed][SyphLow] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphLow];
		Refer_WomenMen[SyphMed][SyphMed] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphMed];
		Refer_WomenMen[SyphMed][SyphHig] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphHig];
		Refer_WomenMen[SyphMed][SyphMSM] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphMSM];
		Refer_WomenMen[SyphMed][SyphNoSex] = new_in_matrix[SyphNumRisksMen + SyphMed][SyphNoSex];

		Refer_WomenMen[SyphHig][SyphLow] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphLow];
		Refer_WomenMen[SyphHig][SyphMed] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphMed];
		Refer_WomenMen[SyphHig][SyphHig] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphHig];
		Refer_WomenMen[SyphHig][SyphMSM] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphMSM];
		Refer_WomenMen[SyphHig][SyphNoSex] = new_in_matrix[SyphNumRisksMen + SyphHig][SyphNoSex];

		/*Referral of partners by men, following symptomatic Syphils*/
		Refer_MenWomen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksWom]);
		Refer_MenWomen[SyphLow][SyphLow] = new_in_matrix[SyphLow][SyphNumRisksMen + SyphLow];
		Refer_MenWomen[SyphLow][SyphMed] = new_in_matrix[SyphLow][SyphNumRisksMen + SyphMed];
		Refer_MenWomen[SyphLow][SyphHig] = new_in_matrix[SyphLow][SyphNumRisksMen + SyphHig];
		Refer_MenWomen[SyphLow][SyphNoSexWom] = new_in_matrix[SyphLow][SyphNumRisksMen + SyphNoSexWom];

		Refer_MenWomen[SyphMed][SyphLow] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphLow];
		Refer_MenWomen[SyphMed][SyphMed] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphMed];
		Refer_MenWomen[SyphMed][SyphHig] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphHig];
		Refer_MenWomen[SyphMed][SyphNoSexWom] = new_in_matrix[SyphMed][SyphNumRisksMen + SyphNoSexWom];

		Refer_MenWomen[SyphHig][SyphLow] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphLow];
		Refer_MenWomen[SyphHig][SyphMed] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphMed];
		Refer_MenWomen[SyphHig][SyphHig] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphHig];
		Refer_MenWomen[SyphHig][SyphNoSexWom] = new_in_matrix[SyphHig][SyphNumRisksMen + SyphNoSexWom];

		Refer_MenWomen[SyphMSM][SyphLow] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphLow];
		Refer_MenWomen[SyphMSM][SyphMed] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphMed];
		Refer_MenWomen[SyphMSM][SyphHig] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphHig];
		Refer_MenWomen[SyphMSM][SyphNoSexWom] = new_in_matrix[SyphMSM][SyphNumRisksMen + SyphNoSexWom];

		/*Choice of Partners of MSM: they can only choose other MSM!*/
		Refer_MenMen.resize(boost::extents[SyphNumRisksMen][SyphNumRisksMen]);
		Refer_MenMen[SyphLow][SyphLow] = 0;
		Refer_MenMen[SyphLow][SyphMed] = 0;
		Refer_MenMen[SyphLow][SyphHig] = 0;
		Refer_MenMen[SyphLow][SyphMSM] = 0;
		Refer_MenMen[SyphLow][SyphNoSex] = 0;

		Refer_MenMen[SyphMed][SyphLow] = 0;
		Refer_MenMen[SyphMed][SyphMed] = 0;
		Refer_MenMen[SyphMed][SyphHig] = 0;
		Refer_MenMen[SyphMed][SyphMSM] = 0;
		Refer_MenMen[SyphMed][SyphNoSex] = 0;

		Refer_MenMen[SyphHig][SyphLow] = 0;
		Refer_MenMen[SyphHig][SyphMed] = 0;
		Refer_MenMen[SyphHig][SyphHig] = 0;
		Refer_MenMen[SyphHig][SyphMSM] = 0;
		Refer_MenMen[SyphHig][SyphNoSex] = 0;

		Refer_MenMen[SyphMSM][SyphLow] = 0;
		Refer_MenMen[SyphMSM][SyphMed] = 0;
		Refer_MenMen[SyphMSM][SyphHig] = 0;
		Refer_MenMen[SyphMSM][SyphMSM] = new_in_matrix[SyphMSM][SyphMSM];
		Refer_MenMen[SyphMSM][SyphNoSex] = 0;

		Refer_MenMen[SyphNoSex][SyphLow] = 0;
		Refer_MenMen[SyphNoSex][SyphMed] = 0;
		Refer_MenMen[SyphNoSex][SyphHig] = 0;
		Refer_MenMen[SyphNoSex][SyphMSM] = 0;
		Refer_MenMen[SyphNoSex][SyphNoSex] = 0;

		/*Referral of women by women: virtually impossible*/
		Refer_WoMenWoMen.resize(boost::extents[SyphNumRisksWom][SyphNumRisksMen]);
		Refer_WoMenWoMen[SyphLow][SyphLow] = 0;
		Refer_WoMenWoMen[SyphLow][SyphMed] = 0;
		Refer_WoMenWoMen[SyphLow][SyphHig] = 0;
		Refer_WoMenWoMen[SyphLow][SyphNoSexWom] = 0;

		Refer_WoMenWoMen[SyphMed][SyphLow] = 0;
		Refer_WoMenWoMen[SyphMed][SyphMed] = 0;
		Refer_WoMenWoMen[SyphMed][SyphHig] = 0;
		Refer_WoMenWoMen[SyphMed][SyphFSW] = 0;
		Refer_WoMenWoMen[SyphMed][SyphNoSexWom] = 0;

		Refer_WoMenWoMen[SyphHig][SyphLow] = 0;
		Refer_WoMenWoMen[SyphHig][SyphMed] = 0;
		Refer_WoMenWoMen[SyphHig][SyphHig] = 0;
		Refer_WoMenWoMen[SyphHig][SyphNoSexWom] = 0;
	}
}
//*******************************End***************************************************//
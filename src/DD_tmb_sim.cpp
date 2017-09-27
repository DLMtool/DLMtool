#include <TMB.hpp>

template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
	using namespace density;

	//-------DATA-------//
	//So_DD
	DATA_SCALAR( So_DD );
	//Alpha_DD
	DATA_SCALAR( Alpha_DD );
	//Rho_DD
	DATA_SCALAR( Rho_DD );
	//ny_DD
	DATA_INTEGER( ny_DD );
	//k_DD
	DATA_INTEGER( k_DD );
	//wa_DD
	DATA_SCALAR( wa_DD );
	//E_hist
	DATA_VECTOR( E_hist );
	//C_hist
	DATA_VECTOR( C_hist );
	//Umsy - prior
	DATA_VECTOR( UMSYprior );
	
	//-------PARAMETERS-------//
	//Umsy
	PARAMETER( log_UMSY_DD );
	//MSY
	PARAMETER( log_MSY_DD );
	//q
	PARAMETER( log_q_DD );

	//--PARAMETER TRANSFORMATION
	Type UMSY_DD = exp(log_UMSY_DD);
	Type MSY_DD = exp(log_MSY_DD);
	Type q_DD = exp(log_q_DD);
	
	//-------FRONT END-------//
	//--DECLARING DERIVED VALUES
	Type SS_DD = So_DD * (1 - UMSY_DD);
	Type Spr_DD = (SS_DD * Alpha_DD/(1 - SS_DD) + wa_DD)/(1 - Rho_DD * SS_DD);
	Type DsprDu_DD = -So_DD * (Rho_DD/(1 - Rho_DD * SS_DD) * Spr_DD + 1/(1 - Rho_DD * SS_DD) * (Alpha_DD/(1 - SS_DD) + SS_DD * Alpha_DD/pow((1 - SS_DD),2)));
	Type Arec_DD = 1/((pow((1 - UMSY_DD),2)) * (Spr_DD + UMSY_DD * DsprDu_DD));
	Type Brec_DD = UMSY_DD * (Arec_DD * Spr_DD - 1/(1 - UMSY_DD))/MSY_DD;
	Type Spr0_DD = (So_DD * Alpha_DD/(1 - So_DD) + wa_DD)/(1 - Rho_DD * So_DD);
	Type Ro_DD = (Arec_DD * Spr0_DD - 1)/(Brec_DD * Spr0_DD);
	Type Bo_DD = Ro_DD * Spr0_DD;
	Type No_DD = Ro_DD/(1 - So_DD);

	//--DECLARING STORAGE VECTORS
	int ny_DDp = ny_DD+1;
	int ny_DDk = ny_DD+k_DD;
	vector<Type> B_DD(ny_DDp);
	vector<Type> N_DD(ny_DDp);
	vector<Type> R_DD(ny_DDk);

	//--INITIALIZE
	B_DD(0) = Bo_DD;
	N_DD(0) = No_DD;
	R_DD.fill(Ro_DD);

	//--STORAGE VECTORS
	vector<Type> Surv_DD(ny_DD);
	vector<Type> Cpred_DD(ny_DD);
	vector<Type> Sp_DD(ny_DD);

	double tiny = 1e-15;

	for(int tt=0; tt<ny_DD; tt++){
		Surv_DD(tt) = So_DD * exp(-q_DD * E_hist(tt));
		Cpred_DD(tt) = B_DD(tt) * (1 - exp(-q_DD * E_hist(tt)));
		//Cpred TRAP
		if(Cpred_DD(tt) < tiny) Cpred_DD(tt) = tiny;
		Sp_DD(tt) = B_DD(tt) - Cpred_DD(tt);
		R_DD[tt + k_DD] = Arec_DD * Sp_DD(tt)/(1 + Brec_DD * Sp_DD(tt));
		B_DD[tt + 1] = Surv_DD(tt) * (Alpha_DD * N_DD(tt) + Rho_DD * B_DD(tt)) + wa_DD * R_DD[tt + 1];
		N_DD[tt + 1] = Surv_DD(tt) * N_DD(tt) + R_DD[tt + 1];
	}

	//--ARGUMENTS FOR NLL
	vector<Type> lCpred_DD = log(Cpred_DD);
	vector<Type> lC_hist = log(C_hist);
	Type lUMSYprior = log(UMSYprior(0));

	// Objective function
	//creates storage for jnll and sets value to 0
	vector<Type> jnll_comp(2); 
	jnll_comp.setZero();
	
	for(int tt=0; tt<ny_DD; tt++){
		jnll_comp(0) -= dnorm(lCpred_DD(tt), lC_hist(tt), Type(0.25), true);
	}
	jnll_comp(1) -= dlognorm(UMSY_DD, lUMSYprior, UMSYprior(1), true);

	//Summing individual jnll
	Type jnll = jnll_comp.sum();
	
	//-------REPORTING-------//
	REPORT( UMSY_DD );
	REPORT( MSY_DD );
	REPORT( q_DD );
	REPORT( jnll_comp );

	ADREPORT( log_UMSY_DD );
	ADREPORT( log_MSY_DD );
	ADREPORT( log_q_DD );
	
	return jnll;
}

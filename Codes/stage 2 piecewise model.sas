PROC IMPORT OUT= work.data_kid 
            DATAFILE= "C:\Users\Jonathan Yu\Dropbox\Aim 1\Data\kidney_cleaned.csv" 
            DBMS=CSV REPLACE;
     * RANGE="data$"; 
     GETNAMES=YES;
    * MIXED=NO;
    * SCANTEXT=YES;
    * USEDATE=YES;
    * SCANTIME=YES;
RUN;

/****************************************************************************************/
/*                                 Patient Survival                                     */
/****************************************************************************************/

/***--- NON-CLUSTERED EVENTS WITHOUT CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.data_kid;
proc nlmixed data=four qpoints=170 noad;
	parms beta0 = 7.7 beta1 = 0.98189 beta2 = 0.98189 beta3 = 0.98189
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1  gamma3 = 0.1  
		  theta=4.7 
		  r01-r10 1e-4;

	bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 theta  >= 0;

	base_haz = r01*fail_r1 + r02*fail_r2 + r03*fail_r3 + r04*fail_r4 + r05*fail_r5
		+ r06*fail_r6 + r07*fail_r7 + r08*fail_r8 + r09*fail_r9 + r10*fail_r10; 
	cum_base_haz = r01*dur_r1 + r02*dur_r2 + r03*dur_r3 + r04*dur_r4 + r05*dur_r5
		+ r06*dur_r6 + r07*dur_r7 + r08*dur_r8 + r09*dur_r9 + r10*dur_r10;
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4) ;
	ll = -cum_base_haz*exp(mu);

	if pstatus=0 then loglik = ll;
	if pstatus=1 then loglik = ll + log(base_haz) + mu;

	model ptime~general(loglik*cbps_weight);

	ods exclude iterhistory parameters;
run;
proc print data=work.fit1 label; 
      format value 8.2; 
run;

/***--- NON-CLUSTERED EVENTS WITH CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.data_kid;
proc nlmixed data=four qpoints=170 noad;
	parms beta0  = 7.7 beta1  = 0.98189 beta2  = 0.98189 beta3  = 0.98189
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1     gamma3 = 0.1  
		  theta=4.7 
		  r01-r10 1e-4;

	bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 theta  >= 0;


	base_haz = r01*fail_r1 + r02*fail_r2 + r03*fail_r3 + r04*fail_r4 + r05*fail_r5
		+ r06*fail_r6 + r07*fail_r7 + r08*fail_r8 + r09*fail_r9 + r10*fail_r10; 
	cum_base_haz = r01*dur_r1 + r02*dur_r2 + r03*dur_r3 + r04*dur_r4 + r05*dur_r5
		+ r06*dur_r6 + r07*dur_r7 + r08*dur_r8 + r09*dur_r9 + r10*dur_r10;

	* cure *;
	check = gamma0 + (gamma1 * (hcvdr_num=3)) +  (gamma2 * (hcvdr_num=2)) + (gamma3 * (hcvdr_num=1));
	expcheck = exp(-check);
    prob= 1/(1+expcheck);
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4);
	ll = -cum_base_haz*exp(mu);

	if pstatus=0 then loglik = log(prob + (1-prob) * exp(ll)); 
	if pstatus=1 then loglik = log(1-prob) + log(base_haz) + ll + mu;
	
	model ptime~general(loglik*cbps_weight);

	ods exclude iterhistory parameters;
run;
proc print data=work.fit2 label; 
      format value 8.2; 
run;

/***--- CLUSTERED EVENTS WITHOUT CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.data_kid;
proc nlmixed data=four qpoints=170 noad;
	parms beta0 = 7.7 beta1 = 0.98189 beta2 = 0.98189 beta3 = 0.98189
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1  gamma3 = 0.1  
		  theta=4.7 
		  r01-r10 1e-4;

	bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 theta  >= 0;


	base_haz = r01*fail_r1 + r02*fail_r2 + r03*fail_r3 + r04*fail_r4 + r05*fail_r5
		+ r06*fail_r6 + r07*fail_r7 + r08*fail_r8 + r09*fail_r9 + r10*fail_r10; 
	cum_base_haz = r01*dur_r1 + r02*dur_r2 + r03*dur_r3 + r04*dur_r4 + r05*dur_r5
		+ r06*dur_r6 + r07*dur_r7 + r08*dur_r8 + r09*dur_r9 + r10*dur_r10;
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4) + nu;
	ll = -cum_base_haz*exp(mu);

	if pstatus=0 then loglik = ll;
	if pstatus=1 then loglik = ll + log(base_haz) + mu;
	model ptime~general(loglik*cbps_weight);
	random nu~normal(0,theta) subject=ctr_code;
	
	ods exclude iterhistory parameters;
run;
proc print data=work.fit3 label; 
      format value 8.2; 
run;

/***--- CLUSTERED EVENTS WITH CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.data_kid;
proc nlmixed data=four qpoints=170 noad;
	parms beta0  = 7.7 beta1  = 0.98189 beta2  = 0.98189 beta3  = 0.98189
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1     gamma3 = 0.1  
		  theta=4.7 lambda = 0.3 
		  r01-r10 1e-4;

	bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 theta  >= 0;


	base_haz = r01*fail_r1 + r02*fail_r2 + r03*fail_r3 + r04*fail_r4 + r05*fail_r5
		+ r06*fail_r6 + r07*fail_r7 + r08*fail_r8 + r09*fail_r9 + r10*fail_r10; 
	cum_base_haz = r01*dur_r1 + r02*dur_r2 + r03*dur_r3 + r04*dur_r4 + r05*dur_r5
		+ r06*dur_r6 + r07*dur_r7 + r08*dur_r8 + r09*dur_r9 + r10*dur_r10;

	* cure *;
	check = gamma0 + (gamma1 * (hcvdr_num=3)) +  (gamma2 * (hcvdr_num=2)) + (gamma3 * (hcvdr_num=1)) + lambda*nu;
	expcheck = exp(-check);
    prob= 1/(1+expcheck);
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4) + nu;
	ll = -cum_base_haz*exp(mu);

	if pstatus=0 then loglik = log(prob + (1-prob) * exp(ll)); 
	if pstatus=1 then loglik = log(1-prob) + log(base_haz) + ll + mu;
	
	model ptime~general(loglik*cbps_weight);
	random nu~normal(0,theta) subject=ctr_code;
	
	ods exclude iterhistory parameters;
run;
proc print data=work.fit4 label; 
      format value 8.2; 
run;

/****************************************************************************************/
/*                                 Kidney graft Survival                                */
/****************************************************************************************/

proc phreg data= work.data_kid;
class hcvdr_num (ref="0")/param=ref;
model gtime_ki*gstatus_ki(0) = hcvdr_num;
run;
/***--- NON-CLUSTERED EVENTS WITHOUT CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.fit5;
proc nlmixed data=work.data_kid qpoints=170 noad;
	parms beta0 = 7.7 beta1 = 0.35952 beta2 = 0.70572 beta3 = 0.545555
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1  gamma3 = 0.1  
		  theta=4.7 
		  g_r01-g_r10 1e-4;

	bounds g_r01 g_r02 g_r03 g_r04 g_r05 g_r06 g_r07 g_r08 g_r09 g_r10 theta  >= 0;

	base_haz = g_r01*g_fail_r1 + g_r02*g_fail_r2 + g_r03*g_fail_r3 + g_r04*g_fail_r4 + g_r05*g_fail_r5
		+ g_r06*g_fail_r6 + g_r07*g_fail_r7 + g_r08*g_fail_r8 + g_r09*g_fail_r9 + g_r10*g_fail_r10; 
	cum_base_haz = g_r01*g_dur_r1 + g_r02*g_dur_r2 + g_r03*g_dur_r3 + g_r04*g_dur_r4 + g_r05*g_dur_r5
		+ g_r06*g_dur_r6 + g_r07*g_dur_r7 + g_r08*g_dur_r8 + g_r09*g_dur_r9 + g_r10*g_dur_r10;
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4) ;
	ll = -cum_base_haz*exp(mu);

	if gstatus_ki=0 then loglik = ll;
	if gstatus_ki=1 then loglik = ll + log(base_haz) + mu;

	model gtime_ki~general(loglik*cbps_weight);

	ods exclude iterhistory parameters;
run;
proc print data=work.fit5 label; 
      format value 8.2; 
run;

/***--- NON-CLUSTERED EVENTS WITH CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.fit6;
proc nlmixed data=work.data_kid qpoints=170 noad;
	parms beta0 = 7.7 beta1 = 0.35952 beta2 = 0.70572 beta3 = 0.545555
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1     gamma3 = 0.1  
		  theta=4.7 
		  g_r01-g_r10 1e-4;

	bounds g_r01 g_r02 g_r03 g_r04 g_r05 g_r06 g_r07 g_r08 g_r09 g_r10 theta  >= 0;


	base_haz = g_r01*g_fail_r1 + g_r02*g_fail_r2 + g_r03*g_fail_r3 + g_r04*g_fail_r4 + g_r05*g_fail_r5
		+ g_r06*g_fail_r6 + g_r07*g_fail_r7 + g_r08*g_fail_r8 + g_r09*g_fail_r9 + g_r10*g_fail_r10; 
	cum_base_haz = g_r01*g_dur_r1 + g_r02*g_dur_r2 + g_r03*g_dur_r3 + g_r04*g_dur_r4 + g_r05*g_dur_r5
		+ g_r06*g_dur_r6 + g_r07*g_dur_r7 + g_r08*g_dur_r8 + g_r09*g_dur_r9 + g_r10*g_dur_r10;

	* cure *;
	check = gamma0 + (gamma1 * (hcvdr_num=3)) +  (gamma2 * (hcvdr_num=2)) + (gamma3 * (hcvdr_num=1));
	expcheck = exp(-check);
    prob= 1/(1+expcheck);
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4);
	ll = -cum_base_haz*exp(mu);

	if gstatus_ki=0 then loglik = log(prob + (1-prob) * exp(ll)); 
	if gstatus_ki=1 then loglik = log(1-prob) + log(base_haz) + ll + mu;
	
	model gtime_ki~general(loglik*cbps_weight);

	ods exclude iterhistory parameters;
run;
proc print data=work.fit6 label; 
      format value 8.2; 
run;

/***--- CLUSTERED EVENTS WITHOUT CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.fit7;
proc nlmixed data=work.data_kid qpoints=170 noad;
	parms beta0 = 7.7 beta1 = 0.35952 beta2 = 0.70572 beta3 = 0.545555
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1  gamma3 = 0.1  
		  theta=4.7 
		  g_r01-g_r10 1e-4;

	bounds g_r01 g_r02 g_r03 g_r04 g_r05 g_r06 g_r07 g_r08 g_r09 g_r10 theta  >= 0;


	base_haz = g_r01*g_fail_r1 + g_r02*g_fail_r2 + g_r03*g_fail_r3 + g_r04*g_fail_r4 + g_r05*g_fail_r5
		+ g_r06*g_fail_r6 + g_r07*g_fail_r7 + g_r08*g_fail_r8 + g_r09*g_fail_r9 + g_r10*g_fail_r10; 
	cum_base_haz = g_r01*g_dur_r1 + g_r02*g_dur_r2 + g_r03*g_dur_r3 + g_r04*g_dur_r4 + g_r05*g_dur_r5
		+ g_r06*g_dur_r6 + g_r07*g_dur_r7 + g_r08*g_dur_r8 + g_r09*g_dur_r9 + g_r10*g_dur_r10;
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4) + nu;
	ll = -cum_base_haz*exp(mu);

	if gstatus_ki=0 then loglik = ll;
	if gstatus_ki=1 then loglik = ll + log(base_haz) + mu;
	model gtime_ki~general(loglik*cbps_weight);
	random nu~normal(0,theta) subject=ctr_code;
	
	ods exclude iterhistory parameters;
run;
proc print data=work.fit7 label; 
      format value 8.2; 
run;

/***--- CLUSTERED EVENTS WITH CURE-RATE PROPORTION ---***/
ods output FitStatistics=work.fit8;
proc nlmixed data=work.data_kid qpoints=190 noad;
	parms beta0 = 7.7 beta1 = 0.35952 beta2 = 0.70572 beta3 = 0.545555
		  gamma0 = 0.1 gamma1 = 0.1 gamma2 = 0.1     gamma3 = 0.1  
		  theta=0.44 lambda = 0.3 
		  g_r01-g_r10 1e-4;

	bounds g_r01 g_r02 g_r03 g_r04 g_r05 g_r06 g_r07 g_r08 g_r09 g_r10 theta  >= 0;


	base_haz = g_r01*g_fail_r1 + g_r02*g_fail_r2 + g_r03*g_fail_r3 + g_r04*g_fail_r4 + g_r05*g_fail_r5
		+ g_r06*g_fail_r6 + g_r07*g_fail_r7 + g_r08*g_fail_r8 + g_r09*g_fail_r9 + g_r10*g_fail_r10; 
	cum_base_haz = g_r01*g_dur_r1 + g_r02*g_dur_r2 + g_r03*g_dur_r3 + g_r04*g_dur_r4 + g_r05*g_dur_r5
		+ g_r06*g_dur_r6 + g_r07*g_dur_r7 + g_r08*g_dur_r8 + g_r09*g_dur_r9 + g_r10*g_dur_r10;

	* cure *;
	check = gamma0 + (gamma1 * (hcvdr_num=3)) +  (gamma2 * (hcvdr_num=2)) + (gamma3 * (hcvdr_num=1)) + lambda*nu;
	expcheck = exp(-check);
    prob= 1/(1+expcheck);
	
	* log-survival prob *;
	mu = beta0*(hcvdr_num=1) + beta1*(hcvdr_num=2)+ beta2*(hcvdr_num=3) + beta3*(hcvdr_num=4) + nu;
	ll = -cum_base_haz*exp(mu);

	if gstatus_ki=0 then loglik = log(prob + (1-prob) * exp(ll)); 
	if gstatus_ki=1 then loglik = log(1-prob) + log(base_haz) + ll + mu;
	
	model gtime_ki~general(loglik*cbps_weight);
	random nu~normal(0,theta) subject=ctr_code;
	
	ods exclude iterhistory parameters;
run;
proc print data=work.fit8 label; 
      format value 8.2; 
run;

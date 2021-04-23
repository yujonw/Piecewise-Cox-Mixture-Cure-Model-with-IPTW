**%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%**;
**                                                                  **;
**                       TWANG for SAS Macro                        **;
**                                                                  **;
**         Tools for Causal Modeling wiht Propensity Scores         **;
**                                                                  **
**                 Version 3.1.2 December 9, 2016                   **;
**                                                                  **;
** This package of macros provide tools for causal modeling using   **;
** propensity scores. They port the functions of the twang package  **;
** in R to SAS via calls to R that are invisible to the users.      **;
** TWANG using Generalized Boosting Models to estimate the pro-     **;
** pensity scores.  The package also provides tools for evaluating  **;
** the utility of the fitted proponsity score by testing covariate  **;
** balance using weighted samples.                                  **;
**                                                                  **;
** To use TWANG for SAS, users must install R and have an active    **;
** SAS license.  Users should check with their local IT support for **;
** information about using SAS. R is a free software environment    **:
** for statistical computing and graphicsThe To install R, point a  **;
** web browers to http://cran.us.r-project.org/ and click on the    **;
** download link that matches your operating platform. Identity     **;
** the folder where R is installed on your machine.  For Windows    **;
** users the defalt location for R is something like:               **;
** C:/Program Files/R/R-3.0.1/bin/x64/R.exe, where 3.0.1 will be    **;
** replaced by the current release version number of R at the time  **;
** it is installed.                                                 **;
**                                                                  **;
** End user macros are:                                             **;
** ps -- estimates propensity scores for two treatments using GBM   **;
** plot -- generate default diagnostics plots for two treatments    **;
** dxwts -- evaluates quality of resulting weights                  **;
** mnps -- fits propensity scores to 3+ tx conditions using GBM     **;
** mnbaltable -- check covariate balance for 3+ treatments          **;
** mnplot -- generate default diagnostics plots for 3+ treatments   **;
** CBPS -- estimates propensity scores using CBPS methods           **;
** update_twang -- Updates the TWANG package in R                   **;
** remove_twang_folder -- removes the twang folder created by the   **;
**                        macros                                    **;
**                                                                  **;
** Help on using the macros can be found in twang_help.txt          **;
**%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%**;

**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**                                                                  **;
** Update Log                                                       **;
** 3.1.1.2                                                          **;
** 12-09-2016: added xsync to the options in all the macros         **;
**             it is the default but it was not on for a user and   **;
**             it created many problems                             **;
** 3.1.1.1                                                          **;
** 09-07-2016: added guessingrow to proc import statements          **;
**                                                                  **;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;

**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**                                                                  **;
** Utility macros                                                   **;
**                                                                  **;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;

* error  handler for untrapped errors from R;
%macro printerr(filename);
* temporarily turn off logging for this macro;
options nosource nonotes;
run;

data raw;
length line $1000;
infile "&filename" lrecl=20000;
input line & $;

data msg;
set raw;
retain flag;
if substr(line,1,5)='Error' then flag=1;

if flag;
run;

data _null_;
set msg;
put line;
run;

options source notes;
run;
%mend;

* macro to trap any warning messages from R output;
%macro printwarn(filename);
* temporarily turn off logging for this macro;
options nosource nonotes;
run;

data raw;
length line $1000;
infile "&filename" lrecl=20000;
input line & $;

data msg;
set raw;
retain flag;
if substr(line,1,7)='Warning' then flag=1;
if flag=1 and substr(line,1,1)='>' then flag=0;
if flag=1 and substr(line,1,25)='Loading required package:' then flag=0;
if flag=1 and substr(line,1,18)='Attaching package:' then flag=0;
if flag;
run;

* print warning message if one exists;
proc sql noprint;
 select count(*) into :nobs from msg;
quit;
%if &nobs>0 %then %do;
	%put WARNING: The following warning message was returned from R:;
    data _null_;
	set msg;
	put line;
	run;
	%put;
%end;
options source notes;
run;
%mend;

* utility macro to take a SAS list and produce a c() list in R notation;
%macro combine(mylist);
  %let k=1;
  %let newlist=;
  %let var = %scan(&mylist, &k, ' ');
     %do %while(("&var" NE ""));
      %if &k=1 %then %let newlist = c("&var";
	  %else %let newlist =&newlist, "&var";
	  %let k = %eval(&k + 1);
      %let var = %scan(&mylist, &k, ' ');
  %end;
  %let newlist = &newlist);
   &newlist
%mend;

* utility macro to take a SAS list of numbers and produce a c() list in R notation;
%macro ncombine(mylist);
  %let k=1;
  %let newlist=;
  %let var = %scan(&mylist, &k, ' ');
     %do %while(("&var" NE ""));
      %if &k=1 %then %let newlist = c(&var;
	  %else %let newlist =&newlist, &var;
	  %let k = %eval(&k + 1);
      %let var = %scan(&mylist, &k, ' ');
  %end;
  %let newlist = &newlist);
   &newlist
%mend;

* utility macro to create a variable list from a formula;
%macro getlist(formula);
  %let ff = %substr(&formula, %index(&formula,~)+1);
  %let varlist = %scan(&ff, 1, '+');
  %let k = 2;
  %let var = %scan(&ff, &k, '+');
  %do %while(("&var" NE ""));
      %let varlist = &varlist &var;
      %let k = %eval(&k + 1);
      %let var = %scan(&ff, &k, '+');
      %end;
  &varlist
%mend;

* macro for checking parameters values are allowable ;
%macro checkmult(paramname,params,allowed);
  %let j=1;
  %let pvar = %scan(&params, &j, ' ');
  %do %while("&pvar" NE "");
    %let k=1;
	%let found=0;
    %let avar = %scan(&allowed, &k, ' ');
    %do %while("&avar" NE "");
      %if &avar eq &pvar %then %do;
		%let found=1;
	  %end;
      %let k = %eval(&k + 1);
      %let avar = %scan(&allowed, &k, ' ');
    %end;
	%if &found=0 %then %do;
      %put ERROR: The term "&pvar" is not a valid value for the &paramname parameter. Please use one of the following: "&allowed".;
      %abort;
    %end;
    %let j = %eval(&j + 1);
    %let pvar = %scan(&params, &j, ' ');
  %end;
%mend;

* macro for checking parameters values are allowable ;
%macro checksingle(paramname,param,allowed);
  %let k=1;
  %let found=0;
  %let var = %scan(&allowed, &k, ' ');
  %do %while(("&var" NE ""));
    %if &var eq &param %then %let found=1;
    %let k = %eval(&k + 1);
    %let var = %scan(&allowed, &k, ' ');
  %end;
	%if &found=0 %then %do;
      %put ERROR: The term "&param" is not a valid value for the &paramname parameter. Please use one of the following: "&allowed".;
      %abort;
  %end;
%mend;

* macro for checking parameters values are allowable -- postive numeric variable ;
%macro checkpnum(paramname,param);
    %let chkval = %eval(%sysfunc(verify(%sysfunc(trim(%sysfunc(left(&param)))),'0123456789.')) +
                        %index(%substr(&param,%index(&param,.)+1), .));	
    %if &chkval > 0 %then %do;
      %put ERROR: The value for the &paramname parameter must be a positive number not in scientific notation. You specified: &param. Please change this to valid value.;
      %abort;
  %end;
%mend;

* macro for checking parameters values are allowable -- postive integer ;
%macro checkpint(paramname,param);
    %let chkval = %eval(%sysfunc(verify(%sysfunc(trim(%sysfunc(left(&param)))),'0123456789')) +
                        %index(%substr(&param,%index(&param,.)+1), .));	
    %if &param = 0 %then %let chkval = 1;
    %if &chkval > 0 %then %do;
      %put ERROR: The value for the &paramname parameter must be a positive integer. You specified: &param. Please change this to valid value.;
      %abort;
  %end;
%mend;

* macro to test for existence of variables in a dataset;
%macro checkexist(varlist,dsn);
  %let k=1;
  %let var = %scan(&varlist, &k, ' ');
  %do %while(("&var" NE ""));
		data _null_;
  	dsid=open("&dsn");
		check=varnum(dsid,"&var");
  	call symput("check",left(put(check, best12.)));
  	run;
		%if &check eq 0 %then %do;
      %put ERROR: The variable &var does not exist in the data set &dsn.;
      %abort;
  	%end;
    %let k = %eval(&k + 1);
    %let var = %scan(&varlist, &k, ' ');
  %end;
%mend;

* check if the weight variables are in a dataset;
%macro checkwgt(stopmethod, estimand, dsn);
  data _null_;
     dsid=open("&dsn");
     call symput("dsid", left(put(dsid, best5.)));
  run;
  %if &dsid > 0 %then %do; 

      %let ii = 1; 
      %let sm1 = %scan(&stopmethod, &ii, ' ');
      %do %while(&sm1 NE %str());
           %let sm1 = %sysfunc(translate(&sm1, '_', '.'))_&estimand;
           data _null_;
              dsid=open("&dsn");
              check = varnum(dsid, "&sm1");
              call symput("check",left(put(check, best12.)));
           run;
           %if &check NE 0 %then %do;
               %put ERROR: The weight variable &sm1 is in the specified output data set &dsn.. To avoid overwriting the variable values, change the output data set, or rename, or delete the variable;
               %abort;
               %end; 
           %let ii = %eval(&ii + 1);
           %let sm1 = %scan(&stopmethod, &ii, ' ');
           %end;   
   %end;  ** ends if dsid > 0;

%mend;

%macro checkps(stopmethod, estimand, dsn);
  data _null_;
     dsid=open("&dsn");
     call symput("dsid", left(put(dsid, best5.)));
  run;
  %if &dsid > 0 %then %do; 

      %let ii = 1; 
      %let sm1 = %scan(&stopmethod, &ii, ' ');
      %do %while(&sm1 NE %str());
           %let sm1 = ps_%sysfunc(translate(&sm1, '_', '.'))_&estimand;
           data _null_;
              dsid=open("&dsn");
              check = varnum(dsid, "&sm1");
              call symput("check",left(put(check, best12.)));
           run;
           %if &check NE 0 %then %do;
               %put ERROR: The propensity score variable &sm1 is in the specified output data set &dsn.. To avoid overwriting the variable values, change the output data set, or rename, or delete the variable;
               %abort;
               %end; 
           %let ii = %eval(&ii + 1);
           %let sm1 = %scan(&stopmethod, &ii, ' ');
           %end;   
   %end;  ** ends if dsid > 0;

%mend;

%macro checkps_mnps(ntxlevs, stopmethod, estimand, dsn);
  data _null_;
     dsid=open("&dsn");
     call symput("dsid", left(put(dsid, best5.)));
  run;
  %if &dsid > 0 %then %do; 

      %let ii = 1; 
      %let sm1 = %scan(&stopmethod, &ii, ' ');
      %do %while(&sm1 NE %str());
           %do jj = 1 %to &ntxlevs;
              %let sm2 = ps_&jj._%sysfunc(translate(&sm1, '_', '.'))_&estimand;
              data _null_;
                 dsid=open("&dsn");
                 check = varnum(dsid, "&sm2");
                 call symput("check",left(put(check, best12.)));
              run;
              %if &check NE 0 %then %do;
                  %put ERROR: The propensity score variable &sm2 is in the specified output data set &dsn.. To avoid overwriting the variable values, change the output data set, or rename, or delete the variable;
                  %abort;
                  %end; 
              %end; ** ends do jj;
           %let ii = %eval(&ii + 1);
           %let sm1 = %scan(&stopmethod, &ii, ' ');
           %end;   
   %end;  ** ends if dsid > 0;

%mend;

* get the directory path part of a fully qualified file name;
%macro getpath(file);
%if %index(&file, /) = 0 %then %let path = %str();
%else %let path=%substr(&file,1,%eval(%index(&file,%scan(&file,-1,'/'))-1));
&path
%mend;

* check if a pathname is valid;
%macro checkdir(dir);

%if %sysfunc(fileexist(&dir)) eq 0  %then %do;
  %put ERROR: The path "&dir" does not exist.;
  %abort;
%end;

%mend;

* utility macro to add "+" between vector of variable names;
%macro addplus(vars);
%let vv = %scan(&vars, 1);
%let form = &vv;
%let i = 2;
%let vv = %scan(&vars, 2);
%do %while(&vv ^= %str());
   %let form = &form + &vv;
   %let i = %eval(&i + 1);
   %let vv = %scan(&vars, &i);
   %end;

   &form

%mend;
 
* utility macro to check that parameter is not blank **;
%macro notblank(var, vname);
%if (%sysevalf(%superq(var)=,boolean)) %then %do; 
     %put ERROR: &vname is missing a value;
     %abort;
     %end;
%mend;

** Utility macro to find the terms in list1 not in list2 **;
%macro setdiff(list1, list2);
%let i = 1;
%let v1 = %scan(&list1, &i);
%let diff = %str();
%do %while(&v1 ^= %str());
  %let loc = %index(&list2, &v1);
  %if &loc = 0 %then %let diff = &diff &v1;
  %let i = %eval(&i + 1);
  %let v1 = %scan(&list1, &i);
  %end;
&diff
%mend setdiff;

** Utility macro to get current path from either SAS log location or user profile **;
%macro getcd;
data _null_;
   call symput("cd", getoption("LOG"));
run;
%if "&cd" ne "" %then %let cd = %sysfunc(translate(&cd,/,\));

%if "&cd" eq "" %then %do;
  %let profile=%sysget(USERPROFILE);
  %let cd =  %sysfunc(translate(&profile,/,\))/;
%end;
%else %let cd = %getpath(&cd);

%mend getcd;

** Utility to check if Rcmd exists or is an executable in the path **;
%macro chkRcmd(Rcmd);

 ** Check that the executable is R.exe or R not Rscript or Rcmd **;

%let rfile = %substr(&Rcmd,%eval(%index(&Rcmd,%scan(&Rcmd,-1,'\'))));
%if (%upcase(&rfile) ^= R and %upcase(&rfile) ^= R.EXE) %then %do;
    %put %str(ERROR: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++);
    %put ERROR: Please use R or R.exe for the executable file in the Rcmd argument. ;
    %put ERROR: The macros do not support Rscript or Rcmd;
    %put %str(ERROR: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++);
    %abort;
    %end;

%if %sysfunc(fileexist(&Rcmd)) eq 0  %then %do;
  %let wdir = %sysfunc(getoption(work));
  %let notfound = TRUE;
  data _null_;
     file "&wdir\wherer.bat";
     put "%str(where %1 > %"&wdir\rchk.txt%")";
  run;
  x "%str(%"&wdir\wherer%" &Rcmd)";
  data _null_;
     infile "&wdir\rchk.txt";
     input chk $ 1 - 10;
     if chk ^= "" then call symput("notfound", "FALSE"); 
  run;
  %if &notfound = TRUE %then %do;
      %put ERROR: Could not find R executable at "&Rcmd".;
      %abort;
      %end;
%end;

%mend chkRcmd;

**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**                                                                  **;
** End User Macros                                                  **;
**                                                                  **;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;


**##################################################################**;
** Macro: ps
** Fits propensity scores using GBM;

%macro ps(treatvar=,
          vars=,
          class=,
          dataset=, 
          ntrees=10000, 
          intdepth=3, 
          shrinkage=0.01, 
          permtestiters=0, 
          stopmethod=ks.mean es.mean, 
          sampw=, 
          estimand=ATE, 
          output_dataset=_inputds,
          Rcmd=, 
          plotname=, 
          objpath=);


options xsync noxwait;

* get cwd;
%getcd;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&treatvar, treatvar)
%notblank(&vars, vars)
%notblank(&dataset, dataset)
%notblank(&stopmethod, stopmethod)
%notblank(&estimand, estimand)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let treatvar =%qsysfunc(compress(&treatvar,%str(%"%')));
%let vars =%qsysfunc(compress(&vars,%str(%"%')));
%let dataset  =%qsysfunc(compress(&dataset,%str(%"%')));
%let ntrees =%qsysfunc(compress(&ntrees,%str(%"%')));
%let intdepth =%qsysfunc(compress(&intdepth,%str(%"%')));
%let shrinkage =%qsysfunc(compress(&shrinkage,%str(%"%')));
%let permtestiters =%qsysfunc(compress(&permtestiters,%str(%"%')));
%let stopmethod =%lowcase(%qsysfunc(compress(&stopmethod,%str(%"%'))));
%let estimand =%upcase(%qsysfunc(compress(&estimand,%str(%"%'))));
%let sampw =%qsysfunc(compress(&sampw,%str(%"%')));
%let plotname =%qsysfunc(compress(&plotname,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate enumerated parameters;
%checkmult(stopmethod,&stopmethod,ks.mean es.mean ks.max es.max ks.max.direct es.max.direct);
%checksingle(estimand,&estimand,ATT ATE);

* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));
%if not(%sysevalf(%superq(plotname)=,boolean)) %then %let plotname = %sysfunc(translate(&plotname,/,\));

* check that the weight variables are not in the specified output_dataset *;
%checkwgt(&stopmethod, &estimand, &output_dataset);

* check that the propensity score variables are not in the specified output_dataset *;
%checkps(&stopmethod, &estimand, &output_dataset);

* check that the covariates are in the dataset;

%checkexist(&vars,&dataset);


* check to make sure treatment variable does not contain missings;
   ** check it is in the dataset **;
%checkexist(&treatvar,&dataset);

   ** check the values **;
proc sql noprint;
select count(*) into :nmiss from &dataset where &treatvar not in (0,1);
%if &nmiss ne 0 %then %do;
  %put ERROR: Treatment variable &treatvar has invalide values (%trim(&nmiss) cases).;
  %put        Treatment variable must take values 0 or 1 only.;
  %abort;
%end;

   ** check the class variables **;
%let chkclass = %sysevalf(&class=, boolean);
%if &chkclass = 0 %then %do;
  %let class =%qsysfunc(compress(&class,%str(%"%')));
  %let notin = %setdiff(&class, &vars);
  %if  &notin ^= %str() %then %do;
     %put ERROR: CLASS variables &notin are not in vars parameter list; 
     %abort;
     %end;
  %end;
   
%if &sampw ^= %str() %then %do;
   ** check it is in the dataset **;
   %checkexist(&sampw,&dataset);

   ** check the values **;
   * check to make sure sampw variable does not contain missing or negative values;
   proc sql noprint;
   select count(*) into :nmiss from &dataset where &sampw < 0;
   %if &nmiss ne 0 %then %do;
     %put ERROR: SAMPW variable &sampw has invalide values (%trim(&nmiss) cases).;
     %put        SAMPW variable must take nonnegative values only.;
     %abort;
   %end;
%end;

** Set all variables used in R script to lower case **;
%let treatvar = %lowcase(&treatvar);
%let vars = %lowcase(&vars);
%if &chkclass = 0 %then %let class = %lowcase(&class);
%if &sampw ^= %str() %then %let sampw = %lowcase(&sampw);

%let formula = &treatvar ~ %addplus(&vars);

filename rdata "&workdir/datafile.csv";
filename rscript "&workdir/ps.R";


* append an internal ID to supplied analysis data set;
data _inputds;
	length tempID 8.;
	set &dataset;
	tempID=_N_;

proc sort;
	by tempID;
run;

* Remove existing versions of files to be created *;
%sysexec(%str(del %"&workdir\datafile.csv%")); 
%sysexec(%str(del %"&workdir\summary.csv%")); 
%sysexec(%str(del %"&workdir\baltab.csv%")); 
%sysexec(%str(del %"&workdir\wts.csv%")); 
%sysexec(%str(del %"&workdir\ps.R%")); 
%sysexec(%str(del %"&workdir\ps.Rout%")); 
%sysexec(%str(del %"&workdir\ps.Rdata%")); 

%if not(%sysevalf(%superq(plotname)=,boolean)) %then %do;
        %let plotpath = %getpath(&plotname);
	* check if path is valid, if one is given;
	%if not(%sysevalf(%superq(plotpath)=, boolean)) %then %checkdir(&plotpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let plotname=%sysfunc(translate(&objpath,/,\))/&plotname;
        %else %let plotname=&cd&plotname;

        %sysexec(%str(del %"&plotname%"));
%end;
	

* export dataset for R;
PROC EXPORT DATA=_inputds OUTFILE=rdata DBMS=csv replace;
run;


* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
*        put "%str(msg <- file(%"&workdir/ps_msg.Rout%", open=%"wt%"))";  ** this was unnecessary so it was not implemented **;
*        put "%str(sink(msg, type=%"message%"))";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
	put "library(twang)";
        put ;
        put ;
        put ;
	put "set.seed(1)";
	put ;
	put "%str(inputds<-read.csv(%"&workdir/datafile.csv%"))";
  put ;
	put "vnames <- names(inputds)          ";
  put "names(inputds) <- tolower(names(inputds))";
  put ;
%if &chkclass = 0 %then
        put %unquote(%str(%'inputds[, %combine(&class)] <- lapply(inputds[,%combine(&class), drop=F], as.factor)%')); ;
        put ;
	put "ps1 <- ps(&formula,";
	put " data = inputds,";
	put " n.trees = &ntrees,";
	put " interaction.depth = &intdepth,";
	put " shrinkage = &shrinkage,";
	put " perm.test.iters = &permtestiters,";
	put " " %unquote(%str(%'stop.method = %combine(&stopmethod)%')) ",";
	put " " %unquote(%str(%'estimand = %"&estimand%"%')) ",";
/*        put " sampw = "; %if &sampw=%str() %then put "NULL"; %else put %"&sampw%";; */
        put " sampw = " %if &sampw=%str() %then "NULL,"; %else "inputds$&sampw,"  ;; 
        put " verbose = FALSE";
	put " )";
	put ;
	put "baltab<-bal.table(ps1)";
	put ;
	put "bnames <- rownames(baltab$unw)                                         ";
	put "bnames1 <- sapply(strsplit(bnames, ':'), function(x){return(x[[1]])})  ";
	put "bnames1 <- vnames[match(bnames1, tolower(vnames))]                     ";
	put "substr(bnames, 1, nchar(bnames1)) <- bnames1                           ";
	put "baltab <- lapply(baltab, function(u){                                  ";
       	put "                 rownames(u) <- bnames                                 ";
      	put "                 return(u)})                                           ";
	put ;
	put "w<-as.data.frame(ps1$w)";
	put "psests<-as.data.frame(ps1$ps)";
	put "%str(names(psests)<-paste(%"ps%",names(psests),sep=%".%"))";
        put "w<-data.frame(w, psests)";
	put "w$tempID<-as.numeric(row.names(w))";
	put "%str(write.table(w,file=%"&workdir/wts.csv%",row.names=FALSE,col.names=TRUE,sep=%',%'))";
	put "baltab <- data.frame(do.call(rbind, baltab), table.name=rep(names(baltab), each=nrow(baltab[[1]])))";
	put "baltab <- data.frame(row_name=row.names(baltab), baltab)";
        put "baltab[baltab==Inf] <- NA";
        put "baltab[baltab==(-Inf)] <- NA";
	put "%str(write.table(baltab,file=%"&workdir/baltab.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
	put "summ<-as.data.frame(rbind(summary(ps1)))";
	put "summ <- data.frame(row_name=row.names(summ), summ)";
	put "%str(write.table(summ,file=%"&workdir/summary.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";

* if plotname parameter is defined, run a plot as well;
%if not(%sysevalf(%superq(plotname)=,boolean)) %then %do;
/*  ** moved the following code above to delete the old versions of the files  **;
        %let plotpath = %getpath(&plotname);
	* check if path is valid, if one is given;
	%if not(%sysevalf(%superq(plotpath)=, boolean)) %then %checkdir(&plotpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let plotname=%sysfunc(translate(&objpath,/,\))/&plotname;
        %else %let plotname=&cd&plotname;
*/
	%put Writing plot to &plotname.;
	put "%str(pdf(%'&plotname%'))";
	put "plot(ps1,plots=1, main='Plot 1 (optimize): GBM Optimization')";
	put "plot(ps1,plots=2, main='Plot 2 (boxplot): Boxplot of Propensity Scores')";
	put "plot(ps1,plots=3, main='Plot 3 (es): Standardized Effect Sizes Pre/Post Weighting')";
	put "plot(ps1,plots=4, main='Plot 4 (t): T-test P-values of Group Means of Covariates')";
	put "plot(ps1,plots=5, main='Plot 5 (ks): K-S P-values of Group Distns of Covariates')";
	*put "plot(ps1,plots=6, main='Plot 6 (histogram): Histogram of Weights by Group')";
	put "dev.off()";
%end;	

* if objpath parameter is defined, save the workspace image;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
	put "%str(save(ps1, file=%'&workdir/ps.RData%'))";
%end;	
	
run;

* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/ps.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/ps.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/ps.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/ps.Rout);
run;

* input weights generated by R;
proc import datafile="&workdir/wts.csv"
     out=_weights
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

proc sort data=_weights;
	by tempID;	

* merge weights back on to analysis data set;
data &output_dataset;
	merge _inputds(in=a) _weights(in=b);
	by tempID;
	drop tempID;
run;

* input summary table generated by R;
proc import datafile="&workdir/summary.csv"
     out=_summ
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

proc print data=_summ;
title "Summary table";
run;

* input balance tables generated by R;
proc import datafile="&workdir/baltab.csv"
     out=_baltab
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

* print balance tables by stop method;
proc sql noprint;
select distinct table_name
into :elements separated by " "
from _baltab
order by table_name desc;

  %let j=1;
  %let key=%scan(&elements.,&j.,%str( ));
  %do %until(&key eq %nrstr( ));
    
    proc print data=_baltab;
		where table_name="&key";
		title "Balance table: &key";
		run;

    %let j=%eval(&j+1);
    %let key=%scan(&elements.,&j.,%str( ));
  %end;

title;
run;

%mend ps;
**##################################################################**;


**##################################################################**;
** Macro: plot
** Fits generates diagnostic plots per user specifications;

%macro plot(inputobj=, 
            plotname=, 
            plotformat=, 
            plots=, 
            subset=, 
            color=TRUE,
            pairwisemax=TRUE,
            figurerows=1,
            Rcmd=, 
            objpath= );

options xsync noxwait;

* get cwd;
%getcd;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&inputobj, inputobj)
%notblank(&plotname, plotname)
%notblank(&plots, plots)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let inputobj =%qsysfunc(compress(&inputobj,%str(%"%')));
%let plotname  =%qsysfunc(compress(&plotname,%str(%"%')));
%let plotformat =%qsysfunc(compress(&plotformat,%str(%"%')));
%let plots =%qsysfunc(compress(&plots,%str(%"%')));
%let color =%qsysfunc(compress(&color,%str(%"%')));
%let subset =%qsysfunc(compress(&subset,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate enumerated parameters;
%checksingle(plots,&plots,1 2 3 4 5 6 optimize boxplot es t ks histogram);
%checksingle(color, &color, T F TRUE FALSE);
%checksingle(pairwisemax, &pairwisemax, T F TRUE FALSE);
%checksingle(figurerows, &figurerows, 1 2 3 4);
%if &plots=optimize or &plots=boxplot or &plots=es or &plots=t or &plots=ks or &plots=histogram %then %do;
  %let plots=%str(%"&plots%");
%end;

* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));
%if not(%sysevalf(%superq(plotname)=,boolean)) %then %let plotname = %sysfunc(translate(&plotname,/,\));
%let inputobj = %sysfunc(translate(&inputobj,/,\));
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let objpath = %sysfunc(translate(&objpath,/,\));



* validate file names;
 ** Check the path for inputobj name if it contains a path **;
%let inputpath = %getpath(&inputobj);
* check if path is valid, if one is given;
%if not(%sysevalf(%superq(inputpath)=, boolean)) %then %checkdir(&inputpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let inputobj=&objpath/&inputobj;
  %else %let inputobj=&cd&inputobj;

* check to make sure inputobj exists;
%if %sysfunc(fileexist(&inputobj)) eq 0  %then %do;
	%put ERROR: The input object "&inputobj" could not be found. Please check your path names or specify an object path.;
  %abort;
%end;



 ** Check the path for plot name if it contains a path **;
%let plotpath = %getpath(&plotname);
* check if path is valid, if one is given;
%if not(%sysevalf(%superq(plotpath)=, boolean)) %then %checkdir(&plotpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let plotname=&objpath/&plotname;
  %else %let plotname=&cd&plotname;
	%put Writing plot to &plotname.;


filename rscript "&workdir/plot.R";

* determine plot format if it exists;
%if %lowcase(&plotformat) eq jpg %then %let fmt = jpeg;
%else %if %lowcase(&plotformat) eq pdf %then %let fmt = pdf;
%else %if %lowcase(&plotformat) eq png %then %let fmt = png;
%else %if %lowcase(&plotformat) eq wmf %then %let fmt = win.metafile;
%else %if %lowcase(&plotformat) eq postscript %then %let fmt = postscript;
%else                                         %let fmt = pdf;



* Remove existing versions of files to be created *;
%sysexec(%str(del %"&plotname%")); 
%sysexec(%str(del %"&workdir\plot.R%")); 
%sysexec(%str(del %"&workdir\plot.Rout%")); 

* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
	put "library(twang)";
        put ;
	put "set.seed(1)";
	put "%str(tmp <- load(%"&inputobj%"))";
	put "%str(&fmt(%'&plotname%'))";
	put "plot(get(tmp),plots=&plots, subset=&subset, color=&color, pairwiseMax=&pairwisemax, figureRows=&figurerows)";
	put "dev.off()";
run;

* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/plot.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/plot.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/plot.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/plot.Rout);

run;
%mend plot;
**##################################################################**;


**##################################################################**;
** Macro: dxwts
** Fits runs diagnostics on weights ;
%macro dxwts(treatvar=, 
             vars=, 
             class=,
             dataset=, 
             weightvars=, 
             estimand=, 
             sampw=, 
             permtestiters=0, 
	     Rcmd=, 
             objpath=);

options xsync noxwait;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&treatvar, treatvar)
%notblank(&vars, vars)
%notblank(&dataset, dataset)
%notblank(&weightvars, weightvars)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let dataset  =%qsysfunc(compress(&dataset,%str(%"%')));
%let weightvars =%qsysfunc(compress(&weightvars,%str(%"%')));
%let estimand =%upcase(%qsysfunc(compress(&estimand,%str(%"%'))));
%let vars =%qsysfunc(compress(&vars,%str(%"%')));
%let treatvar =%qsysfunc(compress(&treatvar,%str(%"%')));
%let sampw =%qsysfunc(compress(&sampw,%str(%"%')));
%let permtestiters =%qsysfunc(compress(&permtestiters,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate enumerated parameters;
%checksingle(estimand,&estimand,ATT ATE);


* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));

filename rdata "&workdir/dxwtsdatafile.csv";
filename rscript "&workdir/dxwts.R";

   ** check the class variables **;
%let chkclass = %sysevalf(&class=, boolean);
%if &chkclass = 0 %then %do;
  %let class =%qsysfunc(compress(&class,%str(%"%')));
  %let notin = %setdiff(&class, &vars);
  %if  &notin ^= %str() %then %do;
     %put ERROR: CLASS variables &notin are not in vars parameter list; 
     %abort;
     %end;
  %end;
 

* check to insure weightvars exist;
%checkexist(&weightvars,&dataset);

* check that weights make sense for given estimand;
proc means noprint nway data=&dataset;
	class &treatvar;
	var &weightvars;
	output out=wmeans;

proc sort data=wmeans;
	by &treatvar _STAT_;

proc transpose data=wmeans out=_wmeans;
	by &treatvar;
	id _STAT_;
	var &weightvars;	

%let warnflag=0;

%if &estimand eq ATE %then %do;  * make sure weights are NOT all 1 for treatment group;
data _null_;
	set _wmeans;
	if &treatvar=1 and mean=1 and std=0 then do;
		call symput('wgtvarname', _NAME_);
		call symput('warnflag', 1);
		stop;
	end;
%end;

%if &estimand eq ATT %then %do;  * make sure weights ARE all 1 for treatment group;
data _null_;
	set _wmeans;
	if &treatvar=1 and not(mean=1 and std=0) then do;
		call symput('wgtvarname', _NAME_);
		call symput('warnflag', 1);
		stop;
	end;
%end;
run;

%if &warnflag eq 1 & &estimand eq ATE %then %do;
  %put WARNING: Estimand = &estimand and weight variable &wgtvarname is identically equal to 1 for treatment group;
  %put ;
%end;

%if &warnflag eq 1 & &estimand eq ATT %then %do;
  %put WARNING: Estimand = &estimand and weight variable &wgtvarname is not identically equal to 1 for treatment group;
  %put ;
%end;

*proc print data=_wmeans;
* end weight check;

* Remove existing versions of files to be created *;
%sysexec(%str(del %"&workdir\dxwtsdatafile.csv%")); 
%sysexec(%str(del %"&workdir\dxsummary.csv%")); 
%sysexec(%str(del %"&workdir\dxwtsbaltab.csv%")); 
%sysexec(%str(del %"&workdir\dxwts.R%")); 
%sysexec(%str(del %"&workdir\dswts.Rout%")); 
	

* export dataset for R;
PROC EXPORT DATA=&dataset OUTFILE=rdata DBMS=csv replace;
run;

** Set all variables used in R script to lower case **;
%let treatvar = %lowcase(&treatvar);
%let vars = %lowcase(&vars);
%let weightvars = %lowcase(&weightvars);
%if &chkclass = 0 %then %let class = %lowcase(&class);
%if &sampw ^= %str() %then %let sampw = %lowcase(&sampw);


* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
	put "library(twang)";
        put ;
	put "set.seed(1)";
	put ;
	put "%str(inputds<-read.csv(%"&workdir/dxwtsdatafile.csv%"))";
        put ;
	put "vnames <- names(inputds)          ";
        put "names(inputds) <- tolower(names(inputds))";
	put ;
%if &chkclass = 0 %then
  put %unquote(%str(%'inputds[, %combine(&class)] <- lapply(inputds[,%combine(&class), drop=F], as.factor)%')); ;
	put ;
	put %unquote(%str(%'x<-as.matrix(subset(inputds,select=%combine(&weightvars)))%'));
	put "dxtmp <- dx.wts(x, ";
	put " data=inputds, ";
	put " " %unquote(%str(%'vars = %combine(&vars)%')) ",";
	put " " %unquote(%str(%'estimand = %"&estimand%"%')) ",";
	put " " %unquote(%str(%'treat.var = %"&treatvar%"%')) ",";
	put " x.as.weights=T, ";
	put " sampw=&sampw, ";
	put " perm.test.iters=&permtestiters";
	put ") ";
	put ;
	put "baltab<-bal.table(dxtmp)";
	put ;
	put "bnames <- rownames(baltab$unw)                                         ";
	put "bnames1 <- sapply(strsplit(bnames, ':'), function(x){return(x[[1]])})  ";
	put "bnames1 <- vnames[match(bnames1, tolower(vnames))]                     ";
	put "substr(bnames, 1, nchar(bnames1)) <- bnames1                           ";
	put "baltab <- lapply(baltab, function(u){                                  ";
        put "                 rownames(u) <- bnames                                 ";
        put "                 return(u)})                                           ";
	put ;
	put "baltab <- data.frame(do.call(rbind, baltab), table.name=rep(names(baltab), each=nrow(baltab[[1]])))";
	put "baltab <- data.frame(row_name=row.names(baltab), baltab)";
        put "baltab[baltab==Inf] <- NA";
        put "baltab[baltab==(-Inf)] <- NA";
	put "%str(write.table(baltab,file=%"&workdir/dxwtsbaltab.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
	put "%str(write.table(dxtmp$summary.tab,file=%"&workdir/dxsummary.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";



	/* ** Code from ps **;
* if objpath parameter is defined, save the workspace image;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
	put "%str(save(ps1, file=%'&workdir/ps.RData%'))";
%end;	
*/
	
run;


* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/dxwts.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/dxwts.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/dxwts.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/dxwts.Rout);


run;

* input summary table generated by R;
proc import datafile="&workdir/dxsummary.csv"
     out=_dxsumm
     dbms=csv
     replace;
     getnames=yes;
     datarow=2;
     guessingrows=MAX;
run;

** drop the iter variable since it is missing for all records **;
data _dxsumm;
   set _dxsumm;
   drop iter;
run;

proc print data=_dxsumm;
title "dxwts Summary table";
run;
	
* input balance tables generated by R;
proc import datafile="&workdir/dxwtsbaltab.csv"
     out=_baltab
     dbms=csv
     replace;
     getnames=yes;
     datarow=2;
     guessingrows=MAX;
run;

* print balance tables by stop method;
proc sql noprint;
select distinct table_name
into :elements separated by " "
from _baltab where(table_name ^= "unw")
order by table_name;

%let elements = unw &elements;

  %let j=1;
  %let key=%scan(&elements.,&j.,%str( ));
  %do %until(&key eq %nrstr( ));
    
    proc print data=_baltab;
		where table_name="&key";
		title "dxwts Balance table: &key";
		run;

    %let j=%eval(&j+1);
    %let key=%scan(&elements.,&j.,%str( ));
  %end;

title;
run;

%mend dxwts;
**##################################################################**;



**##################################################################**;
** Macro: mnps
** Fits propensity scores to 3+ tx conditions using GBM;

%macro mnps(treatvar=,
            vars=,
            class=,
            dataset=, 
            ntrees=10000, 
            intdepth=3, 
            shrinkage=0.01, 
            permtestiters=0, 
            stopmethod=ks.mean es.mean, 
            sampw=, 
            estimand=ATE, 
            treatatt = NULL,
            collapseto = pair,
            output_dataset=_inputds,
            return_ps=FALSE,
            Rcmd=, 
            plotname=, 
            objpath=);


options xsync noxwait;

* get cwd;
%getcd;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&treatvar, treatvar)
%notblank(&vars, vars)
%notblank(&dataset, dataset)
%notblank(&stopmethod, stopmethod)
%notblank(&estimand, estimand)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let treatvar =%qsysfunc(compress(&treatvar,%str(%"%')));
%let vars =%qsysfunc(compress(&vars,%str(%"%')));
%let dataset  =%qsysfunc(compress(&dataset,%str(%"%')));
%let ntrees =%qsysfunc(compress(&ntrees,%str(%"%')));
%let intdepth =%qsysfunc(compress(&intdepth,%str(%"%')));
%let shrinkage =%qsysfunc(compress(&shrinkage,%str(%"%')));
%let permtestiters =%qsysfunc(compress(&permtestiters,%str(%"%')));
%let stopmethod =%lowcase(%qsysfunc(compress(&stopmethod,%str(%"%'))));
%let sampw =%qsysfunc(compress(&sampw,%str(%"%')));
%let estimand =%upcase(%qsysfunc(compress(&estimand,%str(%"%'))));
%let treatatt =%qsysfunc(compress(&treatatt,%str(%"%')));
%let collapseto =%qsysfunc(compress(&collapseto,%str(%"%')));
%let plotname =%qsysfunc(compress(&plotname,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate enumerated parameters;
%checkmult(stopmethod,&stopmethod,ks.mean es.mean ks.max es.max);
%checksingle(estimand,&estimand,ATT ATE);
%checksingle(collapseto,&collapseto,none pair covariate stop.method);

* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));
%if not(%sysevalf(%superq(plotname)=,boolean)) %then %let plotname = %sysfunc(translate(&plotname,/,\));

* Set return_ps to TRUE if equal to t, T, true, or TRUE
%let return_ps = %upcase(&return_ps);
%if %upcase(&return_ps) = T %then %let return_ps = TRUE;

* check that the covariates are in the dataset;

%checkexist(&vars,&dataset);


* check to make sure treatment variable does not contain missings;
   ** check it is in the dataset **;
%checkexist(&treatvar,&dataset);

data _00chk;
   set &dataset(keep=&treatvar);
   mtreat = missing(&treatvar);
   drop &treatvar;
run;

   ** check the values **;
proc sql noprint;
select count(*) into :nmiss from _00chk where mtreat;
quit;
%if &nmiss ne 0 %then %do;
  %put ERROR: Treatment variable &treatvar has missing values (%trim(&nmiss) cases).;
  %abort;
%end;

   ** check that tx has 3+ unique values **;
proc sql noprint;
select count(distinct &treatvar) into :ntx_lev from &dataset;
quit;
%if &ntx_lev < 3 %then %do;
  %put ERROR: Treatment variable &treatvar must have 3+ unique values. It has %trim(&ntx_lev) values.;
  %abort;
%end;

   ** check the class variables **;
%let chkclass = %sysevalf(&class=, boolean);
%if &chkclass = 0 %then %do;
  %let class =%qsysfunc(compress(&class,%str(%"%')));
  %let notin = %setdiff(&class, &vars);
  %if  &notin ^= %str() %then %do;
     %put ERROR: CLASS variables &notin are not in vars parameter list; 
     %abort;
     %end;
  %end;
   
  ** Check the sample weight variable is specified **;
%if &sampw ^= %str() %then %do;
   ** check it is in the dataset **;
   %checkexist(&sampw,&dataset);

   ** check the values **;
   * check to make sure sampw variable does not contain missing or negative values;
   proc sql noprint;
   select count(*) into :nmiss from &dataset where &sampw < 0;
   %if &nmiss ne 0 %then %do;
     %put ERROR: SAMPW variable &sampw has invalide values (%trim(&nmiss) cases).;
     %put        SAMPW variable must take nonnegative values only.;
     %abort;
   %end;
%end;

  ** Check the treatatt var if estimand is ATT **;

     ** Determine variable type (numeric or character of treat variable **;

%let treatvar_type = %sysfunc(vartype(%sysfunc(open(&dataset,i)),%sysfunc(varnum(%sysfunc(open(&dataset,i)), &treatvar))));

%if %upcase(&estimand) = ATT %then %do;
   
     ** check treatatt is specified **;
   %if &treatatt = %str() %then %do;
     %put ERROR: TREATATT must be specified as a value of the treatment variable when estimand = "ATT".;
     %abort;
   %end;

     ** check treatatt equals a value of the treatment variable **;
   proc sql noprint;
   create table _01chk as
   select distinct &treatvar from &dataset;
   select count(*) into :nmatch from _01chk where
      %if &treatvar_type = C %then &treatvar = "&treatatt"; %else &treatvar = &treatatt; ;
   quit;
   %if &nmatch = 0 %then %do;
     %put ERROR: TREATATT must be specified as a value of the treatment variable when estimand = "ATT".;
     %abort;
   %end;
%end;

* check that the weight variables are not in the specified output_dataset *;
%checkwgt(&stopmethod, &estimand, &output_dataset);

* check that the propensity score variables are not in the specified output_data, if return_ps=TRUE *;
%if &return_ps = TRUE %then %checkps_mnps(&ntx_lev, &stopmethod, &estimand, &output_dataset);

** Set collapseto to the collapseto if not none and pair if none
%if &collapseto = none %then %let collapseto = pair;
%else %let collapseto = &collapseto;

** Set all variables used in R script to lower case **;
%let treatvar = %lowcase(&treatvar);
%let vars = %lowcase(&vars);
%if &chkclass = 0 %then %let class = %lowcase(&class);
%if &sampw ^= %str() %then %let sampw = %lowcase(&sampw);

%let formula = &treatvar ~ %addplus(&vars);

filename rdata "&workdir/datafile.csv";
filename rscript "&workdir/mnps.R";


* append an internal ID to supplied analysis data set;
data _inputds;
	length tempID 8.;
	set &dataset;
	tempID=_N_;

proc sort;
	by tempID;
run;

* Remove existing versions of files to be created *;
%sysexec(%str(del %"&workdir\datafile.csv%")); 
%if %upcase(&estimand)=ATT %then %sysexec(%str(del %"&workdir\summary.csv%")); 
%else %if %upcase(&estimand)=ATE %then %do;
    %sysexec(%str(del %"&workdir\summary1.csv%")); 
    %sysexec(%str(del %"&workdir\summary2.csv%"));
    %end; 
%sysexec(%str(del %"&workdir\baltab.csv%")); 
%sysexec(%str(del %"&workdir\wts.csv%")); 
%sysexec(%str(del %"&workdir\mnps.R%")); 
%sysexec(%str(del %"&workdir\mnps.Rout%")); 
%sysexec(%str(del %"&workdir\mnps.Rdata%")); 

%if not(%sysevalf(%superq(plotname)=,boolean)) %then %do;
        %let plotpath = %getpath(&plotname);
	* check if path is valid, if one is given;
	%if not(%sysevalf(%superq(plotpath)=, boolean)) %then %checkdir(&plotpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let plotname=%sysfunc(translate(&objpath,/,\))/&plotname;
        %else %let plotname=&cd&plotname;

        %sysexec(%str(del %"&plotname%"));
%end;
	

* export dataset for R;
PROC EXPORT DATA=_inputds OUTFILE=rdata DBMS=csv replace;
run;


* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
        put "library(twang)";
*        put "library(twang, lib.loc=%str(%"&lib_loc%"))";
        put ;
        put ;
	put "set.seed(1)";
	put ;
	put "%str(inputds<-read.csv(%"&workdir/datafile.csv%"))";
  put ;
	put "vnames <- names(inputds)          ";
  put "names(inputds) <- tolower(names(inputds))";
  put ;
%if &chkclass = 0 %then
        put %unquote(%str(%'inputds[, %combine(&class)] <- lapply(inputds[,%combine(&class), drop=F], as.factor)%')); ;
        put %unquote(%str(%'inputds[,"&treatvar"] <- as.factor(inputds[, "&treatvar"])%')); ;
        put ;
	put "mnps1 <- mnps(&formula,";
	put " data = inputds,";
	put " n.trees = &ntrees,";
	put " interaction.depth = &intdepth,";
	put " shrinkage = &shrinkage,";
	put " perm.test.iters = &permtestiters,";
	put " " %unquote(%str(%'stop.method = %combine(&stopmethod)%')) ",";
	put " " %unquote(%str(%'estimand = %"&estimand%"%')) ",";
/*        put " sampw = "; %if &sampw=%str() %then put "NULL"; %else put %"&sampw%";; */
        put " sampw = " %if &sampw=%str() %then "NULL,"; %else "inputds$&sampw,"  ;; 
	put " treatATT = " %if %upcase(&estimand) = ATE %then "NULL,";
              %else %if &treatvar_type = C %then "%str(%"&treatatt%"),"; %else "&treatatt,";;
        put " verbose = FALSE";
	put " )";
	put ;
        put "if(%str(%"&collapseto%" != %"none%")){                                 ";
	put "%str(baltab<-bal.table(mnps1, collapse.to=%"&collapseto%"))            ";
  ** when collapseto = stop.method there are no variable names **;
  %if &collapseto ^= stop.method %then %do;
	put "bnames <- as.character(baltab$var)                                     ";
	put "bnames1 <- sapply(strsplit(bnames, ':'), function(x){return(x[[1]])})  ";
	put "bnames1 <- vnames[match(bnames1, tolower(vnames))]                     ";
	put "substr(bnames, 1, nchar(bnames1)) <- bnames1                           ";
	put "baltab$var <- bnames                                                   ";
  %end;  ** ends if collapseto ^= stop.method;
        put "baltab[baltab==Inf] <- NA";
        put "baltab[baltab==(-Inf)] <- NA";
	put "%str(write.table(baltab,file=%"&workdir/baltab.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
        put "}                                                                      ";
	put ;
        put "w <- sapply(mnps1$stopMethods, get.weights, ps1=mnps1)                 ";
	put "w<-as.data.frame(w)                                                    ";
        put "names(w) <- paste(mnps1$stopMethods, mnps1$estimand, sep='_')          ";
	put "w$tempid<- inputds$tempid";
	put ;
	put "%str(write.table(w,file=%"&workdir/wts.csv%",row.names=FALSE,col.names=TRUE,sep=%',%'))";
	put "if(%str(%"&estimand%" == %"ATE%")){                                           ";
        put "   summ <- summary(mnps1)                                              ";
        put "   summ1 <- summ[[1]]                                                  "; 
        put "   summ2 <- summ[[2]]                                                  ";
        put "   %str(write.table(summ1,file=%"&workdir/summary1.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
        put "   %str(write.table(summ2,file=%"&workdir/summary2.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
        put "}else{                                                                 ";
        put "   summ<-summary(mnps1)                                                ";
        put "   ctx <- nrow(summ$summaryList[[1]])                                  ";
        put "   ctx <- rep(summ$levExceptTreatATT, each=ctx)                        ";
        put "   summ <- do.call(rbind, summ$summaryList)                            ";
        put "   tmp <- row.names(summ)                                              ";
        put "   rownames(summ) <- NULL                                              ";
        put "   summ <- data.frame(comp_treat=ctx, row_name=tmp, summ)              ";
        put "   %str(write.table(summ,file=%"&workdir/summary.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
        put "}                                                                      ";
%if &return_ps = TRUE %then %do;
        put "txlevids <- data.frame(txlev=mnps1$treatLev, id=rank(mnps1$treatLev))  ";
        put "tempID <- rownames(mnps1$data)                                         ";
        put "psests <- data.frame(tempID=tempID)                                    ";
        put "txxwalk <- NULL                                                        ";
        put "for(i in 1:length(mnps1$psList)){                                      ";
        put "    u <- names(mnps1$psList)[i]                                        ";
        put "    ps <- mnps1$psList[[u]]$ps                                         ";
        put "%str(    txid <- txlevids[which(txlevids$txlev==u),%"id%"]            )";
        put "%str(    names(ps) <- paste(%"ps%", txid, gsub(%"\\.%", %"_%", names(ps)), sep=%"_%")       )";
        put "    txxwalk <- rbind(txxwalk, cbind(u, names(ps)))                     ";
    %if &estimand = ATT %then %do;
        put "%quote(    ttID <- tempID[which(mnps1$data[,mnps1$treat.var] %in% c(u, mnps1$treatATT))]      )";
        put "    ps$tempID <- ttID                                                  ";
        put "%str(    psests <- merge(psests, ps, by.x=%"tempID%", by.y=%"tempID%", all.x=TRUE, sort=FALSE)  )";
    %end; * if ATT;   
    %if &estimand = ATE %then put "    psests <- data.frame(psests, ps)             ";;
        put "    }                                                                  ";
	put "psests$tempID<- inputds$tempid                                         ";
        put "txxwalk <- as.data.frame(txxwalk)                                      ";
        put "%str(names(txxwalk) <- c(%"Tx_Level%", %"Prop_Score_Variable_Name%")  )";
        put "%str(write.table(psests,file=%"&workdir/psests.csv%", row.names=FALSE, col.names=TRUE, na=%'.%', sep=%',%'))";
        put "%str(write.table(txxwalk,file=%"&workdir/txxwalk.csv%", row.names=FALSE, col.names=TRUE, sep=%',%'))";
%end; ** ends if return_ps = TRUE;

* if plotname parameter is defined, run a plot as well;
%if not(%sysevalf(%superq(plotname)=,boolean)) %then %do;
/* ** moved to up to the delete old files code **;
        %let plotpath = %getpath(&plotname);
	* check if path is valid, if one is given;
	%if not(%sysevalf(%superq(plotpath)=, boolean)) %then %checkdir(&plotpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let plotname=%sysfunc(translate(&objpath,/,\))/&plotname;
        %else %let plotname=&cd&plotname;
*/
	%put Writing plot to &plotname.;
	put "%str(pdf(%'&plotname%'))";
	put "plot(mnps1,plots=1, multiPage=TRUE)";
	put "plot(mnps1,plots=2, multiPage=TRUE)";
	put "plot(mnps1,plots=3, multiPage=TRUE)";
	put "plot(mnps1,plots=4, multiPage=TRUE)";
	put "plot(mnps1,plots=5, multiPage=TRUE)";
	put "dev.off()";
%end;	

* if objpath parameter is defined, save the workspace image;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
	put "%str(save(mnps1, file=%'&workdir/mnps.RData%'))";
%end;	
	
run;


* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/mnps.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/mnps.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/mnps.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/mnps.Rout);
run;

* input weights generated by R;
proc import datafile="&workdir/wts.csv"
     out=_weights
     dbms=csv
     replace;
     getnames=yes;
     datarow=2;
     guessingrows=MAX;

proc sort data=_weights;
	by tempID;	

* merge weights back on to analysis data set;
data &output_dataset;
	merge _inputds(in=a) _weights(in=b);
	by tempID;
	%if ^(&return_ps = TRUE) %then drop tempID;;
run;


%if &return_ps = TRUE %then %do;
 
  ** Import the crosswalk from treatment levels to propensity score variables **;
proc import file="&workdir/txxwalk.csv"
     out=_txxwalk
     dbms=csv
     replace;
     getnames=yes;
     datarow=2;
     guessingrows=MAX;
run;

  ** Print crosswalk to the log **;
data _null_;
   set _txxwalk end=end;;
   if _n_ = 1 then do;
      put "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      put "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
      put "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      put "     CROSSWALK FROM TREATMENT LEVELS TO PROPENSITY SCORE VARIABLES     ";
      put "  Tx Level                               Propensity Score Variable"     ;
      end;
   put    @2Tx_Level @50 Prop_Score_Variable_Name ;
   if end then do; 
      put "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      put "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
      put "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
      end;
run;

* input propensity scores generated by R;
proc import datafile="&workdir/psests.csv"
     out=_psests
     dbms=csv
     replace;
     getnames=yes;
     datarow=2;
     guessingrows=MAX;

proc sort data=_psests;
	by tempID;	

* merge propensity scores back on to output_dataset data set;
data &output_dataset;
	merge &output_dataset(in=a) _psests(in=b);
	by tempID;
	drop tempID;
run;
%end;  ** if ps_return=TRUE ;

* input summary table generated by R;
%if(&estimand = ATT) %then %do;
proc import datafile="&workdir/summary.csv"
     out=_summ
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;
run;

proc print data=_summ;
title "Summary table";
title2 "Summary of observations receiving each other treatment";
title3 "weighted to match the observations receiving treatment &treatatt.";
run;
%end;

%if(&estimand = ATE) %then %do;
proc import datafile="&workdir/summary1.csv"
     out=_summ1
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;
run; 

proc print data=_summ1;
title "Summary table";
title2 "Summary of pairwise comparisons";
run;

proc import datafile="&workdir/summary2.csv"
     out=_summ2
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;
run; 

proc print data=_summ2;
title "Summary table";
title2 "Sample sizes and effective sample sizes";
run;
%end;

* If collapseto ^= none input balance tables generated by R;
%if &collapseto ^= none %then %do;
proc import datafile="&workdir/baltab.csv"
     out=_baltab
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

 ** Different Collapses Require Different Pagination for Printing bal.table **;

%if (&collapseto = covariate) or (&collapseto = pair) %then %do;
* print balance tables by stop method;
proc sql noprint;
select distinct stop_method
into :elements separated by " "
from _baltab
order by stop_method desc;

  %let j=1;
  %let key=%scan(&elements.,&j.,%str( ));
  %do %until(&key eq %nrstr( ));
    
    proc print data=_baltab;
		where stop_method="&key";
		title "Balance table: &key";
		run;

    %let j=%eval(&j+1);
    %let key=%scan(&elements.,&j.,%str( ));
  %end;

%end;  * ends if collapseto = pair or collapseto = covariate;

%if &collapseto = stop.method %then %do;

proc print data=_baltab;
title "Balance table";
run;

%end; * ends if collapseto = stop.method;

title;
run;
%end;  ** Ends if collapseto ^= none **;

%mend mnps;
**##################################################################**;



**##################################################################**;
** Macro: mnbaltable
** Fits generates balance tables for multiple treatments following mnps;


* macro for checking parameters values are allowable -- postive numeric variable ;

%macro mnbaltable(inputobj=, 
            collapseto=pair,
            subset_var=,
            subset_treat=, 
            subset_stop_method=, 
            es_cutoff=, 
            ks_cutoff=, 
            p_cutoff=, 
            ks_p_cutoff=, 
            Rcmd=, 
            objpath= );

options xsync noxwait;

* get cwd;
%getcd;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&inputobj, inputobj)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let inputobj =%qsysfunc(compress(&inputobj,%str(%"%')));
%let collapseto =%qsysfunc(compress(&collapseto,%str(%"%')));
%let subset_var =%qsysfunc(compress(&subset_var,%str(%"%')));
%let subset_treat =%qsysfunc(compress(&subset_treat,%str(%"%')));
%let subset_stop_method =%qsysfunc(compress(&subset_stop_method,%str(%"%')));
%let es_cutoff =%qsysfunc(compress(&es_cutoff,%str(%"%')));
%let ks_cutoff =%qsysfunc(compress(&ks_cutoff,%str(%"%')));
%let p_cutoff =%qsysfunc(compress(&p_cutoff,%str(%"%')));
%let ks_p_cutoff =%qsysfunc(compress(&ks_p_cutoff,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));
%let inputobj = %sysfunc(translate(&inputobj,/,\));
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let objpath = %sysfunc(translate(&objpath,/,\));



* validate file names;
 ** Check the path for inputobj name if it contains a path **;
%let inputpath = %getpath(&inputobj);
* check if path is valid, if one is given;
%if not(%sysevalf(%superq(inputpath)=, boolean)) %then %checkdir(&inputpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let inputobj=&objpath/&inputobj;
  %else %let inputobj=&cd&inputobj;

* check to make sure inputobj exists;
%if %sysfunc(fileexist(&inputobj)) eq 0  %then %do;
	%put ERROR: The input object "&inputobj" could not be found. Please check your path names or specify an object path.;
  %abort;
%end;

** validate parameter values ;
   ** Check collapseto **;
%checksingle(collapseto,&collapseto, pair covariate stop.method);
   ** Check es_cutoff **;
%if &es_cutoff ^= %str() %then %checkpnum(es_cutoff, &es_cutoff);
   ** Check ks_cutoff **;
%if &ks_cutoff ^= %str() %then %checkpnum(ks_cutoff, &ks_cutoff);
   ** Check p_cutoff **;
%if &p_cutoff ^= %str() %then %checkpnum(p_cutoff, &p_cutoff);
   ** Check ks_p_cutoff **;
%if &ks_p_cutoff ^= %str() %then %checkpnum(ks_p_cutoff, &ks_p_cutoff);



* Remove existing versions of files to be created *;
%sysexec(%str(del %"&workdir\baltab.csv%")); 
%sysexec(%str(del %"&workdir\mnbaltable.R%")); 
%sysexec(%str(del %"&workdir\mnbaltable.Rout%")); 

filename rscript "&workdir/mnbaltable.R";

/******************************************************************************************************************
** case issues with subset_var  **
** With ps and mnps users can specify variables without respect to case and we turn all the names to lower case.
** If the user is running mnbaltable with an object we created the variables will be in lower cases. But the
** user could use an mnps object from some other source. They might know R is case sensitive and list variables
** using the proper case sensitive names. Or they might not.  So we check if the variable names in the mnps object
** are all lower case. If so we change our variables list to lower case and work with lower case.  The user
** specified variables do not need to match case. They do need to match when changed to lower case.  If the 
** varnames in the mnps object are not all lower case then the user must match them exactly by name and case.
** If the specifies variables with the SUBSET_VAR parameter, we make the output match their case even if the 
** variable in the mnps object are all lower case. We do not know the users preference for variables names 
** if the do not this parameter so we don't do any name changing. 
******************************************************************************************************************/

* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
        put "library(twang)";
*        put "library(twang, lib.loc=%str(%"&lib_loc%"))";
        put ;
	put "set.seed(1)";
	put "%str(tmp <- load(%"&inputobj%"))";
        put ".tmnts <- NULL                  ";
        put ".vars <- NULL                   ";
        put ".smeths <- NULL                 ";
  ** Add a check that the variables specified in the subset_var argument are valid **;
%if &subset_var ^= %str() %then %do;
        put "vvar.names <- get(tmp)$psList[[1]]$gbm.obj$var.names ";
     ** we need to handle mixed case and lower case var.names differently **;
        put "islower <- !any(tolower(vvar.names) != vvar.names)";
        put %unquote(%str(%'.vars <- %combine(&subset_var) %'))                      ;     
        put "vnames <- .vars                 ";
        put "if(islower){ .vars <- tolower(.vars)}                                 ";
        put "if(sum(.vars %nrbquote(%)in%nrbquote(%) vvar.names) != length(.vars)){ ";
        put "stop(paste(%str(%"One or more values of the SUBSET_VARS argument is not a valid name of the covariates used in the model.)"; 
        put "%str(             SUBSET_VARS=%", .vars, %"variables used in modeling were:%", .vvar.names, %"R is case-sensitive make )";
        put "%str(             sure the case of the variable names match.%")))}";
%end;
  ** Add a check that the treatments specified in the subset_treat argument are valid **;
%if &subset_treat ^= %str() %then %do;
        put %unquote(%str(%' .tmnts <- %combine(&subset_treat) %'));
        put "if(sum(.tmnts %nrbquote(%)in%nrbquote(%) levels(get(tmp)$data[,get(tmp)$treat.var])) != length(.tmnts)){";
        put "%str(stop(%"One or more values of the SUBSET_TREAT argument is not a valid value of the treatment variable%")})";
%end;
  ** Add a check that the stop.methods specified in the subset_stop_method argument are valid **;
%if &subset_stop_method ^= %str() %then %do;
        put %unquote(%str(%' .smeths <- tolower(%combine(&subset_stop_method)) %'));
        put "if(sum(.smeths %nrbquote(%)in%nrbquote(%) get(tmp)$stopMethods) != length(.smeths)){ ";
        put "%str(stop(%"One or more values of the SUBSET_STOP_METHOD argument is not a valid value of the stop methods used in modeling%")})";
%end;
  ** write tbe balance table statement **;
	put "baltab <- bal.table(get(tmp),%str(collapse.to=%"&collapseto%"), subset.var=.vars, subset.treat=.tmnts, subset.stop.method=.smeths";
%if &es_cutoff ^= %str() %then
        put "     , es.cutoff = &es_cutoff"; ;
%if &ks_cutoff ^= %str() %then
        put "     , ks.cutoff = &ks_cutoff"; ;
%if &p_cutoff ^= %str() %then
        put "     , p.cutoff = &p_cutoff,"; ;
%if &ks_p_cutoff ^= %str() %then
        put "     , ks.p.cutoff = &ks_p_cutoff"; ;
        put" ) ";
   ** proc import in SAS will not import any empty file so check that baltab is not empty **;
        put "if(nrow(baltab) == 0){stop('The restrictions on the balance table returned no record. Try alternative selection criteria')}";

** subset.var only works with collapse.to = pair -- so variable names only match subset_var when collapsto=pair **;
%if &collapseto = pair %then %do;
%if &subset_var ^= %str() %then %do;
	put "bnames <- as.character(baltab$var)                                     ";
	put "bnames1 <- sapply(strsplit(bnames, ':'), function(x){return(x[[1]])})  ";
        put "if(islower){bnames1 <- vnames[match(bnames1, tolower(vnames))]}else{   ";
	put "            bnames1 <- vnames[match(bnames1, vnames)]         }        ";
	put "substr(bnames, 1, nchar(bnames1)) <- bnames1                           ";
	put "baltab$var <- bnames                                                   ";
%end;
%end;
        put "baltab[baltab==Inf] <- NA";
        put "baltab[baltab==(-Inf)] <- NA";
	put "%str(write.table(baltab,file=%"&workdir/baltab.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";

run;

* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/mnbaltable.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/mnbaltable.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/mnbaltable.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/mnbaltable.Rout);

run;

* input balance tables generated by R;
proc import datafile="&workdir/baltab.csv"
     out=_baltab
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

 ** Different Collapses Require Different Pagination for Printing bal.table **;

%if (&collapseto = covariate) or (&collapseto = pair) %then %do;
* print balance tables by stop method;
proc sql noprint;
select distinct stop_method
into :elements separated by " "
from _baltab
order by stop_method desc;

  %let j=1;
  %let key=%scan(&elements.,&j.,%str( ));
  %do %until(&key eq %nrstr( ));
    
    proc print data=_baltab;
		where stop_method="&key";
		title "Balance table: &key";
		run;

    %let j=%eval(&j+1);
    %let key=%scan(&elements.,&j.,%str( ));
  %end;

%end;  * ends if collapseto = pair or collapseto = covariate;

%if &collapseto = stop.method %then %do;

proc print data=_baltab;
title "Balance table";
run;

%end; * ends if collapseto = stop.method;

title;
run;
%mend mnbaltable;
**##################################################################**;


**##################################################################**;
** Macro: mnplot
** Fits generates diagnostic plots for 3+ treatmens per user specifications;

%macro mnplot(inputobj=, 
            plotname=, 
            plotformat=, 
            plots=, 
            subset=, 
            color=TRUE,
            pairwisemax=TRUE,
            treatments=,
            figurerows=1,
            singleplot=,
            multipage=FALSE,
            Rcmd=, 
            objpath= );

options xsync noxwait;

* get cwd;
%getcd;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&inputobj, inputobj)
%notblank(&plotname, plotname)
%notblank(&plots, plots)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let inputobj =%qsysfunc(compress(&inputobj,%str(%"%')));
%let plotname  =%qsysfunc(compress(&plotname,%str(%"%')));
%let plotformat =%qsysfunc(compress(&plotformat,%str(%"%')));
%let plots =%qsysfunc(compress(&plots,%str(%"%')));
%let subset =%qsysfunc(compress(&subset,%str(%"%')));
%let color =%qsysfunc(compress(&color,%str(%"%')));
%let pairwisemax =%qsysfunc(compress(&pairwisemax,%str(%"%')));
%let treatments =%qsysfunc(compress(&treatments,%str(%"%')));
%let singleplot =%qsysfunc(compress(&singleplot,%str(%"%')));
%let multipage =%qsysfunc(compress(&multipage,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate enumerated parameters;
%checksingle(plots,&plots,1 2 3 4 5 optimize boxplot es t ks);
%checksingle(color, &color, T F TRUE FALSE);
%checksingle(pairwisemax, &pairwisemax, T F TRUE FALSE);
%checksingle(figurerows, &figurerows, 1 2 3 4);
%checksingle(multipage, &multipage, T F TRUE FALSE);
%if &plots=optimize or &plots=boxplot or &plots=es or &plots=t or &plots=ks or &plots=histogram %then %do;
  %let plots=%str(%"&plots%");
%end;

* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));
%if not(%sysevalf(%superq(plotname)=,boolean)) %then %let plotname = %sysfunc(translate(&plotname,/,\));
%let inputobj = %sysfunc(translate(&inputobj,/,\));
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let objpath = %sysfunc(translate(&objpath,/,\));



* validate file names;
 ** Check the path for inputobj name if it contains a path **;
%let inputpath = %getpath(&inputobj);
* check if path is valid, if one is given;
%if not(%sysevalf(%superq(inputpath)=, boolean)) %then %checkdir(&inputpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let inputobj=&objpath/&inputobj;
  %else %let inputobj=&cd&inputobj;

* check to make sure inputobj exists;
%if %sysfunc(fileexist(&inputobj)) eq 0  %then %do;
	%put ERROR: The input object "&inputobj" could not be found. Please check your path names or specify an object path.;
  %abort;
%end;



 ** Check the path for plot name if it contains a path **;
%let plotpath = %getpath(&plotname);
* check if path is valid, if one is given;
%if not(%sysevalf(%superq(plotpath)=, boolean)) %then %checkdir(&plotpath);
	%else %if not(%sysevalf(%superq(objpath)=,boolean)) %then %let plotname=&objpath/&plotname;
  %else %let plotname=&cd&plotname;
	%put Writing plot to &plotname.;

  ** Check the values if treatments is specified **;
%let treatments_quote = NULL;

%if &treatments ^= %str() %then %do;

     ** Count the number of specified treatment is must be one or two **;
   %let treatments_cnt=%sysfunc(countw(&treatments));
   %if &treatments_cnt > 2 %then %do; 
      %put ERROR: 1 or 2 treatments must be specified by the TREATMENTS argument;  
      %abort;
   %end;
   
   ** Next statements only run if 1 or 2 treatments specified by TREATMENTS argument **;
   %let tmp = %scan(&treatments,1);
   %let treatments_quote = %str(%"&tmp%");
   %if &treatments_cnt = 2 %then %do;
       %let tmp = %scan(&treatments,2) ;
       %let treatments_quote = c(&treatments_quote, %str(%"&tmp%")) ;
       %end;

%end; ** ends if treatments not blank ;

** Check the value of singleplot it must be an integer if not blank **;

%if &singleplot = %str() %then %let singleplot = NULL;
%else %checkpint(singleplot, &singleplot) ;

** Check the values of subset **;

%if &subset = %str() %then %let subset = NULL;
%if &subset ^= NULL %then %do;
  ** check if it is an integer **;
  %if &subset = 0 %then %do;
      %put ERROR: The value 0 is not a valid value for the subset parameter. Please use positive integers or stopping rule.;
      %abort;
      %end;
  %let chkval = %eval(%sysfunc(verify(%sysfunc(trim(%sysfunc(compress(&subset)))),'0123456789')) +
                        %index(%substr(&subset,%index(&subset,.)+1), .));	
  %if &chkval = 0 %then %let subset = %ncombine(&subset);
  ** If not an integer check for valid values **;
  %if &chkval > 0 %then %do ;
      %checkmult(subset,&subset,ks.mean es.mean ks.max es.max ks.max.direct es.max.direct);
      %let subset = %combine(&subset);
      %end;
%end;  

filename rscript "&workdir/mnplot.R";

* determine plot format if it exists;
%if %lowcase(&plotformat) eq jpg %then %let fmt = jpeg;
%else %if %lowcase(&plotformat) eq pdf %then %let fmt = pdf;
%else %if %lowcase(&plotformat) eq png %then %let fmt = png;
%else %if %lowcase(&plotformat) eq wmf %then %let fmt = win.metafile;
%else %if %lowcase(&plotformat) eq postscript %then %let fmt = postscript;
%else                                         %let fmt = pdf;



* Remove existing versions of files to be created *;
%sysexec(%str(del %"&plotname%")); 
%sysexec(%str(del %"&workdir\mnplot.R%")); 
%sysexec(%str(del %"&workdir\mnplot.Rout%")); 

* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
        put "library(twang)";
*        put "library(twang, lib.loc=%str(%"&lib_loc%"))";
        put ;
	put "set.seed(1)";
	put "%str(tmp <- load(%"&inputobj%"))";
  ** Add a check that the treatments specified in the TREATMENTS argument are valid **;
%if &treatments ^= %str() %then %do;
        put %unquote(%str(%'.tmnts <- &treatments_quote %'));
        put "if(sum(.tmnts %nrbquote(%)in%nrbquote(%) levels(get(tmp)$data[,get(tmp)$treat.var])) != length(.tmnts)){";
        put "%str(stop(%"One or more values of TREATMENTS argument is not a valid value of the treatment variable%")})";
%end;
        put %unquote(%str(%'.subs <- &subset %'));
        put "if (is.numeric(.subs) & any(.subs > length(get(tmp)$stopMethods))){";
        put "%str(stop(%"One or more values of SUBSET argument is greater than the number of stop methods%")})";
        put "if (is.character(.subs) & !all(.subs %nrbquote(%)in%nrbquote(%) get(tmp)$stopMethods)){";
        put "%str(stop(%"One or more values of SUBSET argument is not a stop method used to fit the model%")})";

	put "%str(&fmt(%'&plotname%'))";
	put "%str(plot(get(tmp),plots=&plots, subset=.subs, color=&color, pairwiseMax=&pairwisemax, treatments=&treatments_quote, figureRows=&figurerows, singlePlot=&singleplot, multiPage=&multipage))";
	put "dev.off()";
run;

* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/mnplot.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/mnplot.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/mnplot.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/mnplot.Rout);

run;
%mend mnplot;
**##################################################################**;



**##################################################################**;
** Macro: CBPS
** Fits propensity scores using CBPS;

%macro CBPS(treatvar=,
          vars=,
          class=,
          dataset=, 
          estimand=ATE, 
          method=over,
          output_dataset=_inputds,
          permtestiters=0, 
          Rcmd=, 
          objpath=);


options xsync noxwait;

/*
* get cwd;
%getcd;
*/
** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

* validate parameters;
*Check that key parameters are not blank;
%notblank(&treatvar, treatvar)
%notblank(&vars, vars)
%notblank(&dataset, dataset)
%notblank(&estimand, estimand)
%notblank(&Rcmd, Rcmd)

* validate parameters;
* remove any extraneous quotes;
%let treatvar =%qsysfunc(compress(&treatvar,%str(%"%')));
%let vars =%qsysfunc(compress(&vars,%str(%"%')));
%let dataset  =%qsysfunc(compress(&dataset,%str(%"%')));
%let estimand =%upcase(%qsysfunc(compress(&estimand,%str(%"%'))));
%let method =%qsysfunc(compress(&method,%str(%"%')));
%let objpath =%qsysfunc(compress(&objpath,%str(%"%')));
%let permtestiters =%qsysfunc(compress(&permtestiters,%str(%"%')));

* make sure R path is correct;
%chkRcmd(&Rcmd);

* validate enumerated parameters;
%checksingle(estimand,&estimand,ATT ATE);
%checksingle(method,&method,over exact);

* validate object path if it exists;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
%checkdir(&objpath);
%end;

* if an object path is given, use it for all output, otherwise use SAStmp as default;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %let workdir = &objpath;
%else %let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

* fix path name separators;
%let workdir = %sysfunc(translate(&workdir,/,\));


* check that the weight variables are not in the specified output_dataset *;
%checkwgt(cbps, &estimand, &output_dataset);

* check that the covariates are in the dataset;

%checkexist(&vars,&dataset);


* check to make sure treatment variable does not contain missings;
   ** check it is in the dataset **;
%checkexist(&treatvar,&dataset);

   ** check the values **;
proc sql noprint;
select count(*) into :nmiss from &dataset where &treatvar not in (0,1);
%if &nmiss ne 0 %then %do;
  %put ERROR: Treatment variable &treatvar has invalide values (%trim(&nmiss) cases).;
  %put        Treatment variable must take values 0 or 1 only.;
  %abort;
%end;

   ** check the class variables **;
%let chkclass = %sysevalf(&class=, boolean);
%if &chkclass = 0 %then %do;
  %let class =%qsysfunc(compress(&class,%str(%"%')));
  %let notin = %setdiff(&class, &vars);
  %if  &notin ^= %str() %then %do;
     %put ERROR: CLASS variables &notin are not in vars parameter list; 
     %abort;
     %end;
  %end;
 

** Set all variables used in R script to lower case **;
%let treatvar = %lowcase(&treatvar);
%let vars = %lowcase(&vars);
%if &chkclass = 0 %then %let class = %lowcase(&class);

%let formula = &treatvar ~ %addplus(&vars);

filename rdata "&workdir/datafile.csv";
filename rscript "&workdir/cbps.R";


* append an internal ID to supplied analysis data set;
data _inputds;
	length tempID 8.;
	set &dataset;
	tempID=_N_;

proc sort;
	by tempID;
run;

* Remove existing versions of files to be created *;
%sysexec(%str(del %"&workdir\datafile.csv%")); 
%sysexec(%str(del %"&workdir\dxsummary.csv%")); 
%sysexec(%str(del %"&workdir\dxwtsbaltab.csv%")); 
%sysexec(%str(del %"&workdir\wts.csv%")); 
%sysexec(%str(del %"&workdir\cbps.R%")); 
%sysexec(%str(del %"&workdir\cbps.Rout%")); 
	

* export dataset for R;
PROC EXPORT DATA=_inputds OUTFILE=rdata DBMS=csv replace;
run;


* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"CBPS%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"CBPS%", repos=%"http://cran.us.r-project.org%"))";
	put "%str(if (!is.element(%"twang%", installed.packages(lib.loc=%"&twangdir%")[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=c(%"twang%",%"CBPS%"))                )";
	put "library(CBPS)";
	put "library(twang)";
        put ;
        put ;
	put "set.seed(1)";
	put ;
	put "%str(inputds<-read.csv(%"&workdir/datafile.csv%"))";
  put ;
	put "vnames <- names(inputds)          ";
  put "names(inputds) <- tolower(names(inputds))";
  put ;
%if &chkclass = 0 %then
        put %unquote(%str(%'inputds[, %combine(&class)] <- lapply(inputds[,%combine(&class), drop=F], as.factor)%')); ;
        put ;
	put "cbps1 <- CBPS(&formula,";
	put " data = inputds,";
%if %lowcase(&estimand) = att %then
	put " ATT = TRUE, ";
%else   put " ATT = FALSE, ";;
        put " type = 'propensity'";
	put " )";
	put ;
	put "w<-data.frame(cbps_&estimand=cbps1$weights, tempID=inputds$tempid)";
        put "%str(tmp <- tapply(w[,1], inputds[,%"&treatvar%"], mean))";
        put "%str(tmp <- ifelse(inputds[,%"&treatvar%"]==min(inputds[,%"&treatvar%"]), tmp[1], tmp[2]))";
        put "w[,1] <- w[,1]/tmp";
	put "%str(write.table(w,file=%"&workdir/wts.csv%",row.names=FALSE,col.names=TRUE,sep=%',%'))";

* if objpath parameter is defined, save the workspace image;
%if not(%sysevalf(%superq(objpath)=,boolean)) %then %do;
	put "%str(save(cbps1, file=%'&workdir/cbps.RData%'))";
%end;	

        put "%str(x<-as.matrix(subset(w ,select=c(%"cbps_&estimand%"))))";
	put "dxtmp <- dx.wts(x, ";
	put " data=inputds, ";
	put " " %unquote(%str(%'vars = %combine(&vars)%')) ",";
	put " " %unquote(%str(%'estimand = %"&estimand%"%')) ",";
	put " " %unquote(%str(%'treat.var = %"&treatvar%"%')) ",";
	put " x.as.weights=T, ";
	put " perm.test.iters=&permtestiters";
	put ") ";
	put ;
	put "baltab<-bal.table(dxtmp)";
	put ;
	put "bnames <- rownames(baltab$unw)                                         ";
	put "bnames1 <- sapply(strsplit(bnames, ':'), function(x){return(x[[1]])})  ";
	put "bnames1 <- vnames[match(bnames1, tolower(vnames))]                     ";
	put "substr(bnames, 1, nchar(bnames1)) <- bnames1                           ";
	put "baltab <- lapply(baltab, function(u){                                  ";
        put "                 rownames(u) <- bnames                                 ";
        put "                 return(u)})                                           ";
	put ;
	put "baltab <- data.frame(do.call(rbind, baltab), table.name=rep(names(baltab), each=nrow(baltab[[1]])))";
	put "baltab <- data.frame(row_name=row.names(baltab), baltab)";
	put "%str(write.table(baltab,file=%"&workdir/dxwtsbaltab.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
	put "%str(write.table(dxtmp$summary.tab,file=%"&workdir/dxsummary.csv%",row.names=FALSE,col.names=TRUE,sep=%',%',na=%'.%'))";
run;


* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/cbps.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/cbps.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/cbps.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/cbps.Rout);
run;

* input weights generated by R;
proc import datafile="&workdir/wts.csv"
     out=_weights
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

proc sort data=_weights;
	by tempID;	

* merge weights back on to analysis data set;
data &output_dataset;
	merge _inputds(in=a) _weights(in=b);
	by tempID;
	drop tempID;
run;

* input summary table generated by R;
proc import datafile="&workdir/dxsummary.csv"
     out=_dxsumm
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

** drop the iter variable since it is missing for all records **;
data _dxsumm;
   set _dxsumm;
   drop iter;
run;

proc print data=_dxsumm;
title "dxwts Summary table";
run;
	
* input balance tables generated by R;
proc import datafile="&workdir/dxwtsbaltab.csv"
     out=_baltab
     dbms=csv
     replace;
     getnames=yes;
	   datarow=2;
     guessingrows=MAX;

* print balance tables by stop method;
proc sql noprint;
select distinct table_name
into :elements separated by " "
from _baltab where(table_name ^= "unw")
order by table_name desc;

%let elements = unw &elements;

  %let j=1;
  %let key=%scan(&elements.,&j.,%str( ));
  %do %until(&key eq %nrstr( ));
    
    proc print data=_baltab;
		where table_name="&key";
		title "dxwts Balance table: &key";
		run;

    %let j=%eval(&j+1);
    %let key=%scan(&elements.,&j.,%str( ));
  %end;

title;
run;

%mend CBPS;
**##################################################################**;



**##################################################################**;
** Macro:  update_twang
** Updates the twang package in R ;

%macro update_twang(Rcmd=); 

options xsync noxwait;

** Check if the c:\user\AppData\Local\TWANG exists and create it if not **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) eq 0  %then x "mkdir &twangdir";;
%let twangdir = %sysfunc(translate(&twangdir,/,\));

%let workdir = %sysfunc(getoption(work)); * Path to SAStmp;

filename rscript "&workdir/update.R";


* generate R script based on supplied parameters;
data _null_;
	file rscript lrecl=20000;
        put "options(warn=1)";
        put "%str(.libPaths(%"&twangdir%"))";
	put "%str(if (!is.element(%"twang%", installed.packages()[,1])) install.packages(%"twang%", repos=%"http://cran.us.r-project.org%"))";
** Make sure the version of twang in &twangdir is always up to date -- it will be the package used by default if it exists **;
        put "%str(update.packages%(lib.loc=%"&twangdir%",                       )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=%"&twangdir%",                        )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
	put "library(twang)";
        put ;
        put "%str(pack <- installed.packages()                                  )";
        put "%str(detach(%"package:twang%", unload=TRUE)                        )";
        put ;
        put "%str(update.packages%(lib.loc=pack[%"twang%",%"LibPath%"],         )";
        put "%str(                repos=%"http://cran.us.r-project.org%",       )"; 
        put "%str(                instlib=pack[%"twang%",%"LibPath%"],          )"; 
        put "%str(                ask=F,                                        )";
        put "%str(                oldPkgs=%"twang%"%)                           )";
run;
 
* Call R using generated script;
%put Calling R with the command: %str(%"&Rcmd%" CMD BATCH --vanilla %"&workdir/update.R%");
%sysexec(%str(echo Starting R, please wait... & %"&Rcmd%" CMD BATCH --vanilla %"&workdir/update.R%"));
%if &sysrc eq 0 %then %put R command completed successfully;
%else %do;
	%put ERROR: R command did not complete successfully.;
	%put  Return message from R is as follows:;
	%put ;
	%printerr(&workdir/update.Rout);
	%put;
	%abort;
%end;

* look for warning messages;
%printwarn(&workdir/update.Rout);
run;

%mend update_twang;
**##################################################################**;


**##################################################################**;
** Macro:  remove_twang_folder
** Removes twang directory ;

%macro remove_twang_folder;
** Check if the c:\user\AppData\Local\TWANG exists and remove it if it does **;
%let user = %sysget(USERPROFILE);
%let twangdir = &user\AppData\Local\TWANG;
%if %sysfunc(fileexist(&twangdir)) neq 0  %then x "rmdir /S &twangdir";;

%mend remove_twang_folder;
**##################################################################**;

**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;
**++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++**;


{smcl}
{* version 1.2.0 29November2023}{...}
{cmd:help xtpb}{right: ({browse "https://doi.org/10.1177/1536867X251322965":SJ25-1: st0768})}
{hline}

{title:Title}

{p2colset 5 13 15 2}{...}
{p2col :{cmd:xtpb} {hline 2}}Pooled Bewley estimator of long-run relationships
in dynamic heterogeneous panels


{title:Syntax}

{p 8 12 2}
{cmd:xtpb} {depvar} [{indepvars}] {ifin} 
[{cmd:,} {it:options}]

{p 4 4 2}
The {cmd:xtpb} command supports balanced and unbalanced panel data.  You must
{cmd:xtset} your data before using {cmd:xtpb}; see {helpb xtset}.  {cmd:xtpb}
uses Mata functions from the {cmd:moremata} package (Jann 2005).  Variable inputs are also
compatible with time-series operators.  The postestimation command
{cmd:predict} cannot be used to obtain residuals of the group-specific error
correcting representations.  These residuals and the demeaned
error-correcting terms (that is, deviations from the long run) can be obtained
using the options below.  The postestimation command predict for fitted values
returns b*x_it.

{synoptset 37}{...}
{synopthdr}
{synoptline}
{synopt:{opt lag:order(#)}}specify the number of lags for the dependent
variable and regressors in the level representation; default is 
{cmd:lagorder(1)} (one lag in level representation, which corresponds to zero
lag order for the first-difference error-correction
representation){p_end}
{synopt:{opt bias:correct(string)}}use small-T bias correction;
{cmd:biascorrect(none)} implements no bias
correction; {cmd:biascorrect(jackknife)} implements half-panel jackknife bias
correction; {cmd:biascorrect(bootstrap)} implements bias correction based on
stochastic simulation (bootstrapping); default is
{cmd:biascorrect(none)}{p_end}
{synopt:{opt full:display}}display full regression output, including
group-specific short-run coefficients from error-correction regressions{p_end}
{synopt:{opt err:orcorrect(string)}}generate a variable named {it:string} that
contains the error-correction terms{p_end}
{synopt:{opt res:iduals(string)}}generate a variable named {it:string} that
contains residuals derived from the error-correction regressions{p_end}
{synopt:{opt boot:strap(#reps [, bootstrap_options])}}compute bootstrapped
confidence intervals (CIs) in addition to the asymptotic CIs; number of
bootstrap replications must be specified;
a large value is recommended (at least 2,000); 
bootstrapped CIs are automatically computed when {cmd:biascorrect(bootstrap)} is specified{p_end}
{synoptline}

{synoptset 20}{...}
{synopthdr:bootstrap_options}
{synoptline}
{synopt:{opt csrobust}}specify the resampling of residuals in the
bootstrapping algorithm for computation of bootstrapped CIs; by default,
it is assumed there is no cross-sectional dependence of errors;
{cmd:csrobust} allows for arbitrary cross-sectional dependence of errors
by resampling column vectors of cross-sectionally stacked residuals{p_end}
{synopt:{opt btx(string)}}specify the bootstrapping algorithm for
resampling regressors in x_it; {cmd:btx(fixed)} or no {cmd:btx()}
specification conducts bootstrapping conditional on x_it; namely, regressors
are fixed across the bootstrap replications; {cmd:btx(varx)} resamples
regressors in bootstrap replications according to the VAR model in first
differences of regressors; {cmd:btx(varxy)} resamples regressors in bootstrap
replications according to the VAR model in first differences of regressors
augmented with lags of the first differences of the dependent variable;
default is {cmd:btx(fixed)}{p_end}
{synopt:{opt btx_lagorder(#)}}specify lag order for the marginal model for
regressors used for resampling in bootstrapping; lag order is specified in the
level representation, similarly to {cmd:lagorder()}; by default, it is set equal to {cmd:lagorder()}{p_end}
{synopt:{opt bcialpha(#)}}set the nominal level for bootstrapped CIs; 
default is {cmd:bcialpha(0.95)}, which computes 95% CIs; any value between 0
and 1 can be chosen{p_end}
{synopt:{opt s:eed(#)}}specify a random seed for reproducibility; default is {cmd:seed(123456)}{p_end}
{synoptline}


{title:Description}

{p 4 4 2}
The pooled Bewley (PB) estimator is for panel-data models with a sufficiently
large time dimension T and a small or large group dimension N.  The PB
estimator is applicable under the same setting as the widely used pooled mean
group estimator of Pesaran, Shin, and Smith (1999).  This setting allows for
group-specific short-run feedbacks from y to x and vice versa, and it assumes
existence of a single homogeneous long-run relationship.  The specification for
the dependent variable is{p_end}

{p 4 4 2}
d(y_it) = c_i - alpha_i*(y_i,t-1 - b*x_i,t-1) + lags of d(y_it) + contemporaneous values and lags of d(x_it) + e_it{p_end}

{p 4 4 2}
{cmd:xtpb} also implements the bias-correction and bootstrapping options
considered by Chudik, Pesaran, and Smith (Forthcoming).


{title:Examples}

{p 4 4 2}PB, no bias correction, asymptotic standard errors{p_end}

{phang2}{cmd:. xtpb y x1 x2, lagorder(1)}

{p 4 4 2}PB with jackknife bias correction{p_end}

{phang2}{cmd:. xtpb y x1 x2, lagorder(1) biascorrect(jackknife)}

{p 4 4 2}PB with jackknife bias correction and cross-sectional robust bootstrapped CIs
(using 3,000 bootstrap replications){p_end}

{phang2}{cmd:. xtpb y x1 x2, lagorder(1) biascorrect(jackknife) bootstrap(3000, btx(varx) btx_lagorder(2) csrobust)}

{p 4 4 2}PB with stochastic simulation bias correction and cross-sectional robust bootstrap
CIs (using 3,000 bootstrap replications){p_end}

{phang2}{cmd:. xtpb y x1 x2, lagorder(1) biascorrect(bootstrap) bootstrap(3000, btx(varx) btx_lagorder(2) csrobust)}

{p 4 4 2}
PB with no bias correction, cross-sectional robust bootstrap 80% CIs (using 3,000
bootstrap replications), generated variable {cmd:error_correct} with error-correction terms, 
generated variable {cmd:residuals} with residuals, 
and display of error-correction regression outputs{p_end}

{phang2}{cmd:. xtpb y x1 x2, lagorder(1) bootstrap(3000, btx(varx) btx_lagorder(2) csrobust bcialpha(0.8)) fulldisplay errorcorrect(error_correct) residuals(residuals)}


{title:Stored results}

{p 4 4 2}
{cmd:xtpb} stores the following in {cmd:e()}:{p_end}

{synoptset 17 tabbed}{...}
{p2col 5 17 21 2: Scalars}{p_end}
{synopt:{cmd: e(N)}}number of groups{p_end}
{synopt:{cmd: e(Tavg)}}average observations per group{p_end}
{synopt:{cmd: e(Tmin)}}minimum observations per group{p_end}
{synopt:{cmd: e(Tmax)}}maximum observations per group{p_end}
{synopt:{cmd: e(bcialpha)}}bootstrap CI nominal level{p_end}

{p2col 5 17 21 2: Macros}{p_end}
{synopt:{cmd: e(cmd)}}{cmd:xtpb}{p_end}
{synopt:{cmd: e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd: e(timevar)}}name of time-period identifier{p_end}
{synopt:{cmd: e(panelvar)}}name of cross-section identifier{p_end}

{p2col 5 17 21 2: Matrices}{p_end}
{synopt:{cmd: e(table)}}matrix containing coefficients with their asymptotic
standard errors, test statistics, {it:p}-values, and CIs{p_end}
{synopt:{cmd: e(V)}}asymptotic variance-covariance matrix of the coefficient vector{p_end}
{synopt:{cmd: e(b)}}coefficient vector{p_end}
{synopt:{cmd: e(BCI)}}bootstrap CIs{p_end}
{synopt:{cmd: e(ec)}}group-specific error-correction coefficients{p_end}
{synopt:{cmd: e(ec_rsq)}}group-specific error-correction regression R-squared
and adjusted R-squared{p_end}
{synopt:{cmd: e(ec_}{it:#}{cmd:)}}group-specific error-correction regression
coefficients with their standard errors, test statistics, {it:p}-values, and
CIs{p_end}

{p2col 5 17 21 2: Functions}{p_end}
{synopt:{cmd: e(sample)}}marks estimation sample{p_end}


{title:References}

{phang}
Chudik, A., M. H. Pesaran, and R. P. Smith. Forthcoming. Pooled Bewley
estimator of long run relationships in dynamic heterogenous panels. 
{it:Econometrics and Statistics}.
{browse "https://doi.org/10.1016/j.ecosta.2023.11.001"}.

{phang}
Jann, B. 2005. moremata: Stata module (Mata) to provide various functions.
Statistical Software Components S455001, Department of Economics,
Boston College. 
{browse "https://ideas.repec.org/c/boc/bocode/s455001.html"}.

{phang}
Pesaran, M. H., Y. Shin, and R. P. Smith. 1999. Pooled mean group estimation
of dynamic heterogeneous panels. 
{it:Journal of the American Statistical Association}
94: 621-634. {browse "https://doi.org/10.2307/2670182"}.


{title:Authors}

{pstd}Priyanka Asnani{p_end}
{pstd}Federal Reserve Bank of Dallas{p_end}
{pstd}Dallas, TX{p_end}
{pstd}{browse "mailto:priyankaasnani13@gmail.com":priyankaasnani13@gmail.com}{p_end}

{pstd}Alexander Chudik{p_end}
{pstd}Federal Reserve Bank of Dallas{p_end}
{pstd}Dallas, TX{p_end}
{pstd}{browse "mailto:alexander.chudik@dal.fed.org":alexander.chudik@dallas.fed.org}{p_end}

{pstd}Braden Strackman{p_end}
{pstd}Federal Reserve Bank of Dallas{p_end}
{pstd}Dallas, TX{p_end}
{pstd}{browse "mailto:braden.strackman@dal.frb.org":braden.strackman@dal.frb.org}{p_end}


{marker see}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 25, number 1: {browse "https://doi.org/10.1177/1536867X251322965":st0768}{p_end}

{p 7 14 2}
Help:  {helpb xtpmg}, {helpb xtpedroni}, {helpb xtcointreg} (if installed){p_end}

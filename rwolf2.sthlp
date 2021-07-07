{smcl}
{* July 07, 2021 @ 10:01:52}{...}
{hline}
help for {hi:rwolf2}
{hline}

{title:Title}

{p 8 20 2}
    {hi:rwolf2} {hline 2} A more flexible syntax to calculate Romano-Wolf stepdown p-values for multiple hypothesis testing

{title:Syntax}

{p 8 14 2}
{cmd:rwolf2}
{cmd:(}{it:method1}
{depvar:1} 
{varlist:1} {ifin} [{it:{help weight:weight}}] [{cmd:,} {options1}]{cmd:)}{p_end}
{p 15 14 2}
{cmd:(}{it:method2}
{depvar:2} 
{varlist:2} {ifin} [{it:{help weight:weight}}] [{cmd:,} {options2}]{cmd:)}{p_end}
{p 15}{it:...}{p_end}
{p 15 14 2}
{cmd:(}{it:methodN}
{depvar:N} 
{varlist:N} {ifin} [{it:{help weight:weight}}] [{cmd:,} {optionsN}]{cmd:)}{cmd:,}
{p_end}
	       {it:indepvars(vars1, vars2, ..., varsN)} [{it:options}]



{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{cmd:indepvar(}{it:varlist}{cmd:)}}This is a required option. Each indepdendent variable should be indicated for which the multiple hypothesis correction is desired.
These refer to variables indicated in each of the independent variable {help varlist}s of models specified in the syntax.  Variables in each model should be separated by commas.
If a correction is desired for a single independent variable from each model this should be specified as (var1, var2, ..., varN).
It is possible to indicate various independent (treatment) variables which should be corrected across models and this quantity can vary for each model, for example
(var1a var1b, var2, ..., varNa VarNb varNc).
{p_end}
{...}
{synopt :{cmd:reps({help bootstrap:#})}}Perform # bootstrap replication; default is {cmd:reps(100)}.
Where possible prefer a considerably larger number of replications for more precise p-values. 
{p_end}
{...}
{synopt :{cmd:seed({help set seed:#})}}Sets seed to indicate the initial value for the pseudo-random number generator.  # can be any integer between 0 and 2^31-1. 
{p_end}
{...}
{synopt :{cmd:holm}}Along with standard output, additionally provide p-values corresponding to the
Holm multiple hypothesis correction. 
{p_end}
{...}
{synopt :{cmd:graph}}Requests that a graph be produced showing the Romano-Wolf null distribution
corresponding to each variable examined.
{p_end}
{...}
{synopt :{cmd:varlabels}}Name panels on the graph of null distributions using their variable labels
rather than their variable names.
{p_end}
{...}
{synopt :{cmd:onesided({help string})}} Indicates that p-values based on one-sided tests should be calculated.
Unless specified, p-values based on two-sided tests are provided, corresponding to the null that
each parameter is equal to 0 (or the values indicated in {cmd:nulls()}). In {cmd:onesided({help string})},
{help string} must be either "positive", in which case the null is that each parameter is greater
than or equal to 0, or "negative" in which case the null is that each parameter is less than or equal to 0.
{p_end}
{...}
{synopt :{cmd:nodots}}Suppress replication dots in bootstrap resamples.
{p_end}
{...}
{synopt :{cmd:noplusone}}Calculate the Resampled and Romano-Wolf adjusted p-values without
adding one to the numerator and denominator.
{p_end}
{...}
{synopt :{cmd:nulls({help numlist})}}Indicates the parameter values of interest used in
each test. If specified, a single scalar value should be indicated for each of the multiple
hypotheses tested, and these should be listed in the same order that variables are
listed as depvars in the command syntax. In the case that multiple {cmd:indepvars}
are specified, null parameters should be specified grouped first by {cmd:indepvars} and
then by {cmd:depvars}. For example, if two independent variables are considered with
four dependent variables, first the four null parameters associated with the first
independent variable should be listed, followed by the four null parameters associated
with the second independent variable. If this option is not used, it is assumed
that each null hypothesis is that the parameter is equal to 0.
{p_end}
{...}
{synopt :{cmd:strata({help varlist})}} specifies the variables identifying strata.  If {cmd:strata()} is specified, bootstrap samples are selected within each stratum when forming the resampled null distributions.
{p_end}
{...}
{synopt :{opth cl:uster(varlist)}} specifies the variables identifying resampling clusters.
If {cmd:cluster()} is specified, the sample drawn when forming the resampled null
distributions is a bootstrap sample of clusters. This option does not cluster standard errors
in each original regression.  If desired, this should be additionally specified using
{cmd:vce(cluster clustvar)}.  It is suggested that these options be used together to ensure that
underlying regression models and bootstrap resampling obey the same clustering schemes.  If
{cmd:vce(cluster clustvar)} is indicated, it is assumed that a clustered bootstrap resample is
desired, and {cmd: cluster()} will cluster on the same {cmd clustvar}.  If this is not desired,
the {cmd:regcluster()} option should be used, which allows for a cluster variable to be passed
only to the underlying regressions, or for different cluster variables to be used for the regression,
and the bootstrap resamples.
{p_end}
{...}
{synopt :{opth idcluster(newvar)}} In cases where clustered resampling is used, the idcluster
option should be used as in {help bsample} to specify a new name for resampled clusters, such
that clustering within each bootstrap is based on newvar rather than the original (repeated)
cluster variable.  If this option is used, any instances of cluster(var) or vce(cluster var)
will be replaced with cluster(newvar) vce(cluster newvar) within each bootstrap resample.
{p_end}
{synopt :{cmd:verbose}} Requests additional output, including display of the initial
(uncorrected) models estimated. XXXXX This will also result in the generation of a summary output
message indicating the number of hypotheses rejected in uncorrected models and
when implementing the Romano-Wolf correction, as well as any dependent variables
for which the null is rejected in the Romano-Wolf procedure.
{p_end}
{...}
{synopt :{cmd:usevalid}} In case of bootstrap replicates where invalid Studentized test
statistics are generated, these will be removed from the count of the total number of
bootstrap replicates, and the calculation of the final p-value.  In cases where such
invalid statistics are generated (for example due to no valid underlying estimate or standard error),
a warning will be issued by rwolf2 that the usevalid option should be specified.  This is highly
recommended.  In cases where bootstrapped models result in errors, the usevalid option
can be used to skip these replicates, and avoid errors with the program.
{p_end}
{...}
{synoptline}
{p2colreset}


{title:Description}

{p 6 6 2}
{hi:rwolf} calculates Romano and Wolf's (2005a,b) step-down adjusted p-values robust to
multiple hypothesis testing. This program follows the resampling algorithm described in
Romano and Wolf (2016), and provides a p-value corresponding to the significance of a
hypothesis test where S tests have been implemented, providing strong control of the
familywise error rate (the probability of committing any Type I error among all	
of the S hypotheses tested).  The {hi:rwolf} algorithm constructs a null
distribution for each of the S hypothesis tests based on Studentized bootstrap replications
of a subset of the tested variables.  Full details of the procedure are described in
Romano and Wolf (2016), and further discussion of this program and its implementation,
plus a full discussion of this ado, is provided in Clarke, Romano and Wolf (2019).

{p 6 6 2}
There are two ways for this command to be used. First, either {cmd:indepvar()}
and {cmd:method()} must be specified if the complete Romano-Wolf procedure should be
implemented including the estimation of bootstrap replications and generation of
adjusted p-values.  Alternatively, the user can provide rwolf with pre-computed
bootstrap or permuted replications of the estimated statistic and standard errors
for each of their multiple hypothesis tests of interest.  In this case, the {cmd:nobootstraps}
and {cmd:pointestimates(numlist)}, {cmd:stderrs(numlist)} and {cmd:stdests(varlist)}
should be indicated, and rwolf calculates the adjusted p-values from the replicates provided. 

{p 6 6 2}
In the former case where {hi:rwolf} takes care of estimating the {help bootstrap} replicates
of each test statistic and its standard error, {hi:rwolf} simply requires that the user
indicates the multiple dependent variables to be tested, the independent variable of
interest, and (optionally) a series of control variables which should be included in
each test.  {hi:rwolf} works with any {help regress:estimation-based regression command}
allowed in Stata, which should be indicated using the {cmd:method()} option. If not
specified, {help regress} is assumed.  In the case that {help ivregress} is specified,
it is assumed that the independent variable is the endogenous variable, and the
instrumental variable(s) should be indicated in the {cmd:iv()} option. If this is not
the case (ie if the treatment variable is an exogenous variable in the IV model), this
should be indicated with the {cmd:indepexog} option. Optionally, regression {help weight}s,
{help if}
or {help in} can be specified.  By default, 100 {help bootstrap} replications are run
for each of the S multiple hypotheses.  Where possible, a larger number of replications
should be preferred given that p-values are computed by comparing estimates to a
bootstrapped null distribution constructed from these replications.  The number of
replications is set using the {cmd:reps({help bootstrap:#})} option, and to replicate
results, the {cmd:seed({help seed:#})} should be set.

{p 6 6 2}
In the case of more complex situations where a user wishes to pre-compute their
test statistics, standard errors, and a large number of {help bootstrap} replicates
of each these, the user can request for only the p-value correction algorithm to
be implemented with the {cmd:bootstrap} option.  This allows for cases where different
estimation methodologies or different independent variables are used in each model
within the family of hypothesis tests, or where more complicated resampling procedures
are used, such as those based on permutation.  

{p 6 6 2}
By default, the re-sampled null distributions are formed using a simple bootstrap
procedure.  However, more complex stratified and/or clustered resampling procedures
can be specified using the {cmd:strata()} and {cmd:cluster()} options.  The
{cmd:cluster()} option refers only to the {help bsample:resampling} procedure, and
not to the standard errors estimated in each original regression model.  If the standard
variance estimator is not desired for regression models, this should be indicated
using the same {help regress:vce()} specification as in the original regression
models, for example {cmd:vce(cluster clustvar)}.  It is suggested that the
{cmd:cluster()} and {cmd:vce(cluster clustvar)} should be used together. If
only {cmd:vce(cluster clustvar)} is indicated, it is assumed that bootstrap resamples
should be conducted over the same {cmd:clustvar}.  If this is not the case, then
the {cmd:regcluster(clustvar)} option should be used, which controls clustering
only in the underlying regressions, or allows for clustering over different variables
in the regressions (using {cmd:clustvar()}) and resamples (using {cmd:cluster()}).

{p 6 6 2}
The command returns the Romano Wolf p-value corresponding to each variable, standard
(bootstrapped) uncorrected p-values, and for reference, the original uncorrected
(analytical) p-value from the initial tests when {hi:rwolf} estimates baseline
regression models.  {hi:rwolf} is an e-class command, and the Romano Wolf p-value for each
variable is returned as a scalar in e(rw_varname).  A matrix is also returned as
e(RW) providing the full set of Romano-Wolf corrected p-values.

{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Use the auto dataset to run multiple regressions of various independent variables on a single dependent variable of interest (weight) controlling for trunk and mpg.  {break}

{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. rwolf headroom turn price rep78, indepvar(weight) controls(trunk mpg) reps(250)}{p_end}

    {hline}

{pstd}Run the same analysis, however using areg to absorb a series of fixed effects {break}

{phang2}{cmd:. rwolf headroom turn price rep78, indepvar(weight) controls(trunk) reps(250) method(areg) abs(mpg)}{p_end}

    {hline}

{pstd}Run an instrumental variables model where the treatment variable (weight) is endogenous and a single instrument (length) is available {break}

{phang2}{cmd:. rwolf headroom turn price rep78, indepvar(weight) controls(trunk)  method(ivregress) iv(length)}{p_end}

{hline}

{pstd}Run multiple hypothesis tests using the National Longitudinal (panel) Survey with an xtreg, fe model.{p_end}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork}{p_end}
{phang2}{cmd:. rwolf wks_ue ln_wage hours tenure, indepvar(nev_mar) controls(i.year age) method(xtreg) seed(51) fe verbose}{p_end}

    {hline}



{marker results}{...}
{title:Saved results}

{pstd}
{cmd:rwolf} saves the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(rw_var1)}}The Romano Wolf p-value associated with variable 1 (var1 will be changed for variable name) {p_end}
{synopt:{cmd:...}} {p_end}

{synopt:{cmd:e(rw_varS)}}The Romano Wolf p-value associated with variable S.  Each of the dependent variables will be returned in this way. {p_end}

{synopt:{cmd:e(rw_depvar1_indepvar1)}} In the case that multiple independent variables are indicated, p-values for each
dependent variable--independent variable pair will be returned using both variables names. {p_end}
{synopt:{cmd:...}} {p_end}

{synoptset 25 tabbed}{...}
{p2col 5 20 24 2: Matrix}{p_end}
{synopt:{cmd:e(RW)}}The full set of Romano-Wolf corrected p-values, as well as the uncorrected p-values estimated by bootstrap and the baseline model (if relevant).{p_end}

{synopt:{cmd:e(RW_indepvar)}}In the case that multiple independent variables are indicated, the full set of Romano-Wolf corrected p-values, as well as the uncorrected p-values estimated by bootstrap and the baseline model (if relevant) are returned corresponding to each {cmd:indepvar}.{p_end}

{p2colreset}{...}


{marker acknowledgements}{...}
{title:Acknowledgements}

{p 6 6 2}
I am grateful to many users of rwolf for providing feedback and suggestions, and to David McKenzie for suggesting the alternative syntax implemented here.

	

{marker references}{...}
{title:References}

{marker RomanoWolf2005a}{...}
{phang}
Romano J.P. and Wolf M., 2005a.
{it:Exact and Approximate Stepdown Methods for Multiple Hypothesis Testing},
Journal of the American Statistical Association 100(469): 94-108.

{marker RomanoWolf2005b}{...}
{phang}
Romano J.P. and Wolf M., 2005b.
{it: Stepwise Multiple Testing as Formalized Data Snooping},
Econometrica 73(4): 1237-1282.

{marker RomanoWolf2016}{...}
{phang}
Romano J.P. and Wolf M., 2016.
{it: Efficient computation of adjusted p-values for resampling-based stepdown multiple testing},
Statistics and Probability Letters 113: 38-40.

{marker Clarketal2019}{...}
{phang}
Clarke, D, Romano J.P. and Wolf M., 2019.
{it: The Romano-Wolf Multiple Hypothesis Correction in Stata}, Forthcoming, Stata Journal.
{p_end}


{title:Author}

{pstd}
Damian Clarke, Department of Economics, University of Chile. {browse "mailto:dclarke@fen.uchile.cl":dclarke@fen.uchile.cl}
{p_end}


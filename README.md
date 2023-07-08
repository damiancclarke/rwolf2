![Stata](https://img.shields.io/badge/stata-2013-green) ![GitHub Starts](https://img.shields.io/github/stars/damiancclarke/rwolf2?style=social) ![GitHub license](https://img.shields.io/github/license/damiancclarke/rwolf2)

# rwolf2 - A more flexible syntax to calculate Romano-Wolf stepdown p-values for multiple hypothesis testing 

`rwolf2` calculates [Romano and Wolf](#references)'s ([2005a](#references),[b](#references))
step-down adjusted *p*-values robust to multiple hypothesis testing.  It provides a more general syntax than that
provided in the [rwolf](https://github.com/damiancclarke/rwolf) command, although the underlying algorithm is the same.
This program follows the resampling algorithm described in Romano and Wolf (2016), and provides a *p*-value corresponding
to the significance of a hypothesis test where *S* tests have been implemented, providing strong control of the
familywise error rate (the
probability of committing any Type I error among all of the *S* hypotheses tested).  The `rwolf2` algorithm constructs a null
distribution for each of the *S* hypothesis tests based on Studentized bootstrap replications of a subset of the tested
variables.  Full details of the procedure are described in [Romano and Wolf (2016)](#references), additional discussion related to the
procedure in Stata is provided in [Clarke, Romano and Wolf (2019)](#references).

This command follows a syntax similar to `sureg` where multiple models of interest should each be specified in a separate set of
parentheses.  These models will be estimated, and bootstrap replicates will be generated, and based on these estimates and
bootstrap replicates, the Romano-Wolf correction will be implemented.  The variables of interest (for which the multiple
hypothesis correction is desired) should be indicated in the required option `indepvars()`.  This syntax can be used to implement
any Stata estimation command that returns a parameter vector and variance-covariance matrix, as well as other non-native
commands including user-written ados such as `reghdfe`, `ivreg2` (and related programs) and `rdrobust`.

In each model, optionally, regression weights, if or in can be specified.  By default, 100 bootstrap replications are run for
each of the *S* multiple hypotheses.  Where possible, a much larger number of replications should be preferred given that
*p*-values are computed by comparing estimates to a bootstrapped null distribution constructed from these replications.  The
number of replications is set using the `reps(#)` option, and to replicate results, the `seed(#)` should be set.

By default, the re-sampled null distributions are formed using a simple bootstrap procedure.  However, more complex stratified
and/or clustered resampling procedures can be specified using the `strata()` and `cluster()` options.  The `cluster()` option refers
only to the resampling procedure, and not to the standard errors estimated in each original regression model.  If the standard
variance estimator is not desired for regression models, this should be indicated in the syntax of each method indicated in the
command syntax.

The command returns the Romano Wolf *p*-value corresponding to each variable, standard (bootstrapped) uncorrected *p*-values, and
for reference, the original uncorrected (analytical) *p*-value from the initial tests when rwolf2 estimates baseline regression
models.  rwolf2 is an e-class command, and a matrix is returned as e(RW) providing the full set of Romano-Wolf corrected
*p*-values (and other *p*-values mentioned above).

To install directly into Stata:
```s
ssc install rwolf2, replace
```
or using ```net install``` command:
```s
net install rwolf2, from("https://raw.githubusercontent.com/damiancclarke/rwolf2/master") replace
```
## Syntax
```s
rwolf2 (method1 depvar1 varlist1 [if] [in] [weight] [, options1])
       (method2 depvar2 varlist2 [if] [in] [weight] [, options2])
        ...
       (methodN depvarN varlistN [if] [in] [weight] [, optionsN]),
       indepvars(vars1, vars2, ..., varsN) [options]
```
+ indepvars(varlist):       This is a required option. Each indepdendent variable should be indicated for which the multiple hypothesis correction is desired.  These refer to variables indicated in each of the independent variable varlists of models
                           specified in the syntax.  Variables in each model should be separated by commas.  If a correction is desired for a single independent variable from each model this should be specified as (var1, var2, ..., varN).  It is
                           possible to indicate various independent (treatment) variables which should be corrected across models and this quantity can vary for each model, for example (var1a var1b, var2, ..., varNa VarNb varNc).
+ reps(*#*):                  Perform # bootstrap replication; default is reps(100).  Where possible prefer a considerably larger number of replications for more precise p-values.
+ seed(*#*):                  Sets seed to indicate the initial value for the pseudo-random number generator.  # can be any integer between 0 and 2^31-1.
+ holm:                     Along with standard output, additionally provide p-values corresponding to the Holm multiple hypothesis correction.
+ graph:                    Requests that a graph be produced showing the Romano-Wolf null distribution corresponding to each hypothesis test considered.
+ varlabels:                Name panels on the graph of null distributions using their variable labels rather than their variable names.
+ onesided(string):         Indicates that p-values based on one-sided tests should be calculated.  Unless specified, p-values based on two-sided tests are provided, corresponding to the null that each parameter is equal to 0 (or the values indicated
                           in nulls()). In onesided(string), string must be either "positive", in which case the null is that each parameter is greater than or equal to 0, or "negative" in which case the null is that each parameter is less than or
                           equal to 0.
+ nodots:                   Suppress replication dots in bootstrap resamples.
+ noplusone:                Calculate the Resampled and Romano-Wolf adjusted p-values without adding one to the numerator and denominator.
+ nulls(numlist):           Indicates the parameter values of interest used in each test. If specified, a single scalar value should be indicated for each of the multiple hypotheses tested, and these should be listed in the same order that variables
                           are listed as indepvars in the command syntax. If this option is not used, it is assumed that each null hypothesis is that the parameter is equal to 0.
+ strata(varlist):          Specifies the variables identifying strata.  If strata() is specified, bootstrap samples are selected within each stratum when forming the resampled null distributions.
+ cluster(varlist):         Specifies the variables identifying resampling clusters.  If cluster() is specified, the sample drawn when forming the resampled null distributions is a bootstrap sample of clusters. This option does not cluster standard
                           errors in each original regression.  If desired, this should be additionally specified as an option in the individual regression models.
+ idcluster(newvar):        In cases where clustered resampling is used, the idcluster option should be used as used in bsample to specify a new name for resampled clusters, such that clustering within each bootstrap is based on newvar rather than
                           the original (repeated) cluster variable.  If this option is used, any instances of cluster(var) or vce(cluster var) will be replaced with cluster(newvar) or vce(cluster newvar) within each bootstrap resample.
+ verbose:                  Requests additional output, including display of the initial (uncorrected) models estimated.
+ usevalid:                 In cases of bootstrap replicates where invalid Studentized test statistics are generated, these will be removed from the multiple hypothesis correction procedure.  In cases where such invalid statistics are generated (for
                           example due to no valid underlying estimate or standard error), a warning will be issued by rwolf2 that the usevalid option should be specified.  This is highly recommended.  In cases where bootstrapped models result in
                           errors causing the command to exit with error, the usevalid option can be used to skip these replicates, and avoid errors with the program.

## Running an example
Below is an example using the system `auto` dataset, estimating two models with three hypothesis tests of interest.
```s
sysuse auto,clear
(1978 automobile data)

rwolf2 (ivregress 2sls price (weight = length)) (regress price weight trunk), indepvars(weight, weight trunk) 
```
The code returns the following results
```s
Bootstrap replications (100). This may take some time.
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5
..................................................     50
..................................................     100


Romano-Wolf step-down adjusted p-values
Number of resamples: 100


------------------------------------------------------------------------------
                | Model p-value    Resample p-value    Romano-Wolf p-value
----------------+-------------------------------------------------------------
price           |     
         weight |     0.0000             0.0099              0.0099
------------------------------------------------------------------------------
price           |     
         weight |     0.0000             0.0099              0.0099
          trunk |     0.5200             0.5644              0.5644
------------------------------------------------------------------------------
```

## Other Links
+ The `rwolf2` command has been written to provide additional functionality and a more flexible syntax than the original `rwolf` command [(Clarke, 2016)](#references). A discussion of the difference between rwolf2 and rwolf as well as documentation of their identical outcomes in cases where identical corrections are implemented is avaialable [here](https://www.damianclarke.net/computation/rwolf2.pdf). 
+ David McKenzie has provided a discussion of multiple hypothesis testing routines in Stata, with some details related to this and other commands.  This is available [here](https://blogs.worldbank.org/impactevaluations/updated-overview-multiple-hypothesis-testing-commands-stata).
+ The Romano-Wolf multiple hypothesis correction is a bootstrap-based resampling correction, and for comparability with regression-based *p*-values, original bootstrap procedures (labelled as Resample p-value in command output) should be very close to regression-based *p*-values (Model p-value) in command output.  If these p-values are signficantly different, alternative bootstrap procedures may be considered (for example block or stratified bootstrapping, as permitted via the `cluster`/`idcluster` or `strata` options).  Some discussion related to this point can be found in Appendix B1 of the following paper by [Leighton (2023)](#https://www.sciencedirect.com/science/article/pii/S0305750X23001535) (open source link). 


## References
D. Clarke. RWOLF: Stata module to calculate Romano-Wolf stepdown pvalues for multiple hypothesis testing. Statistical Software Components, Boston College Department of Economics, Dec. 2016. URL https://ideas.repec.org/c/boc/bocode/s458276.html.

D. Clarke, J. P. Romano, and M. Wolf. **[The Romano–Wolf Multiple­-hypothesis Correction in Stata](https://journals.sagepub.com/doi/abs/10.1177/1536867X20976314)**. *The Stata Journal*, 20(4):812–843, 2020.

M. Leighton, A Martine and J Massaga. **[Fostering early childhood development in low-resource communities: Evidence from a group-based parenting intervention in Tanzania](https://www.sciencedirect.com/science/article/pii/S0305750X23001535)**. *World Development*, 170: 106335, 2023. 

J. P. Romano and M. Wolf. **[Stepwise Multiple Testing as Formalized Data Snooping](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1468-0262.2005.00615.x)**. Econometrica, 73(4): 1237–1282, 2005a.

J. P. Romano and M. Wolf. **[Exact and Approximate Stepdown Methods for Multiple Hypothesis Testing](https://www.tandfonline.com/doi/abs/10.1198/016214504000000539)**. Journal of the American Statistical Association, 100(469):94–108, 2005b.

J. P. Romano and M. Wolf. **[Efficient computation of adjusted p-­values for resampling-­based stepdown multiple testing](https://www.sciencedirect.com/science/article/abs/pii/S0167715216000389)**. *Statistics and Probability Letters*, 113:38–40, 2016.


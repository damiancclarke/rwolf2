*! rwolf2: Romano Wolf stepdown hypothesis testing algorithm, alternative syntax
*! Version 0.0.1 april 20, 2021 @ 22:47:41
*! Author: Damian Clarke
*! Department of Economics
*! University of Chile
*! dclarke@fen.uchile.cl

/*
version highlights:
0.0.1:[03/07/2021]: 
*/

cap program drop rwolf2
program rwolf2, eclass
vers 11.0

*-------------------------------------------------------------------------------
*--- (0) Parsing of equations 
*-------------------------------------------------------------------------------    
local stop 0
local i = 0
while !`stop' {
    gettoken eqn 0 : 0,  parse(" ,[") match(paren)    
    if `"`paren'"'=="(" {
        local ++i
        local eqn`i' `eqn'
    }
    else {
        local stop 1
    }
}
local 0 ", `0'"
local numeqns = `i'

*-------------------------------------------------------------------------------
*--- (1) Syntax
*-------------------------------------------------------------------------------
#delimit ;
syntax , indepvars(string)
[
  reps(integer 100)
  Verbose
  seed(numlist integer >0 max=1)
  holm
  graph
  onesided(name)
  varlabels
  nodots
  holm
  nulls(numlist)
  strata(varlist)
  CLuster(varlist)
];
#delimit cr

*-------------------------------------------------------------------------------
*--- (2) Unpacking baseline options
*-------------------------------------------------------------------------------
if length(`"`verbose'"')==0 local q qui
cap set seed `seed'

if length(`"`nulls'"')!=0 {
    local j = 0
    foreach null of numlist `nulls' {
        local ++j
        local beta0`j' = `null'
    }
    local nnulls = `j'
}

tempname nullvals
tempfile nullfile
file open `nullvals' using "`nullfile'", write all


local diswidth = 15
foreach equation of numlist 1(1)`numeqns' {
    tokenize `eqn`equation''
    ds `2'
    local y`equation' `r(varlist)'
    
    local y`equation'length = length("`y`equation''")
    local diswidth = max(`diswidth',`y`equation'length')
}

*-------------------------------------------------------------------------------
*--- (3) Parsing of independent variables and estimating models
*-------------------------------------------------------------------------------
dis `"`indepvars'"'
tokenize `"`indepvars'"', parse(",")

local inft = 0
local j    = 0

foreach equation of numlist 1(1)`numeqns' {
    `q' dis "Equation `equation' is `eqn`equation''"
    `q' dis "Variables to correct are ``equation''"
    `eqn`equation''
    local xvar`equation'
    foreach var of varlist ``equation'' {
        dis "`var'"
        local xvar`equation' `xvar`equation'' `var'

        local ++j
        if length(`"`nulls'"')==0 local beta0`j' = 0

        local t`j' = abs((_b[`var']-`beta0`j'')/_se[`var'])
        local inft = min(`inft',`t`j'')

        local beta`j' = _b[`var']
        test _b[`var'] = `beta0`j''
        local pv`j' = string(r(p), "%6.4f")
        local pv`j's= r(p)
        dis `pv`j's'

        local cand `cand' `j'
        *file write `nullvals' "b`j'_`var'; se`j'_`var';"
        file write `nullvals' "b`j'; se`j';"
    }
    macro shift

}



*-------------------------------------------------------------------------------
*--- (4) Bootstrap replications
*-------------------------------------------------------------------------------
dis "Bootstrap replications (`reps'). This may take some time."
if length(`"`dots'"')==0 {
    dis "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5"
}
forvalues i=1/`reps' {
    local bmess1 "There has been an issue with a bootstrap replicate"
    local bmess2 "in bootstrap `i'.  Taking next bootstrap sample..."
    local j=0
    preserve
    bsample 
    if length(`"`dots'"')==0 {
        display in smcl "." _continue
        if mod(`i',50)==0 dis "     `i'"
    }

    foreach equation of numlist 1(1)`numeqns' {
        qui `eqn`equation''

        local k=1
        foreach var of varlist `xvar`equation'' {
            local ++j
            if `j'==1&`k'==1 file write `nullvals' _n "`= _b[`var']';`= _se[`var']'"
            else file write `nullvals' ";`= _b[`var']';`= _se[`var']'"
            local ++k
        }
    }
    restore
}
    
preserve
file close `nullvals'
qui insheet using `nullfile', delim(";") names clear case

matrix pvalues = J(`j',3,.)
if length(`"`holm'"')!=0 matrix pvalues = J(`j',4,.)


*-------------------------------------------------------------------------------
*--- (5) Create null t-distribution
*-------------------------------------------------------------------------------
foreach num of numlist 1(1)`j' {    
    qui gen     t`num'=(b`num'-`beta`num'')/se`num'
    qui replace b`num'=abs((b`num'-`beta`num'')/se`num')
}
dis ""
dis "`cand'"

*-------------------------------------------------------------------------------
*--- (6) Create stepdown value in descending order based on t-stats
*-------------------------------------------------------------------------------
local maxt = `inft'-10
local pval = 0
local rank
local Holm = `j'
local prmsm1


**local cand_`var' `cand_`var'' `j'
*******tokenize `varlist'
local ii=0
while length("`cand'")!=0 {
    local ++ii
    local donor_tvals
    
    if length(`"`onesided'"')==0|`"`onesided'"'=="negative" {
        foreach var of local cand {
            if `t`var''>`maxt' {
                dis "D"
                local maxt = `t`var''
                local maxv `var'
                if length(`"`onesided'"')==0  local ovar  b`var'
                if `"`onesided'"'=="negative" local ovar  t`var'
            }
            dis "Maximum t among remaining candidates is `maxt' (variable `maxv')"
            if length(`"`onesided'"')==0  {
                local donor_tvals `donor_tvals' b`var'
            }
            if `"`onesided'"'=="negative" {
                local donor_tvals `donor_tvals' t`var'
            }
        }
        qui egen empiricalDist = rowmax(`donor_tvals')    
        qui count if empiricalDist>=`maxt'  & empiricalDist != .
        local cnum = r(N)
        if length(`"`plusone'"')!=0      local pval = (`cnum')/(`reps')
        else if length(`"`plusone'"')==0 local pval = (`cnum'+1)/(`reps'+1)
    
    
        ***24/02/2020: CHANGED BELOW LINE
        *qui count if `ovar'>=`maxt' & `maxt'!=.
        qui count if `ovar'>=`maxt' & `ovar'!=.
        local cnum = r(N)
        if length(`"`plusone'"')!=0      local pvalBS = `cnum'/`reps'
        else if length(`"`plusone'"')==0 local pvalBS = (`cnum'+1)/(`reps'+1)
    
        local pbs`maxv's= `pvalBS'
        local prm`maxv's= `pval'
        local prh`maxv's= min(`pvalBS'*`Holm',1)
        
        if length(`"`prmsm1'"')!=0 {
            local prm`maxv's=max(`prm`maxv's',`prmsm1')
        }
    
        local prmsm1 = `prm`maxv's'
        if length(`"`graph'"')!=0 {
            gen tDist`ii'=empiricalDist
            local tt`ii' = `t`maxv''
            local ll`ii' "`label`maxv''"
            local yy`ii' ``maxv''
            local pp`ii' = string(`prm`maxv's',"%6.4f")
        }         
    }
    drop empiricalDist 
    local rank `rank' `maxv' `minv'
    local candnew
    foreach c of local cand {
        local match = 0
        foreach r of local rank {
            if `r'==`c' local match = 1
        }
        if `match'==0 local candnew `candnew' `c'
    }
    local cand `candnew'
    local maxt = `inft'-10
    local maxv = 0
    local --Holm    
}



*-------------------------------------------------------------------------------
*--- (7) Graphing null distributions
*-------------------------------------------------------------------------------
if length(`"`graph'"')!= 0 {
    local graphs
    if length(`"`onesided'"')==0 {
        foreach num of numlist 1(1)`ii' {
            local title "Variable `yy`num''"
            if length(`"`varlabels'"')!=0 local title "`ll`num''"
            
            #delimit ;
            twoway hist tDist`num', bcolor(gs12) ||
            function y=sqrt(2)/sqrt(c(pi))*exp(-x^2/2), range(0 6)
            xline(`tt`num'', lcolor(red)) name(g`num', replace)
            scheme(sj) nodraw legend(off) title(`title')
            note("p-value = `pp`num''");
            #delimit cr
            local graphs `graphs' g`num'
        }
    }
    else {
        foreach num of numlist 1(1)`ii' {
            local title "Variable `yy`num''" 
            if length(`"`varlabels'"')!=0 local title "`ll`num''"
            #delimit ;
            twoway hist tDist`num', bcolor(gs12) ||
            function y=1/sqrt(2*c(pi))*exp(-x^2/2), range(-4 4)
            xline(`tt`num'', lcolor(red)) name(g`num', replace)
            scheme(sj) nodraw legend(off) title(`title')
            note("p-value = `pp`num''");
            #delimit cr
            local graphs `graphs' g`num'
        }
    }
    graph combine `graphs', scheme(sj)
}
restore


*-------------------------------------------------------------------------------
*--- (8) Report and export p-values
*-------------------------------------------------------------------------------
ereturn clear
local vardisplay
local linelength=0

foreach equation of numlist 1(1)`numeqns' {
    foreach var of varlist `xvar`equation'' {
        local linelength=`linelength'+length("`var'")
        if `linelength'>57 {
            local vardisplay "`vardisplay'" _newline _col(22) "`var'"
            local linelength=length("`var'")
        }
        else local vardisplay "`vardisplay' `var'"
    }
}
local crit = (100-c(level))/100
local lev  = c(level)

dis _newline
dis _newline
dis "Romano-Wolf step-down adjusted p-values"
dis "Number of resamples: `reps'"
dis _newline
dis "{hline 78}"
dis "{dup `diswidth': } | Model p-value    Resample p-value    Romano-Wolf p-value"
dis "{hline `=`diswidth'+1'}+{hline `=76-`diswidth''}"

local j = 0
foreach equation of numlist 1(1)`numeqns' {
    local DUP = `diswidth' - `y`equation'length'
    display as text "`y`equation'' {dup `DUP': }{c |}     "
    if length(`"`holm'"')==0 {
        foreach var of varlist `xvar`equation'' {
            local ++j
            display as text %`diswidth's abbrev("`var'",19) " {c |}     "   /*
            */  as result %6.4f `pv`j's' "             "   /*
            */  as result %6.4f `pbs`j's' "              "   /*
            */  as result %6.4f `prm`j's'
            ereturn scalar rw_`var'=`prm`j's'
            matrix pvalues[`j',1]=`pv`j's'
            matrix pvalues[`j',2]=`pbs`j's'
            matrix pvalues[`j',3]=`prm`j's'
            
        }
    }    
    dis "{hline 78}"
}
dis _newline
ereturn matrix RW=pvalues

end

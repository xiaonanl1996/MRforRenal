library(stringr)
options("scipen"=100)

# Run chunks of R code in arbitrary working directory
# From Hadley Wickham
# https://github.com/yihui/knitr/issues/38
in_dir <- function(dir, code) {
  cur <- getwd()
  setwd(dir)
  on.exit(setwd(cur))
  
  force(code)
}

pretty_dp <- function(x, dp=0, pct=FALSE, comma=FALSE){
  if(pct){x <- 100*x}
  if(comma){
    format(round(x, dp), digits=dp, nsmall=dp, big.mark=",") %>% trimws
  } else {
    format(round(x, dp), digits=dp, nsmall=dp) %>% trimws
  }
}

pretty_confint <- function(lci, uci, dp, pct=FALSE){
  paste0("(", pretty_dp(x=lci, dp=dp, pct=pct), ", ", pretty_dp(x=uci, dp=dp, pct=pct), ")")
}

pretty_pval <- function(p, cutoff=0.001, string="<0.001", dp=3){
  ifelse(p<cutoff, string, pretty_dp(p, dp))
}

lower <- function(x){
  paste0(tolower(substring(x, 1,1)), substring(x, 2))
}

upper <- function(x){
  paste0(toupper(substring(x, 1,1)), substring(x, 2))
}

prettyfunc <- function(x, pnames=list(), upper=FALSE, bold=FALSE, flist=c()){
  out <- x
  if(x %in% names(pnames)){
    out <- pnames[[x]]
    if(upper){
      out <- upper(out)
    }
    if(bold){
      out <- paste0("**", out, "**")
    }
  }
  if(x %in% flist){
    if(!exists("footnote_no")){
      footnote_no <<- 1
    }
    out <- paste0(out, "^", footnote_no, "^")
    footnote_no <<- footnote_no + 1
  }
  return(out)
}

diagcollist <- function(colstring, sep="", ncols) {
  x <- 0:ncols
  colstr <- paste0("`",paste0(colstring, sep, x, collapse="`, `"),"`")
  return(colstr)
}

table1_standardised <- function(data, varlist, adj, stratify, strata=NULL, 
                                pretty_names=list(), singlecol=TRUE, dp=1, show_crude=FALSE){
  #' Create age-standardised Table 1
  #'
  #' @param data The data
  #' @param varlist The list of variables to include in the table
  #' @param adj The variable to adjust by, as a factor. Eg age in single-year brackets.
  #' @param stratify The variable to stratify by, as a factor
  #' @param strata If you don't want to use all levels of the stratifying variable, specify desired levels here
  #' @param pretty_names List of human readable names corresponding to variable names
  #' @param singlecol Whether to stack variable names and levels into one column (TRUE - default) or spread across two columns (FALSE)
  #' @param dp Number of decimal places to display in the table (default 1)
  #' @param show_crude Include crude proportions in the form "crude (adjusted)". Advisable for sanity check but not for final presentation.
  #' 
  #' @return A dataframe formatted appropriately to output as Table 1
  #' @export
  #'
  #' @examples
  #' 
  if(is.null(strata)) { strata <- levels(data[[stratify]]) }
  
  table <- c()
  colnames <- c()
  
  for(s in strata) {
    colnames <- c(colnames, paste0(s, " (N=", nrow(data[data[[stratify]]==s,]), ")"))
    col <- c()
    
    for(var in varlist){
      if(is.factor(data[[var]])){
        if(singlecol){ col <- c(col, "") }
        for(l in levels(data[[var]])){
          
          count <- table(data[[adj]][data[[var]]==l & data[[stratify]]==s], useNA='ifany')
          pop <- table(data[[adj]][data[[stratify]]==s], useNA='ifany')
          stdpop <- table(data[[adj]], useNA='ifany')
          
          proportions <- ageadjust.direct(count=count, pop=pop, stdpop=stdpop)
          crude.prop <- pretty_dp(100*proportions[["crude.rate"]],dp)
          adj.prop <- pretty_dp(100*proportions[["adj.rate"]],dp)
          
          if(show_crude){
            col <- c(col, paste0(crude.prop, " (", adj.prop, ")"))
          } else {
            col <- c(col, adj.prop)
          }
        }
      } else {
        col <- c(col, paste0(pretty_dp(mean(data[[var]][data[[stratify]]==s]), dp), 
                             " (", pretty_dp(sd(data[[var]][data[[stratify]]==s]), dp), ")"))
      }
    }
    table <- cbind(table, col)
  }
  coefflist <- list()
  # Prepare the list of coefficients - variables and levels for factors or blanks for continuous
  for(var in varlist){
    coefflist[[var]] <- preparecoefflist(df=data, varname=var, pretty_names=pretty_names, onecol=singlecol)
  }
  coeffnames <- do.call(rbind, coefflist)
  
  colnames(table) <- colnames
  table <- cbind(coeffnames%>%select(-IDcol), table)
  rownames(table) <- NULL
  
  return(table)
}

# Normal Table 1
# Print numbers and proportions for factors, median and IQR or mean and 95% CI for continuous variables
# Optionally provide p-values from chi-squared (categorical) and t-test (continuous)
descriptivetable <- function(df, varlist, contavg='mean', assocvar=NULL, singlecol=FALSE, 
                             pretty_names=list(), footnote_list=c()){
  if(!exists("footnote_no")){footnote_no <<- 1} # Note use of <<- instead of <- to assign this globally
  outtable <- c()
  for(var in varlist){
    if(is.factor(df[[var]])){ # Categorical variables (factors) need a row per level, n's and %'s
      n <- table(df[[var]], useNA='ifany')
      pct <- pretty_dp(prop.table(n), dp=1, pct=TRUE)
      npct <- paste0(n, " (", pct, "%)")
      variable <- c(prettyfunc(var, pnames=pretty_names, upper=TRUE, flist=footnote_list))
      levels <- names(n)
      if(!is.null(assocvar)){
        tab <- table(df[[assocvar]], df[[var]])
        chi <- chisq.test(tab)
        pval <- c(ifelse(chi$p.value<0.001, "<0.001", round(chi$p.value,3)))
        outtable <- rbind(outtable, 
                          c(var, paste0("**", variable, "**"), "", "", "", pval), 
                          cbind(paste0(var, levels), levels, n, pct, npct, ""))
      } else{
        outtable<- rbind(outtable, 
                         c(var, paste0("**", variable, "**"), "", "", ""), 
                         cbind(paste0(var, levels), levels, n, pct, npct))
      }
    } else { # Continuous variables need the mean (and SD) or median (and IQR)
      if(contavg=="mean"){
        n <- pretty_dp(mean(df[[var]], na.rm=TRUE), dp=1, comma=TRUE)
        pct <- pretty_dp(sd(df[[var]], na.rm=TRUE), dp=1, comma=TRUE)
        npct <- paste0(n, " (", pct, ")")
        variable <- paste0("Mean ", prettyfunc(var, pretty_names, upper=FALSE, flist=footnote_list), " (SD)")
      } else if (contavg=="median"){
        n <- pretty_dp(median(df[[var]], na.rm=TRUE), dp=1, comma=TRUE)
        IQR <- pretty_dp(quantile(df[[var]], na.rm=TRUE), dp=1, comma=TRUE)
        pct <- paste0("(", IQR[2], "-", IQR[4], ")")
        npct <- paste0(n, " ", pct)
        variable <- paste0("Median ", prettyfunc(var, pnames=pretty_names, upper=FALSE, flist=footnote_list), " (IQR)")
      } else if(contavg=="n"){
        n <- nrow(df[!is.na(df[[var]]),])
        pct <- NA
        npct <- NA
        variable <- prettyfunc(var, pnames=pretty_names, upper=TRUE, flist=footnote_list)
      }
      if(!is.null(assocvar)){
        tt <- t.test(df[[var]][df[[assocvar]]==TRUE], df[[var]][df[[assocvar]]==FALSE])
        p <- pretty_pval(tt$p.value)
        outtable <- rbind(outtable, cbind(var, paste0("**", variable, "**"), n, pct, npct, p))
      } else {
        outtable<- rbind(outtable, cbind(var, paste0("**", variable, "**"), n, pct, npct))
      }
    }
  }
  rownames(outtable) <- c()

  outdf <- as.data.frame(outtable, stringsAsFactors=FALSE)
  if(singlecol){
    outdf <- outdf %>% select(-c(n, pct))
  } else {
    outdf <- outdf %>% select(-npct)
  }
  return(outdf)
}

printMIresults <- function(df, varlist, modeloutput, pretty_names=list(), onecol=FALSE, IDcol=FALSE){
  require(dplyr)
  
  coefflist <- list()
  # Prepare the list of coefficients - variables and levels for factors or blanks for continuous
  for(var in varlist){
    coefflist[[var]] <- preparecoefflist(df=df, varname=var, pretty_names=pretty_names, onecol=onecol)
  }
  coeffnames <- do.call(rbind, coefflist)
  
  regression <- data.frame(
    IDcol=modeloutput$term,
    HR=pretty_dp(exp(modeloutput$estimate),dp=2),
    CI=pretty_confint(exp(modeloutput$estimate-1.96*modeloutput$std.error),
                            exp(modeloutput$estimate+1.96*modeloutput$std.error),
                            dp=2),
    p=pretty_pval(modeloutput$p.value),
    stringsAsFactors=FALSE
  )
  
  results <- left_join(coeffnames, regression, by="IDcol")
  if(onecol){
    results$HR[is.na(results$HR) & (results$IDcol != results$Coefficient & !is.na(results$Coefficient))] <- "1"
  } else {
    results$HR[is.na(results$HR) & !is.na(results$Coefficient)] <- "1"
  }
  
  coeffcols <- colnames(coeffnames)
  if(IDcol==FALSE){
    coeffcols <- coeffcols[coeffcols != "IDcol"]
  }
  results <- results[,c(coeffcols, "HR", "CI", "p")]
  names(results) <- c(coeffcols, "HR", "95% CI", "p")
  
  
  rownames(results) <- NULL
  # https://www.r-bloggers.com/regression-on-categorical-variables/
  return(results)
}

# Prettyprint the results from a Cox model
# To use this, 
# model <- coxph(Surv(time_to_dementia, dementia_status) ~ age, data=data)
# kable(printcoxresults(model), caption="")
printresults <- function(df, varlist, modeloutput, pretty_names=list(), onecol=FALSE, IDcol=FALSE,forplot=FALSE){
  require(dplyr)
  
  coefflist <- list()
  # Prepare the list of coefficients - variables and levels for factors or blanks for continuous
  for(var in varlist){
    coefflist[[var]] <- preparecoefflist(df=df, varname=var, pretty_names=pretty_names, onecol=onecol)
  }
  coeffnames <- do.call(rbind, coefflist)
  
  if(class(modeloutput)[1]=="coxph"){
    summ <- summary(modeloutput)
    coeff <- summ$coefficients
    conf <- summ$conf.int
    
    regression <- data.frame(
      IDcol=rownames(coeff),
      HR_num=coeff[,2],
      LCI=conf[,3],
      UCI=conf[,4],
      HR=pretty_dp(coeff[,2], dp=2), # HR
      CI=pretty_confint(conf[,3], conf[,4], dp=2), # 95% CI
      p=pretty_pval(coeff[,5]), # p-value
      stringsAsFactors=FALSE
    )
  }else if (class(modeloutput)[1]=="glm"){
    summ <- summary(modeloutput)
    coeff <- summ$coefficients
    # NOMVAR <- rownames(coeff)
    regression <- data.frame(
      IDcol=(rownames(coeff)),
      HR=pretty_dp(exp(coeff[,1]), dp=2), # OR
      HR_num=exp(coeff[,1]),
      LCI=exp(coeff[,1]-(1.96*coeff[,2])),
      UCI=exp(coeff[,1]+(1.96*coeff[,2])),
      CI=pretty_confint(exp(coeff[,1]-(1.96*coeff[,2])),
                        exp(coeff[,1]+(1.96*coeff[,2])),
                        dp=2), # 95% CI
      p=pretty_pval(coeff[,4]), # p-value
      stringsAsFactors=FALSE
    )
  }
  
  results <- left_join(coeffnames, regression, by="IDcol")
  if(onecol){
    results$HR[is.na(results$HR) & (results$IDcol != results$Coefficient & !is.na(results$Coefficient))] <- "1"
    results$HR_num[is.na(results$HR_num) & (results$IDcol != results$Coefficient & !is.na(results$Coefficient))] <- 1
  } else {
    results$HR[is.na(results$HR) & !is.na(results$Coefficient)] <- "1"
    results$HR_num[is.na(results$HR_num) & !is.na(results$Coefficient)] <- 1
  }
  
  coeffcols <- colnames(coeffnames)
  if(IDcol==FALSE){
    coeffcols <- coeffcols[coeffcols != "IDcol"]
  }
  
  if(forplot==FALSE){
    # If not for plot
    results <- results[,c(coeffcols, "HR", "CI", "p")]
    
    if(class(modeloutput)[1]=="coxph"){
      
      names(results) <- c(coeffcols, "HR", "95% CI", "p")
      
    }else if (class(modeloutput)[1]=="glm"){
      names(results) <- c(coeffcols, "OR", "95% CI", "p")
    }
  } else {
    # If for plot
    results <- results[,c(coeffcols, "HR_num", "LCI", "UCI")]
    
    if(class(modeloutput)[1]=="coxph"){
      
      names(results) <- c(coeffcols, "HR_num", "LCI", "UCI")
      
    }else if (class(modeloutput)[1]=="glm"){
      names(results) <- c(coeffcols, "OR_num", "LCI", "UCI")
    }
  }
  
  
  rownames(results) <- NULL
  # https://www.r-bloggers.com/regression-on-categorical-variables/
  return(results)
}

# Pretty print the results from a logistic regression model
printlogresults <- function(model, coeffnames=NULL, IDcol=FALSE){
  summ <- summary(model)
  coeff <- summ$coefficients
  # NOMVAR <- rownames(coeff)
  regression <- data.frame(
    IDcol=(rownames(coeff)),
    OR=pretty_dp(exp(coeff[,1]), dp=2), # OR
    CI=pretty_confint(exp(coeff[,1]-(1.96*coeff[,2])),
                      exp(coeff[,1]+(1.96*coeff[,2])),
                      dp=2), # 95% CI
    p=pretty_pval(coeff[,4]), # p-value
    stringsAsFactors=FALSE
  )
  if(is.null(coeffnames)){
    results <- regression
    names(results) <- c("Coefficient", "OR", "95% CI", "p")
  } else {
    results <- merge(coeffnames, regression, all.x=TRUE)
    results$OR[is.na(results$OR)] <- "1"
    results <- results[match(coeffnames$IDcol, results$IDcol),]
    
    coeffcols <- colnames(coeffnames)
    if(IDcol==FALSE){
      coeffcols <- coeffcols[coeffcols != "IDcol"]
    }
    results <- results[,c(coeffcols, "OR", "CI", "p")]
    names(results) <- c(coeffcols, "OR", "95% CI", "p")
  }
  rownames(results) <- NULL
  # https://www.r-bloggers.com/regression-on-categorical-variables/
  return(results)
}


# Round to the nearest m
mround <- function(x, base){
  base*round(x/base)
}


# Make a pretty proportion table
# To use this,
# tab <- table(data$VIhyp, data$prevHBP, useNA='ifany')
# kable(propped(tab), caption="")
propped <- function(table, margin=NULL) {
  prop <- round(100*prop.table(table, margin=margin),2)
  tabsums <- addmargins(table)
  dimt <- dim(table)
  for(i in c(1:dimt[1])) {
    if(length(dimt)>1){
      for(j in c(1:dimt[2])) {
        tabsums[i,j] <- paste0(table[i,j], " (", prop[i,j], "%)")
      }
    }
    else{
      tabsums[i] <- paste0(table[i], " (", prop[i], "%)") 
    }
  }
  return(tabsums)
}

# Check correlation among all pairs of covariates in a given dataframe 
# (this will need extending to be more robust if I want to use it for non-categorical covariates)
corr_mat <- function(data){
  corr_matrix <- matrix(0L, nrow=ncol(data), ncol=ncol(data))
  colnames(corr_matrix) <- colnames(data)
  rownames(corr_matrix) <- colnames(data)
  for(i in c(1:(ncol(data)-1))){
    x <- data[[i]]
    for(j in c((i+1):ncol(data))){
      y <- data[[j]]
      if(is.factor(x) & is.factor(y)){
        corr_matrix[i, j] <- (cramerV(table(x, y)))
      } else if (is.numeric(x) & is.numeric(y)) {
        corr_matrix[i, j] <- cor(x, y, method="pearson")
      } else {
        print("Unanticipated type")
      }
    }
  }
  return(corr_matrix)
}

preparecoefflist_1col <- function(df, varname, pretty_names=list()){
  pretty_varname <- prettyfunc(varname, pnames=pretty_names, bold=TRUE, upper=TRUE)
  if(is.factor(df[[varname]])){
    if(is.ordered(df[[varname]])){
      poly <- length(levels(df[[varname]]))
      levels <- c("Ref", ".L", ".Q", ".C")
      if(poly > 4){
        powers <- c(4:(poly-1))
        levels <- c(levels, paste0("^", powers))
      }
      levels <- levels[1:poly]
    } else {
      levels <- levels(df[[varname]])
    }
    variable <- c(pretty_varname, levels)
    coeffname <- c(pretty_varname, paste0(varname, levels))
  } else {
    variable <- pretty_varname
    coeffname <- varname
  }
  output <- data.frame(coeffname, variable, stringsAsFactors=FALSE)
  colnames(output) <- c("IDcol", "Coefficient")
  rownames(output) <- NULL
  return(output)
}

preparecoefflist_2col <- function(df, varname, pretty_names=list()){
  pretty_varname <- prettyfunc(varname, pnames=pretty_names, upper=TRUE)
  if(is.factor(df[[varname]])){
    if(is.ordered(df[[varname]])){
      poly <- length(levels(df[[varname]]))
      levels <- c("Ref", ".L", ".Q", ".C")
      if(poly > 4){
        powers <- c(4:(poly-1))
        levels <- c(levels, paste0("^", powers))
      }
      levels <- levels[1:poly]
    } else {
      levels <- levels(df[[varname]])
    }
    variable <- c(pretty_varname,
                  rep(NA, length(levels)-1))
    coeffname <- paste0(varname,levels)
  } else {
    levels <- NA
    variable <- pretty_varname
    coeffname <- varname
  }
  output <- data.frame(coeffname, variable, levels, stringsAsFactors=FALSE)
  colnames(output) <- c("IDcol", "Coefficient", "Levels")
  rownames(output) <- NULL
  return(output)
}

preparecoefflist <- function(onecol=FALSE, ...){
  if(onecol) {
    preparecoefflist_1col(...)
  } else {
    preparecoefflist_2col(...)
  }
}

regressiontable <- function(df, outcome, varlist, regresstype, adjvarlist=c("agegrp", "gender"), pretty_names=list(), IDcol=TRUE){
  coefflist <- list()
  # Prepare the list of coefficients - variables and levels for factors or blanks for continuous
  for(var in varlist){
    coefflist[[var]] <- preparecoefflist(df=df, varname=var, pretty_names=pretty_names)
  }
  
  if(regresstype=="univariable"){
    modellist <- list()
    for(var in varlist){
      coeffnames <- coefflist[[var]]
      
      # Prepare the formula and pass it to the model
      formula <- paste0(outcome, " ~ ", var)
      model <- glm(formula, data=df, family="binomial")
      
      # Add the pretty-formatted outputs to the list
      modellist[[var]] <- printlogresults(model, coeffnames, IDcol=IDcol)
    }
    # Vertically concatenate all the pretty outputs into one output table
    outdf <- do.call(rbind, modellist)
    
  } else if (regresstype=="adjusted"){
    modellist <- list()
    
    # Run the regressions for age and gender separately to go on top of the table
    for(adjvar in adjvarlist){
      coeffnames <- preparecoefflist(df=df, varname=adjvar)
      
      formula <- paste0(outcome, " ~ ", paste(adjvarlist, collapse="+"))
      model <- glm(formula, data=df, family="binomial")
      
      # Add the pretty-formatted outputs to the list
      modellist[[adjvar]] <- printlogresults(model, coeffnames, IDcol=IDcol)
    }
    
    # Putting age or gender in the regression twice would confuse it, so make sure they're not in the varlist
    varlist <- varlist[!varlist %in% adjvarlist]
    
    for(var in varlist){
      coeffnames <- coefflist[[var]]
      
      # Prepare the formula and pass it to the model
      formula <- paste0(outcome, " ~ ", paste(adjvarlist, collapse="+"), "+", var)
      model <- glm(formula, data=df, family="binomial")
      
      # Add the pretty-formatted outputs to the list
      modellist[[var]] <- printlogresults(model, coeffnames, IDcol=IDcol)
    }
    outdf <- do.call(rbind, modellist)
    
  } else if (regresstype=="multivariable"){
    coeffnames <- do.call(rbind, coefflist)
    formula <- paste0(outcome, " ~ ", paste(varlist, collapse=" + "))
    model <- glm(formula, data=df, family="binomial")
    outdf <- printlogresults(model, coeffnames, IDcol=IDcol)
  }
  
  rownames(outdf) <- NULL
  return(outdf)
}


output1_2<-function(outcome="TEU_paren_status",RF,return="table"){
  
  #' @param outcome Outcome indicator
  #' @param RF vector of overlapped RFs fitted in both models
  #' @param return raw model output or merged table or VIF 
  #' @return if return = "table", then Merged table output ; if return = "model", then list of 2 raw model outputs; 
  
  # Without BMI
  model1<-glm(formula=as.formula(paste0(outcome,"~",paste(RF, collapse="+"))),data=analysis_data,na.action = na.exclude,family = "binomial")
  
  model1_output<-printresults(analysis_data,varlist = RF,model1,pretty_names = pretty_names,onecol = TRUE)
  
  # With BMI
  
  model2<-glm(formula=as.formula(paste0(outcome,"~",paste(c(RF,'TEU_BSM_BMIcat'), collapse="+"))),data=analysis_data,na.action = na.exclude,family = "binomial")
  
  model2_output<-printresults(analysis_data,varlist = c(RF,'TEU_BSM_BMIcat'),model2,pretty_names = pretty_names,onecol = TRUE)
  
  if(return=="table"){
    # Merge model outputs together
    tot_output1<-left_join(model2_output,model1_output,by=c("Coefficient"))%>%
      `colnames<-`(c("Variable", "OR", "95% CI" , "p","OR" ,"95% CI", "p"))
    
    return(tot_output1)
    
  } else if(return=="model"){
    # Append raw model outputs as list
    models<-list(model1,model2)
    names(models)<-c("model1","model2")
    return(models)
    
  }
  
}

# Create function for flexible switching between outcomes and PRS.

MR_outputs <- function(outcome="Surv(TEU_paren_time,TEU_paren_status)",PRS="TEU_SBP_PRS",
                       exposure="TEU_BlP_SBP.avg",
                       scale_unit=5,
                       data=analysis_data,RF,
                       return,
                       forestplot=FALSE){
  
  # First build linear regression (BP ~ BP PRS) for relevance check and scale for interpretation
  lm <- lm(formula = as.formula(paste0(exposure,"~",PRS,"+",paste(RF, collapse="+"))),data)
  
  if (return=="relevance"){
    # Then return summary tab of lm and quit
    return(summary(lm))
  } else if (return=="MR"){
    # Keep running the following
    
    # Scale PRS for easy interpretation 
    # (e.g. if scale_unit=5, then OR of BP PRS would mean OR per 5mmHg increase in BP from BP PRS)
    
    data[[PRS]]=data[[PRS]]*lm$coefficients[PRS]/scale_unit
    
    # Build Cox regression
    model <- coxph(formula = as.formula(paste0(outcome,"~",PRS,"+",paste(RF, collapse="+"))),data)
    
    # In fact, the beta of PRS from model above actually corresponds to betaZY/betaXY (when scale_unit=1)
    
    tab<-printresults(df=data,varlist = c(PRS,RF),modeloutput = model,pretty_names = pretty_names)
    
    if(forestplot==FALSE){
      
      return(tab)
      
    }else{
      # If one wants to create a forest plot, need numeric values
      
      summ <- summary(model)
      conf <- as.data.frame(summ$coeff)
      conf[,'lci'] <- exp(conf[,1] - (1.96*conf[,3]))
      conf[,'uci'] <- exp(conf[,1] + (1.96*conf[,3]))
      #conf[,'Estimate'] <- exp(conf[,1])
      plot_df <- data.frame(hr = c(conf[rownames(conf) %in% c(PRS),2]),
                           lci = c(conf[rownames(conf) %in% c(PRS),6]),
                           uci = c(conf[rownames(conf) %in% c(PRS),7]),
                           pvalue = c(conf[rownames(conf) %in% c(PRS),5])
      )
      
      return(list(tab=tab,plot=plot_df))
    
    }
    
    
  }
  
}

combine_forestplot<-function(PRS1,PRS2,PRS3,which_exposure="SBP",which_PRS="SBP PRS",which_outcome="Renal parenchyma cancer"){
  
  # Combine forest plots
  tabletext <- cbind(c("Versions of PRS",paste0("Full ",which_PRS),paste0(which_PRS," (SNP-confounder p<5e-8)"),paste0(which_PRS," (SNP-confounder p<5e-6)")),
                     c("HR",PRS1$tab$HR[1],PRS2$tab$HR[1],PRS3$tab$HR[1]),
                     c("95% CI",PRS1$tab$`95% CI`[1],PRS2$tab$`95% CI`[1],PRS3$tab$`95% CI`[1]),
                     c("p-value",PRS1$tab$p[1],PRS2$tab$p[1],PRS3$tab$p[1])
  )
  
  p <- forestplot(tabletext,
                  graph.pos=4,
                  #legend = c("PRS"),
                  #legend_args = fpLegend(title="Quintiles of BP PRS",
                  # gp = gpar(col = "#CCCCCC", fill = "#F9F9F9")),
                  mean = t(cbind(NA,PRS1$plot$hr, PRS2$plot$hr,PRS3$plot$hr)),
                  lower = t(cbind(NA,PRS1$plot$lci, PRS2$plot$lci,PRS3$plot$lci)),
                  upper = t(cbind(NA,PRS1$plot$uci, PRS2$plot$uci,PRS3$plot$uci)),
                  #col = fpColors(box = c("green", "blue", "purple", "red", "orange")),
                  xlab = paste0("HR per 5mmHg increase in ",which_exposure," from ",which_PRS),
                  vertices=TRUE,
                  xlog=TRUE,
                  xticks=c(1,1.1,1.2,1.3),
                  col = fpColors(box = "royalblue",
                                 line = "darkblue"
                  ),
                  title=paste0("Forest plot of ",which_PRS," with ",which_outcome)
  )
  return(p)
  
  
}



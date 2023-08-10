library(tidyverse);library(mediation);library(pbapply)
select <- dplyr::select

# datasets were loaded as `dat.all.clean` for NACC data and `mmselast2` for cognition data

# variables to filter complete cases
inclu_criteria <- c("NACCAPOE4", "NPADNC_dic", "NACCLEWY_dic", "NACCAMY_dic") 
npi_items <- c("DEL_last", "HALL_last", "AGIT_last", "DEPD_last", "ANX_last", 
               "ELAT_last", "APA_last", "DISN_last", "IRR_last", "MOT_last", 
               "NITE_last", "APP_last")
sim_time <- 5000

# Function mediation analysis for repetition
med_analysis <- function(out.var, indep.var, med.var, 
                         med.family, out.family, polr = F,
                         adjust, adjust.med = "", simulation, dat){
    if (polr == T){
        require(MASS)
        # proportional odds logistic regression
        med.mod <- polr(as.formula(paste(med.var, 
                                         "~", "SEX+AGE_at_visit+EDUC+last_to_death_year+",
                                         indep.var, adjust.med)), 
                        data = dat, method = med.family, Hess = T) 
    }
    else if (polr == F){
        med.mod <- glm(as.formula(paste(med.var, 
                                        "~", "SEX+AGE_at_visit+EDUC+last_to_death_year+",
                                        indep.var, adjust.med)), 
                       data = dat, family = med.family) 
    }
    else {
        warning("polr arguement is logical")
    }
    y.mod <- glm(as.formula(paste(out.var,
                                  "~", adjust, "+",
                                  indep.var, "+",
                                  med.var)),
                 data = dat, family = out.family)
    set.seed(1234)
    med1 <- mediate(med.mod, y.mod, 
                    treat = indep.var, mediator = med.var,
                    robustSE = T, sims = simulation)
    return(list(y.mod, med.mod, med1))
}

# Function for mediation analysis with parallel mediators
source("mediate.new.R")

# Function for extract mediation output
xtract_mediation_summary <- function(x){ 
    
    clp <- 100 * x$conf.level
    isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                       (inherits(x$model.y, "glm") && x$model.y$family$family == 
                            "gaussian" && x$model.y$family$link == "identity") || 
                       (inherits(x$model.y, "survreg") && x$model.y$dist == 
                            "gaussian"))
    
    printone <- !x$INT && isLinear.y
    
    if (printone) {
        
        smat <- c(x$d1, x$d1.ci, x$d1.p)
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        
        rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
        
    } else {
        smat <- c(x$d0, x$d0.ci, x$d0.p)
        smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
        smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
        smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
        smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
        
        rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                            "ADE (control)", "ADE (treated)", "Total Effect", 
                            "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                            "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
        
    }
    
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                        paste(clp, "% CI Upper", sep = ""), "p-value")
    smat
    
}

# Data preprocessing
dat.model <- left_join(dat.all.clean, mmselast2, by = "NACCID") %>% 
    dplyr::select(NACCID, AGE_at_visit, EDUC, SEX, del_recall_harmo,
                  last_to_death_year, zero_path, NACCAPOE, 
                  NPTHAL, NACCBRAA, NACCNEUR,
                  all_of(inclu_criteria), all_of(npi_items)) %>%
    mutate_at(names(.)[-1:-12], as.factor) 

# Dataframe with logical memory test delay recall completed
dat.med.delay <- dat.model %>% 
    mutate_at(inclu_criteria, function(x){as.numeric(x)-1}) %>% 
    drop_na(all_of(inclu_criteria), all_of(npi_items), del_recall_harmo, EDUC) 

skimr::skim(dat.med.delay)

# Preliminary Association Analyses
## E4 -> AD/LBD/CAA (ADRD pathology) (adjusted)
mod1 <- lapply(inclu_criteria[-1], 
               function(x){
                   glm(as.formula(
                       paste(x,"~SEX+AGE_at_visit+EDUC+last_to_death_year+NACCAPOE4")),
                       data = dat.med.delay, 
                       family = binomial(link = "logit"))
               }
)

sjPlot::tab_model(mod1)

# Primary Mediation Analyses
## Mediation 1: E4 -> neuropathologies (AD,LBD,CAA) -> Logical Memory (delay)
med.ad <- glm(NPADNC_dic~SEX+AGE_at_visit+EDUC+last_to_death_year+NACCAPOE4, 
              data = dat.med.delay, family = "binomial") 
med.lewy <- glm(NACCLEWY_dic~SEX+AGE_at_visit+EDUC+last_to_death_year+NACCAPOE4, 
                data = dat.med.delay, family = "binomial") 
med.amy <- glm(NACCAMY_dic~SEX+AGE_at_visit+EDUC+last_to_death_year+NACCAPOE4, 
               data = dat.med.delay, family = "binomial")

y.mod <- glm(del_recall_harmo~SEX+AGE_at_visit+EDUC+last_to_death_year+NACCAPOE4
             +NACCLEWY_dic+NPADNC_dic+NACCAMY_dic,
             data = dat.med.delay, family = "gaussian")

med1 <- mediate_multi(med.lewy, y.mod, 
                      treat = "NACCAPOE4", mediator = "NACCLEWY_dic",
                      robustSE = T, sims = sim_time, seed = 1234)

med2 <- mediate_multi(med.ad, y.mod, 
                      treat = "NACCAPOE4", mediator = "NPADNC_dic",
                      robustSE = T, sims = sim_time, seed = 1234)

med3 <- mediate_multi(med.amy, y.mod, 
                      treat = "NACCAPOE4", mediator = "NACCAMY_dic",
                      robustSE = T, sims = sim_time, seed = 1234)

mod1b.sumtb <- integrate_multimed(med1, med2, med3)

## Mediation 2: AD/LBD/CAA -> Logical Memory (delay) -> NPI items
adrd <- c(inclu_criteria[-1],
          "NPADNC_dic + NACCLEWY_dic + NACCAMY_dic") 
# grid for mapply
grid.adjust.adrd <- expand.grid(adrd[-4], npi_items) %>% 
    mutate(Var3 = rep(c("NACCLEWY_dic + NACCAMY_dic",
                        "NPADNC_dic + NACCAMY_dic",
                        "NPADNC_dic + NACCLEWY_dic"), 12)) 

med2b.1.out <- pbmapply(function(x, y, z){
    med_analysis(out.var = paste(x), indep.var = paste(y), 
                 med.var = "del_recall_harmo", 
                 med.family = "gaussian", out.family = "binomial", 
                 adjust = paste("SEX+AGE_at_visit+EDUC+last_to_death_year+",
                                z, sep = ""),
                 adjust.med = paste("+", z),
                 simulation = sim_time, dat = dat.med.delay)
}, grid.adjust.adrd$Var2, grid.adjust.adrd$Var1, grid.adjust.adrd$Var3)

## Mediation 3: E4 -> Logical Memory (delay) -> NPI items controlled for ADRD
med4b.out <- pblapply(npi_items, med_analysis,
                      med.var = "del_recall_harmo", indep.var = "NACCAPOE4", 
                      med.family = "gaussian", out.family = "binomial", 
                      adjust = paste(
                          "SEX+AGE_at_visit+EDUC+last_to_death_year+
                     NPADNC_dic+NACCLEWY_dic+NACCAMY_dic"),
                      adjust.med = paste("+NPADNC_dic+NACCLEWY_dic+NACCAMY_dic"),
                      dat = dat.med.delay, simulation = sim_time,
                      cl = 5)

## Mediation 4: E4 -> AD/LBD/CAA -> NPI items controlled for Logical Memory (delay)
med5b.out <- pbmapply(function(x, y, z){
    med_analysis(out.var = paste(x), indep.var = "NACCAPOE4", 
                 med.var = paste(y), 
                 med.family = "binomial", out.family = "binomial", 
                 adjust = paste("SEX+AGE_at_visit+EDUC+last_to_death_year+del_recall_harmo+",
                                z, sep = ""),
                 adjust.med = paste("+", z),
                 simulation = sim_time, dat = dat.med.delay)
}, grid.adjust.adrd$Var2, grid.adjust.adrd$Var1, grid.adjust.adrd$Var3)



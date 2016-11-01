
#' @title birdegg_ssd.8 data
#' @description Data set on bird egg sexual size dimorphism meta-analysis. Data was taken from Rutkowska et al. (2014). In this study more than one egg came from the same study which introduces non-independence in the data (i.e. correlated structure). As a sensitivity analysis Rutkowska et al. (2014) used an effective sample size proposed by Higgins and Green (2009), which requires an estimate of the intra-class correlation coefficient (ICC) in the calculation of effect sizes. In this case, an ICC estimate of 0.8 was used. 
#' @format A raw csv file containing the following columns and descriptors
#' \enumerate{
#' \item animal:	species Latin name
#' \item Spp:	species common name
#' \item Lndim:	logaritm of adult mass dimorphism
#' \item Measure:	method of eggs measurments
#' \item Neggs:	number of eggs
#' \item Nclutches:	number of clutches
#' \item ESr:	effect size, correlation
#' \item Type:	how the effect size was calculated: row data/inferential statistics
#' \item StudyID:	reference
#' \item Year:	year of publication
#' \item D:	design effect
#' \item EN:	effective sample size
#' \item Zr:	Fisher's Z
#' \item VZr:	sampling variance
#'}	
#' @references Rutkowska, J., Dubiec, A. and Nakgawa, S. (2014). All eggs are made equal: a meta-analysis of egg sexual size dimorphism. Journal of Evolutionary Biology, 27:153-160
#' @references Higgins, J.P.T. and Green, S. (2009). Cochran handbook for systematic reviews and interventions. Wiley-Blackwell. pgs 412-414.
#' @name birdegg_ssd.8
#' @docType data
NULL

#' @title birdegg_ssd.5 data
#' @description Data set on bird egg sexual size dimorphism meta-analysis. Data was taken from Rutkowska et al. (2014). In this study more than one egg came from the same study which introduces non-independence in the data (i.e. correlated structure). As a sensitivity analysis Rutkowska et al. (2014) used an effective sample size proposed by Higgins and Green (2009), which requires an estimate of the intra-class correlation coefficient (ICC) in the calculation of effect sizes. In this case, an ICC estimate of 0.5 was used.
#' @format A raw csv file containing the following columns and descriptors
#' \enumerate{
#' \item animal:	species Latin name
#' \item Spp:	species common name
#' \item Lndim:	logaritm of adult mass dimorphism
#' \item Measure:	method of eggs measurments
#' \item Neggs:	number of eggs
#' \item Nclutches:	number of clutches
#' \item ESr:	effect size, correlation
#' \item Type:	how the effect size was calculated: row data/inferential statistics
#' \item StudyID:	reference
#' \item Year:	year of publication
#' \item D:	design effect
#' \item EN:	effective sample size
#' \item Zr:	Fisher's Z
#' \item VZr:	sampling variance
#'}	
#' @references Rutkowska, J., Dubiec, A. and Nakgawa, S. (2014). All eggs are made equal: a meta-analysis of egg sexual size dimorphism. Journal of Evolutionary Biology, 27:153-160
#' @references Higgins, J.P.T. and Green, S. (2009). Cochran handbook for systematic reviews and interventions. Wiley-Blackwell. pgs 412-414.
#' @name birdegg_ssd.5
#' @docType data
NULL

#' @title matdiet
#' @description Data set on how maternal diet impacts copying styles in rodents. Data was taken from Besson et al. (2016). 
#' @format A raw csv file containing the following columns and descriptors
#' @references Besson, A. A., Lagisz, M., Senior, A.M., Hector, K.L. and Nakagawa, S. (2016). Effect of maternal diet on offspring coping styles in rodents: a systematic review and meta-analysis. Biological Reviews, 91:1065-1080
#' @name matdiet
#' @docType data
NULL

#' @title moSizetradeoff
#' @description Data set investing the correlations between maternal size and offspring size (dataset 1: MO.size). Data was taken from Lim et al. (2014). 
#' @format A raw csv file containing the following columns and descriptors
#' @references Lim, J. N., Senior, A. M. and Nakagawa, S. (2014) Heterogeneity in individual quality and reproductive trade-offs within species. Evolution, 68:2306-2318.
#' @name moSizetradeoff
#' @docType data
NULL

#' @title moNumtradeoff
#' @description Data set investing the correlations between maternal size and offspring number (dataset 2: MO.Number). Data was taken from Lim et al. (2014). 
#' @format A raw csv file containing the following columns and descriptors
#' @references Lim, J. N., Senior, A. M. and Nakagawa, S. (2014) Heterogeneity in individual quality and reproductive trade-offs within species. Evolution, 68:2306-2318.
#' @name moNumtradeoff
#' @docType data
NULL

#' @title ooAdjtradeoff
#' @description Data set investing the correlations between offspring size and offspring number after adjusting for maternal size(dataset 3: OO.Adjusted). Data was taken from Lim et al. (2014). 
#' @format A raw csv file containing the following columns and descriptors	
#' @references Lim, J. N., Senior, A. M. and Nakagawa, S. (2014) Heterogeneity in individual quality and reproductive trade-offs within species. Evolution, 68:2306-2318.
#' @name ooAdjtradeoff
#' @docType data
NULL

#' @title ooUnAdjtradeoff
#' @description Data set investing the correlations between offspring size and offspring number not adjusting for maternal size (dataset 4: OO.Unadjusted). Data was taken from Lim et al. (2014). 
#' @format A raw csv file containing the following columns and descriptors	
#' @references Lim, J. N., Senior, A. M. and Nakagawa, S. (2014) Heterogeneity in individual quality and reproductive trade-offs within species. Evolution, 68:2306-2318.
#' @name ooUnAdjtradeoff
#' @docType data
NULL

#' @title mtor data
#' @description Data set investing whether changes in mTOR (rapamycin) signaling extend lifespan. Data was taken from Garratt et al. (2016). 
#' @format A raw csv file containing the following columns and descriptors	
#' @references Garratt, M., Nakagawa, S., Simons, M.J.P. (2016) Comparative idiosyncrasies in life extension by reduced mTOR signaling and its distinctiveness from dietary restriction. Aging Cell, 15: 737-743.
#' @name mtor
#' @docType data
NULL


#' @title heat shock and life extension data
#' @description Data set investing the impact of heat shock exposure on longevity. Data was taken from Lagisz et al. (2013).
#' \enumerate{
#' \item dataID: Unique identifier for each comparison of a pair of survival curves between the control group (con) and the heat-shocked group (hs).
#' \item papID: Unique identifier for each study (i.e. paper) in the dataset.
#' \item expID: Unique identifier for separate experiments within each study (usually determined by having separate control group).
#' \item author: Study first author nam.
#' \item pubyear: Year of publication of the study.
#' \item journal: Name of journals in which the study was published.
#' \item IF: Journal Impact Factor (2011, JCR).
#' \item name: Common name of the species.
#' \item genus: Genus name of the species.
#' \item species: Latin name of the species.
#' \item strain: Unique strain name for particular.
#' \item class: Taxonomic classe of the species.
#' \item family:  Taxonomic family of the species.
#' \item animal: Unique identifier for  the species.
#' \item sex: Sex category (M = males, F = females, H = none/hermaphrodite or B = males and females).
#' \item age: Age of animals at the time of the first heat shock (days from egg) [day].
#' \item relage: Age of animals individuals at the time of the first heat shock (days from egg) divided by maximum life span observed in the control group (con_maxLS).
#' \item conN: Number of individuals in the control group.
#' \item hsN: Number of individuals in the heat-shocked group.
#' \item conTemp: Standard (control) temperature in which animals were kept [C].
#' \item hsTemp: Heat-shock temperature [C].
#' \item tempDiff: Difference between control andheat-shock temperature [C].
#' \item duration: Single heat-shock duration [h].
#' \item repeats: How many times heat-shock was repeated.
#' \item spacing: Time interval between heat shocks [h].
#' \item maxCon: Time from the start of the experiment when the last animal in the control group died (maximum observed longevity) [day].
#' \item maxCon_LS: Time from the start of the experiment when the last animal in the control group died, corrected for the age at the start of the experiment (maximum lifespan) [day].
#' \item p10con: Time since the start of the experiment when 10% of control animals were still alive [day].
#' \item p10con_LS: Time since the start of the experiment when 10% of control animals were still alive, corrected for the age at the start of the experiment [day].
#' \item p90,p80,p70,p60,p50,p40,p30,p20,p10: Proportion of heat-shocked group animals that were alive when 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20, 10%, respectively, of control animals was still alive [proportion].
#' \item lnHR: Overall ln(HR) for 9 time intervals.
#' \item varHR: Variance of overall ln(HR) for the 9 time intervals.
#' \item conN_uniq: Number of individuals in the control group, unique within each experiment.
#' \item hsN_uniq: Number of individuals in the heat-shocked group, unique within each experiment.
#' \item data_source: Name of the figure or table from which the survival data was extracted (or raw data, tables with additional information).
#'} 
#' @format A raw csv file containing the following columns and descriptors	
#' @references Lagisz, M., Hector, K.L. and Nakagawa, S. (2013) Life exstension after heat shock exposure: Assessing meta-analytic evidence for hormesis. Ageing Research Reviews, 12:653-660.
#' @name heatshock
#' @docType data
NULL


#' @title sparrow data
#' @description  Correlation between 
#' @format A raw csv file containing the following columns and descriptors	
#' \enumerate{
#' \item StudyID: Study number
#' \item Place: Location of sparrow population
#' \item Correlation: correlation coefficient
#' \item SampleSize: Sample size for study
#'}
#' @references Nakgawa, S., Ockendon, N., Gillespie, D.O.S., Hatchwell, B.J. and Burke, T. (2007) Assessing the function of house sparrows' bib size using a flexible meta-analysis method. Behavioural Ecology, 18:831-840.
#' @name sparrows
#' @docType data
NULL

#' @title Host manipulation and parasite data
#' @description  Correlation between 
#' @format A raw csv file containing the following columns and descriptors. UI = Uninfected; U = Infected	
#' \enumerate{
#' \item Parasite.species: Species of parasite
#' \item Phylum: Phylum
#' \item Class: Class
#' \item Family: Family
#'}
#' @references Poulin, R. (2000) Manipulation of host behaviour by parasites: a weakening paradigm? Proceedings of the Royal Society B: Biological Sciences, 267: 1471-2954.
#' @name parasites
#' @docType data
NULL


#' @title Within pair paternity data and age
#' @description  Proportion of within and extra-pair paternity for different ages. Data taken from Cleasby and Nakagawa (2012).
#' @format A raw csv file containing the following columns and descriptors. 
#' @references Cleasby, I.R. and Nakagawa, S. (2012) The influence of male age on within-pair and extra-pair paternity in passerines. Ibis, 154: 318-324.
#' @name age_epp_wpp
#' @docType data
NULL

#' @title Multivariate wpp and epp
#' @description  Proportion of within and extra-pair paternity for different ages. Data taken from Cleasby and Nakagawa (2012).
#' @format A raw csv file containing the following columns and descriptors. 
#' @references Cleasby, I.R. and Nakagawa, S. (2012) The influence of male age on within-pair and extra-pair paternity in passerines. Ibis, 154: 318-324.
#' @name multi_epp_wpp
#' @docType data
NULL

#' @title plumage and dominance
#' @description  A meta-analysis of the relationship between plumage traits and dominance in birds. Data taken from Santos et al. (2011).
#' @format A raw csv file containing the following columns and descriptors. 
#' @references Santos, E.S.A., Scheck, D. and Nakagawa, S. (2011) Dominance and plumage traits: meta-analysis and meta-regression analysis. Animal Behaviour, 82: 3-19.
#' @name plumage
#' @docType data
NULL
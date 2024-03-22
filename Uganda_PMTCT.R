#Uganda PMTCT 
# Contact: Sahal Thahir sahal.thahir@unchealth.unc.edu
library(tidyverse)
library(table1)
library(REDCapR)
install.packages('DescTools')
install.packages('manipulate')

#RedCAP API export
data1 <- redcap_read_oneshot(redcap_uri = "https://global.redcap.unc.edu/api/",
                             token = "")$data

#Hep B tested subset
HepB <- subset(data1, !is.na(v1_hbv_res))

# Study specific (read prior to running)
  # Use for prelim PMTCT study: screening prior to 21 January 2022
    HepB <- subset(HepB, screen_date < '2022-01-21')
  

#Factorizing Existing columns---- 
  #HBV factor testing
  HepB$v1_hbv_res = factor(HepB$v1_hbv_res, levels = 0:1, labels = c("HBV-", "HBV+"))

  #Education
  HepB$edu_level = factor(HepB$edu_level, levels = 0:3, labels = c("None", "Primary","Secondary", "Tertiary"))
 
  #Marital status
  HepB$marital_status = factor(HepB$marital_status, levels = 0:5, labels = c("Not married", "Married", "Domestic", "Separated", "Divorced","Widowed"))
  #Religion
  HepB$religion = factor(HepB$religion, levels = 0:4, labels = c("Protestant", "Catholic", "Pentecostal", "Muslim", "Other"))

  #Subcounty
  HepB$sub_county = factor(HepB$sub_county, levels= 0:1, labels = c("Bugoye", "Other"))
 
  #Maternity center
  HepB$v1_site = factor(HepB$v1_site, levels= 0:4, labels= c("Bugoye","Katooke","Kibirizi","Mukhati","Other"))
  
  #HIV
  HepB$v1_hiv_res = factor(HepB$v1_hiv_res, levels= 0:1, labels= c("Negative", "Positive"))
  
  #Malaria
  HepB$v1_mrdt_res = factor(HepB$v1_mrdt_res, levels=0:1, labels= c("Negative", "Positive"))
 
  #Syphilis
  HepB$v1_sphy_res = factor(HepB$v1_sphy_res, levels= 0:1, labels= c("Negative", "Positive"))
  
  #Delivery location
  HepB$deliv_where = factor(HepB$deliv_where, levels= 0:4, labels= c("Home","Home with provider","Health Center","Hospital","Other"))
  
  #Any children?
    HepB$kids <- factor(HepB$kids, levels= 0:1, labels= c("No", "Yes"))
    
# New Columns----
    #HBV viremia (Y/N): HBV_VL
      HepB$HBV_VL <- ifelse(HepB$v1_hbv_vl > 0, TRUE, FALSE)
      HepB$HBV_VL[is.na(HepB$HBV_VL)] <- FALSE
    #Gestational age at presenetation (1st, 2nd, 3rd trimester)
      HepB$trim_v1 <- cut(HepB$v1_weeks, breaks=c(0,12, 24, Inf)
          
    #Gestational age difference: Diff_weeks and deliv_weeks
      HepB$Diff_weeks <- round(as.numeric(difftime(HepB$deliv_date, HepB$v1_date, units= "weeks")), 1)
      HepB$deliv_weeks <- HepB$v1_weeks + HepB$Diff_weeks
    #Term vs Prematurity: deliv_term
      HepB$deliv_term <- ifelse(HepB$deliv_weeks > 37, "Term", "Premature")
    #Normal birthweight vs. LBW
      HepB$deliv_wt<- ifelse(HepB$deliv_weight> 2.5, "Normal birthweight", "Low Birthweight")
        #Binomial (for regression analysis)
          HepB$deliv_wt_bi <- ifelse(HepB$deliv_wt == "Low Birthweight", 1, 0)
          table(HepB$deliv_wt_bi)
      #APGARs (scores >/=7)
      HepB$deliv_apgar1 <- ifelse(HepB$deliv_apgar1 >6, ">7", "<7")
      HepB$deliv_apgar5 <- ifelse(HepB$deliv_apgar5 >6, ">7", "<7")
      
      #Documented delivery
        HepB$delivery <- ifelse(is.na(HepB$deliv_date), FALSE, TRUE)
      
      #C-section: "sect", "oper", and "blade
      HepB$deliv_Csection <- grepl("sect", HepB$deliv_comp, fixed= TRUE) | grepl("oper", HepB$deliv_comp, fixed= TRUE) | grepl("blade", HepB$deliv_comp, fixed= TRUE)
      #Hemorrhage: "rhage"
      HepB$deliv_hemorrhage <- grepl("rhage", HepB$deliv_comp, fixed= TRUE) 
      #Stillbirth: "died" or "dead
      HepB$deliv_stillbirth <- grepl("dead", HepB$deliv_comp, fixed= TRUE) | grepl("died", HepB$deliv_comp, fixed= TRUE) 
      
      #Infant gge at Epi1, Epi2, Epi3
      HepB$epi1_age<- round(as.numeric(difftime(HepB$epi1_date, HepB$deliv_date , units = "weeks")))
      HepB$epi2_age<- round(as.numeric(difftime(HepB$epi2_date, HepB$deliv_date , units = "weeks")))
      HepB$epi3_age<- round(as.numeric(difftime(HepB$epi3_date, HepB$deliv_date , units = "weeks")))
      HepB$Vac1_age<- round(as.numeric(difftime(HepB$penta1_date, HepB$deliv_date , units = "weeks")))
      HepB$Vac2_age<- round(as.numeric(difftime(HepB$penta2_date, HepB$deliv_date , units = "weeks")))
      HepB$Vac3_age<- round(as.numeric(difftime(HepB$penta3_date, HepB$deliv_date , units = "weeks")))
      
      #Vaccination age per EPI/Weyer protocol
      HepB$EPI1 <- ifelse(HepB$Vac1_age<6, "<6 wo", ">6 wo")
      HepB$EPI2 <- ifelse(HepB$Vac2_age <14, "<14 wo", ">14 wo")
      HepB$EPI3 <- ifelse(HepB$Vac3_age <14, "<14 wo", ">14 wo")
#Subsets----  
      #Demographics (HepB_melt)
      HepB_melt <- subset(HepB, select= c(study_id, v1_site, v1_hbv_res, hbv_vl_copies, HBV_VL, age, trim_v1,  v1_weeks, kids, kids_num, marital_status, edu_level, religion, v1_hiv_res, v1_sphy_res, v1_mrdt_res, sub_county, epi1_vacc) )
     #Clinical characteristics (Clin + ClinPos)
       Clin <- subset(HepB, select= c(study_id, v1_site, v1_hbv_res, hbv_vl_copies, age, v1_hiv_res, v1_sphy_res, v1_mrdt_res, sub_county, v1_sbp, v1_dbp, hgb:alt, deliv_where, deliv_weight, deliv_apgar1, deliv_apgar5) )
      Clin_pos <- subset(Clin, v1_hbv_res=="HBV+") 
      #Hep B positive only (HepBPos)
      HepBPos <- subset(HepB, v1_hbv_res =="HBV+")
      #High Risk HBV limited
      HR_HepB <- subset(HepB, hbv_vl_copies > 200000)
      #Vaccination data (Vaccine) 
      Vaccine <- subset(HepB, select= c(study_id, v1_hbv_res, delivery, EPI1, EPI2, EPI3)) 
      #Delivered infants only
      Delivery <- subset(HepB, delivery== TRUE)
      #Outcome dataframe
      Outcome <- subset(Delivery, select= c(study_id, v1_hbv_res, HBV_VL, age, trim_v1, deliv_weeks, deliv_wt, deliv_term, deliv_hemorrhage, deliv_Csection, deliv_stillbirth ))
      Outcome <- na.omit(Outcome) 
        #Required special column: delivery complications cleaned up 
        Outcome$deliv_comp <- case_when(Outcome$deliv_hemorrhage == TRUE ~"Hemorrhage",
                                      Outcome$deliv_Csection == TRUE ~"C-section",
                                      Outcome$deliv_stillbirth == TRUE ~"Stillbirth")
            
      
      
      
#Tables (Demographics/Testing) ----

#Table 1 Demographics (Age, gestation, children, demographics based on HBV status)
  label(HepB_melt$age) <- "Age"
  label(HepB_melt$v1_weeks) <- "Gestational age at screening "
  label(HepB_melt$trim_v1) <- "Gestational trimester at screening"
  label(HepB_melt$kids) <- "Living children"
  label(HepB_melt$sub_county) <- "Subcounty of residence"
  label(HepB_melt$marital_status) <- "Marital status"
  label(HepB_melt$edu_level) <- "Education"
  label(HepB_melt$religion) <- "Religion"
  label(HepB_melt$HBV_VL) <- "HBV viremia"
  units(HepB_melt$age) <- "years"
  units(HepB_melt$v1_weeks) <- "weeks"
  
  table1(~ age + v1_weeks + trim_v1 + kids +  sub_county + marital_status + edu_level + religion | v1_hbv_res, data=HepB_melt)
 # Table 1: HBV viremia
  table1(~ age + v1_weeks + kids +  sub_county + marital_status + edu_level + religion | HBV_VL, data=HepB_melt, overall= F)
  
# Table: HBV to other tested diseases
  label(HepB_melt$v1_hiv_res) <- "Rapid HIV"
  label(HepB_melt$v1_mrdt_res) <- "Malaria RDT"
  label(HepB_melt$v1_sphy_res) <- "Syphilis RPR"
  
  table1(~ v1_hiv_res + v1_sphy_res + v1_mrdt_res | v1_hbv_res, data= HepB_melt)

# Table: Testing site, screening results
  label(HepB_melt$v1_site) <- "Testing Center"
  label(HepB_melt$v1_hiv_res) <- "Rapid HIV"
  label(HepB_melt$v1_mrdt_res) <- "Malaria RDT"
  label(HepB_melt$v1_sphy_res) <- "Syphilis RPR"
  label(HepB_melt$v1_hbv_res) <- "HBV Antigen"
  
  table1(~ v1_hiv_res + v1_sphy_res + v1_mrdt_res + v1_hbv_res | v1_site, data= HepB_melt)
  
  
#Hep B Pos limited
  HepBPos$Risk <- ifelse(HepBPos$hbv_vl_copies > 200000, "High", "Low")
  table(data$hbv_vl_copies)
 
# Time between viral testing and result
   HepBPos$hbv_vl_time <- difftime(HepBPos$hbv_vl_date, HepBPos$v1_date)
  quantile(na.omit(HepBPos$hbv_vl_time))


 


#Outcomes dataframes and analysis---- 


                                  

  #Outcome analysis
    #HBV status
      label(Outcome$age) <- "Maternal age"
      label(Outcome$trim_v1) <- "Gestational trimester at initial visit "
      label(Outcome$deliv_weeks) <- "Gestational age at delivery"
      label(Outcome$deliv_wt) <- "Birthweight"
      label(Outcome$deliv_term) <- "Term or Preterm delivery"
      label(Outcome$deliv_comp) <- "Delivery complications"
      label(Outcome$HBV_VL) <- "HBV Viremia"
      
      units(Outcome$age) <- "years"
      units(Outcome$deliv_weeks) <- "weeks"
      units(Outcome$deliv_wt) <- "kg"
      table1(~ age + +trim_v1+ deliv_weeks + deliv_wt + deliv_term + deliv_comp| v1_hbv_res, data=Outcome)
    #HBV viremia (Y/N)
      table1(~ age + deliv_weeks + deliv_wt + deliv_term + deliv_comp| HBV_VL, data=Outcome, overall= F)
  
      #Vaccinations Cascade
  label(Vaccine$v1_hbv_res) <- "HBV status"
  label(Vaccine$delivery) <- "Documented delivery"
  label(Vaccine$EPI1) <- "Immunization at 1-6 weeks"
  label(Vaccine$EPI2) <- "Immunization #2 at 2-14 weeks"
  label(Vaccine$EPI3) <- "Immunization #3 at 2-14 weeks"

  table1(~delivery + EPI1 +EPI2 + EPI3 | v1_hbv_res, data= Vaccine, overall= F)
  
  #Binomial Regression for LBW (1) v Normal BW (0) using deliv_wt_bi
    # Continuous: maternal age +gestational age
    # Categorical: HBsAg, educational status
 
       Fit_con <- glm(deliv_wt_bi ~ age + deliv_weeks+ v1_hbv_res + edu_level, data = Delivery, family= "binomial")
        summary(Fit_con)    
 
  
#Old/Misc ----
  # Clinical Data Analysis-
  
  #hbv_risk (titer level risk stratification)
  Clin_pos$hbv_risk <- case_when(Clin_pos$hbv_vl_copies >199999 ~ "high", Clin_pos$hbv_vl_copies > 10000 ~ "medium", is.na(Clin_pos$hbv_vl_copies) ~ "Missing", TRUE ~ "low")
  #APRI (higher level AST set to 34 per McClendon guidelines
  # source: https://www.uncmedicalcenter.org/mclendon-clinical-laboratories/available-tests/aspartate-aminotransferase-ast/) )
  Clin_pos$APRI <- (Clin_pos$ast/34*100)/Clin_pos$plts
  
  #APRI_risk (<0.5, 0-1.5, >/=1.5)c
  Clin_pos$APRI_risk <- case_when(Clin_pos$APRI<= 0.5 ~"<= 0.5", Clin_pos$APRI >0.5 & Clin_pos$APRI<1.5 ~"0.5-1.5", Clin_pos$APRI >= 1.5 ~ ">=1.5", Clin_pos$APRI==NA ~ "Missing")
  
  # LFTs (AST/ALT cutoffs of 34/40 respectively )        
  Clin_pos$LFT <- if_else(Clin_pos$ast >34 | Clin_pos$alt >40 , "Abnormal", "Normal", "Missing")
  
  #Table: HBV status and clinical factors 
  label(Clin$age) <- "Age"
  label(Clin$v1_sbp) <- "Systolic pressure"
  label(Clin$v1_dbp) <- "Diastolic pressure"
  label(Clin$deliv_weight) <- "Birthweight"
  label(Clin$deliv_apgar1) <- "1 minute APGAR"
  label(Clin$deliv_apgar5) <- "5 minute APGAR"
  label(Clin$deliv_where) <- "Delivery location"
  
  units(Clin$age) <- "years"
  units(Clin$v1_sbp) <- "mmHg"
  units(Clin$v1_dbp) <- "mmHg"
  units(Clin$deliv_weight) <- "kg"
  
  table1(~ age + v1_sbp + v1_dbp + deliv_weight + deliv_apgar1 + deliv_apgar5 + deliv_where | v1_hbv_res, data= Clin)
  
  #HBV viral titer analysis
  
  # Table: HBV risk, clinical characteristics
  label(Clin_pos$age) <- "Age"
  label(Clin_pos$v1_sbp) <- "Systolic pressure"
  label(Clin_pos$v1_dbp) <- "Diastolic pressure"
  label(Clin_pos$hgb) <- "Hemoglobin"
  label(Clin_pos$plts) <- "Platelets"
  label(Clin_pos$creat) <- "Creatinine"
  label(Clin_pos$ast) <- "AST"
  label(Clin_pos$alt) <- "ALT"
  
  units(Clin_pos$age) <- "years"
  units(Clin_pos$v1_sbp) <- "mmHg"
  units(Clin_pos$v1_dbp) <- "mmHg"
  units(Clin_pos$hgb) <- "mg/dL"
  units(Clin_pos$plts) <- "10^9 mL"
  units(Clin_pos$creat) <- "mg/dL"
  units(Clin_pos$ast) <- "u/L"
  units(Clin_pos$alt) <- "u/L"
  
  table1(~ age + v1_sbp + v1_dbp + hgb + plts + creat + ast + alt | hbv_risk, data= Clin_pos)
  
  #Table: Abnormal LFTs, APRI with HBV 
  label(Clin_pos$hbv_risk) <- "HBV Risk"
  label(Clin_pos$LFT) <- "Liver Function Test"
  label(Clin_pos$APRI_risk) <- "APRI`"
  
  table1(~ LFT + APRI_risk | hbv_risk, data= Clin_pos)
  
  #Location clinical characteristics
  label(Clin_pos$deliv_where) <- "Delivery site"
  label(Clin_pos$v1_site) <- "Maternity site"
  label (Clin_pos$hbv_risk) <- "HBV titers"
  
  units(Clin_pos$hbv_risk) <- "IU/µL"
  table1(~ hbv_risk + deliv_where | v1_site, data= Clin_pos)
  
  
  #Mother cascade of care----
  #Mothers
  #underwent viral titer sampling
  table(HepBPos$v1_hbv_vl)
  # Titer resulted note the "0" from the above vector
  table(HepBPos$hbv_vl_yn)
  # High risk titer (>200k)  
  view(HepBPos$hbv_vl_copies)
  # TDF +/-
  table(HepBPos$tdf_start)
  # post-partum TDF +/-
  table(HepBPos$epi1_tdf)
  
  #Viral load 
  ViralLoad <- subset(HR_HepB, select= c(study_id, hbv_vl_copies, epi1_tdf, epi1_vl_copies, epi2_vl_copies, epi3_vl_copies))
  
  
  
  #Infant Cascade of care: HBV + subset) "Infant"----
  
  #Infant Age in weeks (UTD)
  HepBPos$IfnAge<- round(as.numeric(difftime(Sys.Date(), HepBPos$deliv_date, units = "weeks")))

  #6 week/14 week eligibility
  HepBPos$Wk6 <- ifelse(HepBPos$IfnAge> 5, "Yes", "No")
  Infant <- subset(HepBPos, select= c(study_id, Risk, IfnAge, Wk6, epi1_age, Vac1_age, epi2_age, Vac2_age, epi3_age, Vac3_age, deliv_date, epi1_date, epi1_hb_rdt, epi1_vacc, epi2_date, epi2_vacc1:penta3_date, epi3_date, epi3_hb_rdt, epi3_hb_rdt_res, epi1_hb_rdt_res) )
  #6 weeksh
  table(Infant$Wk6)
  table(Infant$epi1_date)

#First off: 
#only outcome.measure.new == 
#"Hazard Ratio, Odds Ratio, Peto Odds Ratio, Risk Ratio, Risk difference, Mean difference and Std. Mean Difference"

#outcome.type == bin: 
#not all(events1 == 0) & all(events2 == 0). 
#additionally: 62301 and 94519 (copas convergence issues)

#outcome.type == cont:
#not all(mean1 == 0) & all(mean2 == 0) and all(sd1 == 0) | all(sd2 == 0). 
#addtionally: - which(meta.id.vector > 42716 & meta.id.vector < 42725) (single patient data)

#outcome.type == surb:
#none



#cor.pearson analysis
#outcome.type != "surv"
#additionally: -157083 (is one study and has total1 = 5 throughout).., 159329 (gives issues with copas)

#cohens'd analysis
#outcome.type != "surv"

#zscore analysis
#outcome.type != "surv"
#additionally: -157083 (is one study and has total1 = 5 throughout).., 159329 (gives issues with copas)


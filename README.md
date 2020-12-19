# Kaplan-Meier-estimate-Cox-model

full code in code.R

Your task is to write a report on the relative effectiveness of two drugs A and B as long term treatments for a critical illness. The data assumes that 100 patients have been given drug A and another 100 have been given drug B. Times are assumed to be in years. The aim of the report is to assess if there is any difference in prolonging life between the two drugs A and B.

The R-script Assignment1.r (available on Vision) generates censored survival data (times to claim) for two groups A and B: the survival times are in Time, the censoring information is in Censor (0 indicates that the observation was censored) and group membership is in Group (1 indicates group A). Group A are those individuals who are taking drug A, while group B are those who are taking drug B. Generate your data by running the script.

You should also fit a Cox proportional hazards model to the data. Denote by &#955;A(t) and &#955;B(t) the instantaneous forces of mortality at time t for groups A and B respectively. You should suppose that the Cox proportional hazards model holds with group A taken as the reference group and &#955;B(t) = &#955;A(t)e&#914;:

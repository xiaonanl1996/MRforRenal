#--------------------------------------------------------------------------------------------------------------
# Xiaonan Liu 29/10/2020
# Update: 23/11/2020
# Create general function for producing descriptive characteristics table (in HTML format)
# Below are 2 packages for producing descriptive table 
#--------------------------------------------------------------------------------------------------------------
library(summarytools)
library(arsenal)


# 1. summarytools package
# dfSummary() function

# 2. arsenal package

# Customise summary statistics table
# (Nmiss: only show missing when there is missing data; 
#  Nmiss2: always show missing number even though there is none)
my_controls <- tableby.control(
    test = F,
    total = T,
    #numeric.test = "kwt", cat.test = "chisq",
    numeric.stats = c("meansd", "medianq1q3", "range", "Nmiss"), 
    cat.stats = c("countpct", "Nmiss"),
    stats.labels = list(
      meansd = "Mean (SD)",
      medianq1q3 = "Median (IQR)",
      range = "Min - Max",
      Nmiss = "Missing",
      Npct="N (Pct)"
    ),
    digits = 2L
  )
  

  

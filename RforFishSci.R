#*********************************************************
#****** Introduction to R for Fisheries Scientists *******
#****************** 28 February 2017 *********************
#**************** TNAFS Annual Meeting *******************
#**************** Knoxville, Tennessee *******************
#****** Author: Daniel Walker, dwalke44@vols.utk.edu *****
#*********************************************************


# Set the working directory:
# 1. Have a folder for every analysis to store copies of data, code sheets, output graphics, etc.
# NOTE: I highly recommend that you have an archival copy of all of your data saved as a MS Excel workbook (.xlsx), and save copies of each worksheet as individual .csv pages for analysis
# 2. Use setwd() at the beginning of every analysis to access data and store output in same location:
# setwd("C:/Path/to_your/analysis/folder")

setwd("C:/Users/danwa/OneDrive/Teaching/TNAFS R workshop")

# NOTE: always use the FORWARD SLASH "/" in pathnames

# The working directory path will be unique to a single computer, must be changed if running analysis on another machine

# The set.seed() command: state after setwd() command
# Many R procedures (e.g., randomly selecting training and test observations) utilize random number generation at some point. R uses an algorithm to generate pseudo-random numbers 
# (numbers that look random but are produced from a repeatable process). For reproducibility, you want to use the same starting point for that pseudo-number generation for your analysis, 
# and be able to recreate your initial results every time you run the full analysis. Use the set.seed() command to ensure the same output every time.
# Try to use a seed with some meaning not related to the analysis to prevent accusations of 'p-hacking'!

sample(1:20, 15, replace = F) #run this line multiple times to see different randomly selected numbers and order

set.seed(865)
sample(1:20, 15, replace = F)

set.seed(42)
sample(1:20, 15, replace = F)

set.seed(865)
sample(1:20, 15, replace = F)

# Use set.seed() right after setwd() so that any RNG in the analysis uses a single seed, and you end with a stable, reproducible outcome. You only need set.seed() and setwd() once for a particular analysis.

# *********************************************************************************************************************
# 3. Data and data structures
# + Numeric, integer, character, factor
# + Integer = whole numbers only
# + Numeric = fractions, decimals, whole numbers

14 * 6
1 + 1 + 3
4 / 2 + 1 #=3
4 / (2 + 1) #=1.333

10^3
log10(100) #base 10
log(100) #natural log

# Variable assignment
X = 4 + 6 - (9*2)
A = B = 7

7 = X = Y #only works in 'right-hand' assignment

# **************************************************************************************************************************

# + Vectors - a collection of elements all of the same type;lacks directionality
# + Operations applied to all elements of a vector automatically = greatly reduces need for loops

 
x = c(1,2,3,4) #numeric
x1 = c(1:4) #same as above, shortened; R is inclusive so range is 1 to 4, not 2 to 4 like Python

y = c("Fort Loudoun", "Watts Bar", "Chickamauga", "Nickajack") #character (aka string)

z = factor(c("Little River", "Nolichucky River")) #factor (= categorical variable, no need for arbitrary coding)

x3 = x * 3 #operation applied to every element of the vector

#What can we learn about the vector?
length(x)

length(y)
nchar(y) #includes spaces as characters

is.numeric(x)
is.character(y)
is.factor(z)

#Create a new vector variable with information from previous vector
counts_y = nchar(y)

# ********************************************************************************************************************************
# + Factors ~= categorical variables
# + levels = unique values of that factor variable
# 
# + data.frames, the R way of handling data
# + data.frame = a spreadsheet-like object in R, can be a mixture of variable/data types
# + data.frames, matrices, lists

#Create a dataframe from several vectors
sampleID = as.vector(c(1:15)) #integer

length = as.vector(rnorm(15, 10, 2)) #numeric

species = as.factor(c("crappie", "crappie", "crappie", "crappie", "crappie", "bluegill", "bluegill", "bluegill", "bluegill", "bluegill", "redbreast", "redbreast", "redbreast", "redbreast", "redbreast")) #factor variable

df1 = data.frame(sampleID, length, species)

View(df1) # each column is a variable and each row is an observation

# Accessing a single column/vector/variable:
View(df1$species)

# Accessing portions of the datatframe
# Accessing a single record

df1[2,3] # returns the value in the cell in the second row, third column

# Accessing multiple records
# column 2, rows 4-10
df1[4:10, 2]

# All rows, column 3 (returns all values of the species column)
df1[ , 3]

# Dropping variables with NA in a single column:
df_wq = df_all[-which(is.na(df_all$pH)), ]

# All columns of row 7
obs7 = df1[7, ]
View(obs7)

# Create a new variable in data frame: > data.frame.name$name_of_new_var = 

df1$logLength = log(df1$length) #Which log is this?
View(df1)

# Matrices - rectangluar, every single element must be same type (usually numeric)
# Matrices are sometimes the required format for functions/procedures (instead of dataframes)

A = matrix(length)
A

B = matrix(df1$logLength)
B

C = matrix(species)
C

# Lists - a container that can store data of any type, other lists, etc.
#NOTE: Where I use lists: create an empty list to store output from another operation

results_of_my_procedures = list() #creates an empty list that results can be output to

# a more computationally efficient, R-like method
better_output_storage_vector = vector(mode = "numeric", length = 100) #creates an empty vector w/ 100 open slots


# *************************************************************************************************************************************


# 4. Operations in R
# + General R syntax: object_name = function.name(xvar, yvar, data = dataset, option1 = TRUE, option2 = "orange", ... )
# + Best  practices - variable_name, fxn.name, "=" v. "<-" , camelCase
# + document every procedure in R doc
# + Packages
# + CRAN
# + Install and call a package
# + Functions vs. packages, how to get help with syntax: ?fxn.name, library(help = "package.name")
# + Error messages after operations

 
# installing a package: install.packages('package1', 'package2')
install.packages('MASS', 'ggplot2')

# call the newly installed package to the current session
library(ggplot2)

# sometimes two packages do not get along, so you need to remove one from the session 
detach(package:ggplot2)

# packages are bundles of functions that perform various operations 
# functions are bundles of code that perform specific operations. syntax: function.name(options, settings, etc.)

# methods of accessing help documentation within R: mostly useful for getting syntax of functions right
# accessing help for entire package:
library(help = "ggplot2")
# provides list of all functions with brief descriptions of each one. you can then access help for specific functions by prefacing them with a question mark:
?geom_bar()

error_message = geom_bar(N)
# How I investigate errors: copy error message in to Google, look for applicable solutions!

# Every time R is opened, the package 'stats' is automatically loaded and ready to access. It contains a wide variety of regular statistical tools like ANOVA, linear models, t-tests, etc. See >library(help = 'stats') for the full list of functions.

```
# **********************************************************************************************************************************************

# + Reading data into R


# clear Environment of previously generated and loaded objects: click on broom in top right corner of Environment window

# always make sure your working directory is set = reduces manually typing pathnames = reduced chances for errors
setwd("C:/Users/danwa/OneDrive/Teaching/TNAFS R Workshop")

# there are as many functions to read data into R as there are file types
# two most common filetypes and the functions to read them into R:

#1. .xlsx (MS Excel Workbook). Can be read in one page at a time with read.excel(){readxl}
install.packages(readxl)
library(readxl)
?read_excel
data_xl = read_excel("RforFishSci.xlsx", sheet = 1, col_names = TRUE)
View(data_xl)

#2. .csv (comma-separated). Readable by native R, MS Excel, SAS, and ArcGIS, among others. 
data_csv = read.csv("electrofishing_data")
data_csv = read.csv("electrofishing_data.csv")

# export dataframes from R as new files in working directory with write.csv(){base}
data_csv$logLength = log(data_csv$length)
View(data_csv)
write.csv(data_csv, "electrofishing_logL.csv")
write.csv(data_csv, "electrofishing_logL.csv", row.names = FALSE)
# to save a file in another location outside of working directory, specify full pathname
write.csv(data_csv, "C:/Users/danwa/Documents/electrofishing_logL.csv", row.names = FALSE)
# the .csv file can then be opened by many other programs
 
# *********************************************************************************************************

# 5. Basic statistical operations
# + Mean, median, summary
# + Correlations - Pearson and Spearman
# + T-tests - one way, two way
# + ANOVAs
 
# clear environment of previously generated and loaded objects: click on broom in top right corner of Environment window
setwd("C:/Users/danwa/OneDrive/Teaching/TNAFS R workshop")
set.seed(865)

data = read.csv("electrofishing_data.csv")
View(data)
data$logLength = log(data$length)
data$logWeight = log(data$weight)

# summary statistics
summary(data) # general summary information for each of the columns in the dataframe. first step in data QC

meanLength = mean(data$length)
# meanLength_noNA = mean(data$Length, na.rm = TRUE)
medianWeight = median(data$weight)
rangeWeight = range(data$weight)
rangeWeight
min(data$logLength)
# variance, standard deviation of numeric variables
var(data$length)
sd(data$weight)

# summarize the factor variable species
modeSpecies = mode(data$species)# mode only works on numeric variables

table(data$species) # counts of each species
prop.table(table(data$species)) #proportion of total sample of each species

# correlations
corrSizeP = cor(data$length, data$weight) # Pearson correlation = default
corrSizeS = cor(data$length, data$weight, method = "spearman") # Spearman rank correlation

# T-tests
# 1. one sample - is the mean of a variable >,<, or = some value?
mu.test1 = 250 #value we want to test lengths against
mu.real = mean(data$length)# = 173.517

greater = t.test(data$length, alternative = "greater", mu = mu.test1)
greater

less = t.test(data$length, alternative = "less", mu = mu.test1) # alternative hypothesis = mu.real is less than mu.test1
less

equal = t.test(data$length, alternative = "two.sided", mu = mu.test1) # alternative hypothesis = mu.real is not equal to mu.test1
equal

# two sample t-test. Compare the means of two samples
# create subset of data to compare bluegill data to yellow perch
# only contains bluegill and yellow perch observations
newdat = data[which(data$species == 'bluegill' | data$species == 'yellowPerch'), ] #only rows containing species = bluegill OR yellowPerch, all columns

# test whether means are significantly different between the two species with a two-sided t-test
two_samp = t.test(length~species, data = newdat, var.equal = TRUE)
two_samp
# p value = 0.5864, fail to reject null hypothesis that the mean lengths are significantly different

# ANOVA - compare means of multiple groups
# use original dataset with all species
lengthANOVA = aov(length~species, data = data)
summary(lengthANOVA)

# p-value is significant, so need to test among groups (species)
# e.g. with Tukey post-hoc test
postHoc = TukeyHSD(lengthANOVA, "species")
postHoc
plot(postHoc)
plot(postHoc, las = 1)
# diff = mean difference in the mean lengths of the species pair
# positive diff = species on left is longer than species on right
# p adj - adjusted p-value for the pairwise comparison
# interpretation: two bass species + channel catfish are significantly larger than the panfish

# *************************************************************************************************************************

# 
# 6. Plotting
# + scatterplot
# + histogram/kernel density plot
# + bar plot
# + ggplot2


dat = read.csv("electrofishing_data.csv")

# scatterplot
# length of catfish vs. weight of catfish
# 1. isolate catfish measurements in new dataframe 

catfish = dat[which(dat$species == "channelCatfish"),]

# 2. plot scatterplot of length vs. weight using plot()
plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight")

# add title
plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight")
title(main = "Channel Catfish Length vs. Weight") 
#if you receive an error about plot.new() missing, run title() and plot() commands simultaneously, not in separate lines

plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight")
title(main = "Channel Catfish", sub = "Length vs. Weight")

# ***
# using return to change code readability
plot(catfish$length, 
catfish$weight, 
xlab = "Channel Catfish Length", 
ylab = "Channel Catfish Weight")

title(main = "Channel Catfish Length vs. Weight")

# ***


# change symbol of points with pch = 
plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight", pch = 16)
title(main = "Channel Catfish Length vs. Weight")

# change size of symbols with cex = 
plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight", pch = 1, cex = 1) # default

plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight", pch = 1,  cex = 3)

# change line weight of symbols with lwd = 
plot(catfish$length, catfish$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight", lwd = 5)
title(main = "Channel Catfish Length vs. Weight")

# simple histogram
# weight distribution of bluegill

bluegill = dat[which(dat$species == "bluegill"),]
bg_hist = hist(bluegill$weight, breaks = 6)
bg_hist

# plot several different bin numbers
bg3 = hist(bluegill$weight, breaks = 3)
bg3
bg5 = hist(bluegill$weight, breaks = 5)
bg5
hist(bluegill$weight, breaks = 5)
hist(bluegill$weight, breaks = 11)
hist(bluegill$weight, breaks = 50)
# hist() will use closest bin value to requested while maintaining equal bin widths in the plot

# add fill color to bins
bg_hist_blue = hist(bluegill$weight, breaks = 11, col = "blue")

# kernel density = preferable to histograms, which are dependent on the number of bins used
# 1. calculate density
bg_dens = density(bluegill$weight)
# 2. plot kernel density values
plot(bg_dens)
plot(bg_dens, main = "Bluegill Kernel Density Plot")

# compare multiple kernel density curves with sm.density.compare()
# install.packages("sm")
library(sm)

# create vector of corrected species names for labels
species_labels = factor(dat$species,
            levels = c(1:5),
            labels = c("Bluegill",
                       "Channel Catfish",
                       "Largemouth Bass", 
                       "Spotted Bass", 
                       "Yellow Perch"))

# plot multiple density curves with corrected label
density_compare = sm.density.compare(dat$length, dat$species)
legend("topright", levels(species_labels), fill = 2+(0:nlevels(species_labels))) #weird fill numbering = how colors are referenced (i.e. white = 0)

density_compare = sm.density.compare(dat$length, dat$species)
legend(400, 0.010, levels(species_labels), fill = 2+(0:nlevels(species_labels))) #use XY coordinates for top left corner of legend to place in specific location

# barplot
species_counts = table(dat$species)
species_counts

barplot(species_counts, main = "Species Counts", xlab = "Number of Individuals")

# add x-axis line for visual appeal
barplot(species_counts, main = "Species Counts", xlab = "Number of Individuals", axis.lty = 1)


# barplot of mean weights by species using the aggregate() function

mean_weights = aggregate(dat$weight, by = list(dat$species), FUN = mean, na.rm = TRUE)# object in by = must be a list
mean_weights
barplot(mean_weights[,2], 
names.arg = mean_weights[ ,1], 
main = "Mean Weight by Species", 
xlab = "Species", 
ylab = "Mean Weight (g)")

barplot(mean_weights[,2], 
names.arg = mean_weights[ ,1], 
main = "Mean Weight by Species", 
xlab = "Species", 
ylab = "Mean Weight (g)",
axis.lty = 1)



# ggplot: the Grammar of Graphics, package written by Hadley Wickham.
# Utilizes a standard language for all graphics function
# Standardized language different from the base plotting functions

# *********************************************************************************************************************
# 7. Regression
# + Simple linear regression
# + Multiple regression
# + Logisitic regression
# + Regression diagnostics
 
# use the lm function to generate a simple linear regression
# predictor = channel catfish length, response = channel catfish weight
# step 1 = extract catfish-specific length and weights from sample data
setwd("C:/Users/danwa/OneDrive/Teaching/TNAFS R workshop")
set.seed(865)
data = read.csv("electrofishing_data.csv")

lrData = data[which(data$species == "channelCatfish"),]
View(lrData)

# simple linear regression syntax: lm(response~predictor, data = dataset_name)
simpleLR = lm(weight~length, data = lrData)
simpleLR # only prints the intercept and slope coefficients

summary(simpleLR)
# both terms significant at alpha = 0.10. Very low adjusted r-squared = 0.05348

# View diagnostic plots of simple LR model
plot(simpleLR)

# Plot LR line over scatterplot of input data using plot() command
#plot(x, y, xlab= "x label text", ylab = "y label text")
# step 1: scatterplot of input
# step 2: add regression line over scatterplot with abline() command, which must be entered immediately following  plot() command
# step 3: add model statement + adjusted R squared value
plot(lrData$length, lrData$weight, xlab = "Channel Catfish Length", ylab = "Channel Catfish Weight")
abline(simpleLR)
text(345, 420, labels = "Weight = 486.045 + 0.434*Length, Adj. R-sq = 0.053")


# Multiple linear regression - one response variable, multiple predictors

# create new variable in lrData: new response variable percent fullness' which we will compare to length and weight
lrData$fullness = runif(length(lrData$species), min = 0, max = 1)
View(lrData)
# fullness = percentage value, reported between 0.0 and 1.0.
# two predictors - length and weight, one response - fullness
# generate multiple linear regression using same lm() function

multipleLR = lm(fullness~length + weight, data = lrData)
summary(multipleLR)

# to retrieve just the coefficients of the model:
coeffs1 = coef(multipleLR)
coeffs1

# Logistic regression - binary response variable - glm()
# step 1 - generate two-species data set from 5 species data set
logreg_dat = data[which(data$species == 'bluegill' | data$species == 'channelCatfish'), ]
View(logreg_dat)

logreg = glm(species~length+weight, data = logreg_dat, family = "binomial")
# glm() with family = "gaussian" is same as simple linear regression e.g.:
# glm_ex = glm(weight~length, data = lrData, family = "gaussian") produces same result as line (simpleLR)
summary(logreg)
logreg_coef = coef(logreg)
exp_logregc = exp(logreg_coef)
exp_logregc
# interpretation: bluegills are reference case, for every cm increase in length, there is a 3% increase in chance the fish are catfish and 44% increase in chance for every added g in weight

# **************************************************************************************************************************************************

# 9. Fisheries analyses
# + Writing your own functions
# Writing a function:
# General syntax with command *function()*:

# function.name = function(input1, input2, input3, or no inputs at all...){
# output1 = interior_function1(input1)
# output2 = interior_function2(input2, input 3)
# return(c(output1, output2))
# }



hello.world = function(){
print("Hello World")
}

hello.world()

twice.and.thrice = function(x){
twice = x*2
thrice = x*3
output = list(twice, thrice)
return(c(twice,thrice))
}

twice.and.thrice()

twice.and.thrice(3)

easy.math = twice.and.thrice(5)
easy.math
# 

# + CPUE
# + Length frequency
# + Relative Weight
# + length-specific weight of sampled fish relative to species standard; input length to regression formula
# + Channel cat: Wr = -5.8+3.294*log10(length)
# + Exponentiate (10^X) to convert to (g), then calc (actual wt/ideal wt)*100. Value>100 = fat & happy
# + Von Bertalanaffy curve
# + Catch-curve and mortality
# 

setwd("C:/Users/danwa/OneDrive/Teaching/TNAFS R workshop")
set.seed(865)

dat1 = read.csv("electrofishing_data.csv") # lengths and weights for 5 species caught during sampling event
spot = read.csv("spot.csv") # length-at-age data for a sample of marine Spot
lmb = read.csv("largemouth.csv") # counts at age for a sample of FL largemouth bass

# CPUE function
cpue.fun = function(catch, effort){
cpue = catch/effort
print(paste0("CPUE is: ", cpue, " fish per hour of effort"))
# return(cpue)
}
# Using CPUE function, calculate bluegill CPUE in catch per hour from dat1 given effort

sp.counts = table(dat1$species)
sample = 600 #seconds of shocking
effort = sample/3600 #seconds per hour
sp.counts
catch = 165 #bluegill caught in 600s shocking

bluegill.cpue = cpue.fun(catch,effort)
bluegill.cpue


# length-frequency histogram

bins = as.integer(sqrt(nrow(spot))-1) # subtract 1 since bins command is number of divisions, not number of bins
length_freq = hist(spot$tl, breaks = bins)
length_freq1 = hist(spot$tl, 
        breaks = bins,
        plot = TRUE, #specifies desired output is a plot and not summary of plot
        xlim = c(0, 15),
        xlab = "Spot Length (TL mm)"
)


# relative weight - Wr, length-specific weight of sampled fish compared to a species standard
# provides information on relative condition of sampled fish
# E.g.: channel catfish species standard Wr = -5.8 + 3.294*log10(length)

# 1. create channel-catfish only data
catfish = dat1[which(dat1$species == "channelCatfish"), ]

# 2. calculate log10 of channel catfish lengths
catfish$logLength = log10(catfish$length)
View(catfish)

#3. calculate relative weight
catfish$Wr = (-5.8 + (3.294*catfish$logLength))# this equation will change for different species
View(catfish)

#4. exponentiate Wr
catfish$Wr.exp = 10^catfish$Wr
View(catfish)

#5. calculate ratio of exponentiated ideal relative weight to actual weight to assess stock
catfish$ratio = ((catfish$weight/ catfish$Wr.exp)*100)
View(catfish)

#6. Assess how many catfish weight ratios are > 100 (indicating fat and happy fish)
catfish$weightcat[catfish$ratio >= 100] = 1
catfish$weightcat[catfish$ratio < 100] = 0
table(catfish$weightcat)
# 39 channel catfish are relatively heavier for their length than species standard (89%)
avg.rel.weight = mean(catfish$ratio)
avg.rel.weight # = 1091.54 >>> 100, population is not prey-limited


# catch curve analysis
# use regression of N at age to determine instantaneous mortality (Z), annual total survial (S), and annual total mortality (A)

View(lmb)

# 1. Calcluate ln(N) for each age class
lmb$lnN = log(lmb$N)

# 2. plot age vs. lnN
plot(lmb$age, lmb$lnN)

# 2.1 Drop age-0 observations, which are not recruited to the gear
lmb_rec = lmb[which(lmb$age>0), ]
View(lmb_rec)

#3. calculate regression of age vs. lnN
cc_reg = lm(lnN~age, data = lmb_rec)
cc_reg # coefficient for age = -0.5958. Z = negative(coefficient)

#4. plot age vs. lnN with regression equation
plot(lmb_rec$age, lmb_rec$lnN)
abline(cc_reg)
text(8,5, "lnN = -0.5958*age + 6.4247")

plot(lmb_rec$age, lmb_rec$lnN)
abline(cc_reg)
text(8,5, "lnN = -0.5958*age + 6.4247\n Z = 0.5958")

S = exp(-0.5958)
A = 1 - S
A
# total annual mortality is around 45%


# Von Bertalanffy growth function (VBGF)
#0.1. describe mean length at age for Spot dataset
meanTL_ages = as.data.frame(aggregate(spot[, 1], list(spot$age), mean))
View(meanTL_ages)
names(meanTL_ages)[1] = "Age"
names(meanTL_ages)[2] = "TL"
plot(meanTL_ages) #visualize mean lengths at age

#0.2 Fit VB models using package FSA
#0.3 install.packages("FSA")
library(FSA)

# Fitting 'typical' VB function:
# 0.1 Initialize the VB function. State which function you want to solve
vb = vbFuns("typical")
# Other available functions to solve in package
vbFuns("BevertonHolt") #synonymous with 'typical'
vbFuns("original")
vbFuns("vonBertalanffy")#synonymous with 'original'

#1. Find starting values for K, t0, Linf
VBstart = vbStarts(tl~age, data = spot, type = "typical", meth0 = "poly", plot = TRUE)
# meth0 = "poly" = fit 2nd degree polynomial fit to mean length-at-age to calculate t0 and L0
# meth0 = "yngAge" = starting values for t0 and L0 found algebraically 
# see ?vbStarts for more information 
VBstart

# manually specify starting values
VBstartFIXED = vbStarts(tl~age, data = spot, type = "typical", meth0 = "poly", plot = TRUE, fixed = list(Linf = 15, K = 0.3, t0 = 0))

#2. Fit typical VB function using generated start values
# uses nonlinear least-squares estimation to calculate parameters
VBfit = nls(tl~vb(age, Linf, K, t0), data = spot, start = VBstart)
summary(VBfit, correlation = TRUE)
VB_coeffs = coef(VBfit)
VB_coeffs

VBcoeffs_confidence = confint(VBfit)
VBcoeffs_confidence

# diagnostic plots
residPlot(VBfit)

#3. Visualize VB fit with data
clr = rgb(0,0,0,0.5) #establish grayscale color gradient
plot(tl~age, data = spot, xlab = "Age", ylab = "Total Length", pch = 16)
curve(vb(x,VB_coeffs), from = 0, to = 4, n=500, lwd = 2, col = "red", add = TRUE)

#4. Make predictions for TL with age using VBGF equation
ages = c(5:10) #ages to predict lengths for
TLpreds = predict(VBfit, data.frame(age = ages))
TLpreds

# *******************************************************************************************************


# 10. Classification
# + logistic regression already covered
# + principal components analysis
# + linear discriminant analysis
# + CART


setwd("C:/Users/danwa/OneDrive/Teaching/TNAFS R workshop")
set.seed(865)

dat = read.csv("BarentsFish.csv")
View(dat)
# x1 and x2 = random, uninformative variables
dat = dat[1:89, ] #remove empty last row
View(dat)
summary(dat)
# Need categorical dependent variable
# We can make a categorical variable based on the quality of the cod catch (in numbers caught) per site
# Low = first quartile or less, < 53 cod caught
# Average = 53 - 295 cod caught
# High = >295 cod caught

dat$codCatch[dat$Cod < 53] <- "Low"
dat$codCatch[dat$Cod >= 53 & dat$Cod <= 295] <- "Average"
dat$codCatch[dat$Cod > 295] <- "High"
View(dat)
table(dat$codCatch)

# Objective: produce infomrative linear combinations of variables relating dependent categorical variable of cod catch quality to environmental variables depth, temperature, x1, x2

# PCA
pca.preds = dat[ , c(4:5, 9:10)] #isolate input variables
View(pca.preds)
pca.y = dat[, 11] #isolate dependent variable

pca1 = prcomp(pca.preds, scale = TRUE, center = TRUE) #recommended always scaling and centering PCA predictors
summary(pca1) #2 Principal components = 70% variation explained (good dimension reduction)
plot(pca1, type = "l") #elbow at either 2 or 3 PC's

# check PC loadings on predictors
print(pca1)

biplot(pca1)
#NOTE: More advanced biplot generation in 'ggbiplot' package


# linear discriminant analysis
# install.packages("MASS")
library(MASS)

lda.fit=lda(pca.y~dat$Depth+dat$Temp+dat$x1+dat$x2,
prior = c(1,1,1)/3, #specify uniform prior probabilities of site membership in cod quality classes
CV = FALSE #  CV=TRUE fits LOOCV)

lda.fit
# proportion of trace = between class variance explained by linear discriminant terms
plot(lda.fit) # plots the first two linear discriminant directions
plot(lda.fit, dimen=1, type="both") # distribution of y-var categories along LD1


# Classification and Regression Trees
# classification tree = categorical y-var
# regression tree = continuous y-var
# install.packages("rpart", "rpart.plot")
library(rpart)
library(rpart.plot)

tree = rpart(codCatch~Depth+Temp+x1+x2, 
             data = dat, 
             method = "class")
tree.full = rpart.plot(tree, main = "Full Classification Tree")
# Each node displays predicted probability for each class in that node (3 values = 3 classes) and percent of total observations in that node
# intensity of color corresponds to purity of the node

# CART algorithms typically overfit models, so it is helpful to 'prune' the trees
# Reduce number of splits to increase interpretability at a slight cost of classification accuracy

plotcp(tree)
# there is a minimum relative error at tree size = 3 splits, complexity parameter = 0.066
# so we can recreate the tree stopping at only 3 splits while minimizing our relative classification error

prune.tree = rpart(codCatch~Depth+Temp+x1+x2,
                   data = dat,
                   method = "class",
                   cp = 0.066)

tree.pruned = rpart.plot(prune.tree, main = "Pruned Classification Tree")
# interpretation: highest cod catches occur deeper than 394 meters at temps greater than 3.4C in Berents Sea

# ********************************************************************************************************

# 
# 11. Final Quiz  
# 1. Start with an empty environment, new R session
# 2. Set working directory
# 3. Load the channel catfish length-weight-age sheet from 'RforFishSci.xlsx' as a dataframe named 'catfish'
# 4. Load the crappie TL sheet from 'RforFishSci.xlsx' as a separate dataframe named 'crappie'
# 5. Load the spot length-at-age data from the 'spot.csv' file as a dataframe named 'spot'
# 6. Generate a length-frequency histogram from the crappie TL data using the default settings
# 7. Make a new length-frequency histogram by changing at least two plot settings
# 8. Calculate the ratio of ideal to actual weights for the channel catfish 
# 9. Classify whether each catfish's weight ratio is greater or less than 100
# 10. Identify the proportion of channel catfish who exceed a weight ratio of 100
# 11. Calculate the average relative weight for the catfish population
# 12. Export channel catfish data frame with intermediate steps to new csv file called 'catfish_new'
# 13. Solve an original VB function for the spot length at age data
# 14. Plot the spot length at age data and VB growth curve
# 
# 
# 
# 12. Additional Resources
# + Google: " 'function/operation' and/or 'error message' r "
# + www.r-tutor.com
# + www.statmethods.net
# + www.stats.stackexchange.com
# + http://derekogle.com/fishR/courses
# + 'R for Everyone' first ed. - Lander 2013, $25 on Amazon. Second ed. releasing April, $40
# + 'Data Analysis and Graphics Using R' - Maindonald and Braun 2003, hardcover, ~$100 on Amazon

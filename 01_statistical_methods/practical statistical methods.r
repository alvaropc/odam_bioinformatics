#Practise Alvaro Ponce Cabrera
#Before start the exercises is necesary to set the work directory and load the data.
setwd("~/Desktop/omicscopia/BIOINFORMATICA/practises/practise 1") #Set work directory
rm(list=ls()) #remove all objects in the global environment
bladder<-read.table("bladder.16.txt", header=T, sep="") #load the dataset
attach(bladder) #Attach the data make easier the access
head(bladder) 

###########################################################
#1. Convert variables “y” and “gender” to factor variables.
#2. Define the labels for variables “y” and “gender”.
###########################################################

bladder.y<- factor(bladder$y, 
                   levels=0:1, labels=c("control", "case")) #Y factor variable, called bladder.y
bladder.y
bladder.gender<- factor(bladder$gender, 
                        levels=1:2, labels=c("male", "female")) #Gender factor variable, called bladder.gender
bladder.gender

##################################################################################################
#3. Build a frequency table for “y” that contains, both, the absolute and the relative frequencies. 
#The same for “gender”.
################################################################

freq.y<-table(bladder.y)   #Frequency table for the variable y
freq.y
freq.rel.y<-prop.table(freq.y) #Relative frequency table for y
freq.rel.y
cbind(freq.rel.y,freq.y)  #Creation of a table with relative and
                          #frequency tables data
freq.gender<-table(bladder.gender)  #Frequency table for gender
freq.gender
freq.rel.gender<-prop.table(freq.gender) #Relative frequency table for gender
freq.rel.gender
cbind(freq.rel.gender,freq.gender)  #Creation of a table with relative and
                                    #frequency tables data

###########################################################
#4. Plot gene1 levels as a function of “gender”.
############################################################

#Relative frequency table for y

plot(bladder.gender, gene1, col=c("red","green") ) 
#Plot of gene1 levels as function of gender. 
#Red is for male, green for female.


####################################################################################
#5. Test the normality of the expression levels of the 20 genes (use function apply). 
#How many genes are not normally distributed and which are these genes?
#####################################################################################

s.test<-(apply(bladder[,4:23],2,shapiro.test)) #Test for the normality of each 20 genes expresion level 
names(s.test)

#To obtain the names of the genes that are positives or negatives in the shapiro.test (positive=p.value<0.05)
positive.test=NULL
negative.test=NULL
for (i in 1:20)
{                                                       #eval() evaluate of gene"i" inside the s.test enviroment
  gen.test<-eval(parse(text=paste0("gene",i)),s.test)   #parse() give an object to evaluate to eval(), and paste0() provide 
  gen.test$data.name<-paste0("gene",i)                  #this name to parse(). Also we save the name of each gene
  if (gen.test$p.value<0.05)
  {
    positive.test<-c(gen.test$data.name,positive.test)    #The names are saved into 2 variables, positives results
                                                          #in positive.test and negatives in negative.test
  }
  else
  {
    negative.test<-c(gen.test$data.name, negative.test)
  }
}

positive.test #These genes are not normaly distributed, because p.value<0.05
negative.test


####################################################################################
#6. Test whether mean expression levels of gene1 and gene2 are equal. 
#Note: this is a test for equality of means with paired samples. 
#This test is performed by computing the difference of the two variables (gene1-gene2) 
#and testing whether the mean of the difference is equal to zero.
#####################################################################################

#T.test for the equality of two means with paired samples
d<- gene1-gene2
test<-t.test(d,mu=0)
test$p.value
#we do not reject H0

#####################################################################################
#7. Test if the mean expression levels of gene1 are equal between cases and controls
####################################################################################
#First we need to test if the variances are equal or not
var.test(gene1~bladder.y)

#Variances are not equal, so we compute a t.test with var.equal=F and show the p.value of the test.
t.test(gene1~bladder.y,var.equal=F)$p.value



###############################################################################
#8. Obtain a 95% bootstrap confidence interval for the 20th percentile of gene1 
###############################################################################

#open the library that contain the function we are going to use.
library(boot)

#An argument of boot() function is a function of one statistic. 
#In this exercise we need to obtain a quantile, so we will use quantile as our statistic
Sq<-function(data,i)
{
  q<-quantile(data[i],0.2)  #Sq() function call a data i times to calculate the 0.2 quantile
  return(q)                 #save this data un q variable and return q.
}

Sq(gene1) #We can check the function works

boot.outq<-boot(gene1, Sq, R=1000)    #We apply the boot function to our data and save it in a variable

boot.ci(boot.outq, conf =0.95, type = "perc") #boot.ci receive a boot object and calculate a confidence interval
                                              #here it computes a percentil confidence interval (type="perc)


########################################################################################################
#9. Contrast using a permutation test whether the 20th percentile of gene1 for male and female are equal
#######################################################################################################

#Generate of the permutation test, first we obtain the observable data.
SFunction<-function(x){
    
  male.q<-quantile(x[bladder.gender=="male","gene1"],0.2)
  female.q<-quantile(x[bladder.gender=="female","gene1"],0.2)

  S<-male.q-female.q
  return(S)
}

SFunction(bladder)

Sobs<-SFunction(bladder)

SNull<-function(x){
  s.gender<- sample(bladder.gender)
  male.q<-quantile(x[s.gender=="male","gene1"],0.2)
  female.q<-quantile(x[s.gender=="female","gene1"],0.2)
  S<-male.q-female.q
  return(S)
}
SNull(bladder)

#Then we generate the null hypothesis curve or null permutational distribution using our data....
#and running 1000 permutation sums
# Our H0 is that q1-q2=0

  S.nullpermutdist <- replicate(1000,SNull(bladder))
  S.nullpermutdist
  # Plot the histogram of the permutational null distribution 
hist(abs(S.nullpermutdist))
#We need the absolute value of variables in order to not ignore values below our Sobs 
#in the calculation of p.value

# we compute the permutation test p.value for the  alternative H1:
#p-value=P(S.nullpermutdist>Sobs)

p.value<-(sum(abs(S.nullpermutdist)>abs(Sobs))+1)/(length(S.nullpermutdist)+1) 
p.value 
#We can reject the null hypothesis, but p.value is not very smaller than 0.05, so maybe we should take
#some care with this result if we work with its.



###########################################################################################
#10. Perform a nonparametric test for association of gender and the risk of disease. 
#Provide the OR (change the levels of “gender” if necessary in order that the given OR is larger than 1)
###########################################################################################

#We open the library epitools in order to use after that the function oddsratio()
library("epitools")
#oddsratio must calculate cases/controls and the final value should be bigger than 1 to be easier 
#to understand, that's why the argument rev="both" is used, to reverse columns and rows. 
odds<-oddsratio(table(bladder.gender,bladder.y),rev="both")
odds
odds$measure
#Oddsratio Value
odds$p.value
#P.value is bigger than 0.05, so the oddsratio value is not signifcant.

##############################################################################
#11. Explore for possible relationship between methylation and gene expression.
##############################################################################
#FIrst, it is necesary to know if methyl is normally distributed or not.
shapiro.test(methyl)
#Apply run cor.test over each gene from bladder.
c.test<-apply(bladder[,4:23],2,function (x) cor.test(x,methyl,method="spearman"))

#Using the same for loop that excercise 5 we obtain the names of the genes which p.value is <0.05 or >0.05
positive.rho=NULL
negative.rho=NULL
for (i in 1:20)
{                                                       
  gen.rho<-eval(parse(text=paste0("gene",i)),c.test)   
  gen.rho$data.name<-paste0("gene",i)                  
  if (gen.rho$p.value<0.05)
  {
    positive.rho<-c(gen.rho$data.name,positive.rho)    
  }
  else
  {
    negative.rho<-c(gen.rho$data.name, negative.rho)
  }
}

positive.rho
#this genes have correlation with methyl, lets check the level of this correlation (rho)
c.test$gene16$estimate #low correlation rho=0.08
c.test$gene4$estimate  #high correlation rho=0.7

negative.rho


########################################################################################
#12. Identify genes that are related to the risk of bladder cancer using a multivariate
#logistic regression model with stepwise variable selection. Denote the selected model as “best.model”. 
#Interpret the obtained model. 
##################################################################################

#Generation of the model using all variables.
model<-glm(y~gene1+gene2+gene3+gene4+gene5+gene6+gene7+gene8+gene9+gene10+gene11+gene12+gene13+
     gene14+gene15+gene16+gene17+gene18+gene19+gene20, family=binomial())
#We compute for the model with the less AIC using step() function
step(model,direction="both")$coefficients

#the best model is this one.
best.model<- glm(y~gene7+gene11+gene15+gene16+gene19)

summary(best.model)
###############################################
#13. Analyze the classification ability of “best.model”
#(ROC curve and AUC) according to the following schemes:
######################################################


#####################################
#a. Apparent validation of “best.model” using the same data that was used for model building.


#Open the library
library(ROCR)

#Save the linear.predictros of best models in the variable lp
lp<-best.model$linear.predictors
hist(lp)
 #prediction is using in every start of ROCR evaluation. 
#prediction standarize the data and clasify it with the labels given un the second argument
pred <- prediction(lp, bladder.y)
pred
#performance representate a prediction object using the axis that we indicate, "tpr" and "fpr" in this case
perf <- performance(pred, "tpr", "fpr" )
plot(perf)
#It give us a curve done by all the standarized linear.predictors points
#but the value is lower than 0.5, so we need to do the complementary
#########complementary

pred <- prediction(1-lp, bladder.y)
pred

perf <- performance(pred, "tpr", "fpr" )
plot(perf)

title("ROC curve for best.model; bladder data")

performance(pred,"auc")

##################################################
## b. Crosvalidation with k = 5 for “best.model”.

K <- 5
n <- nrow(bladder)  
bp <- c(0, ceiling((1:K)*n/K))   

bp
# index permutation
index <- sample(rownames(bladder))

# CV
pred <- NULL
for(i in 1:K){
  indTest <- index[(bp[i]+1):bp[i+1]] 
  indTrain <- index[!(index %in% indTest)]  
  model.i <- glm(y~gene7+gene11+gene15+gene16+gene19, data=bladder[indTrain,], family=binomial())
  pred.i <- predict(model.i, newdata=bladder[indTest, ])
  pred <- c(pred,pred.i)
}  


# pred contains the concatenated linear predictor of the permuted individuals
# pred.y contains the ordered linear predictions
pred.y <- pred[rownames(bladder)]
#ROC curve
prediction <- prediction(pred.y, bladder.y)
perf <- performance(prediction, "tpr", "fpr" )
plot(perf)
performance(prediction,"auc")

###complementary

prediction <- prediction(1-pred.y, bladder.y)
perf <- performance(prediction, "tpr", "fpr" )
plot(perf)
title("ROC curve for best.model; crosvalidation")

performance(prediction,"auc")

#################################################################################################
#c. Explain why both the "a" and "b" overestimate the classification accuracy of “best.model”. 
#Suggest a scheme for determining an unbiased estimate of the classification accuracy of “best.model”.

#Both validations overestimate the classification accuracy of “best.model” because they use 
#the same data for build the model and test it.
#To obtain a better validation we need our training being different to the test data.

########################################################################################
#14. For each variable “genei”, i=1, …, 20, define a new factor variable 
#“levels.genei” with values “high” if gene levels are positive or zero and “low” if
#gene levels are negative.
#########################################################################################

#Creation of the new factor variables running apply inside bladder[,4:23]. It generates a matrix object
levels.gene.matrix<-apply(bladder[,4:23],2, function(x) factor(x>=0,levels=c(TRUE,FALSE),labels=c("high","low")))
class(levels.gene.matrix)
df.levels.gene<-data.frame(levels.gene.matrix)
#The matrix is converted into a data frame.
names(df.levels.gene)
names(df.levels.gene)<-gsub("gene", "levels.gene", names(df.levels.gene))
#Change of the names into the data frame, gsub() change the first argument to the second inside an object
#given in the last argument, in this case, our data frame.
names(df.levels.gene)
#So now we have our variables "levels.genei" but inside a data frame.
attach(df.levels.gene)

levels.gene20
#Checking the data frame
names(df.levels.gene)

####################################################################################
#15. Perform a nonpametric test of association between each variable “levels.genei” 
#and the risk of disease. Adjust the p-values for multiple testing according to an fdr 
#threshold equal to 0.1. Interpret the results.
####################################################################################

#Apply of the chisq.test. Chisq.test test the the association between two categorical variables
chi.test<-apply(df.levels.gene,2,function (x) chisq.test(table(x,bladder.y)))
names(chi.test)
#For loop to create a vector with the p.values of each leves.gene
p.vector=NULL
for (i in 1:20)
{
  gen.chi<-eval(parse(text=paste0("levels.gene",i)),chi.test) 
  p.vector<-c(gen.chi$p.value, p.vector)
}
p.vector
#with p.vector now we can adjust by the fdr method the p.values.
p.vector.adjust<-p.adjust(p.vector,method="fdr")
which(p.vector.adjust<0.1)
#We can obtain the positons of the levels.genes that are true for this logical condition (p.vector.adjus<1)
#The levels.genes are ordered from 1 to 20 so this positions give us the numbers of levels.genei 
#which have a p.value adjusted lower tha 0.1. In this case the levels.gene are 5,6,10 and 14.
#All of them reject H0 and are related with the risk of disease (badder.y)

###################################################################################################
#16. Using the last 30 variables corresponding to gene expression levels in three different pathways, 
#perform a clustering analysis (hierarchical and k-means) and explore groups of genes that have similar 
#expression. Explain the results. 
###################################################################################################


#Hclust generate clusters from the distances between the variables given.
#The trasposition was done in order to generate a better plot that really provide us the relation
#between each pathway.
hcgenes<-hclust(dist(t(bladder[,25:54]),method="euclidian"),method="average")
plot(hcgenes)

#We can see 3 mains clusters path1.i path2.i and path3.i but some path are not related in their group.
#also, path3.i is very relate with path1.i they are in one same big cluster.



#Using Kmeans we create function other way to see the relationship between the pathways by 
#coercing the data to be distributed in just k clusters, in this case, 3.
cl<-kmeans(t(bladder[,25:54]),3,nstart=10)
cl
#inside this obeject you can access to the cluster distribution.
cl$cluster
#sort is just for organize the data.
sort(cl$cluster)
#Again we can see the same distribution of clusters that we explained before.
#The explicacion is more extended in the word document.

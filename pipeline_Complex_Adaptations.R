###############################################################################
##########              Pipeline Growth dynamics and rates          ###########
##########                        Philipp Hummer                    ###########
##########                          31/08/2022                      ###########
###############################################################################
############################### OBJECTIVES ####################################
# 
############################# INSTRUCTIONS ####################################
# 
########################## SAFETY INSTRUCTIONS ################################
# /!\  ######   /!\
############################### LIBRARIES #####################################

# Packages requested 
list.of.packages <- c("tidyverse","ggpubr","MASS")
# Checks if any requested package is installed 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
## Checks for package updates (optional)
if (is.null(old.packages())) {
  update.packages(ask = FALSE)
}
# Loads requested packages
for (i in list.of.packages){
  library(i, character.only = TRUE)
}

#####Alternative section for compiling multiple CSV files into one table#######
############################## PARAMETERS #####################################

### Specific parameters ###

Path <- "<C:/Path/>" # Path to your working directory
FilesPath <- "<C:/Path/>" 
# Path to the folder containing all the cleaned CSV files AND the template files

########################## Bring files together ################################

Files <- list.files(FilesPath) # Take the names of the files in the folder
Files <- Files[which(substring(Files, 1, 1) != "~")] # Remove the hidden files
Files <- Files[which(substr(Files, nchar(Files)-4+1, nchar(Files)) == ".csv")] # Keep only CSV files

#these template files must be in your directory
template_all<-read.csv("~/all cultures.csv",sep=";",header=TRUE)
template_red<-read.csv("~/reduced cultures.csv",sep=";",header=TRUE)

Total <- c()

for (f in 1:length(Files)) { # For each files in the folder
  
  CSV <- read.csv(file.path(FilesPath, Files[f]), skip = 7, header = FALSE, sep= ",") # Read the CSV file
  
  d <- as.numeric(substr(Files[f], 4,6)) # Get the experiment day from the file name
  oldornew <- substr(Files[f], 8,8)# Get the experiment day from the file name
  w <- substr(Files[f], 11,12)# Get the experiment week from the file name

  stock <- c()
  
  #depending on the day, one of the two templates will be used
  if (d<=28) {reduction<-0} else {reduction<-1}
  if (reduction==TRUE) {cultures<-template_red} else {cultures<-template_all}
  if (oldornew=="c") {cultures$Continuation<-0} else {cultures$Continuation<-1}
  if (reduction==1) {rows<-56} else {rows<-80}
  
  #fill up the template with your the data from the CSV files
  cultures$DunaAlive[1:rows]<-as.numeric(CSV[,10])
  cultures$DunaDead[1:rows]<-CSV[,13]
  cultures$Day<-d
  cultures$Week<-w
  if (reduction==TRUE) {treatments<<-as.list(cultures$Treatment.ID[1:56])} else 
    {treatments<<-as.list(cultures$Treatment.ID[1:80])}
  for (z in treatments){
    # print(treatments)
    cultures$DunaAlive[cultures$Sample.ID==z]<-mean(cultures$DunaAlive[cultures$Treatment.ID==z],na.rm=TRUE) 
    #the standard error plots can most easily be produced by replacing the mean() with sd()/sqrt(3)
    cultures$DunaDead[cultures$Sample.ID==z]<-mean(cultures$DunaDead[cultures$Treatment.ID==z],na.rm=TRUE)
  }
  
  complexad<-rbind(Total, cultures) # Add all this informations in a big dataset with all the informations
  
}


#####################shortcut using the compiled CSV############################

complexad<-read.csv("<C:/Path/>Complexad_report26.8.csv",sep=",",header=TRUE)
template_all<-read.csv("<C:/Path/>all cultures.csv",sep=";",header=TRUE)
template_red<-read.csv("<C:/Path/>reduced cultures.csv",sep=";",header=TRUE)

################################################################################
#restrict the data to a subset that only contains the means of the relevant measurements for the analysis

complexad1<-subset(complexad,Week=="01"|Week==1)
l<-length(template_all[,1])
day_week<-c(rep(0,l),rep(1,l),rep(2,l),rep(3,l),rep(4,l),rep(7,l),rep(8,l),rep(9,l),rep(10,l),rep(11,l))
dim(day_week)<-c(length(day_week),1)
complexad1$day_week<-day_week

complexad2<-subset(complexad,Week=="02"|Week==2)
day_week<-c(rep(0,l),rep(1,l),rep(2,l),rep(3,l),rep(4,l),rep(7,l),rep(8,l),rep(10,l),rep(11,l))
dim(day_week)<-c(length(day_week),1)
complexad2$day_week<-day_week

complexad6<-subset(complexad,Week=="06"|Week==6)
lr<-length(template_red[,1])
day_weekr<-c(rep(0,lr),rep(1,lr),rep(2,lr),rep(3,lr),rep(4,lr),rep(7,lr),rep(8,lr),rep(9,lr),rep(10,lr),rep(11,lr))
dim(day_weekr)<-c(length(day_weekr),1)
complexad6$day_week<-day_weekr

complexad10<-subset(complexad,Week=="10"|Week==10)
complexad10$day_week<-day_weekr

mweeks<-rbind(complexad1,complexad2,complexad6,complexad10)#mweeks is later used again

complexad<-mweeks%>% 
  mutate(Week=as.numeric(Week))
complexad<-subset(complexad,Replicate=="pooled") #this means we only keep the means
complexad<-subset(complexad,Salinity!="3,5") 
#depending on what kind of plots or analysis you want, you may not want to execute this line


############################PLOTTING############################################
#####################growth dynamics plots######################################

l <- length(unique(complexad$Temperature))*length(unique(complexad$Culture))
list.plot<-c()
list_plot <- vector(mode = "list", length = l) 
#the plots will be stored in this list and can then by addressed with list_plot$
i <- 1
for(t in unique(complexad$Temperature)) {
  for (s in unique(complexad$Culture)) {
    df <- subset(complexad, Culture == s & Temperature == t)
    if (s=="mix") lab<-paste("Population dynamics of the ",s," at ",t,"?",sep="") else 
      lab<-paste("Population dynamics of strain S",s," at ",t,"?",sep="")
    
    plt<-ggplot(df, aes(day_week,DunaAlive,color=factor(Week,ordered=TRUE)))+
      geom_point(size=2, shape=16)+
      geom_line(aes(group=interaction(Sample.ID,Week),linetype=Salinity))+
      theme_bw(base_size=14)+
      scale_colour_manual( values=c("red","orange","green","cyan"))+#taste the rainbow
      scale_y_log10(limits = c(1e2,1.5e6))+
      theme(aspect.ratio=1)+
      labs(title=lab,x="Days",y="se of the concentration of living cells/mL",color="Week")+ 
      theme(plot.title = element_text(size=8))
    plt
    list_plot[[i]] <- plt
    names(list_plot)[i] <- paste("complexad", t, s, sep = "_")
    i=i+1
  }
}

figurenew_ad<-ggarrange(list_plot$complexad_24_1C,list_plot$complexad_36_1C,list_plot$complexad_24_3,
                        list_plot$complexad_36_3,labels=c("a","b","c","d"),ncol=2,nrow=2)
figurenew_eh<-ggarrange(list_plot$complexad_24_13,list_plot$complexad_36_13,list_plot$complexad_24_mix,
                        list_plot$complexad_36_mix,labels=c("e","f","g","h"),ncol=2,nrow=2)
#I plot these two parts separately because that makes it easier to handle them but it could be one big ggarrange too
figurenew_s11<-ggarrange(list_plot$complexad_24_11,list_plot$complexad_36_11,labels=c("a","b"),ncol=2,nrow=2) 

#########################growth rate analysis###################################
complexad<-mweeks%>% 
  mutate(Week=as.numeric(Week))

complexad<-subset(complexad,Replicate==1  |  Replicate==2  |  Replicate==3)#gets rid of the means
complexad<-subset(complexad,Culture !="11" & Salinity !="3,5")
#depending on your analysis you might want to adapt this

#########Use this code if you want to analyse 3.5M for S3#######################
#complexad24<-subset(complexad,Temperature==24 & Culture !="11" )
#complexad36<-subset(complexad,Temperature==36 & Culture !="11" )
#complexad36_2M<-subset(complexad36,Salinity==2)
#complexad36_3.5M<-subset(complexad36,Salinity=="3,5" & Culture=="3")
#complexad36_4M<-subset(complexad36,Salinity==4 & Culture!="3")
#complexad<-rbind(complexad24,complexad36_2M,complexad36_3.5M,complexad36_4M)


#####################get the growth rates and standard errors###################
predictors_all<-c()#these empty variables will be filled by the loop
predictors_all_factor<-c()
for(t in unique(complexad$Temperature)) {
  for (s in unique(complexad$Culture)) {
    df_1 <- subset(complexad, Culture == s & Temperature == t)
    for (r in unique(df_1$Salinity)) {
      df <- subset(df_1,Salinity == r) %>% #creates the subset and changes the formats of its elements by row
        mutate(Sample.ID = as.factor(Sample.ID), Replicate = as.factor(Replicate), Culture = as.factor(Culture),
               Salinity= as.factor(Salinity), Temperature  = as.factor(Temperature),
               Treatment.ID = as.factor(Treatment.ID))
      
      #create empty data frames to fill
      number_weeks<-4 #the number of the weeks that we want to calculate the growth rate of, not the total number!
      Daycol<-rep(rep(c(0,1,2,3,4,7,8,9,10,11),number_weeks))
      Weeks<-rep(c("1","2","6","10"),each=10) #3 replicates x 10 days = 30
      predictors2<-cbind(Daycol,Weeks)
      predictors2<-as.data.frame(predictors2)
      names(predictors2) <-c("day_week","Week")
      predictors<-c(1,2,6,10)
      dim(predictors)<-c(number_weeks,1)
      empty<-rep(NA,number_weeks*8)
      dim(empty)<-c(number_weeks,8)
      predictors<-cbind(predictors,empty)
      colnames(predictors)<-c("Week","r0","se_r0","Condition","Temperature","Culture","Salinity","Start_exp_growht","End_exp_growth")
      predictors_l<-as.data.frame(predictors)
      
      #create the model and fill its predicts into the template
      model<-glm.nb(DunaAlive ~ as.numeric(day_week)*as.numeric(Week) ,data=df)  #Week either as numeric or factor
      response_predict2<- predict(model, newdata= predictors2, type="link", se.fit = TRUE)
      response_predict2<- as.data.frame(response_predict2)      
      names(response_predict2) <- c("fit2", "se_fit2", "residual_scale2")
      prediction2 <-cbind(predictors2,response_predict2)
      
      #define the boundaries of the exponential growht phase, this can be adapted with if conditions for the treatent conditions
      start_phase<-1
      if (unique(df$Temperature) == 24 & unique(df$Culture) != "1C") end_phase<-4 else end_phase<-7
      duration_phase<-end_phase-start_phase
      
      #fill our templates with all the values we want, including the growth rate
      predictors_l$r0  <- (prediction2[prediction2$day_week == end_phase,]$fit2 - prediction2[prediction2$day_week == start_phase,]$fit2)/duration_phase  
      predictors_l$se_r0 <- sqrt((prediction2[prediction2$day_week == end_phase,]$se_fit2)**2 + (prediction2[prediction2$day_week == start_phase,]$se_fit2)**2)/duration_phase
      predictors_l$Condition<-paste(s,t,r,sep="_")
      predictors_l$Temperature<-t
      predictors_l$Culture<-s
      predictors_l$Salinity<-r
      predictors_l$Start_exp_growht<-start_phase
      predictors_l$End_exp_growth<-end_phase
      predictors_l$aic<-model$aic
      predictors_all<-rbind(predictors_all,predictors_l)
      
      #now the same thing again with a model that uses as.factor(Week)
      predictors<-c("1","2","6","10")
      dim(predictors)<-c(number_weeks,1)
      predictors<-cbind(predictors,empty)
      colnames(predictors)<-c("Week","r0","se_r0","Condition","Temperature","Culture","Salinity","Start_exp_growht","End_exp_growth")
      predictors_l<-as.data.frame(predictors)
      model_f<-glm.nb(DunaAlive ~ as.numeric(day_week)*as.factor(Week) ,data=df)  #Week either as numeric or factor
      response_predict2<- predict(model_f, newdata= predictors2, type="link", se.fit = TRUE)
      response_predict2<- as.data.frame(response_predict2)      
      names(response_predict2) <- c("fit2", "se_fit2", "residual_scale2")
      prediction2 <-cbind(predictors2,response_predict2)
      predictors_l$r0  <- (prediction2[prediction2$day_week == end_phase,]$fit2 - prediction2[prediction2$day_week == start_phase,]$fit2)/duration_phase  
      predictors_l$se_r0 <- sqrt((prediction2[prediction2$day_week == end_phase,]$se_fit2)**2 + (prediction2[prediction2$day_week == start_phase,]$se_fit2)**2)/duration_phase
      predictors_l$Condition<-paste(s,t,r,sep="_")
      predictors_l$Temperature<-t
      predictors_l$Culture<-s
      predictors_l$Salinity<-r
      predictors_l$Start_exp_growht<-start_phase
      predictors_l$End_exp_growth<-end_phase
      predictors_l$aic<-model_f$aic
      predictors_all_factor<-rbind(predictors_all_factor,predictors_l)
      
    }
  }
}

#################plot the growth rates##########################################
list.grplot<-c()
list_grplot <- vector(mode = "list", length = l)
i <- 1
for (x in unique(predictors_all$Culture)){
  #x<-"mix"
  plt_subset<-subset(predictors_all,Culture==x)
  plt_subset_factor<-subset(predictors_all_factor,Culture==x)
  if (x=="mix") lab<-paste("The change of the growth rates of the ",x,sep="") else lab<-paste("The change of the growth rates of strain ",x,sep="")
  
  plt<-ggplot(plt_subset_factor, aes(as.numeric(Week),r0,color=as.character(Temperature)))+ #growth is stored into the global environment
    geom_point(aes(color=as.factor(Temperature),shape=Salinity))+
    geom_line(data=plt_subset,aes(group=Condition,linetype=Salinity))+
    theme_bw(base_size=14)+
    scale_colour_manual( values=c("blue", "red"))+
    theme(aspect.ratio=1)+
    labs(title=lab,x="Week",y="Growth rate",color="Temperature")+
    theme(plot.title = element_text(size=12)) #you might need to change the limits to fit your data
  plt
  list_grplot[[i]] <- plt
  names(list_grplot)[i] <- paste("r0", x, sep = "_")
  i=i+1
}

figure_r0<-ggarrange(list_grplot$r0_1C,list_grplot$r0_3,list_grplot$r0_13,list_grplot$r0_mix, labels=c("a","b","c","d"),ncol=2,nrow=2)
figure_r0

############get the AIC values from your models#################################
#bring the AIC values of all treatments from both models into one data.frame
aic_table_n<-subset(predictors_all,Week==1)
aic_comp_n<-cbind(aic_table_n$Condition,aic_table_n$aic)
mean_aic_n<-c("mean",mean(predictors_all$aic))
sd_aic_n<-c("sd",sd(predictors_all$aic))
aic_comp_n<-rbind(aic_comp_n,mean_aic_n,sd_aic_n)
num<-rep("numeric",length(aic_comp_n[,1]))
dim(num)<-c(length(num),1)
aic_comp_n<-cbind(aic_comp_n,num)

aic_table_f<-subset(predictors_all_factor,Week==1)
aic_comp_f<-cbind(aic_table_f$Condition,aic_table_f$aic)
mean_aic_f<-c("mean",mean(predictors_all$aic))
sd_aic_f<-c("sd",sd(predictors_all$aic))
aic_comp_f<-rbind(aic_comp_f,mean_aic_f,sd_aic_f)
num<-rep("factor",length(aic_comp_n[,1]))
dim(num)<-c(length(num),1)
aic_comp_f<-cbind(aic_comp_f,num)

aic_comp<-rbind(aic_comp_n,aic_comp_f)
aic_comp<-as.data.frame(aic_comp)
names(aic_comp)<-c("Condition","AIC","Week_in_model")

##########################plot the AIC values###################################
#remove the sd and mean rows
aic_comp<-aic_comp[-(35:36),]
aic_comp<-aic_comp[-(17:18),]

aic_plot<-ggplot(aic_comp,aes(factor(Condition,ordered=TRUE),as.numeric(AIC),color=Week_in_model))+
  geom_point()+
  theme_bw(base_size=14)+
  scale_colour_manual( values=c("green", "orange"))+
  theme(aspect.ratio=1)+
  labs(title="AIC values compared between the two models for each condition",x="Conditions",y="AIC",color="Model")+
  theme(plot.title = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aic_plot

########################## THE END #############################################
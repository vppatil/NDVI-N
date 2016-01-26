#this script contains two functions:
#julian() takes three arguments (month, day, and year), and returns the corresponding julian day

#julian.df() is a convenience function that produces a vector of julian days from a data frame containing columns named 'Month', 'Day', and 'Year'

julian<-function(day,month,year){
	julian.vector<-vector()
	#first, I created two vectors containing the number of days in each month for regular and leap years
	
	monthmatrix1<-c(0,31,28,31,30,31,30,31,31,30,31,30)
	monthmatrix2<-c(0,31,29,31,30,31,30,31,31,30,31,30)
	
	#the first step in the function is to determine whether the year is a leap year by checking if it has a remainder when divided by four
	#the appropriate vector is then assigned to the object called 'julian.vector'
	if (year/4-round(year/4)==0){julian.vector<-monthmatrix2}else{
		julian.vector<-monthmatrix1}
  		
	#finally, the julian date is calculated by adding up the number of days in all previous months during the year, 
	#plus the days that have elapsed in the current month
	
	jdate<-day+sum(julian.vector[1:month])
	return(jdate)
	}

julian.r<-function(jdate,year){
  julian.vector<-vector()
  #first, I created two vectors containing the number of days in each month for regular and leap years
  
  monthmatrix1<-c(0,31,28,31,30,31,30,31,31,30,31,30)
  monthmatrix2<-c(0,31,29,31,30,31,30,31,31,30,31,30)
  
  #the first step in the function is to determine whether the year is a leap year by checking if it has a remainder when divided by four
  #the appropriate vector is then assigned to the object called 'julian.vector'
  if (year/4-round(year/4)==0){julian.vector<-monthmatrix2}else{
    julian.vector<-monthmatrix1}
  
  j.cumsum=cumsum(julian.vector)
  for(i in 1:(length(j.cumsum)-1))
  {
    if(jdate>j.cumsum[i] & jdate<j.cumsum[i+1]){month = i}
  }
  
  #finally, the julian date is calculated by adding up the number of days in all previous months during the year, 
  #plus the days that have elapsed in the current month
  return(month)
}

julian.df<-function(data){

#this function doesn't work if a row of data doesn't contain values for all three arguments
options(na.action='na.omit')
jdate<-vector()

#basically, the julian() function is run sequentially for each row of data
for (i in 1:length(data$Year)){
	jdate[i]<-ifelse(is.na(data$Month[i]),NA,julian(data$Day[i],data$Month[i],data$Year[i]))}
	return(jdate)
	}

########################
#methods-CuffClusterSet.R
#
#Author: Loyal A. Goff
#
#Date created: 5-17-2011
#
#Description:
#########################

#################
#Initialize		#
#################


#################
#Validate		#
#################
#TODO: Add validity constraints
setValidity("CuffClusterSet",function(object){
			TRUE
			#Add test for genes with no expression
		}
)		



#################
#Class Methods	#
#################

#################
#Subsetting		#
#################
#TODO: Add subset methods to return a CuffFeature object
#setMethod("[","CuffFeatureSet",function(object,featureID){
#			
#		}
#)

#################
#Accessors
##################
#clustering

#################
#Plotting		#
#################
#plot cluster images


#################
#Misc			#
#################
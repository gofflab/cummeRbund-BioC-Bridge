# TODO: Add comment
# 
# Author: lgoff
###############################################################################


#####################
#File Archetype parsing
#####################
#Genes
loadGenes<-function(fpkmFile,
		diffFile,
		promoterFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		promoterArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {

	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)

	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########

	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate genes table
	######
	genesTable<-full[,c(1:3,5,7:9)]
	write("Writing genes table",stderr())
	#dbWriteTable(dbConn,'genes',genesTable,row.names=F,append=T)
	insert_SQL<-'INSERT INTO genes VALUES(:tracking_id, :class_code, :nearest_ref_id, :gene_short_name, :locus, :length, :coverage)'
	bulk_insert(dbConn,insert_SQL,genesTable)
	
	######
	#Populate geneData table
	######
	write("Reshaping geneData table",stderr())
	genemelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
	colnames(genemelt)[colnames(genemelt)=='variable']<-'sample_name'
	#Clean up and normalize data
	genemelt$measurement = ""
	
	genemelt$measurement[grepl("_FPKM$",genemelt$sample_name)] = "fpkm"
	genemelt$measurement[grepl("_conf_lo$",genemelt$sample_name)] = "conf_lo"
	genemelt$measurement[grepl("_conf_hi$",genemelt$sample_name)] = "conf_hi"
	genemelt$measurement[grepl("_status$",genemelt$sample_name)] = "status"

	genemelt$sample_name<-gsub("_FPKM$","",genemelt$sample_name)
	genemelt$sample_name<-gsub("_conf_lo$","",genemelt$sample_name)
	genemelt$sample_name<-gsub("_conf_hi$","",genemelt$sample_name)
	genemelt$sample_name<-gsub("_status$","",genemelt$sample_name)
	
	#Adjust sample names with make.db.names
	genemelt$sample_name <- make.db.names(dbConn,as.vector(genemelt$sample_name),unique=FALSE)
	
	#Recast
	write("Recasting",stderr())
	genemelt<-as.data.frame(cast(genemelt,...~measurement))
	
	#debugging
	#write(colnames(genemelt),stderr())
	
	#Write geneData table
	write("Writing geneData table",stderr())
	#dbWriteTable(dbConn,'geneData',as.data.frame(genemelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
	insert_SQL<-'INSERT INTO geneData VALUES(:tracking_id,:sample_name,:fpkm,:conf_hi,:conf_lo,:status)'
	bulk_insert(dbConn,insert_SQL,genemelt[,c(1:2,5,3,4,6)])
	
	#######
	#Handle gene_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		#Something like this to make sure sample names are treated as character values and not numeric, logical, etc.
		#diffArgs$colClasses<-c(rep('character',7),rep('numeric',6),'character')
		diff<-as.data.frame(do.call(read.table,diffArgs))
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)
			
			write("Writing geneExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			
			#debugging
			#write(colnames(diff[,diffCols]),stderr())
			
			#dbWriteTable(dbConn,'geneExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO geneExpDiffData VALUES(:test_id,:sample_1,:sample_2,:status,:value_1,:value_2,?,:test_stat,:p_value,:q_value,:significant)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in", diffFile),stderr())
		}
	
	}
	
	########
	#TODO: Handle promoters.diff
	########
	if(file.exists(promoterFile)){
		#Read promoterFile
		write(paste("Reading ",promoterFile,sep=""),stderr())
		promoterArgs$file = promoterFile
		promoter<-as.data.frame(do.call(read.table,promoterArgs))
		
		write("Writing promoterDiffData table",stderr())
		promoterCols<-c(2,5:14)
		if(dim(promoter)[1]>0){
			#dbWriteTable(dbConn,'promoterDiffData',promoter[,promoterCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO promoterDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,promoter[,promoterCols])
		}else{
			write(paste("No records found in", promoterFile),stderr())
		}
	}
	
	#########
	#Handle Feature Data (this will actually be done on CuffData objects instead...but I may include something here as well)
	#########
	
}
	
#Isoforms
loadIsoforms<-function(fpkmFile,
		diffFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {
	
	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)
	
	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########
	
	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			write(samples,stderr())
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate isoforms table
	######
	isoformCols<-c(1,4,5,6,2,3,7:9)
	isoformsTable<-full[,isoformCols]
	
	#This is a temporary fix until p_id is added to the 'isoforms.fpkm_tracking' file
	isoformsTable<-cbind(isoformsTable[,1:2],data.frame(CDS_id=rep("NA",dim(isoformsTable)[1])),isoformsTable[,-c(1:2)])
	#print (head(isoformsTable))
	write("Writing isoforms table",stderr())
	#dbWriteTable(dbConn,'isoforms',as.data.frame(isoformsTable),row.names=F,append=T)
	insert_SQL<-'INSERT INTO isoforms VALUES(?,?,?,?,?,?,?,?,?,?)'
	bulk_insert(dbConn,insert_SQL,isoformsTable)
	
	######
	#Populate isoformData table
	######
	write("Reshaping isoformData table",stderr())
	isoformmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
	colnames(isoformmelt)[colnames(isoformmelt)=='variable']<-'sample_name'
	#Clean up and normalize data
	isoformmelt$measurement = ""
	
	isoformmelt$measurement[grepl("_FPKM$",isoformmelt$sample_name)] = "fpkm"
	isoformmelt$measurement[grepl("_conf_lo$",isoformmelt$sample_name)] = "conf_lo"
	isoformmelt$measurement[grepl("_conf_hi$",isoformmelt$sample_name)] = "conf_hi"
	isoformmelt$measurement[grepl("_status$",isoformmelt$sample_name)] = "status"
	
	isoformmelt$sample_name<-gsub("_FPKM$","",isoformmelt$sample_name)
	isoformmelt$sample_name<-gsub("_conf_lo$","",isoformmelt$sample_name)
	isoformmelt$sample_name<-gsub("_conf_hi$","",isoformmelt$sample_name)
	isoformmelt$sample_name<-gsub("_status$","",isoformmelt$sample_name)
	
	#Adjust sample names with make.db.names
	isoformmelt$sample_name <- make.db.names(dbConn,as.vector(isoformmelt$sample_name),unique=FALSE)
	
	#Recast
	write("Recasting",stderr())
	isoformmelt<-as.data.frame(cast(isoformmelt,...~measurement))
	
	#Write geneData table
	write("Writing isoformData table",stderr())
	#dbWriteTable(dbConn,'isoformData',as.data.frame(isoformmelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
	insert_SQL<-"INSERT INTO isoformData VALUES(?,?,?,?,?,?)"
	bulk_insert(dbConn,insert_SQL,isoformmelt[,c(1:2,5,3,4,6)])
	
	#######
	#Handle isoform_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		diff<-as.data.frame(do.call(read.table,diffArgs))
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)
		
			write("Writing isoformExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			#dbWriteTable(dbConn,'isoformExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO isoformExpDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in",diffFile),stderr())
		}
	}
	
}

#TSS groups
loadTSS<-function(fpkmFile,
		diffFile,
		splicingFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		splicingArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {
	
	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)
	
	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########
	
	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate genes table
	######
	tssTable<-full[,c(1:5,7:9)]
	write("Writing TSS table",stderr())
	#dbWriteTable(dbConn,'TSS',tssTable,row.names=F,append=T)
	if (nrow(tssTable)>0){
		insert_SQL<-"INSERT INTO TSS VALUES(?,?,?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,tssTable)
		
		######
		#Populate geneData table
		######
		write("Reshaping TSSData table",stderr())
		tssmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
		colnames(tssmelt)[colnames(tssmelt)=='variable']<-'sample_name'
		#Clean up and normalize data
		tssmelt$measurement = ""
		
		tssmelt$measurement[grepl("_FPKM$",tssmelt$sample_name)] = "fpkm"
		tssmelt$measurement[grepl("_conf_lo$",tssmelt$sample_name)] = "conf_lo"
		tssmelt$measurement[grepl("_conf_hi$",tssmelt$sample_name)] = "conf_hi"
		tssmelt$measurement[grepl("_status$",tssmelt$sample_name)] = "status"
		
		tssmelt$sample_name<-gsub("_FPKM$","",tssmelt$sample_name)
		tssmelt$sample_name<-gsub("_conf_lo$","",tssmelt$sample_name)
		tssmelt$sample_name<-gsub("_conf_hi$","",tssmelt$sample_name)
		tssmelt$sample_name<-gsub("_status$","",tssmelt$sample_name)
		
		#Adjust sample names with make.db.names
		tssmelt$sample_name <- make.db.names(dbConn,as.vector(tssmelt$sample_name),unique=FALSE)
		
		#Recast
		write("Recasting",stderr())
		tssmelt<-as.data.frame(cast(tssmelt,...~measurement))
		
		#Write geneData table
		write("Writing TSSData table",stderr())
		#dbWriteTable(dbConn,'TSSData',as.data.frame(tssmelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)

		insert_SQL<-"INSERT INTO TSSData VALUES(?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,tssmelt[,c(1:2,5,3,4,6)])
	}else{
		write(paste("No records found in",fpkmFile),stderr())
		write("TSS FPKM tracking file was empty.",stderr())
	}
	#######
	#Handle tss_groups_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		diff<-as.data.frame(do.call(read.table,diffArgs))
		
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)

			write("Writing TSSExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			#dbWriteTable(dbConn,'TSSExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO TSSExpDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in",diffFile),stderr())
		}
	}
	
	#########
	#TODO: Handle splicing.diff
	########
	if(file.exists(splicingFile)){
		#Read promoterFile
		write(paste("Reading ",splicingFile,sep=""),stderr())
		splicingArgs$file = splicingFile
		splicing<-as.data.frame(do.call(read.table,splicingArgs))
		
		if(dim(splicing)[1]>0){
			write("Writing splicingDiffData table",stderr())
			splicingCols<-c(1:2,5:14)
			#dbWriteTable(dbConn,'splicingDiffData',splicing[,splicingCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO splicingDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,splicing[,splicingCols])
		}else{
			write(paste("No records found in",splicingFile),stderr())
		}
	}
	
}

#CDS
loadCDS<-function(fpkmFile,
		diffFile,
		CDSDiff,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		CDSDiffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {
	
	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)
	
	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########
	
	
	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate genes table
	######
	cdsTable<-full[,c(1:5,6:9)]
	write("Writing CDS table",stderr())
	#dbWriteTable(dbConn,'CDS',cdsTable,row.names=F,append=T)
	if (nrow(cdsTable)>0){
		insert_SQL<-"INSERT INTO CDS VALUES(?,?,?,?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,cdsTable)
		
		######
		#Populate geneData table
		######
		write("Reshaping CDSData table",stderr())
		cdsmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
		colnames(cdsmelt)[colnames(cdsmelt)=='variable']<-'sample_name'
		#Clean up and normalize data
		cdsmelt$measurement = ""
		
		cdsmelt$measurement[grepl("_FPKM$",cdsmelt$sample_name)] = "fpkm"
		cdsmelt$measurement[grepl("_conf_lo$",cdsmelt$sample_name)] = "conf_lo"
		cdsmelt$measurement[grepl("_conf_hi$",cdsmelt$sample_name)] = "conf_hi"
		cdsmelt$measurement[grepl("_status$",cdsmelt$sample_name)] = "status"
		
		cdsmelt$sample_name<-gsub("_FPKM$","",cdsmelt$sample_name)
		cdsmelt$sample_name<-gsub("_conf_lo$","",cdsmelt$sample_name)
		cdsmelt$sample_name<-gsub("_conf_hi$","",cdsmelt$sample_name)
		cdsmelt$sample_name<-gsub("_status$","",cdsmelt$sample_name)
		
		#Adjust sample names with make.db.names
		cdsmelt$sample_name <- make.db.names(dbConn,as.vector(cdsmelt$sample_name),unique=FALSE)
		
		#Recast
		write("Recasting",stderr())
		cdsmelt<-as.data.frame(cast(cdsmelt,...~measurement))
		
		#Write geneData table
		write("Writing CDSData table",stderr())
		#dbWriteTable(dbConn,'CDSData',as.data.frame(cdsmelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
		insert_SQL<-"INSERT INTO CDSData VALUES(?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,cdsmelt[,c(1:2,5,3,4,6)])
	
	}else {
		write(paste("No records found in",fpkmFile),stderr())
		write("CDS FPKM tracking file was empty.",stderr())
	}
	
	
	#######
	#Handle cds_groups_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		diff<-as.data.frame(do.call(read.table,diffArgs))
		
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)
			
			write("Writing CDSExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			#dbWriteTable(dbConn,'CDSExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO CDSExpDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in",diffFile),stderr())
		}
	}
	
	#########
	#TODO: Handle CDS.diff
	########
	if(file.exists(CDSDiff)){
		#Read promoterFile
		write(paste("Reading ",CDSDiff,sep=""),stderr())
		CDSDiffArgs$file = CDSDiff
		CDS<-as.data.frame(do.call(read.table,CDSDiffArgs))
		if(dim(CDS)[1]>0){
			write("Writing CDSDiffData table",stderr())
			CDSCols<-c(2,5:14)
			#dbWriteTable(dbConn,'CDSDiffData',CDS[,CDSCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO CDSDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,CDS[,CDSCols])
		}else{
			write(paste("No records found in",CDSDiff),stderr())
		}
	}
	
}

########################
#Add FeatureData
########################


#####################
#Database Setup Functions
#####################

createDB<-function(dbFname="cuffData.db",driver="SQLite") {
	#Builds sqlite db at 'dbFname' and returns a dbConnect object pointing to newly created database.
	#May be deprecated if I can load first and index later...
	
	drv<-dbDriver(driver)
	db <- dbConnect(drv,dbname=dbFname)
	
	schema.text<-'
-- Creator:       MySQL Workbench 5.2.33/ExportSQLite plugin 2009.12.02
-- Author:        Loyal Goff
-- Caption:       CuffData.db Model
-- Project:       cummeRbund
-- Changed:       2011-08-02 14:03
-- Created:       2011-05-02 12:52
PRAGMA foreign_keys = OFF;
PRAGMA synchronous = OFF;
PRAGMA journal_mode = MEMORY;
-- Schema: cuffData
BEGIN;
DROP TABLE IF EXISTS "genes";
CREATE TABLE "genes"(
  "gene_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT
);
DROP TABLE IF EXISTS "biasData";
CREATE TABLE "biasData"(
  "biasData_id" INTEGER PRIMARY KEY NOT NULL
);
DROP TABLE IF EXISTS "samples";
CREATE TABLE "samples"(
  "sample_index" INTEGER PRIMARY KEY NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL
);
DROP TABLE IF EXISTS "TSS";
CREATE TABLE "TSS"(
  "TSS_group_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45) NOT NULL,
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_TSS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "TSS.fk_TSS_genes1" ON "TSS"("gene_id");
CREATE INDEX "TSS.fk_TSS_genes2" ON "TSS"("gene_short_name");
DROP TABLE IF EXISTS "TSSData";
CREATE TABLE "TSSData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_TSSData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "TSSData.fk_TSSData_TSS1" ON "TSSData"("TSS_group_id");
CREATE INDEX "TSSData.fk_TSSData_samples1" ON "TSSData"("sample_name");
DROP TABLE IF EXISTS "CDS";
CREATE TABLE "CDS"(
  "CDS_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_CDS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_CDS_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
CREATE INDEX "CDS.fk_CDS_genes1" ON "CDS"("gene_id");
CREATE INDEX "CDS.fk_CDS_genes2" ON "CDS"("gene_short_name");
CREATE INDEX "CDS.fk_CDS_TSS1" ON "CDS"("TSS_group_id");
DROP TABLE IF EXISTS "CDSData";
CREATE TABLE "CDSData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_CDSData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "CDSData.fk_CDSData_CDS1" ON "CDSData"("CDS_id");
CREATE INDEX "CDSData.fk_CDSData_samples1" ON "CDSData"("sample_name");
DROP TABLE IF EXISTS "splicingDiffData";
CREATE TABLE "splicingDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_splicingDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_splicingDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples1" ON "splicingDiffData"("sample_1");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples2" ON "splicingDiffData"("sample_2");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_TSS1" ON "splicingDiffData"("TSS_group_id");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_genes1" ON "splicingDiffData"("gene_id");
DROP TABLE IF EXISTS "TSSExpDiffData";
CREATE TABLE "TSSExpDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_TSSExpDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_TSSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_TSS1" ON "TSSExpDiffData"("TSS_group_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples1" ON "TSSExpDiffData"("sample_1");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples2" ON "TSSExpDiffData"("sample_2");
DROP TABLE IF EXISTS "CDSDiffData";
CREATE TABLE "CDSDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples1" ON "CDSDiffData"("sample_1");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples2" ON "CDSDiffData"("sample_2");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_genes1" ON "CDSDiffData"("gene_id");
DROP TABLE IF EXISTS "CDSExpDiffData";
CREATE TABLE "CDSExpDiffData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSExpDiffData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_CDS1" ON "CDSExpDiffData"("CDS_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples1" ON "CDSExpDiffData"("sample_1");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples2" ON "CDSExpDiffData"("sample_2");
DROP TABLE IF EXISTS "promoterDiffData";
CREATE TABLE "promoterDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_promoterDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_promoterDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_promoterDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "promoterDiffData.fk_promoterDiffData_genes1" ON "promoterDiffData"("gene_id");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples1" ON "promoterDiffData"("sample_1");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples2" ON "promoterDiffData"("sample_2");
DROP TABLE IF EXISTS "geneFeatures";
CREATE TABLE "geneFeatures"(
  "gene_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_geneFeatures_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "geneFeatures.fk_geneFeatures_genes1" ON "geneFeatures"("gene_id");
DROP TABLE IF EXISTS "TSSFeatures";
CREATE TABLE "TSSFeatures"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_TSSFeatures_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
CREATE INDEX "TSSFeatures.fk_TSSFeatures_TSS1" ON "TSSFeatures"("TSS_group_id");
DROP TABLE IF EXISTS "CDSFeatures";
CREATE TABLE "CDSFeatures"(
  "CDS_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_CDSFeatures_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
CREATE INDEX "CDSFeatures.fk_CDSFeatures_CDS1" ON "CDSFeatures"("CDS_id");
DROP TABLE IF EXISTS "geneData";
CREATE TABLE "geneData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "geneData.fk_geneData_genes1" ON "geneData"("gene_id");
CREATE INDEX "geneData.fk_geneData_samples1" ON "geneData"("sample_name");
DROP TABLE IF EXISTS "phenoData";
CREATE TABLE "phenoData"(
  "sample_name" VARCHAR(45) NOT NULL,
  "parameter" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45),
  CONSTRAINT "fk_phenoData_samples"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "phenoData.fk_phenoData_samples" ON "phenoData"("sample_name");
DROP TABLE IF EXISTS "geneExpDiffData";
CREATE TABLE "geneExpDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_geneExpDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_genes1" ON "geneExpDiffData"("gene_id");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples1" ON "geneExpDiffData"("sample_1");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples2" ON "geneExpDiffData"("sample_2");
DROP TABLE IF EXISTS "isoforms";
CREATE TABLE "isoforms"(
  "isoform_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "gene_id" VARCHAR(45),
  "CDS_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_isoforms_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_isoforms_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_isoforms_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "isoforms.fk_isoforms_TSS1" ON "isoforms"("TSS_group_id");
CREATE INDEX "isoforms.fk_isoforms_CDS1" ON "isoforms"("CDS_id");
CREATE INDEX "isoforms.fk_isoforms_genes1" ON "isoforms"("gene_id");
CREATE INDEX "isoforms.fk_isoforms_genes2" ON "isoforms"("gene_short_name");
DROP TABLE IF EXISTS "isoformData";
CREATE TABLE "isoformData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT NOT NULL,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_isoformData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
CREATE INDEX "isoformData.fk_isoformData_samples1" ON "isoformData"("sample_name");
CREATE INDEX "isoformData.fk_isoformData_isoforms1" ON "isoformData"("isoform_id");
DROP TABLE IF EXISTS "isoformExpDiffData";
CREATE TABLE "isoformExpDiffData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_isoformExpDiffData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_isoforms1" ON "isoformExpDiffData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples1" ON "isoformExpDiffData"("sample_1");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples2" ON "isoformExpDiffData"("sample_2");
DROP TABLE IF EXISTS "isoformFeatures";
CREATE TABLE "isoformFeatures"(
  "isoform_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_isoformFeatures_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
CREATE INDEX "isoformFeatures.fk_isoformFeatures_isoforms1" ON "isoformFeatures"("isoform_id");
COMMIT;

'
		create.sql <- strsplit(schema.text, "\n")[[1]]
		create.sql <- paste(collapse="\n", create.sql)
		create.sql <- strsplit(create.sql, ";")[[1]]
		create.sql <- create.sql[-length(create.sql)] #nothing to run here
				
		tmp <- sapply(create.sql,function(x) sqliteQuickSQL(db,x))
		db
}

createDB_noIndex<-function(dbFname="cuffData.db",driver="SQLite") {
	#Builds sqlite db at 'dbFname' and returns a dbConnect object pointing to newly created database.
	#No indexes are present
	
	drv<-dbDriver(driver)
	db <- dbConnect(drv,dbname=dbFname)
	
	schema.text<-'
-- Creator:       MySQL Workbench 5.2.33/ExportSQLite plugin 2009.12.02
-- Author:        Loyal Goff
-- Caption:       CuffData.db Model
-- Project:       cummeRbund
-- Changed:       2011-08-02 14:03
-- Created:       2011-05-02 12:52
PRAGMA foreign_keys = OFF;
PRAGMA synchronous = OFF;
-- Schema: cuffData
BEGIN;
DROP TABLE IF EXISTS "genes";
CREATE TABLE "genes"(
  "gene_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT
);
DROP TABLE IF EXISTS "biasData";
CREATE TABLE "biasData"(
  "biasData_id" INTEGER PRIMARY KEY NOT NULL
);
DROP TABLE IF EXISTS "samples";
CREATE TABLE "samples"(
  "sample_index" INTEGER PRIMARY KEY NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL
);
DROP TABLE IF EXISTS "TSS";
CREATE TABLE "TSS"(
  "TSS_group_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45) NOT NULL,
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_TSS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSData";
CREATE TABLE "TSSData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_TSSData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDS";
CREATE TABLE "CDS"(
  "CDS_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_CDS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_CDS_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "CDSData";
CREATE TABLE "CDSData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_CDSData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "splicingDiffData";
CREATE TABLE "splicingDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_splicingDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_splicingDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSExpDiffData";
CREATE TABLE "TSSExpDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_TSSExpDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_TSSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDSDiffData";
CREATE TABLE "CDSDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "CDSExpDiffData";
CREATE TABLE "CDSExpDiffData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSExpDiffData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "promoterDiffData";
CREATE TABLE "promoterDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_promoterDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_promoterDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_promoterDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneFeatures";
CREATE TABLE "geneFeatures"(
  "gene_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_geneFeatures_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSFeatures";
CREATE TABLE "TSSFeatures"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_TSSFeatures_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "CDSFeatures";
CREATE TABLE "CDSFeatures"(
  "CDS_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_CDSFeatures_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
DROP TABLE IF EXISTS "geneData";
CREATE TABLE "geneData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "phenoData";
CREATE TABLE "phenoData"(
  "sample_name" VARCHAR(45) NOT NULL,
  "parameter" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45),
  CONSTRAINT "fk_phenoData_samples"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneExpDiffData";
CREATE TABLE "geneExpDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_geneExpDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoforms";
CREATE TABLE "isoforms"(
  "isoform_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "gene_id" VARCHAR(45),
  "CDS_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_isoforms_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_isoforms_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_isoforms_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "isoformData";
CREATE TABLE "isoformData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT NOT NULL,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_isoformData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "isoformExpDiffData";
CREATE TABLE "isoformExpDiffData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_isoformExpDiffData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoformFeatures";
CREATE TABLE "isoformFeatures"(
  "isoform_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_isoformFeatures_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
COMMIT;

			'
	create.sql <- strsplit(schema.text, "\n")[[1]]
	create.sql <- paste(collapse="\n", create.sql)
	create.sql <- strsplit(create.sql, ";")[[1]]
	create.sql <- create.sql[-length(create.sql)] #nothing to run here
	
	tmp <- sapply(create.sql,function(x) sqliteQuickSQL(db,x))
	db
}


createIndices<-function(dbFname="cuffData.db",driver="SQLite",verbose=F){
	
	drv<-dbDriver(driver)
	db <- dbConnect(drv,dbname=dbFname)
	
	index.text<-
'CREATE INDEX "genes.gene_id_index" ON "genes"("gene_id");
CREATE INDEX "genes.gsn_index" ON "genes"("gene_short_name");
CREATE INDEX "genes.cc_index" ON "genes"("class_code");
CREATE INDEX "TSS.fk_TSS_genes1" ON "TSS"("gene_id");
CREATE INDEX "TSS.fk_TSS_genes2" ON "TSS"("gene_short_name");
CREATE INDEX "TSSData.fk_TSSData_TSS1" ON "TSSData"("TSS_group_id");
CREATE INDEX "TSSData.fk_TSSData_samples1" ON "TSSData"("sample_name");
CREATE INDEX "TSS.PRIMARY" ON "TSS"("TSS_group_id");
CREATE INDEX "CDS.PRIMARY" ON "CDS"("CDS_id");
CREATE INDEX "CDS.fk_CDS_genes1" ON "CDS"("gene_id");
CREATE INDEX "CDS.fk_CDS_genes2" ON "CDS"("gene_short_name");
CREATE INDEX "CDS.fk_CDS_TSS1" ON "CDS"("TSS_group_id");
CREATE INDEX "CDSData.fk_CDSData_CDS1" ON "CDSData"("CDS_id");
CREATE INDEX "CDSData.fk_CDSData_samples1" ON "CDSData"("sample_name");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples1" ON "splicingDiffData"("sample_1");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples2" ON "splicingDiffData"("sample_2");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_TSS1" ON "splicingDiffData"("TSS_group_id");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_genes1" ON "splicingDiffData"("gene_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_TSS1" ON "TSSExpDiffData"("TSS_group_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples1" ON "TSSExpDiffData"("sample_1");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples2" ON "TSSExpDiffData"("sample_2");
CREATE INDEX "TSSExpDiffData.TSSExpDiffData_sig_index" ON "TSSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples1" ON "CDSDiffData"("sample_1");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples2" ON "CDSDiffData"("sample_2");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_genes1" ON "CDSDiffData"("gene_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_CDS1" ON "CDSExpDiffData"("CDS_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples1" ON "CDSExpDiffData"("sample_1");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples2" ON "CDSExpDiffData"("sample_2");
CREATE INDEX "CDSExpDiffData.CDSExpDiffData_sig_index" ON "CDSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_genes1" ON "promoterDiffData"("gene_id");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples1" ON "promoterDiffData"("sample_1");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples2" ON "promoterDiffData"("sample_2");
CREATE INDEX "geneFeatures.fk_geneFeatures_genes1" ON "geneFeatures"("gene_id");
CREATE INDEX "TSSFeatures.fk_TSSFeatures_TSS1" ON "TSSFeatures"("TSS_group_id");
CREATE INDEX "CDSFeatures.fk_CDSFeatures_CDS1" ON "CDSFeatures"("CDS_id");
CREATE INDEX "geneData.fk_geneData_genes1" ON "geneData"("gene_id");
CREATE INDEX "geneData.fk_geneData_samples1" ON "geneData"("sample_name");
CREATE INDEX "phenoData.fk_phenoData_samples" ON "phenoData"("sample_name");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_genes1" ON "geneExpDiffData"("gene_id");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples1" ON "geneExpDiffData"("sample_1");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples2" ON "geneExpDiffData"("sample_2");
CREATE INDEX "geneExpDiffData.geneExpDiff_status_index" ON "geneExpDiffData"("status");
CREATE INDEX "geneExpDiffData.geneExpDiff_sig_index" ON "geneExpDiffData"("significant","p_value","q_value","test_stat");
CREATE INDEX "isoforms.PRIMARY" ON "isoforms"("isoform_id");
CREATE INDEX "isoforms.fk_isoforms_TSS1" ON "isoforms"("TSS_group_id");
CREATE INDEX "isoforms.fk_isoforms_CDS1" ON "isoforms"("CDS_id");
CREATE INDEX "isoforms.fk_isoforms_genes1" ON "isoforms"("gene_id");
CREATE INDEX "isoforms.fk_isoforms_genes2" ON "isoforms"("gene_short_name");
CREATE INDEX "isoformData.fk_isoformData_samples1" ON "isoformData"("sample_name");
CREATE INDEX "isoformData.fk_isoformData_isoforms1" ON "isoformData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_isoforms1" ON "isoformExpDiffData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples1" ON "isoformExpDiffData"("sample_1");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples2" ON "isoformExpDiffData"("sample_2");
CREATE INDEX "isoformExpDiffData.isoformExpDiffData_sig_index" ON "isoformExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "isoformFeatures.fk_isoformFeatures_isoforms1" ON "isoformFeatures"("isoform_id");
'
	create.sql <- strsplit(index.text,"\n")[[1]]
	
	tmp <- sapply(create.sql,function(x){
			if (verbose){
						write(paste(x,sep=""),stderr())
					}
			sqliteQuickSQL(db,x)
	})
}


getSamples<-function(fpkmDF){
	sample_name<-unique(fpkmDF$sample)
	#sample_name<-as.data.frame(sample_name)
}

getSamplesFromColnames<-function(fpkmDF){
	samples<-gsub("_FPKM$","",colnames(fpkmDF)[grepl("_FPKM$",colnames(fpkmDF))])
}

populateSampleTable<-function(samples,dbConn){
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	samples<-data.frame(index=c(1:length(samples)),sample_name=samples)
	dbWriteTable(dbConn,'samples',samples,row.names=F,append=T)
}

bulk_insert <- function(dbConn,sql,bound.data)
{
	dbBeginTransaction(dbConn)
	dbGetPreparedQuery(dbConn, sql, bind.data = bound.data)
	dbCommit(dbConn)
}

#############
#readCufflinks
#############
#TODO: Add directory pointer
readCufflinks<-function(dir = getwd(),
						dbFile="cuffData.db",
						geneFPKM="genes.fpkm_tracking",
						geneDiff="gene_exp.diff",
						isoformFPKM="isoforms.fpkm_tracking",
						isoformDiff="isoform_exp.diff",
						TSSFPKM="tss_groups.fpkm_tracking",
						TSSDiff="tss_group_exp.diff",
						CDSFPKM="cds.fpkm_tracking",
						CDSExpDiff="cds_exp.diff",
						CDSDiff="cds.diff",
						promoterFile="promoters.diff",
						splicingFile="splicing.diff",
						driver = "SQLite",
						rebuild = FALSE,
						...){
	
	#Set file locations with directory
	dbFile=file.path(dir,dbFile)
	geneFPKM=file.path(dir,geneFPKM)
	geneDiff=file.path(dir,geneDiff)
	isoformFPKM=file.path(dir,isoformFPKM)
	isoformDiff=file.path(dir,isoformDiff)
	TSSFPKM=file.path(dir,TSSFPKM)
	TSSDiff=file.path(dir,TSSDiff)
	CDSFPKM=file.path(dir,CDSFPKM)
	CDSExpDiff=file.path(dir,CDSExpDiff)
	CDSDiff=file.path(dir,CDSDiff)
	promoterFile=file.path(dir,promoterFile)
	splicingFile=file.path(dir,splicingFile)
					
					
	#Check to see whether dbFile exists
	if (!file.exists(dbFile) || rebuild == TRUE){
		#if not, create it
		write(paste("Creating database ",dbFile,sep=""),stderr())
		dbConn<-createDB_noIndex(dbFile)
		
		#populate DB
				
		loadGenes(geneFPKM,geneDiff,promoterFile,dbConn)
		loadIsoforms(isoformFPKM,isoformDiff,dbConn)
		loadTSS(TSSFPKM,TSSDiff,splicingFile,dbConn)
		loadCDS(CDSFPKM,CDSExpDiff,CDSDiff,dbConn)
		
		#Create Indexes on DB
		write("Indexing Tables...",stderr())
		createIndices(dbFile)
		
		#load Distribution Tests
		#loadDistTests(promoterFile,splicingFile,CDSDiff)
		
	}
	dbConn<-dbConnect(dbDriver(driver),dbFile)
	return (
			new("CuffSet",DB = dbConn,
					genes = new("CuffData",DB = dbConn, tables = list(mainTable = "genes",dataTable = "geneData",expDiffTable = "geneExpDiffData",featureTable = "geneFeatures"), filters = list(),type = "genes",idField = "gene_id"),
					isoforms = new("CuffData", DB = dbConn, tables = list(mainTable = "isoforms",dataTable = "isoformData",expDiffTable = "isoformExpDiffData",featureTable = "isoformFeatures"), filters = list(),type="isoforms",idField = "isoform_id"),
					TSS = new("CuffData", DB = dbConn, tables = list(mainTable = "TSS",dataTable = "TSSData",expDiffTable = "TSSExpDiffData",featureTable = "TSSFeatures"), filters = list(),type = "TSS",idField = "TSS_group_id"),
					CDS = new("CuffData", DB = dbConn, tables = list(mainTable = "CDS",dataTable = "CDSData",expDiffTable = "CDSExpDiffData",featureTable = "CDSFeatures"), filters = list(),type = "CDS",idField = "CDS_id"),
					promoters = new("CuffDist", DB = dbConn, table = "promoterDiffData",type="promoter",idField="gene_id"),
					splicing = new("CuffDist", DB = dbConn, table = "splicingDiffData",type="splicing",idField="TSS_group_id"),
					relCDS = new("CuffDist", DB = dbConn, table = "CDSDiffData",type="relCDS",idField="gene_id")
			)
	)	
							
}

############
# Handle GTF file
############
loadGTF<-function(gtfFile,dbConn) {
	
	#Error Trapping
	if (missing(gtfFile))
		stop("GTF file cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	gtf<-read.table(gtfFile,sep="\t",header=F)
	
	
	attributes<-melt(strsplit(as.character(gtf$V9),"; "))
	colnames(attributes)<-c("attribute","featureID")
	attributes<-paste(attributes$attribute,attributes$featureID)
	attributes<-strsplit(as.character(attributes)," ")
	attributes<-as.data.frame(do.call("rbind",attributes))
	
	colnames(attributes)<-c("attribute","value","featureID")
	attributes<-attributes[,c(3,1,2)]
	
	#Grab only gene_ID and transcript_ID to add to features table
	id.attributes<-attributes[attributes$attribute %in% c("gene_id","transcript_id"),]
	id.attributes$featureID<-as.numeric(as.character(id.attributes$featureID))
	id.attributes<-cast(id.attributes,...~attribute)
	
	#Main features table
	features<-gtf[,c(1:8)]
	colnames(features)<-c("seqname","source","type","start","end","score","strand","frame")
	features$featureID<-as.numeric(as.character(rownames(features)))
	
	#Merge features and id.attributes
	features<-merge(features,id.attributes,by.x='featureID',by.y='featureID')
	features<-features[,c(1,10:11,2:9)]
	
	#strip gene_id and transcript_id from attributes
	attributes<-attributes[!(attributes$attribute %in% c("gene_id","transcript_id")),]
	
	#Write features table
	write("Writing features table",stderr())
	#dbWriteTable(dbConn,'geneData',as.data.frame(genemelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
	dbWriteTable(dbConn,'features',as.data.frame(features),append=F)
	
	#Write features table
	write("Writing feature attributes table",stderr())
	dbWriteTable(dbConn,'attributes',as.data.frame(attributes),append=F)
	
}
	

#######
#Unit Test
#######

#dbConn<-createDB()
#date()
#loadGenes("genes.fpkm_tracking","gene_exp.diff",dbConn)
#loadIsoforms("isoforms.fpkm_tracking","isoform_exp.diff",dbConn)
#loadTSS("tss_groups.fpkm_tracking","tss_group_exp.diff",dbConn)
#loadCDS("cds.fpkm_tracking","cds_exp.diff",dbConn)
#date()
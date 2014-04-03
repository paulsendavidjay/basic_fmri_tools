#.libPaths("~/R-packages")
# 06/12/13
# 	changed (sum(cope_image.2D[,i]==0) > 0) to (count(cope_image.2D[,i] == 0)[2,2] > 0) to exclude from analysis voxels for which 
#	even just 1 subject was missing data
# install_github("robustlmm", "kollerma",ref="911bb686ba31f8ef2924b23eaebbbe94976bc05d") # use the commit reference on github that's taged with 1.1 instead of 1.2 -- after inspecting github log https://github.com/kollerma/robustlmm/commits/master
#
#	DEPENDENCIES: robustlmm MUST BE VERSION 1.1 TO ALLOW USING lme4 VERSION 0.999999-2
#					version 1.2 forces lme4 version to be 1.1 or greater, which isn't well documented and is not backwards compatible
#				  lme4 MUST BE VERSION 0.999999-2 IN ORDER TO COLLECT AIC AND TO USE mcmcsample()


#number_of_processors_to_use <- 7
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)

# current_cope = 1; model_name="age_beh3"

current_cope <- args[1]
model_name <- args[2]


if (length(args) < 2) {
	stop("       Must specify current_cope, model name, in that order")
}



#####################################################################################
# LOAD LIBRARIES AFTER DATA HAS BEEN CHECKED FOR

#.libPaths("~/R-packages")

suppressMessages(library(boot)) # for robust mixed effects regression
suppressMessages(library(AnalyzeFMRI)) # for loading manipulating fmri data
suppressMessages(library(R.utils)) # for gunzip
suppressMessages(library(nlme)) # for mixed effects regression amenable to MCMC
suppressMessages(library(reshape))
suppressMessages(library(doMC)) # for parallel processing


#####################################################################################
#							FUNCTIONS

# RANDOM STRING FOR MAKING .NII FILES IN CASE OF MULTIPLE OVERLAPPING FILE CREATIONS/DELETIONS
makeRandomString <- function(n=1, length=12) {
    randomString <- c(1:n) # initialize vector
    for (i in 1:n) {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS), length, replace=TRUE), collapse="")
    }
    return(randomString)
}

# LOADING FMRI DATA
load_gz_volumes <- function(file_name) { # returns a list of [[1]] header [[2]] data
		if (! substr(file_name, nchar(file_name) - 1, nchar(file_name)) == "gz") {
			print("not a gz file")
		} else {
			nii_file_name <- paste(c("/Volumes/Phillips/bars/analysis/temp/", makeRandomString(length=15), ".nii"), sep="", collapse="")
			gunzip(file_name, remove=FALSE, destname= nii_file_name)
			outfile_hdr <- f.read.nifti.header(nii_file_name)
			outfile_data <- f.read.nifti.volume(nii_file_name)
			unlink(nii_file_name)
			list(outfile_hdr, outfile_data)
		}
}

# CHANGE NAMES FROM LMER OUTPUT
change_stat_value_names <- function(matrix) {
	matrix <- gsub('\\(Intercept\\)', 'intcpt', matrix)
	matrix <- gsub('I\\(Age\\^2\\)', 'age2', matrix)
	matrix <- sub('Age\\:beh', 'ageXbeh', matrix)
	matrix <- gsub('Age', 'age', matrix)
	matrix <- sub('Estimate', 'pe', matrix)
	matrix <- sub('Std. Error', 'se', matrix)
	matrix <- sub('t value', 'tval', matrix)
	matrix <- sub('I\\(log\\(Age\\)\\)', 'logAge', matrix)
	matrix <- sub('I\\(1\\/Age\\)', 'invAge', matrix)

	matrix <- sub('I\\(age_mean_adj\\^2\\)', 'age_mean_adj2', matrix)

	matrix <- sub('Variance', 'var', matrix)
	matrix <- sub('Std.Dev.', 'sd', matrix)
	matrix
}

# COLLECT P-VALUES FROM MCMC FIXED EFFECTS DISTRIBUTION
mcpvals <- function(fixed_effects_mat) {
	nr <- nrow(fixed_effects_mat)
	prop <- colSums(fixed_effects_mat > 0)/nr
	ans <- t(matrix(2 * pmax(0.5/nr, pmin(prop, 1 - prop))))
	colnames(ans) <- colnames(fixed_effects_mat)
	rownames(ans) <- ""
	ans
}

se <- function(x) {sd(x)/sqrt(length(x))}

#####################################################################################




main_dir <- "/Volumes/Phillips/bars"
data_file_identifier <- "_MODEL2_whole_trial_n197_masked.nii.gz" # portion of data file (after cope number) indicating how many model & subjects
data_dir <- paste(c(main_dir, "/analysis/0_level2_volumes/"), sep="", collapse="")
output_dir <- paste(c(main_dir, "/analysis/2_level2_volumes_mcmc/MODEL2/boot_lmm/"), sep="", collapse="")

cat("\tChecking existence of cope file\n")
cope_name <- paste(c("cope", current_cope, data_file_identifier), sep="", collapse="")
cope_full_name <- paste(c(data_dir, cope_name), sep="", collapse="")

if (! file.exists(cope_full_name)) {
	print(cope_name)
	stop()
} else {
	cat("\t and it's a go\n")
}


# # DEFINE MODEL TO RUN



# RANDOM ERRORS ARE SPECIFIED DIRECTLY IN THE CALL TO lme()
 if (model_name == "age2") {	
	model <- "current_data ~ Age_mean_adj + I(Age_mean_adj^2)"
	
}else if (model_name == "age") {
	model <- "current_data ~ Age_mean_adj"

} else if (model_name == "intercept") {
	model <- "current_data ~ 1"

} else if (model_name == "invage") {
	model <- "current_data ~ I(1/Age)"
	
} else if (model_name == "age_beh") {
	model <- "current_data ~ Age + beh + Age:beh"

} else if (model_name == "age_sex") {
	model <- "current_data ~ Age + sex + Age:sex"

} else if (model_name == "age_beh3") {
	model <- "current_data ~ Age + beh + Age:beh"

} else if (model_name == "age2_beh3") {
	model <- "current_data ~ Age + I(Age^2) + beh + Age:beh"
}



model <- as.formula(model)




#####################################################################################







# LOAD CURRENT COPE VOLUME
vol <- load_gz_volumes(cope_full_name)
output_txt <- paste(c("\tLoading", cope_full_name), sep="", collapse=" ") 
write(output_txt, stdout())


#n_voxels <- (91*109*91)
x_dim <- 91
y_dim <- 109
z_dim <- 91
t_dim <- 197


# PULL TOGETHER SUBJECT LIST AND ACCURACY INFORMATION
info_list <- read.table(paste(c(data_dir, "_info_list_n197.txt"), sep="", collapse=""), header=T)
info_list$Age_mean_adj <- info_list$Age - mean(info_list$Age)

accuracy_name <- paste( c(main_dir, "/analysis/accuracy_good_fmri_subs.txt"), sep="", collapse="")
accuracy = read.table(accuracy_name, header=T)
accuracy$total <- rowMeans(accuracy[,c("rew_acc","neut_acc", "pun_acc")])

# ARCSINE TRANSFORM ACCURACY DATA
accuracy[,c("rew_acc", "neut_acc", "pun_acc", "total")] <- asin(accuracy[,c("rew_acc", "neut_acc", "pun_acc", "total")])

info_list <- merge(info_list, accuracy)

if (model_name == "age_beh") {
	if (current_cope %in% c(1,4)) {
		beh <- info_list[,"rew_acc"]
	} else if (current_cope %in% c(2)) {
		beh <- info_list[,"neut_acc"] - mean(info_list[,"neut_acc"])
	} else if (current_cope %in% c(3,5)) {
		beh <- info_list[,"pun_acc"]
	}
	info_list$beh <- beh
}
if (model_name == "age_beh2") {
	if (current_cope %in% c(1,4)) {
		beh <- info_list[,"rew_acc"]
	} else if (current_cope %in% c(2)) {
		beh <- info_list[,"neut_acc"] - mean(info_list[,"neut_acc"])
	} else if (current_cope %in% c(3,5)) {
		beh <- info_list[,"pun_acc"]
	}
	info_list$beh <- beh
}
if (model_name == "age_beh3") {
	if (current_cope %in% c(1,4)) {
		beh <- info_list[,"rew_acc"]
	} else if (current_cope %in% c(2)) {
		beh <- info_list[,"neut_acc"]
	} else if (current_cope %in% c(3,5)) {
		beh <- info_list[,"pun_acc"]
	}
	
	lm.1 <- lm(beh ~ Age, data=info_list)
	info_list$beh <- lm.1$resid
}
if (model_name == "age2_beh3") {
	if (current_cope %in% c(1,4)) {
		beh <- info_list[,"rew_acc"]
	} else if (current_cope %in% c(2)) {
		beh <- info_list[,"neut_acc"]
	} else if (current_cope %in% c(3,5)) {
		beh <- info_list[,"pun_acc"]
	}
	
	lm.1 <- lm(beh ~ Age + I(Age^2), data=info_list)
	info_list$beh <- lm.1$resid
}
#####################################################################################
# CHECK Z-SLICES FOR NON-ZERO VOXELS



# FOR TESTING ON WALLACE (ALSO ADDED TO SAVE NIFTI FILES AT EACH ITERATION OF Z)
subdir <- paste(c(output_dir, model_name), sep="", collapse="")
if (! file.exists(subdir)) {
	dir.create(file.path(subdir), recursive=T)
}
# LOAD NIFTI HEADER INFORMATION (hdr)
output_name <- paste(c(subdir, "/", model_name, "_cope", current_cope), sep="", collapse="")

if (file.exists(output_name)) {
	load(output_name)
	if (length(GLM_OUTPUT_VOLUMES) == 0) {
		z_start <- 1		
	} else {
		finished_vox <- which(GLM_OUTPUT_VOLUMES[["intcpt_pe"]]>0, arr.ind=T)
		finished_slice <- finished_vox[nrow(finished_vox),3]
		z_start <- finished_slice + 1
	}
} else {
	GLM_OUTPUT_VOLUMES <- list()
	z_start <- 1
}

# usually starts at z=23
for (z in z_start:z_dim) {
	
	# CHECK TO SEE IF THE Z-DIMENSION HAS ANY DATA
	if (sum(vol[[2]][,,z,] != 0)) { # IF THERE IS ANY DATA AT ALL IN THE SLICE, PROCEED

		
		current_array <- matrix(vol[[2]][,,z,], nrow=t_dim, ncol=x_dim*y_dim, byrow=T)				
		current_array <- round(current_array, digits = 5)
		non_zero_cols <- which(apply(current_array, 2, function(x) { sum(x != 0) == 197 }) )
		
		# IF THERE IS ONLY PARTIAL DATA, MOVE ON TO NEXT SLICE
		if (length(non_zero_cols) == 0) {
			next
		}
		

		# IF PROCESSING ON SKYNET, USE MAX_NUM_JOBS FILE TO DETERMINE N_PROCESSORS
		if (Sys.info()["nodename"] == "skynet") {
			number_of_processors_to_use <- as.numeric(readLines(file("/Volumes/Phillips/bars/code_t2/max_num_jobs.txt")))
		} else {
			number_of_processors_to_use <- 7
		}
		registerDoMC(number_of_processors_to_use) # number of processors to use


		txt <- paste(c("\t", model_name, ": cope ", current_cope, " slice ", z), sep="", collapse="")
		write(txt, stdout()) # there is data
		
		# doMC LOOP
		
		#counter <- 1
		#list_output <- list()
		list_output <- foreach(i = non_zero_cols) %dopar% { 
		#for (i in non_zero_cols) {
			
			
			info_list$current_data <- current_array[,i]

			
			# USE LINEAR MIXED EFFECTS MODEL
			# lmer(model, data=info_list) == lme(current_data ~ Age + beh + Age:beh, random = ~ 1|subjectID, data=info_list)
			#current_output <- lmer(model, data=info_list) # use for lmer object
			
			
			bootcoef <- function(data,index){
				dat<-data[index,]
				mod<-lme(model, random = ~ 1|subjectID, data=dat)
				fixef(mod)
			}
			
			#current_output <- lme(current_data ~ Age + beh + Age:beh, random = ~ 1|subjectID, data=info_list) # use for nlme object
			lme.1 <- lme(model, random = ~ 1|subjectID, data=info_list)

			DFs <- summary(lme.1)$fixDF$X			
			coefficients.names <- names(lme.1$coefficients$fixed)
			coefficients.names <- change_stat_value_names(coefficients.names)

						
			boot.out <- boot(info_list, bootcoef, 500)

			pe <- apply(boot.out$t, 2, mean)
			bias <- pe - boot.out$t0
			se <- apply(boot.out$t, 2, sd)
			tval <- pe/se
			pval <- (1-pt(abs(tval), df=DFs))*2
			zval <- qnorm(1 - (pval/2)) * sign(pe) # returns pos or neg zvals
			
			
			
			current_output.df <- as.data.frame(rbind(pe, bias, se, tval, pval, zval))
			names(current_output.df) <- coefficients.names
			current_output.df$names <- row.names(current_output.df)
			
			current_output.df <- melt(current_output.df, id=c("names"))
			
			# NOT AS LENGTHY PROCESS FOR GETTING COEFFICIENT INFO INTO READABLE DATAFRAME FORMAT				
			
			current_output.df $names <- apply(current_output.df, 1, function(x) { paste(c(x["variable"], x["names"]), sep="", collapse="_") })
			current_output.df <- current_output.df[,c("names", "value")]

			# ADD AIC TO DATAFRAME
			aic <- summary(lme.1)$AIC
			current_output.df <- rbind(current_output.df, data.frame(names=as.factor("aic"), value=aic))
			
		
			# UNCOMMENT FOR DOMC	
			current_output.df
			
			# UNCOMMENT FOR TRADITIONAL LOOP
			#list_output[[counter]] <- current_output.df
			#counter <- counter + 1
			
			
		} # END DO LOOP
		
		# formula to identify the current vector of values given the flattened matrix with dim(t_dim, x_dim*y_dim)
		# your_flattened_matrix[,(y(in FSL vox))*x_dim + x(FSL vox) + 1]



		stat_names = list_output[[1]]$names
		n_stats = length(stat_names)

		# SET UP EMPTY LIST OF MATRICES FOR STATS IF NECESSARY
		if (length(GLM_OUTPUT_VOLUMES) == 0) { # list is empty and we need to build empty volumes
			for (i in 1:n_stats) {
				GLM_OUTPUT_VOLUMES[[ stat_names[i] ]] <- array(data=0, dim=c(x_dim, y_dim, z_dim))
			}
		}
		
		# CREATE A BLANK ARRAY OF ZEROS, THEN PASTE OUTPUT INTO THIS ARRAY
		blank <- array(data=0, dim=c(n_stats, x_dim*y_dim))
		for (i in 1:length(non_zero_cols)) {
			blank[,non_zero_cols[i]] <- list_output[[i]]$value
		}
		
		for (i in 1:n_stats) {
			GLM_OUTPUT_VOLUMES[[ stat_names[i] ]] [,,z] <- array(data=blank[i,], dim=c(x_dim,y_dim))
		}

	
				
	} # END IF Z IS NOT EMPTY

	# FOR TESTING ON WALLACE (ALSO ADDED TO SAVE NIFTI FILES AT EACH ITERATION OF Z)
	subdir <- paste(c(output_dir, model_name), sep="", collapse="")
	if (! file.exists(subdir)) {
			dir.create(file.path(subdir), recursive=T)
	}
	# LOAD NIFTI HEADER INFORMATION (hdr)
	output_name <- paste(c(subdir, "/", model_name, "_cope", current_cope), sep="", collapse="")
	save(GLM_OUTPUT_VOLUMES, file=output_name)

	
} # END CYCLE THROUGH Z SLICES


#############

subdir <- paste(c(output_dir, model_name), sep="", collapse="")

if (! file.exists(subdir)) {
		dir.create(file.path(subdir), recursive=T)
}

# LOAD NIFTI HEADER INFORMATION (hdr)
load("/Volumes/Phillips/bars/code_t2/AnalyzeFMRI_nifti_hdr.R")
for (stat_name in names(GLM_OUTPUT_VOLUMES)) {
	output_name <- paste(c(subdir, "/", model_name, "_cope", current_cope, "_", stat_name), sep="", collapse="")
	hdr$file.name <- stat_name
	f.write.nifti(GLM_OUTPUT_VOLUMES[[stat_name]], output_name, size="float", hdr, nii=TRUE)
}

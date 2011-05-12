diffDataNew<-function(storage.mode = c("lockedEnvironment","environment","list"),...){
	storage.mode<-match.arg(storage.mode) ##defaults to "lockedEnvironment"
	diffData<-switch(storage.mode,
					lockedEnvironment =,
					environment = new.env(parent=emptyenv()),
					list = list())
	arglist <-list(...)
	if((length(arglist)>0L) && ((is.null(names(arglist))) || any(names(arglist)=="")))
		stop("One or more arguments of this function are not named - please provide named arguments!")
	
	for (nm in names(arglist)) diffData[[nm]] <- arglist[[nm]]
	if (storage.mode == "lockedEnvironment") diffDataEnvLock(diffData)
	msg <- diffDataValidMembers(diffData)
	if(!is.logical(msg)) stop(msg)
	diffData
}

diffDataValidMembers <- function(diffData, required) {
	msg <- NULL
	eltNames <-
			if ("list" == diffDataStorageMode(diffData)) names(diffData)
			else ls(diffData)
	if (!missing(required)) {
		absent <- required[!required %in% eltNames]
		if (length(absent) != 0)
			msg <- paste(msg,
					paste("'DiffData' missing'",absent,"'",sep="",collapse="\n\t"),
					sep="\n")
	}
	dimsOk <-
			sapply(eltNames, function(elt)
						tryCatch(length(dim(diffData[[elt]]))>1,
								error=function(err) FALSE))
	if(!all(dimsOk))
		msg <- c(msg, paste("'DiffData' elements with invalid dimensions: '",
						paste(eltNames[!dimsOk], collapse="' '"), "'", sep=""))
	if (length(diffData)>1) {
		eltRowNames <- rownames(diffData[[eltNames[[1]]]])
		rowNamesOk <-
				all(sapply(eltNames[-1], function(elt)
									all(eltRowNames == rownames(diffData[[elt]]))))
		if(!rowNamesOk)
			msg <- c(msg, "'diffData' elements with different rowNames")
	}
	if(is.null(msg)) TRUE else msg
								
}

diffDataStorageMode <- function(object) {
	if (is(object, "list"))
		"list"
	else if (environmentIsLocked(object))
		"lockedEnvironment"
	else
		"environment"
}

setMethod("storageMode", "DiffData", diffDataStorageMode)


diffDataStorageModeReplace <- function(object, value) {
	storageMode <- diffDataStorageMode(object)
	if (storageMode == value) return(object)
	names <- if (storageMode == "list") names(object) else ls(object)
	switch(value,
			lockedEnvironment = {
				diffData <- new.env(parent=emptyenv())
				for (nm in names) diffData[[nm]] <- object[[nm]]
				diffDataEnvLock(diffData)
				diffData
			}, environment = {
				diffData <- new.env(parent=emptyenv())
				for (nm in names) diffData[[nm]] <- object[[nm]]
				diffData
			}, list = as.list(object))
}

setReplaceMethod("storageMode",
		signature=c(object="DiffData", value="character"),
		diffDataStorageModeReplace)

diffDataEnvLock <- function(diffData)
	lockEnvironment(diffData, bindings=TRUE)

diffDataSubsetElements <-
		function(object, elts, storageMode = diffDataStorageMode(object)) {
	if (any(duplicated(elts)))
		stop("'DiffData' element names must be unique")
	names <-
			if (storageMode(object)=="list") names(object)
			else ls(object)
	if (!all(elts %in% names))
		stop("'DiffData' missing elements: '",
				paste(elts[!elts %in% names], collapse="', '", sep=""), "'")
	switch(storageMode,
			lockedEnvironment = {
				diffData <- new.env(parent = emptyenv())
				for (nm in elts) diffData[[nm]] <- object[[nm]]
				diffDataEnvLock(diffData)
				diffData
			},
			environment = {
				diffData <- new.env(parent = emptyenv())
				for (nm in elts) diffData[[nm]] <- object[[nm]]
				diffData
			},
			list = {
				object[elts]
			})
}


setMethod("diffData",
		signature=signature(object="DiffData"),
		function(object) object)

.diffDataDimnames <- function(diffData) {
	switch(storageMode(diffData),
			lockedEnvironment=,
			environment=eapply(diffData, dimnames),
			list=lapply(diffData, dimnames))
}

setMethod("sampleNames", signature(object="DiffData"),
		function(object) {
			if (!length(object))
				return(character(0))
			safe.colnames <-
					function(x) if (ncol(x) == 0) character(0) else colnames(x)
			switch(diffDataStorageMode(object),
					list=safe.colnames(object[[1]]),
					safe.colnames(object[[ls(object)[1]]]))
		})

setReplaceMethod("sampleNames", c("DiffData", "ANY"), function(object, value) {
			dims <- 
					switch(diffDataStorageMode(object),
							lockedEnvironment=,
							environment = eapply(object, ncol),
							list = lapply(object, ncol))
			if (length(dims)==0 && length(value) !=0)
				return(object)                    # early exit; no samples to name
			if (!all(dims==length(value)))
				stop("'value' length (", length(value),
						") must equal sample number in DiffData (",dims[[1]], ")")
			switch(diffDataStorageMode(object),
					lockedEnvironment = {
						object <- copyEnv(object)
						for (nm in ls(object)) colnames(object[[nm]]) <- value
						diffDataEnvLock(object)
					},
					environment = for (nm in ls(object)) colnames(object[[nm]]) <- value,
					list = for (nm in names(object)) colnames(object[[nm]]) <- value
			)
			object
		})

setReplaceMethod("sampleNames",
		signature=signature(
				object="DiffData",
				value="list"),
		function(object, value) 
		{
			.names_found_unique <- function(names, table)
			{
				ok <- !is.null(names) && all(names %in% table) && 
						!any(duplicated(names))
				if (!ok)
					stop("'sampleNames' replacement list must have unique named elements corresponding to diffData element names")
			}
			switch(diffDataStorageMode(object),
					lockedEnvironment = {
						.names_found_unique(names(value), ls(object))
						object <- copyEnv(object)
						for (nm in names(value))
							colnames(object[[nm]]) <- value[[nm]]
						diffDataEnvLock(object)
					}, environment = {
						.names_found_unique(names(value), ls(object))
						for (nm in names(value))
							colnames(object[[nm]]) <- value[[nm]]
					}, list= {
						.names_found_unique(names(value), names(object))
						for (nm in names(value))
							colnames(object[[nm]]) <- value[[nm]]
					})
			object
})

setMethod("featureNames", signature(object="DiffData"),
		function(object) {
			if (!length(object))
				return(character(0))
			safe.rownames <-
					function(x) if (nrow(x) == 0) character(0) else rownames(x)
			switch(diffDataStorageMode(object),
					list=safe.rownames(object[[1]]),
					safe.rownames(object[[ls(object)[1]]]))
		})


setReplaceMethod("featureNames", signature(object="DiffData", value="ANY"),
		function(object, value) {
	dims <- 
			switch(diffDataStorageMode(object),
					lockedEnvironment=,
					environment = eapply(object, nrow),
					list = lapply(object, nrow))
	if (length(dims)==0 && length(value) !=0)
		return(object)                    # early exit; no features to name
	if (!all(dims==length(value)))
		stop("'value' length (", length(value),
				") must equal feature number in DiffData (",dims[[1]], ")")
	switch(diffDataStorageMode(object),
			lockedEnvironment = {
				object <- copyEnv(object)
				for (nm in ls(object)) rownames(object[[nm]]) <- value
				diffDataEnvLock(object)
			},
			environment = for (nm in ls(object)) rownames(object[[nm]]) <- value,
			list = for (nm in names(object)) rownames(object[[nm]]) <- value,
	)
	object
})

setMethod("combine", c("DiffData", "DiffData"), function(x, y, ...) {
	storage.mode <- diffDataStorageMode(x)
	nmfunc <- if ("environment"==class(x)) ls else names
	
	if (diffDataStorageMode(y) != storage.mode)
		stop(paste("diffData must have same storage, but are ",
						storage.mode, ", ", diffDataStorageMode(y), sep=""))
	if (length(nmfunc(x)) != length(nmfunc(y)))
		stop("diffData have different numbers of elements:\n\t",
				paste(nmfunc(x), collapse=" "), "\n\t",
				paste(nmfunc(y), collapse=" "))
	if (!all(nmfunc(x) == nmfunc(y)))
		stop(paste("diffData have different element names:",
						paste(nmfunc(x), collapse=" "),
						paste(nmfunc(y), collapse=" "), sep="\n\t"))
	if ("list" == storage.mode) {
		dData <- lapply(names(x), function(nm) combine(x[[nm]],y[[nm]]))
		names(dData) <- names(x)
	} else {
		dData <- new.env(parent=emptyenv())
		for (nm in ls(x)) dData[[nm]] <- combine(x[[nm]], y[[nm]])
		if ("lockedEnvironment" == storage.mode) diffDataEnvLock(dData)
	}
	dData
})

diffDataDim <- function(object) {
	nms <- if (diffDataStorageMode(object) == "list") names(object) else ls(object)
	if ( length( nms ) == 0 ) return( NA )
	d <- dim( object[[ nms[[1]] ]])
	names(d) <- c( "Features", "Samples", rep("...", max(length(d)-2, 0)))
	d
}

#Not sure what this is for yet.
diffDataDims <- function(object) {
	nms <- if (diffDataStorageMode(object) == "list") names(object) else ls(object)
	if (length(nms) == 0)
		return(matrix(integer(0), nrow = 2, ncol = 0, 
						dimnames = list(c("Features", "Samples"), character(0))))
	d <- sapply(nms, function(i) dim(object[[i]]))
	rownames(d) <- c("Features", "Samples", rep("...", nrow(d)-2))
	colnames(d) <- nms
	d[,order(colnames(d)), drop=FALSE]
}
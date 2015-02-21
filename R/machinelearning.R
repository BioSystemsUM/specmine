#summary function for multi class with ROC metric
multiClassSummary <- function (data, lev = NULL, model = NULL){
    
    #Load Libraries
    require(Metrics)
    require(caret)
    
    #Check data
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
        stop("levels of observed and predicted data do not match")
    
    #Calculate custom one-vs-all stats for each class
    prob_stats <- lapply(levels(data[, "pred"]), function(class){
        
        #Grab one-vs-all data for the class
        pred <- ifelse(data[, "pred"] == class, 1, 0)
        obs <- ifelse(data[, "obs"] == class, 1, 0)
        prob <- data[,class]
        
        #Calculate one-vs-all AUC and logLoss and return
        cap_prob <- pmin(pmax(prob, .000001), .999999)
        prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
        names(prob_stats) <- c('ROC', 'logLoss')
        return(prob_stats)
    })
    prob_stats <- do.call(rbind, prob_stats)
    rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
    
    #Calculate confusion matrix-based statistics
    CM <- confusionMatrix(data[, "pred"], data[, "obs"])
    
    #Aggregate and average class-wise stats
    #Todo: add weights
    class_stats <- cbind(CM$byClass, prob_stats)
    class_stats <- colMeans(class_stats)
    
    #Aggregate overall stats
    overall_stats <- c(CM$overall)
    
    #Combine overall with class-wise stats and remove some stats we don't want
    stats <- c(overall_stats, class_stats)
    stats <- stats[! names(stats) %in% c('AccuracyNull',
                                         'Prevalence', 'Detection Prevalence')]
    
    #Clean names and return
    names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
    return(stats)
    
}


#MODEL SELECTION AND PREDICTION

#train a classifier and predict new samples
# dataser: data and metadata
# new.samples: dataframe with new samples to predict the class label
# column.class: metadata column class
# model: model to be used in training
# validation: validation method (boot, boot632, cv, repeatedcv, LOOCV, LGOCV, oob(only for random forests))
# num.folds: number of folds in validation
# num.repeats: number of repeats
# tunelength: number of levels for each tuning parameters
# tunegrid: dataframe with possible tuning values
train.and.predict = function(dataset, new.samples, column.class, model, validation, num.folds = 10, 
                             num.repeats = 10, tunelength = 10, tunegrid = NULL, metric = NULL, 
                             summary.function = defaultSummary) {
	train.result = train.classifier(dataset, column.class, model, validation, num.folds, num.repeats, 
                                  tunelength, tunegrid, metric, summary.function)
	predict.result = predict.samples(train.result, new.samples)
	result = list(train.result = train.result, predictions.result = predict.result)
	result
}

# train classifier
train.classifier = function(dataset, column.class, model, validation, num.folds = 10, 
                            num.repeats = 10, tunelength = 10, tunegrid = NULL, metric = NULL, summary.function = defaultSummary, class.in.metadata = T) {
  if(class.in.metadata)
	  train.result = trainClassifier(dataset$data, dataset$metadata[,column.class], model, validation, 
                                 num.folds, num.repeats, tunelength, tunegrid, metric, summary.function)
  else
    train.result = trainClassifier(dataset$data, dataset$data[column.class,], model, validation, 
                                   num.folds, num.repeats, tunelength, tunegrid, metric, summary.function)
	train.result
}


#model: caret classifier
#validation: caret resampling method -> boot, boot632, cv, repeatedcv, LOOCV, LGOCV, oob(only for random forests)
trainClassifier <- function(datamat, sampleclass, model, validation, num.folds = 10, num.repeats = 10, 
                            tunelength = 10, tunegrid = NULL, metric = NULL, summary.function = defaultSummary, class.in.metadata = T)
{
	require(caret)
	samples.df.ml = data.frame(t(datamat))
	rnames = gsub('[-\ ]','_',rownames(datamat))
	colnames(samples.df.ml) = paste("X",rnames,sep="")
	#names(samples.df.ml) = paste("X",rownames(datamat),sep="")
	rownames(samples.df.ml) = colnames(datamat)
	train.metric = metric
	if (class.in.metadata){
		samples.df.ml$class = sampleclass
		if (is.null(metric) && is.factor(samples.df.ml$class)){
			train.metric = "Accuracy"
		} else if (is.null(metric) && !is.factor(samples.df.ml$class)){
			train.metric = "RMSE"
		}
	} else {
		if (is.null(metric) && is.factor(samples.df.ml$sampleclass)){
			train.metric = "Accuracy"
		} else if (is.null(metric) && !is.factor(samples.df.ml$sampleclass)){
			train.metric = "RMSE"
		}
	}
	if (train.metric == "ROC") class.probs = T
	else class.probs = F
	train.control = trainControl(method=validation, number = num.folds, repeats=num.repeats, classProbs = class.probs, 
					summaryFunction= summary.function)
	if (class.in.metadata) 
		result.train = train(class ~., data = samples.df.ml, method=model, tuneLength = tunelength, metric = train.metric,
						trControl = train.control, tuneGrid = tunegrid)
	else
		result.train = train(sampleclass ~., data = samples.df.ml, method=model, tuneLength = tunelength, metric = train.metric, 
						trControl = train.control, tuneGrid = tunegrid)
	result.train
} 

#predict new samples
#train.result: result object of training a classifier
#new.samples: dataframe with new samples
predict.samples = function(train.result, new.samples){
	new.samples.df = data.frame(t(new.samples))
	rnames = gsub('[-\ ]','_',rownames(new.samples))
	colnames(new.samples.df) = paste("X",rnames, sep="")
	predict.result = predict(train.result, newdata = new.samples.df)
	result = data.frame(sample = rownames(new.samples.df), predicted.class = predict.result)
	result$sample = as.character(result$sample)
	result
}


train.models.performance = function(dataset, models, column.class, validation, num.folds = 10, 
                                    num.repeats = 10, tunelength = 10, tunegrid = NULL, metric = NULL, summary.function = "default", class.in.metadata = T){
	result.df = NULL
	classification.flag = FALSE
	vars.imp = list()
	final.result = list()
  full.results = list()
  if ( is.factor(dataset$metadata[,column.class])){
	classification.flag = TRUE
  	confusion.matrices = list()
  }
  if (is.character(summary.function)){
	  if (!is.null(metric) && metric == "ROC" && summary.function == "default"){
		summary.function = multiClassSummary
	  } else if (summary.function == "default"){
		summary.function = defaultSummary
	  }
  }
  
  best.tunes = list()
  final.models= list()
	for (i in 1:length(models)){
		train.result = train.classifier(dataset, column.class, models[i], validation, num.folds, 
                                    num.repeats, tunelength, tunegrid, metric, summary.function, class.in.metadata = class.in.metadata)
		vips = var.importance(train.result)
		rownames(vips) = substring(rownames(vips), 2, nchar(rownames(vips)))
		vips$Mean = apply(vips, 1, mean) 
		bestTune = train.result$bestTune
		result.df = rbind(result.df, train.result$result[rownames(bestTune),-1])
		vars.imp[[i]] = vips[order(vips$Mean, decreasing=T),]
		vips = NULL
    full.results[[i]] = train.result$results
    if (classification.flag) confusion.matrices[[i]] = try(confusionMatrix(train.result), TRUE)
    best.tunes[[i]] = train.result$bestTune
    final.models[[i]] = train.result$finalModel
	}
	rownames(result.df) = models
	names(vars.imp) = models
	names(full.results) = models
  if (classification.flag) names(confusion.matrices) = models
  names(best.tunes) = models
  names(final.models) = models
  result.df = result.df[,colnames(result.df) %in% c("RMSE", "Rsquared", "RMSESD","RsquaredSD","Accuracy","AccuracySD","Kappa","KappaSD",
													"ROC","Sensitivity","Specificity","SensitivitySD","SpecificitySD","ROCSD")]
  final.result$performance = result.df
	final.result$vips = vars.imp
  final.result$full.results = full.results
  final.result$best.tunes = best.tunes
  if (classification.flag) final.result$confusion.matrices = confusion.matrices
  final.result$final.models = final.models
	final.result
}

# VARIABLE IMPORTANCE

var.importance = function(train.result){
	require(caret)
	vip = varImp(train.result)
	vip$importance
}

summary.var.importance = function(performances, number.rows){
	for(i in 1:length(performances$vips)){
		performances$vips[[i]] = performances$vips[[i]][1:number.rows,]
	}
	performances$vips
}

# PCA PLOTS

"pca.plot.3d" = function(dataset, model, var.class, pcas = 1:3, pt.cex = 1, cex = 1, ...) {
  require(scatterplot3d)
  if (length(pcas) != 3) stop("Wrong dimension in parameter pcas")
  classes = dataset$metadata[,var.class]
  scatterplot3d(model$scores[,pcas], color=as.integer(classes), pch=17, 
                xlab = "Component 1", ylab = "Component 2", zlab = "Component 3")
  legend(-1.5, 2.5, levels(classes), col = 1:length(classes), pt.cex = pt.cex, pch= 17, cex = cex, ...)
}

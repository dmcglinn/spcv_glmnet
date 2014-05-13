#' Get the row indices for the training data sets
#' 
#' For each sample in the dataset a set of indices for the training
#' data to compare to that sample is provided based on a threshold
#' distance. If the threshold distance is zero then all other samples
#' are included in the training set. 
#' 
#' @param coords 2 dimensional matrix of spatial coordinates
#' @param dist_thres the threshold distance, if longlat = TRUE
#' this is in units of kilometers
#' @longlat boolean, if TRUE great circle distances are computed 
#' rather than euclidean distances
#' @export
#' @examples
#' coords = expand.grid(1:4, 1:4)
#' training_rows = get_training_row(coords, 3)
#' coords[training_rows[[i]]
get_training_rows = function(coords, dist_thres, longlat=FALSE) {
    ## Computing the distance matrix from the data:
    if (longlat) {
        require(sp)
        dist_matrix = spDists(coords, longlat=TRUE)
    }
    else 
        dist_matrix = as.matrix(dist(coords))
    ## Initializing the row indices of the dataset to be used in training 
    training_rows = list()
    ## Creating the sets of training indices
    for (i in 1:nrow(dist_matrix)) {
        # Keeping only the observations far enough of the i-st observation by using the threshold distance
        num_cell = which(dist_matrix[i, ] > dist_thres)
        training_rows[[i]] = num_cell
    }
    return(training_rows)
}

#' Spatial cross-validation for glmnet
#'
#' Does leave-one-out cross-validation for glmnet,  produces a plot,  and returns
#' a value for lambda
#'
#' @param x x matrix as in glmnet 
#' @param y response y as in glmnet.
#' @param coords the spatial coordinates of the samples
#' @param dist_thres the spatial distance threshold
#' @param longlat boolean, if TRUE great circle distances are computed 
#' rather than euclidean distances
#' @param ... other variables that can be supplied to cv.glmnet
#' @export
#' @examples
#' x = matrix(rnorm(100*20),100,20)
#' y = rnorm(100)
#' coords = expand.grid(1:10, 1:10)
#' nsp_cv = cv.glmnet(x,y)
#' sp_cv = spcv.glmnet(x, y, coords, 5, FALSE)
#' par(mfrow=c(1,2))
#' plot(nsp_cv)
#' plot(sp_cv)
spcv.glmnet = function(x, y, coords, dist_thres, longlat, weights,
                       offset=NULL, lambda=NULL,
                       type.measure=c("mse", "deviance", "class", "auc", "mae"),
                       nfolds=10, training_rows, grouped=TRUE, keep=FALSE,
                       parallel=FALSE, ...){
    require(glmnet)
    if (missing(type.measure))
        type.measure = "default"
    else 
        type.measure = match.arg(type.measure)
    if (!is.null(lambda) && length(lambda) < 2)
        stop("Need more than one value of lambda for cv.glmnet")
    N = nrow(x)
    if (missing(weights))
        weights = rep(1.0, N)
    else 
        weights = as.double(weights)
    ## Fit the model once to get dimensions etc of output
    y = drop(y) # we dont like matrix responses unless we need them
    ## Next we construct a call,  that could recreate a glmnet object - tricky
    ## This if for predict,  exact=TRUE
    glmnet.call = match.call(expand.dots = TRUE)
    which = match(c("type.measure", "nfolds", "training_rows", "grouped", "keep"),
                  names(glmnet.call), FALSE)
    if (any(which))
        glmnet.call = glmnet.call[-which]
    glmnet.call[[1]] = as.name("glmnet") 
    glmnet.object = glmnet(x, y, weights=weights, offset=offset,
                           lambda=lambda, ...)
    glmnet.object$call = glmnet.call
    is.offset = glmnet.object$offset
    lambda = glmnet.object$lambda
    if (inherits(glmnet.object, "multnet")){
        nz = predict(glmnet.object, type="nonzero")
        nz = sapply(nz, function(x) sapply(x, length))
        nz = ceiling(apply(nz, 1, median))
    }
    else  
        nz = sapply(predict(glmnet.object, type="nonzero"), length)
    #if (missing(foldid)) 
    #    foldid = sample(rep(seq(nfolds), length=N))
    #else 
    #    nfolds = max(foldid)
    training_rows = get_training_rows(coords, dist_thres, longlat)
    ## drop empty traning sets
    training_rows = training_rows[sapply(training_rows, length) > 0]
    nfolds = length(training_rows)
    #if (nfolds < 3)
    #    stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))
    ## Now fit the nfold models and store them
    ## First try and do it using foreach if parallel is TRUE
    if  (parallel && require(foreach)) {
        outlist = foreach (i=seq(nfolds), .packages=c("glmnet")) %dopar% {
            which = foldid == i
            if (is.matrix(y))
                y_sub=y[!which, ]
            else 
                y_sub = y[!which]
            if (is.offset)
                offset_sub = as.matrix(offset)[!which, ]
            else 
                offset_sub=NULL
            glmnet(x[!which, , drop=FALSE], y_sub, lambda=lambda,
                   offset=offset_sub, weights=weights[!which], ...)
        }
    }
    else {
        for (i in seq(nfolds)) {
            rows = training_rows[[i]]
            if (is.matrix(y))
                y_sub = y[rows, ]
            else 
                y_sub = y[rows]
            if (is.offset)
                offset_sub = as.matrix(offset)[rows[[i]], ]
            else 
                offset_sub=NULL
            outlist[[i]] = glmnet(x[rows, , drop=FALSE], y_sub,
                                  lambda=lambda, offset=offset_sub,
                                  weights=weights[rows], ...)
        }
    }
    ## What to do depends on the type.measure and the model fit
    fun = paste("spcv", class(glmnet.object)[[1]], sep=".")
    cvstuff = do.call(fun, list(outlist, lambda, x, y, coords, dist_thres,
                                longlat, weights, offset, training_rows, 
                                type.measure, grouped, keep))
    cvm = cvstuff$cvm
    cvsd = cvstuff$cvsd
    cvname = cvstuff$name
    out = list(lambda=lambda, cvm=cvm, cvsd=cvsd, cvup=cvm + cvsd,
               cvlo=cvm - cvsd, nzero=nz, name=cvname, glmnet.fit=glmnet.object)
    if (keep)
        out = c(out, list(fit.preval=cvstuff$fit.preval, foldid=foldid))
    if (type.measure=="auc") 
        lamin = getmin(lambda, -cvm, cvsd)
    else 
        lamin = getmin(lambda, cvm, cvsd)
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    return(obj)
}

spcv.elnet = function (outlist, lambda, x, y, coords, dist_thres, longlat,
                       weights, offset, training_rows, 
                       type.measure, grouped, keep = FALSE) {
    typenames = c(deviance = "Mean-Squared Error", mse = "Mean-Squared Error", 
                  mae = "Mean Absolute Error")
    if (type.measure == "default") 
        type.measure = "mse"
    if (!match(type.measure, c("mse", "mae", "deviance"), FALSE)) {
        warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
        type.measure = "mse"
    }
    if (!is.null(offset)) 
        y = y - drop(offset)
    predmat = matrix(NA, length(y), length(lambda))
    nfolds = length(training_rows)
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
        rows = training_rows[[i]]
        fitobj = outlist[[i]]
        fitobj$offset = FALSE
        preds = predict(fitobj, x[-rows, , drop = FALSE])
        nlami = length(outlist[[i]]$lambda)
        predmat[-rows, seq(nlami)] = preds
        nlams[i] = nlami
    }
    N = length(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure, mse = (y - predmat)^2, 
                   deviance = (y - predmat)^2, mae = abs(y - predmat))
    if ((length(y)/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", 
                call. = FALSE)
        grouped = FALSE
    }
    if (grouped) {
        cvob = cvcompute(cvraw, weights, foldid, nlams)
        cvraw = cvob$cvraw
        weights = cvob$weights
        N = cvob$N
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
                      w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
    if (keep) 
        out$fit.preval = predmat
    out
}

spcv.fishnet = function(outlist, lambda, x, y, coords, dist_thres, longlat,
                        weights, offset, training_rows, type.measure,
                        grouped, keep = FALSE) {
    typenames = c(mse = "Mean-Squared Error", mae = "Mean Absolute Error", 
                  deviance = "Poisson Deviance")
    if (type.measure == "default") 
        type.measure = "deviance"
    if (!match(type.measure, c("mse", "mae", "deviance"), FALSE)) {
        warning("Only 'deviance', 'mse' or 'mae'  available for Poisson models; 'deviance' used")
        type.measure = "deviance"
    }
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else 
        is.offset = FALSE
    devi = function(y, eta) {
        deveta = y * eta - exp(eta)
        devy = y * log(y) - y
        devy[y == 0] = 0
        2 * (devy - deveta)
    }
    predmat = matrix(NA, length(y), length(lambda))
    nfolds = length(training_rows)
    nlams = double(nfolds)
    for (i in seq(nfolds)) {
        rows = training_rows[[i]]
        fitobj = outlist[[i]]
        if (is.offset) 
            off_sub = offset[-rows]
        preds = predict(fitobj, x[-rows, , drop = FALSE], offset = off_sub)
        nlami = length(outlist[[i]]$lambda)
        predmat[-rows, seq(nlami)] = preds
        nlams[i] = nlami
    }
    N = length(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure, mse = (y - exp(predmat))^2, 
                   mae = abs(y - exp(predmat)), deviance = devi(y, predmat))
    if ((length(y)/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", 
                call. = FALSE)
        grouped = FALSE
    }
    if (grouped) {
        cvob = cvcompute(cvraw, weights, foldid, nlams)
        cvraw = cvob$cvraw
        weights = cvob$weights
        N = cvob$N
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean, 
                      w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
    if (keep) 
        out$fit.preval = predmat
    out
}

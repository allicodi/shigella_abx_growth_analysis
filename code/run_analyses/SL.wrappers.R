# ------------------------------------------------------------------------------
# SuperLearner wrappers for abx-growth meta-analysis
# ------------------------------------------------------------------------------

# Custom outcome models:
# - GLM with spline for age, spline for bl HAZ
# - Step Forward all + spline for age, spline for HAZ
# - NOTE these assume names are age = age, enrollment HAZ = enr_haz
# - also force all_abx into models for outcome

# GLM with spline for age, HAZ
SL.glm.spline.age.haz <- function(Y, X, newX, family, obsWeights, ...){
  
  X_var_names <- colnames(X)
  
  # Remove age and enr_haz
  X_var_names <- setdiff(X_var_names, c("age", "enr_haz"))
  
  # Add splines for age and enr_haz (assuming the variables are named age and haz)
  formula <- as.formula(paste("Y ~ splines::ns(age, df = 3) + splines::ns(enr_haz, df = 3) +", 
                              paste(X_var_names, collapse = " + ")))
  
  # Fit model
  sl.glm.spline.age.haz_fit <- glm(formula,
                                       data = X, 
                                       family = family,
                                       weights = obsWeights)
  
  # Make predictions
  pred <- predict(sl.glm.spline.age.haz_fit, newdata = newX, type = 'response')
  
  # Store fitted model
  fit <- list(fitted_model = sl.glm.spline.age.haz_fit)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.glm.spline.age.haz"
  
  return(out)
}
predict.SL.glm.spline.age.haz <- function(object, newdata, ...){
  pred <- predict(object$fitted_model, newdata = newdata, type="response")
  return(pred)
}

# Step forward with spline for age, HAZ
SL.step.forward.spline.age.haz <- function (Y, X, newX, family, direction = "forward", trace = 0, k = 2, ...) {
  # Get var names
  X_var_names <- colnames(X)
  
  # Remove age and enr_haz
  X_var_names <- setdiff(X_var_names, c("age", "enr_haz"))
  
  # Add splines for age and enr_haz (assuming the variables are named age and haz)
  formula <- as.formula(paste("Y ~ splines::ns(age, df = 3) + splines::ns(enr_haz, df = 3) +", 
                              paste(X_var_names, collapse = " + ")))
  
  fit.glm <- glm(formula, data = X, family = family)
  fit.step <- step(glm(as.formula("Y ~ splines::ns(age, df = 3) + splines::ns(enr_haz, df = 3)"), data = X, family = family), scope = formula(fit.glm), 
      direction = direction, trace = trace, k = k)
  pred <- predict(fit.step, newdata = newX, type = "response")
  fit <- list(object = fit.step)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.step.forward.spline.age.haz")
  return(out)
}
predict.SL.step.forward.spline.age.haz <- function(object, newdata, ...){
    pred <- predict(object = object$object, newdata = newdata, type = "response")
    return(pred)
}

# Step forward with spline for age, HAZ - force abx in too
SL.step.forward.spline.age.haz.abx <- function (Y, X, newX, family, direction = "forward", trace = 0, k = 2, ...) {
  # Get var names
  X_var_names <- colnames(X)
  
  # Remove age and enr_haz
  X_var_names <- setdiff(X_var_names, c("age", "enr_haz"))
  
  # Add splines for age and enr_haz (assuming the variables are named age and haz)
  formula <- as.formula(paste("Y ~ splines::ns(age, df = 3) + splines::ns(enr_haz, df = 3) + ", 
                              paste(X_var_names, collapse = " + ")))
  
  fit.glm <- glm(formula, data = X, family = family)
  fit.step <- step(glm(as.formula("Y ~ splines::ns(age, df = 3) + splines::ns(enr_haz, df = 3) + all_abx"), data = X, family = family), scope = formula(fit.glm), 
                   direction = direction, trace = trace, k = k)
  pred <- predict(fit.step, newdata = newX, type = "response")
  fit <- list(object = fit.step)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.step.forward.spline.age.haz.abx")
  return(out)
}
predict.SL.step.forward.spline.age.haz.abx <- function(object, newdata, ...){
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  return(pred)
}

# Custom glmnet screener to handle abx as a factor
SL.screen.abx.glmnet <- function (
    Y, X, family, alpha = 1, 
    minscreen = 1, nfolds = 10, 
    nlambda = 100, ...
){
  SuperLearner:::.SL.require("glmnet")
  p_original_X <- ncol(X)
  
  abx_variable_label <- "all_abx"
  abx_variables_original_X <- which(grepl(abx_variable_label, colnames(X)))
  other_variables_original_X <- seq_len(p_original_X)[-abx_variables_original_X]
  
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
  }
  
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance", 
                             nfolds = nfolds, family = family$family, alpha = alpha, 
                             nlambda = nlambda)
  
  abx_variables <- which(grepl(abx_variable_label, colnames(X)))
  n_abx_variables <- length(abx_variables)
  other_variables <- seq_len(ncol(X))[-abx_variables]
  
  all_coefs <- as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1]
  whichVariable <- rep(FALSE, p_original_X)
  whichVariable[other_variables_original_X] <- all_coefs[other_variables] != 0
  whichVariable[abx_variables_original_X] <- any(abs(all_coefs[abx_variables]) > 1e-5)
  
  if (sum(whichVariable) < minscreen) {
    # just return abx if nothing passes screen
    whichVariable[abx_variables_original_X] <- TRUE
  }
  return(whichVariable)
}


# Custom glmnet screener to handle abx as a factor
SL.screen.abx.lt.min_prop <- function (
    Y, X, family, obsWeights, id, pathogen_suffix = "_new", min_prop = 0.05, ...
){
  
  # Get the indices of pathogen quantity variable (all end with pathogen_suffix)
  X_colnames <- colnames(X)
  idx_pathogens <- grep(pathogen_suffix, X_colnames)
  
  # create vector of whichVariables
  whichVariable <- rep(TRUE, ncol(X))
  
  # for pathogens, see if proportion detected < min_prop
  # if < min_prop, screen out
  screen_paths <- lapply(idx_pathogens, function(i){
    path_col <- X[,i]
    path_ind <- ifelse(path_col > 0, 1, 0)
    if(mean(path_ind, na.rm = TRUE) < min_prop){
      return(FALSE)
    } else{
      return(TRUE)
    }
  })
  
  whichVariable[idx_pathogens] <- unlist(screen_paths)
  
  return(whichVariable)
}



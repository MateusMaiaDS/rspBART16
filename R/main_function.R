# source("R/debugging_rspBART.R")
# rm(list=ls())
# source("R/other_functions.R")
# source("R/sim_functions.R")
# source("R/tree_functions.R")
# source("inst/debugging_rspBART.R")
# devtools::load_all()
set.seed(42)

# Creating the main function from the rspBART
rspBART <- function(x_train,
                    y_train,
                    x_test,
                    n_tree = 10,
                    node_min_size = 15,
                    n_mcmc = 2000,
                    n_burn = 500,
                    alpha = 0.95,
                    beta = 2,
                    df = 3,
                    sigquant = 0.9,
                    kappa = 2,
                    # Splines parameters
                    nIknots = 3,
                    dif_order = 1,
                    tau = 100,
                    scale_bool = TRUE,
                    stump = FALSE,
                    numcut = 100L, # Defining the grid of split rules
                    usequants = FALSE,
                    motrbart_bool = FALSE,
                    use_bs = FALSE,
                    plot_preview = FALSE,
                    all_var = TRUE,
                    scale_init = FALSE,
                    update_tau_beta = FALSE,
                    main_effects_pred = FALSE,
                    interaction_term =  FALSE,
                    interaction_list = NULL,
                    store_tree_fit = FALSE,
                    mle_prior = FALSE,
                    linero_sampler = TRUE,
                    varimportance_bool = FALSE
) {


  # # Another constrain
  # if(interaction_term & is.null(interaction_list)){
  #   stop("Define the interaction list.")
  # }


  # Making interaction FALSE for univariate
  if(NCOL(x_train)==1){
    interaction_term <- FALSE
  }

  if(isFALSE(interaction_term) & !is.null(interaction_list)){
    stop("Is a model without the interaction but a interaction list is defined")
  }

  # Verifying if x_train and x_test are matrices
  if(!is.data.frame(x_train) || !is.data.frame(x_test)){
    stop("Insert valid data.frame for both data and xnew.")
  }


  # Getting the valid
  dummy_x <- base_dummyVars(x_train)

  # Create a list
  if(length(dummy_x$facVars)!=0){
    for(i in 1:length(dummy_x$facVars)){
      # See if the levels of the test and train matches
      if(!all(levels(x_train[[dummy_x$facVars[i]]])==levels(x_test[[dummy_x$facVars[i]]]))){
        levels(x_test[[dummy_x$facVars[[i]]]]) <- levels(x_train[[dummy_x$facVars[[i]]]])
      }
      df_aux <- data.frame( x = x_train[,dummy_x$facVars[i]],y)
      formula_aux <- stats::aggregate(y~x,df_aux,mean)
      formula_aux$y <- rank(formula_aux$y)
      x_train[[dummy_x$facVars[i]]] <- as.numeric(factor(x_train[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

      # Doing the same for the test set
      x_test[[dummy_x$facVars[i]]] <- as.numeric(factor(x_test[[dummy_x$facVars[[i]]]], labels = c(formula_aux$y)))-1

    }
  }


  # Getting the train and test set
  x_train_scale <- as.matrix(x_train)
  x_test_scale <- as.matrix(x_test)

  # Scaling x
  x_min <- apply(as.matrix(x_train_scale),2,min)
  x_max <- apply(as.matrix(x_train_scale),2,max)

  # Storing the original
  x_train_original <- x_train
  x_test_original <- x_test


  # Normalising all the columns
  for(i in 1:ncol(x_train)){
    x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
    x_test_scale[,i] <- normalize_covariates_bart(y = x_test_scale[,i],a = x_min[i], b = x_max[i])
  }



  # Creating the numcuts matrix of splitting rules
  xcut_m <- matrix(NA,nrow = numcut,ncol = ncol(x_train_scale))
  for(i in 1:ncol(x_train_scale)){

    if(nrow(x_train_scale)<numcut){
      xcut_m[,i] <- sort(x_train_scale[,i])
    } else {
      xcut_m[,i] <- seq(min(x_train_scale[,i]),
                        max(x_train_scale[,i]),
                        length.out = numcut+2)[-c(1,numcut+2)]
    }
  }




  # =========================================================================================================
  # Getting the Splines Basis functions
  # =========================================================================================================
  # Creating a list of basis functions
  B_train_obj <- B_test_obj <- vector("list",length = length(dummy_x$continuousVars) + NCOL(interaction_list))

  # Setting new parameters for the spline
  ndx <- nIknots+1
  ord_ <- 4
  degree_ <- 3
  x_min_sp <- apply(x_train_scale,2,min)
  x_max_sp <- apply(x_train_scale,2,max)
  dx <- (x_max_sp-x_min_sp)/ndx

  # New_knots
  new_knots <- matrix()
  new_knots <- matrix(mapply(x_min_sp,x_max_sp,dx, FUN = function(MIN,MAX,DX){seq(from = MIN-(ord_-1)*DX, to = MAX+(ord_-1)*DX, by = DX)}), ncol = length(dummy_x$continuousVars)) # MIN and MAX are 0 and 1 respectively, because of the scale
  colnames(new_knots) <- dummy_x$continuousVars



  # Selecting the basis size.
  if(use_bs){
    basis_size <- (nIknots+degree_)
  } else {
    if(interaction_term){ # Checking the interaction
      # Checking all possible interactions
      basis_size <- (nrow(new_knots)-ord_)     # Change this value to the desired size of each sublist
    } else {
      basis_size <- (nrow(new_knots)-ord_)     # Change this value to the desired size of each sublist
    }
  }

  # Setting the indexes for the betas
  D_seq <- 1:(basis_size*length(dummy_x$continuousVars))  # Replace this with the columns of D

  # Creating a vector
  basis_subindex_main <- split(D_seq, rep(1:((length(D_seq) %/% basis_size)), each = basis_size, length.out = length(D_seq)))

  # Getting the interactions subindex only
  if(interaction_term){
    interaction_list <- combn(1:NCOL(x_train),2,simplify = TRUE)
  }

  D_int_seq <- (length(D_seq)+1):(basis_size*length(dummy_x$continuousVars)+(basis_size^2)*NCOL(interaction_list))
  interactions_subindex <- split(D_int_seq, rep(1:((length(D_int_seq) %/% basis_size)), each = basis_size^2, length.out = length(D_int_seq)))

  # Getting the final basis subindex
  if(interaction_term){
      basis_subindex <- append(basis_subindex_main,interactions_subindex)
  } else {
      basis_subindex <- basis_subindex_main
  }


  # Creating the natural B-spline for each predictor
  for(i in 1:length(dummy_x$continuousVars)){


    # Modify the basis only with respect to the main effects at the moment
    if(use_bs){
      B_train_obj[[i]] <- splines::bs(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                      df = nIknots+3,
                                      degree = 3,intercept = FALSE,
                                      Boundary.knots = c(-2,2)*range(x_train_scale[,dummy_x$continuousVars[i]]))

      B_train_obj[[i]] <- splines::bs(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE])
    } else {
      B_train_obj[[i]] <- splines::spline.des(x = x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                              knots = new_knots[,dummy_x$continuousVars[i]],
                                              ord = ord_,
                                              derivs = 0*x_train_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design
    }



    # For the test setting
    if(use_bs){
      B_test_obj[[i]] <- predict(object = B_train_obj[[i]],newx = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE])
    } else {
      B_test_obj[[i]] <- splines::spline.des(x = x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],
                                             knots = new_knots[,dummy_x$continuousVars[i]],
                                             ord = ord_,
                                             derivs = 0*x_test_scale[,dummy_x$continuousVars[i], drop = FALSE],outer.ok = TRUE)$design
    }


  }

  # Interaction matrix list to be used in the penalised
  interaction_matrix_list <- vector("list",length = NCOL(interaction_list))

  # Adding the interaction basis
  if(interaction_term){
    jj_ = length(dummy_x$continuousVars)
    for (jj in 1:NCOL(interaction_list)) {
      jj_ = jj_ +1
      B_train_obj[[jj_]] <- multiply_matrices_general(A = B_train_obj[[interaction_list[1,jj]]],B = B_train_obj[[interaction_list[2,jj]]])
      B_test_obj[[jj_]] <- multiply_matrices_general(A = B_test_obj[[interaction_list[1,jj]]],B = B_test_obj[[interaction_list[2,jj]]])
    }
  }

  # Renaming the columns of inthe basis subindex
  main_effects_names <- paste0(1:length(dummy_x$continuousVars))

  # Getting a string of all pairwise possibilities
  concatenate_columns <- function(matrix) {
    result <- apply(matrix, 2, function(col) paste(col, collapse = ""))
    return(result)
  }

  # Renaming only the interactions
  if(interaction_term){
    interaction_names <- concatenate_columns(interaction_list)
  }


  # Renaming only the interactions
  if(interaction_term){
    names(basis_subindex)[(length(dummy_x$continuousVars)+1):length(basis_subindex)] <- concatenate_columns(interaction_list)
    names(B_train_obj) <- names(B_test_obj) <- c(main_effects_names,interaction_names)
  } else {
    names(B_train_obj) <- names(B_test_obj) <- main_effects_names
  }

  # And the main effects


  if(motrbart_bool){

    # for(i)
    stop("MOTR-BART is not available anymore since we are not setting up intercepts")

    B_train <- x_train_scale
    D_test <- x_test_scale

    basis_size <- 1 # Change this value to the desired size of each sublist
    D_seq <- 1:ncol(D_train)  # Replace this with the columns of D

    # Creating a vector
    basis_subindex <- split(D_seq, rep(1:(length(D_seq) %/% basis_size), each = basis_size, length.out = length(D_seq)))
  }


  # Scaling the y
  min_y <- min(y_train)
  max_y <- max(y_train)

  # Getting the min and max for each column
  min_x <- apply(x_train_scale,2,min)
  max_x <- apply(x_train_scale, 2, max)

  # Scaling "y"
  if(scale_bool){
    y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)

    # New update
    tau_mu <- 4*n_tree*(kappa^2)
  } else {
    y_scale <- y_train

    # New parameter update
    tau_mu <- (4*n_tree*(kappa^2))/((max_y-min_y)^2)


  }




  # Getting the naive sigma value
  nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

  # Calculating tau hyperparam
  a_tau <- df/2

  # Calculating lambda
  qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
  lambda <- (nsigma*nsigma*qchi)/df
  d_tau <- (lambda*df)/2


  # Getting hyperparameters for \tau_beta_j
  tau_beta_df <- 3
  a_tau_beta_j <- tau_beta_df/2
  sigquant_beta <- 0.9
  nsigma_beta <- (tau_mu)^(-1/2)
  # nsigma_beta <- (tau_mu)^(-1/2)

  # nsigma_beta <- 1

  # Calculating lambda
  qchi_beta <- stats::qchisq(p = 1-sigquant_beta,df = df,lower.tail = 1,ncp = 0)
  lambda_beta <- (nsigma_beta*nsigma_beta*qchi_beta)/tau_beta_df
  d_tau_beta_j <- (lambda_beta*tau_beta_df)/2



  # Visualising the prior
  if(plot_preview){
    rgamma(n = 1000,shape = a_tau_beta_j,rate = d_tau_beta_j) |> density() |> plot(main = "density prior \tau_beta_j")
    mean(rgamma(n = 1000,shape = a_tau_beta_j,rate = d_tau_beta_j))
  }

  # Call the bart function
  tau_init <- nsigma^(-2)

  mu_init <- mean(y_scale)

  # Creating the vector that stores all trees
  all_tree_post <- vector("list",length = round(n_mcmc-n_burn))


  # =====================================================================
  # ========= From here I gonna initialise the BART function itself =====
  # =====================================================================
  n_post <- (n_mcmc-n_burn)
  all_trees <- vector("list", n_mcmc)
  all_betas <- vector("list",n_mcmc)
  if(interaction_term){
    tau_beta <- rep(tau_mu,length(dummy_x$continuousVars)+NCOL(interaction_list))
    # tau_beta <- rep(tau_mu, NCOL(x_train_scale))
  } else {
    tau_beta <- rep(tau_mu,length(dummy_x$continuousVars))
  }

  # all_df <- cbind(x_train,y_train)
  # mod <- mgcv::gam(y_train ~ bs(x.10), data = all_df)
  # mod$coefficients

  # In this first scenario we are going to work with a single value of \tau
  if(interaction_term){
    all_tau_beta <- matrix(NA, nrow = (n_mcmc), ncol = NCOL(x_train_scale)+NCOL(interaction_list))
  } else {
    all_tau_beta <- matrix(NA, nrow = (n_mcmc), ncol = NCOL(x_train_scale))
  }
  # all_delta <- numeric(n_mcmc)
  all_tau <- numeric(n_mcmc)

  all_y_hat <- matrix(NA,nrow = n_mcmc,ncol = nrow(x_train_scale))
  all_y_hat_test <- matrix(NA, nrow = n_mcmc, ncol = nrow(x_test_scale))
  all_trees_fit <- vector("list",n_mcmc)
  all_trees <- vector("list",n_mcmc)
  forest <- vector("list",n_tree)

  # Partial component pieces
  # (OPTIONAL)
  if(store_tree_fit){
    partial_train_fits <- vector("list", n_tree)
  }

  proposal_outcomes <- setNames(data.frame(matrix(nrow = 0, ncol = 6)),
                                c("tree_number" , "proposal", "status","mcmc_iter", "new_tree_loglike", "old_tree_loglike"))
  all_train_indexes <- data.frame(matrix(data = NA,nrow = nrow(xcut_m),ncol = ncol(xcut_m)))

  # Creating a list to store all the main effects for the sum of trees
  if(main_effects_pred){

    # Calculating the main effects
    if(interaction_term){

      # Interacting
      main_effects_train_list_norm <- main_effects_test_list_norm <- main_effects_train_list <- main_effects_test_list <- vector("list", length = length(dummy_x$continuousVars)+NCOL(interaction_list))

      for(list_size in 1:length(main_effects_train_list)){
        main_effects_train_list[[list_size]] <- main_effects_train_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_train))
        main_effects_test_list[[list_size]] <- main_effects_test_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_test))
      }

      # Renaming the main list
      if(interaction_term){
        internames_ <- apply(interaction_list,2,function(vars_){paste0("x",paste0(vars_,collapse = ":"))})
        names(main_effects_train_list_norm) <- names(main_effects_test_list_norm) <- names(main_effects_train_list) <- names(main_effects_test_list) <- c(dummy_x$continuousVars, internames_)
      } else {
        names(main_effects_train_list_norm) <- names(main_effects_test_list_norm) <- names(main_effects_train_list) <- names(main_effects_test_list) <-dummy_x$continuousVars
      }


    } else { # In case there are no interactions
      main_effects_train_list_norm <- main_effects_test_list_norm <- main_effects_train_list <- main_effects_test_list <- vector("list", length = length(dummy_x$continuousVars))

      for(list_size in 1:length(dummy_x$continuousVars)){
        main_effects_train_list[[list_size]] <- main_effects_train_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_train))
        main_effects_test_list[[list_size]] <- main_effects_train_list_norm[[list_size]] <- matrix(0,nrow = n_mcmc,ncol = nrow(x_test))
      }

      names(main_effects_train_list_norm) <- names(main_effects_test_list_norm) <-  names(main_effects_train_list) <- names(main_effects_test_list) <- dummy_x$continuousVars
    }
  } else {
    main_effects_train_list_norm <- main_effects_test_list_norm <- main_effects_train_list <- main_effects_test_list <- NULL
  }

  # Gonna create a list of lists to store all the indexes for all split rules and cutpoints
  all_var_splits <- vector("list",ncol(x_train_scale))
  names(all_var_splits) <- colnames(x_train_scale)



  # Iterating over all possible x.columns
  for(i in 1:length(all_var_splits)){

    # Creating the dummy for a list of index to store all numeric split values
    all_cut_points <- vector("list", nrow(xcut_m))


    for(j in 1:length(all_cut_points)){

      # Getting the node indexes object
      left_train_list <- vector("list",length = 1L)
      names(left_train_list) <- "left_train"
      right_train_list <- vector("list",length = 1L)
      names(right_train_list) <- "right_train"
      left_test_list <- vector("list",length = 1L)
      names(left_test_list) <- "left_test"
      right_test_list <- vector("list",length = 1L)
      names(right_test_list) <- "right_test"

      node_index <- append(left_train_list, right_train_list) |>
        append(left_test_list) |> append(right_test_list)

      all_cut_points[[j]]$left_train <-  which(x_train_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_train <-  which(x_train_scale[,i] >= xcut_m[j,i])
      all_cut_points[[j]]$left_test <-  which(x_test_scale[,i] < xcut_m[j,i])
      all_cut_points[[j]]$right_test <-  which(x_test_scale[,i] >= xcut_m[j,i])

    }

    all_var_splits[[i]] <- all_cut_points

  }

  # Creating the penalty matrix
  if(dif_order==0){
    P_train_main <- diag(nrow = NCOL(B_train_obj[[1]]))
  } else {
    P_train_main <- P_gen(D_train_ = B_train_obj[[1]],dif_order_ = dif_order,eta = 1)
  }


  # Creating the penalty matrix for the interaction
  if(interaction_term){
    if(dif_order == 0){
      P_interaction <- diag(nrow = NCOL(B_train_obj[[1]])^2)
    } else {
      P_interaction <-  kronecker(P_train_main,P_train_main)
    }

  } else {
    P_interaction <- NULL
  }


  # Here define a new vector of priors for \tau_beta_j
  # Can we fit the bigger model
  if(mle_prior){

      all_p <- length(basis_subindex)

      nll2 <- function(pars) {
        sigma <- pars[1]
        sigma_b <- pars[2:(all_p+1)]
        Big_var <- matrix(0, nrow(x_train), nrow(x_train))

        for(pp in 1:(all_p)) {

          if(pp <= NCOL(x_train)){
            Big_var <- Big_var +  (sigma_b[pp]^(2))*B_train_obj[[pp]]%*%t(B_train_obj[[pp]])
          } else {
            Big_var <- Big_var +  (sigma_b[pp]^(2))*B_train_obj[[pp]]%*%t(B_train_obj[[pp]])
          }

        }
        return(-mvnfast::dmvn(y_train, rep(0, nrow(x_train)),
                              diag(sigma^2, nrow(x_train)) + Big_var,
                              log = TRUE))
      }

      ans2 <- nlminb(rep(1, (length(basis_subindex)+1)), nll2, lower = rep(1e-5, (length(basis_subindex)+1)))

      # Defining tau beta and removing the first element
      tau_beta <-  (ans2$par^(-2))[-1]

      # Rescaling based on the scaled
      tau_beta <- tau_beta*n_tree*((max_y-min_y)^2)
      d_tau_beta_j <- numeric(all_p)

      for(pp in 1:all_p){
          sigquant_beta <- 0.9
          nsigma_beta <- (tau_beta[pp]*n_tree)^(-1/2)
          # nsigma_beta <- 1

          # Calculating lambda
          qchi_beta <- stats::qchisq(p = 1-sigquant_beta,df = df,lower.tail = 1,ncp = 0)
          lambda_beta <- (nsigma_beta*nsigma_beta*qchi_beta)/df
          d_tau_beta_j[pp] <- (lambda_beta*df)/2
      }

  }

  # Tau_beta_priors


  #most of the functions
  data <- list(x_train = x_train_scale,
               x_test = x_test_scale,
               y_train = y_scale,
               xcut_m = xcut_m,
               B_train = B_train_obj,
               B_test = B_test_obj,
               dif_order = dif_order,
               alpha = alpha,
               beta = beta,
               all_var_splits = all_var_splits,
               n_tree = n_tree,
               tau_mu = tau_mu,
               tau = nsigma^(-2),
               a_tau = a_tau,
               d_tau = d_tau,
               tau_beta = tau_beta,
               basis_subindex = basis_subindex,
               P = P_train_main,
               P_interaction  = P_interaction,
               node_min_size = node_min_size,
               all_var = all_var,
               stump = stump,
               a_tau_beta_j = a_tau_beta_j,
               d_tau_beta_j = d_tau_beta_j,
               interaction_term = interaction_term,
               interaction_list = interaction_list,
               dummy_x = dummy_x,
               linero_sampler = linero_sampler)

  #   So to simply interepret the element all_var_splits each element correspond
  #to each variable. Afterwards each element corresponds to a cutpoint; Finally,
  #inside that level we would have the index for the the left and right nodes;



  # Create a list of arrays to store the contribution of each tree to the main effects
  if(main_effects_pred){
      tree_main_effects <- vector("list",n_mcmc)
      for( iii_ in 1:n_mcmc){
        tree_main_effects[[iii_]] <- array(0,dim = c(250,length(data$basis_subindex),data$n_tree))
      }
  }
  # Initialing for storing post samples
  post <- 0

  # Initialising all the stumps
  for(k in 1:data$n_tree){
      forest[[k]] <- stump(data = data)
  }

  # and tree predictions
  trees_fit <- matrix(0,nrow = n_tree,ncol = nrow(x_train_scale))
  trees_fit_test <- matrix(0,nrow = n_tree, ncol  = nrow(x_test_scale))

  # For cases where the tree is greater than one;
  if(scale_init){
    if(n_tree>1){
      # Initial prediction
      for(i in 1:n_tree){
        trees_fit[i,] <- y_scale/n_tree
      }
    }
  }


  # Setting a matrix to store the frequency that a variable appear within terminal nodes
  if(varimportance_bool){
    variable_importance_matrix <- matrix(0,nrow =n_mcmc, ncol = length(data$basis_subindex))
  }
  # Initialsing the loop
  for(i in 1:n_mcmc){

    # Initialising the partial train tree fits
    if(store_tree_fit){
      partial_train_fits[[i]] <- vector("list",data$n_tree)
      names(partial_train_fits[[i]]) <- paste0("tree",1:data$n_tree)
    }

    # Initialising orogress bar
    progress <- i / n_mcmc * 100

    x1_pred <- numeric(nrow(x_train))

    for(t in 1:data$n_tree){

      # Calculating the partial residuals
      if(n_tree>1){
        partial_residuals <- y_scale-colSums(trees_fit[-t,,drop = FALSE])
      } else {
        partial_residuals <- y_scale
      }


      # Sample a verb
      verb <- sample(c("grow","prune", "change"), prob = c(0.3,0.3,0.4),size = 1)
      # verb <- sample(c("grow","prune"),size = 1)

      # Verifuy this line here, if does really make change
      # if(i<500){
      #   verb <- "change"
      # }

      # Checking the trees variables

      # Sampling a verb
      if(verb == "grow"){
          forest[[t]] <- grow(tree = forest[[t]],
                              curr_part_res = partial_residuals,
                              data = data)


      } else if (verb == "prune"){
          forest[[t]] <- prune(tree = forest[[t]],
                               curr_part_res = partial_residuals,
                               data = data)

      } else if (verb == "change"){
          forest[[t]] <- change(tree = forest[[t]],
                                curr_part_res = partial_residuals,
                                data = data)

      }

      # Updating the betas
      update_betas_aux <- updateBetas(tree = forest[[t]],
                                      curr_part_res = partial_residuals,
                                      data = data,trees_fit = trees_fit,tree_number =  t)

      forest[[t]] <- update_betas_aux$tree

      # Getting the predictions
      # tree_predictions <- getPredictions(tree = forest[[t]],
      #                                    data = data)
      if(main_effects_pred){
          tree_main_effects[[i]][,,t] <- update_betas_aux$y_hat_train
      }
      trees_fit[t,] <- rowSums(update_betas_aux$y_hat_train)
      trees_fit_test[t,] <- rowSums(update_betas_aux$y_hat_test)
      if(store_tree_fit){
        partial_train_fits[[t]] <- update_betas_aux$y_hat_train
      }

      # Adding up the contribution for each tree with respect to the covariate (i)
      if(main_effects_pred){
        for(ii in 1:length(main_effects_train_list)){
          main_effects_train_list[[ii]][i,] <- main_effects_train_list[[ii]][i,] + update_betas_aux$y_hat_train[,ii]
          main_effects_test_list[[ii]][i,] <- main_effects_test_list[[ii]][i,] + update_betas_aux$y_hat_test[,ii]
        }

      }

    }


    # Counting the number of times that a vriable was selected
    if(varimportance_bool){
        variable_importance_matrix[i, ] <- varimportance(forest = forest,data = data)
    }

    # Getting final predcition
    y_hat <- colSums(trees_fit)
    y_hat_test <- colSums(trees_fit_test)


    # Seeing the results for the unidimensional cases.
    if(NCOL(data$x_train)==1){
      if(plot_preview){
        plot(x_train_scale,y_scale)
        for(plot_i in 1:n_tree){
          points(x_train_scale,trees_fit[plot_i,],pch=20,col = ggplot2::alpha(plot_i,0.2))
        }
        points(x_train_scale,y_hat,col = "blue",pch=20)
      }
    }

    # Updating all other parameters
    if(update_tau_beta){
      data$tau_beta <- update_tau_betas_j(forest = forest,data = data)
    }


    # Updating delta
    # data$delta <- update_delta(data = data)


    # Getting tau
    data$tau <- update_tau(y_train_hat = y_hat,
                           data = data)


    # Storing all predictions and (TREES depending on the boolean)
    if(store_tree_fit){
      all_trees[[i]] <- forest
      all_trees_fit[[i]] <- partial_train_fits
    }

    all_tau[[i]] <- data$tau
    all_y_hat[i,] <- y_hat
    all_y_hat_test[i,] <- y_hat_test
    all_tau_beta[i,] <- data$tau_beta
    # all_delta[i] <- data$delta


    # Print progress bar
    cat("\rProgress: [", paste(rep("=", floor(progress / 5)), collapse = ""),
        paste(rep(" ", floor((100 - progress) / 5)), collapse = ""),
        "] ", sprintf("%.1f%%", progress))

    # Flush the output
    flush.console()

    # Simulate some work
    Sys.sleep(0.1)

  }



  # Normalising elements
  all_tau_norm <- numeric(n_mcmc)

  if(interaction_term){
    all_tau_beta_norm <- matrix(NA,nrow = n_mcmc,ncol = NCOL(x_train_scale)+NCOL(interaction_list))
  } else {
    all_tau_beta_norm <- matrix(NA,nrow = n_mcmc,ncol = NCOL(x_train_scale))
  }

  if(store_tree_fit){
    all_trees_fit_norm <- vector("list",n_mcmc)
  }

  all_y_hat_norm <- matrix(NA,nrow = nrow(all_y_hat),ncol = ncol(all_y_hat))
  all_y_hat_test_norm <- matrix(NA,nrow = nrow(all_y_hat_test),ncol = ncol(all_y_hat_test))

  # Returning to the original scale
  if(scale_bool){

    all_tau_norm <- all_tau/((max_y-min_y)^2)
    all_tau_beta_norm <- all_tau_beta/((max_y-min_y)^2)

    for(post_iter in 1:n_mcmc){

      if(store_tree_fit){
        for(tree_number in 1:n_tree){
          all_trees_fit_norm[[post_iter]][[tree_number]] <- unnormalize_bart(z = all_trees_fit[[post_iter]][[tree_number]],a = min_y,b = max_y)
        }
      }
      all_y_hat_norm[post_iter,] <- unnormalize_bart(z = all_y_hat[post_iter,],a = min_y,b = max_y)
      all_y_hat_test_norm[post_iter, ] <- unnormalize_bart(z = all_y_hat_test[post_iter,],a = min_y,b = max_y)

      if(main_effects_pred){
        for(ii in 1:length(main_effects_train_list)){
          main_effects_train_list_norm[[ii]][post_iter,] <- unnormalize_bart_me(z = main_effects_train_list[[ii]][post_iter,],a = min_y,b = max_y)
          main_effects_test_list_norm[[ii]][post_iter,] <- unnormalize_bart_me(z = main_effects_test_list[[ii]][post_iter,],a = min_y,b = max_y)
        }
      }
    }
  } else {
    all_tau_norm <- all_tau
    all_tau_beta_norm <- all_tau_beta
    all_trees_fit_norm <-all_trees_fit

    all_y_hat_norm <- all_y_hat
    all_y_hat_test_norm <- all_y_hat_test
  }

  # Visualizing some plots
  if(plot_preview){
    if(interaction_term){

      # Main effect range difference
      range_basis_j_predictions <- numeric(length = length(main_effects_train_list_norm))
      range_tree_bais_j_predictions <- matrix(NA, ncol = n_tree, nrow = length(main_effects_train_list_norm))

      par(mfrow = c(2,floor(NCOL(data$x_train)/2)))
      for(jj in 1:length(main_effects_test_list_norm)){

        if(jj <= NCOL(data$x_train)){
            plot(x_train[,jj],colMeans(main_effects_train_list_norm[[jj]][1:i,, drop = FALSE]),main = paste0('X',jj),
                 ylab = paste0('G(X',jj,')'),pch=20,ylim = c(-4,4),xlab = paste0('x.',jj), col = alpha("black",1.0))
        } else if(jj == 5){
            plot(x_train[,1], x_train[,2], cex = colMeans(main_effects_train_list_norm[[jj]][1:i,,drop = FALSE]))
        }


        range_basis_j_predictions[jj] <- abs(diff(range(colMeans(main_effects_train_list_norm[[jj]][1:i,, drop = FALSE]))))

        for(tree_number in 1:data$n_tree){
            # Create an auxiliar matrix for the main effect
            aux_main_effect_matrix <- matrix(0, nrow = n_mcmc,ncol = nrow(data$x_train))


            if(main_effects_pred){
                for(mcmc_aux in 1:n_mcmc){
                  aux_main_effect_matrix[mcmc_aux,] <- unnormalize_bart_me(tree_main_effects[[mcmc_aux]][,jj,tree_number],a = min_y,b = max_y)
                }
            }

            if(jj <= NCOL(data$x_train)){
                  points(x_train[,jj],colMeans(aux_main_effect_matrix),main = paste0('X',jj),
                       col = ggplot2::alpha(tree_number,0.1), pch = 20)
            }
            range_tree_bais_j_predictions[tree_number, jj] <- diff(range(colMeans(aux_main_effect_matrix)))

        }

        # points(x_train[,jj],colMeans(main_effects_train_list[[jj]][1:i,, drop = FALSE]),main = paste0('X',jj),
        #      ylab = paste0('G(X',jj,')'),ylim = c(-0.5,0.5),pch=20,xlab = paste0('x.',jj))
        # if( jj ==3 ){
        #   points(x_train[,jj],20*(x_train[,jj]-0.5)^2,pch = 20, col = "blue")
        # } else if ( jj == 4){
        #   points(x_train[,jj],10*x_train[,jj],pch = 20, col = "blue")
        # } else if ( jj == 5){
        #   points(x_train[,jj],5*x_train[,jj],pch = 20, col = "blue")
        # }

      }
      par(mfrow = c(1,1))
    }

  }

  # Quick check over the RMSE
  if(plot_preview){
    rmse(all_y_hat_test_norm[(n_burn+1):n_mcmc,] %>% colMeans(), y = sim_test$y)
    rmse(sim_test$y, aux$yhat.test.mean)
  }

  # currr_ <- numeric(n_mcmc)
  # get_which <- 0
  # for(test_mcmc_ in 1:n_mcmc){
  #     for(test_ in 1:data$n_tree){
  #       if(any(tree_main_effects[[test_mcmc_]][,11,test_]>0)){
  #         currr_[test_mcmc_] <- currr_[test_mcmc_] + 1
  #         get_which <- test_
  #       }
  #     }
  # }
  # plot(currr_)
  # # Some extra analysis
  # par(mfrow=c(2,2))
  # plot(all_tau_norm, type = 'l', ylab = expression(tau), xlab = "MCMC_iter",main = expression(tau))
  # abline(v = 1318, lty = 'dashed', col = 'blue')
  #
  # plot(all_tau_beta[,11], type = 'l', ylab = expression(tau[11]), xlab = "MCMC_iter",main = expression(tau[11]))
  # abline(v = 1318, lty = 'dashed', col = 'blue')
  #
  # plot(all_tau_beta[,1], type = 'l', ylab = expression(tau[1]), xlab = "MCMC_iter",main = expression(tau[1]))
  # abline(v = 1318, lty = 'dashed', col = 'blue')
  #
  # plot(all_tau_beta[,2], type = 'l', ylab = expression(tau[2]), xlab = "MCMC_iter",main = expression(tau[2]))
  # abline(v = 1318, lty = 'dashed', col = 'blue')


  # plot(colMeans(all_y_hat_norm),y_train)
  # sqrt(crossprod((colMeans(all_y_hat_norm)-y_train))/n_)

  # plot(10*sin(pi*x_train_scale[,1]*x_train_scale[,2]),main_effects_train_list[[11]][1:i,] %>% colMeans(), ylab = "y")

  # ====== Few analyses from the results ======
  #           (Uncomment to run those)
  # ===========================================
  #

  # plot(all_tau,type = "l")
  # plot(all_tau,type = "l")
  #
  # y1_hat <- matrix(0,nrow = n_post,ncol = nrow(data$x_train))
  #
  # for(i in 1:nrow(y1_hat)){
  #   for(t in 1:n_tree){
  #     y1_hat[i,] <- all_trees_fit[[i]][[t]][,1]
  #   }
  # }
  #
  # plot(x_train_scale[,1],colMeans(y1_hat[1801:3000,,drop = FALSE]))
  # plot(colMeans(y1_hat[1801:3000,,drop = FALSE]),y_train)
  # # plot(x_train[,1],10 * sin(pi * x_train[, 1])) # For x1
  # # plot(x_train[,2],20 * (x_train[, 2] - 0.5)^2) # For x1
  #
  # plot(all_tau_beta, type = "l")
  # plot(all_delta, type = "l")
  # plot(x_train_scale,y_scale)
  # points(x_train_scale,colMeans(all_y_hat[1501:2000,]),pch= 20, col = "blue")
  # # ============================================
  #
  # curr <- 0
  # tree_lengths <- numeric()
  #
  # for(i in 1:length(all_trees)){
  #
  #     curr <- curr + 1
  #
  #     for(j in 1:n_tree){
  #         tree_lengths[curr] <- length(all_trees[[i]][[j]])
  #     }
  #
  # }
  # rmse(x = sim_test$y,y = all_y_hat_test_norm %>% colMeans())

  # Return the list with all objects and parameters
  return(list(y_train_hat = all_y_hat_norm,
              y_test_hat = all_y_hat_test_norm,
              all_tau = all_tau_norm,
              all_tau_beta = all_tau_beta_norm,
              prior = list(n_tree = n_tree,
                           alpha = alpha,
                           beta = beta,
                           tau_mu = tau_mu,
                           a_tau = a_tau,
                           d_tau = d_tau),
              mcmc = list(n_mcmc = n_mcmc,
                          n_burn = n_burn,
                          all_trees = all_trees,
                          main_effects_train = main_effects_train_list_norm,
                          main_effects_test = main_effects_test_list_norm),
              data = list(x_train = x_train,
                          y_train = y_train,
                          B_train = B_train_obj,
                          x_test = x_test,
                          B_test = B_test_obj,
                          basis_subindex = basis_subindex)))

}





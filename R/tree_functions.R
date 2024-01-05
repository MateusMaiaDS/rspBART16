# A function to get which variables are used within each terminal node
varimportance <- function(forest, data){

  # Getting the auxilinar vector
  var_counter <- var_counter_aux <- numeric(length = length(data$basis_subindex))

  for(ii_ in 1:data$n_tree){
    tree_ <- forest[[ii_]]
    terminal_names <- get_terminals(tree = tree_)
    n_terminal <- length(terminal_names)
    var_counter_aux <- numeric(length = length(data$basis_subindex))

    for(jj_ in 1:length(terminal_names)){
      var_counter_aux[tree_[[terminal_names[jj_]]]$pred_vars] <- var_counter_aux[tree_[[terminal_names[jj_]]]$pred_vars] + 1
    }
    var_counter <- var_counter + var_counter_aux/n_terminal
  }

  return(var_counter/data$n_tree)
}
# Creating a stump for a tree
stump <- function(data){

  # Creating the base node
  node <- list()

    node[["node0"]] <- list(
      # Creating the node number
      node_number = 0,
      pred_vars = 1:length(data$basis_subindex),
      inter = NA,
      isRoot = TRUE,
      # Creating a vector with the tranining index
      train_index = 1:nrow(data$x_train),
      test_index = 1:nrow(data$x_test),
      depth_node = 0,
      node_var = NA,
      node_cutpoint_index = NA,
      left = NA,
      right = NA,
      parent_node = NA,
      ancestors = NA,
      terminal = TRUE,
      betas_vec = rep(0, length(unlist(data$basis_subindex)))
  )



  # Returning the node
  return(node)

}

# Get all the terminal nodes
get_terminals <- function(tree){

  # Return the name of the termianl nodes
  return(names(tree)[unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)])
}

# Get nog terminal nodes
get_nogs <- function(tree){

  # Return the name of the termianl nodes
  non_terminal <- names(tree)[!unlist(lapply(tree, function(x){x$terminal}),use.names =  TRUE)]

  # In case there are non nonterminal nondes
  if(length(non_terminal)==0){
    return(non_terminal)
  }

  bool_nog <- vector("logical",length = length(non_terminal))
  for(i in 1:length(bool_nog)){
    # Checking if both children are terminal
    if( tree[[tree[[non_terminal[i]]]$left]]$terminal & tree[[tree[[non_terminal[i]]]$right]]$terminal) {
      bool_nog[i] <- TRUE
    }
  }

  return(  non_terminal[bool_nog])
}

# Getting the maximum node index number
get_max_node <- function(tree){

  # Return the name of the termianl nodes
  return(max(unlist(lapply(tree, function(x){x$node_number}),use.names =  TRUE)))
}




# A function to calculate the loglikelihood
nodeLogLike <- function(curr_part_res,
                        j_,
                        index_node,
                        data){

  # Subsetting the residuals
  curr_part_res_leaf <- curr_part_res[index_node]

  # Getting the number of observationsin the terminal node
  n_leaf <- length(index_node)
  d_basis <- length(j_)
  ones <- matrix(1,nrow = n_leaf)


  if(length(j_)==0){
    stop(" Node Log-likelihood: No variables")
  }

  # Using the Andrew's approach I would have
  mean_aux <- rep(0,length(curr_part_res_leaf))
  cov_aux <- matrix(0,nrow = length(index_node),ncol = length(index_node))
  # diag_tau_beta_inv <- diag(x = 1/unique(data$tau_beta), nrow = )

  for(jj in 1:length(j_)){
      # Adding the quantities with respect to the interaction
      if(j_[jj] <= length(data$dummy_x$continuousVars)){
        cov_aux <- cov_aux + data$B_train[[j_[jj]]][index_node,,drop = FALSE]%*%solve((data$tau_beta[j_[jj]]*data$P),t(data$B_train[[j_[jj]]][index_node,,drop = FALSE]))
      } else {
        cov_aux <- cov_aux + data$B_train[[j_[jj]]][index_node,,drop = FALSE]%*%solve((data$tau_beta[j_[jj]]*data$P_interaction),t(data$B_train[[j_[jj]]][index_node,,drop = FALSE]))
      }
  }
  cov_aux <- diag(x = (data$tau^(-1)),nrow = n_leaf) + cov_aux

  result <- mvnfast::dmvn(X = curr_part_res_leaf,mu = mean_aux,
                          sigma = cov_aux ,log = TRUE)


  return(c(result))

}



# Grow a tree
grow <- function(tree,
                 curr_part_res,
                 data){

  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  g_node_name <- sample(terminal_nodes,size = 1)
  g_node <- tree[[g_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0

  # acceptance_grid <- numeric(100)
  # for(kk in 1:100){
  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2

    # Sample a split var
    # ===== Uncomment this line below after ========
    p_var <- sample(1:NCOL(data$x_train),size = 1)
    # ==============================================
    # p_var <- 7

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[g_node$train_index,p_var])

    # Case of invalid range
    if(length(valid_range_grow)==0){
      return(tree)
    }

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)
    # sample_cutpoint <- valid_cutpoint[kk]

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% g_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% g_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% g_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% g_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(g_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(g_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }


    # === Uncomment those lines after
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }

  # Calculating loglikelihood for the grown node, the left and the right node
  # Recover the g_node index

  node_index_var <- g_node$pred_vars

  g_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           j_ = node_index_var,
                           index_node = g_node$train_index,
                           data = data)


  left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                               j_ = node_index_var,
                               index_node = left_index,
                               data = data)

  right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                j_ = node_index_var,
                                index_node = right_index,
                                data = data)

  # Calculating the prior
  prior_loglike <- log(data$alpha*(1+g_node$depth_node)^(-data$beta)) + # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+g_node$depth_node+1)^(-data$beta)) - # plus the prior of the two following nodes being terminal
    log(1-data$alpha*(1+g_node$depth_node)^(-data$beta)) # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_nog_nodes+1))-log(0.3/n_t_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(-g_loglike+left_loglike+right_loglike+prior_loglike+log_trasition_prob)
  # acceptance_grid[kk] <- acceptance



  # par(mfrow=c(1,2))
  # plot(data$xcut_m[,2],(acceptance_grid), main = "Acceptance to split on X2", xlab = "X2", ylab = "Prob. Acceptance")

  if(data$stump) {
    acceptance <- acceptance*(-1)
  }

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    if(any(is.na(g_node$ancestors))){
      new_ancestors <- p_var
    } else {
      new_ancestors <- c(g_node$ancestors,p_var)
    }

    left_node <- list(node_number = max_index+1,
                      j = g_node$j,
                      pred_vars = g_node$pred_vars,
                      inter = g_node$inter,
                      isRoot = FALSE,
                      train_index = left_index,
                      test_index = left_test_index,
                      depth_node = g_node$depth_node+1,
                      node_var = p_var,
                      node_cutpoint_index = sample_cutpoint,
                      left = NA,
                      right = NA,
                      parent_node = g_node_name,
                      ancestors = new_ancestors,
                      terminal = TRUE,
                      betas_vec = g_node$betas_vec)

    right_node <- list(node_number = max_index+2,
                       j = g_node$j,
                       pred_vars = g_node$pred_vars,
                       inter = g_node$inter,
                       isRoot = FALSE,
                       train_index = right_index,
                       test_index = right_test_index,
                       depth_node = g_node$depth_node+1,
                       node_var = p_var,
                       node_cutpoint_index = sample_cutpoint,
                       left = NA,
                       right = NA,
                       parent_node = g_node_name,
                       ancestors = new_ancestors,
                       terminal = TRUE,
                       betas_vec = g_node$betas_vec)

    # Modifying the current node
    tree[[g_node_name]]$left = paste0("node",max_index+1)
    tree[[g_node_name]]$right = paste0("node",max_index+2)
    tree[[g_node_name]]$terminal = FALSE

    tree[[paste0("node",max_index+1)]] <- left_node
    tree[[paste0("node",max_index+2)]] <- right_node


  } else {

    # Do nothing

  }

  # Return the new tree
  return(tree)
}

# Pruning a tree
prune <- function(tree,
                  curr_part_res,
                  data){


  # Getting the maximum index number
  max_index <- get_max_node(tree)

  # Sampling a terminal node
  terminal_nodes <- get_terminals(tree)
  n_t_nodes <- length(terminal_nodes)
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)

  # Just in case to avoid errors
  if(n_nog_nodes==0){
    return(tree)
  }

  # Selecting a node to be pruned
  p_node_name <- sample(nog_nodes,size = 1)
  p_node <- tree[[p_node_name]]

  # Getting the indexes from the left and right children from the pruned node
  children_left_index <- tree[[p_node$left]]$train_index
  children_right_index <- tree[[p_node$right]]$train_index
  children_left_ancestors <- tree[[p_node$left]]$ancestors
  children_right_ancestors <- tree[[p_node$right]]$ancestors

  # Calculating loglikelihood for the grown node, the left and the right node

  if(!any(is.na(p_node$inter))){
    # node_index_var <- c(p_node$j,which( names(data$basis_subindex) %in% paste0(p_node$j,sort(p_node$inter))))
    node_index_var <- p_node$pred_vars
  } else {
    # node_index_var <- p_node$j
    node_index_var <- p_node$pred_vars
  }

  p_loglike <- nodeLogLike(curr_part_res = curr_part_res,
                           index_node = p_node$train_index,
                           j_ = node_index_var,
                           data = data)


  p_left_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                 index_node =  children_left_index,
                                 j_ = node_index_var,
                                 data = data)

  p_right_loglike <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = children_right_index,
                                  j_ = node_index_var,
                                  data = data)

  # Calculating the prior
  prior_loglike <- log(1-data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the new terminal node
    log(data$alpha*(1+p_node$depth_node)^(-data$beta)) - # Prior of the grown node becoming nonterminal
    2*log(1-data$alpha*(1+p_node$depth_node+1)^(-data$beta))  # plus the prior of the two following nodes being terminal
  # minus the probability of the grown node being terminal

  # Transition prob
  log_trasition_prob  = log(0.3/(n_t_nodes))-log(0.3/n_nog_nodes)

  # Calculating the acceptance probability
  acceptance <- exp(p_loglike-p_left_loglike-p_right_loglike+prior_loglike+log_trasition_prob)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1)<acceptance){

    # Erasing the terminal nodes
    tree[[p_node$left]] <- NULL
    tree[[p_node$right]] <- NULL

    # Modifying back the pruned node
    tree[[p_node_name]]$left <- NA
    tree[[p_node_name]]$right <- NA
    tree[[p_node_name]]$terminal <- TRUE

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}


# Change a tree
change <- function(tree,
                   curr_part_res,
                   data){

  # Changing the stump
  if(length(tree)==1 & (!data$all_var)){
    change_stump_obj <- change_stump(tree = tree,
                                     curr_part_res = curr_part_res,
                                     data = data)
    return(change_stump_obj)
  }


  # For the seocnd case
  if(length(tree)==1){
    return(tree)
  }

  # Sampling a terminal node
  nog_nodes <- get_nogs(tree)
  n_nog_nodes <- length(nog_nodes)
  c_node_name <- sample(nog_nodes,size = 1)
  c_node <- tree[[c_node_name]]


  valid_terminal_node <- TRUE
  valid_count <- 0


  while(valid_terminal_node){
    # Convinience while to avoid terminal nodes of 2
    # Sample a split var
    p_var <- sample(1:ncol(data$x_train),size = 1)

    # Selecting an available cutpoint from this terminal node
    valid_range_grow <- range(data$x_train[c_node$train_index,p_var])

    # Subsetting the indexes of
    valid_cutpoint <- which(data$xcut_m[,p_var]>valid_range_grow[1] & data$xcut_m[,p_var]<valid_range_grow[2])

    # When there's no valid cutpoint on the sampled terminal node
    if(length(valid_cutpoint)==0){
      return(tree)
    }

    # Getting which cutpoints are valid and sample onde index
    sample_cutpoint <- sample(valid_cutpoint,
                              size = 1)

    # Getting the left & right index
    left_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_train %in% c_node$train_index]
    right_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_train %in% c_node$train_index]

    left_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$left_test %in% c_node$test_index]
    right_test_index  <- data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test[data$all_var_splits[[p_var]][[sample_cutpoint]]$right_test %in% c_node$test_index]



    # Verifying that the correct number was used
    if((length(left_index)+length(right_index))!=length(c_node$train_index)){
      stop("Something went wrong here --- train grown index doest match")
    }

    if((length(left_test_index)+length(right_test_index))!=length(c_node$test_index)){
      stop("Something went wrong here --- test grown index doest match")
    }

    # Avoiding having terminal nodes with just one observation
    if( (length(left_index) > data$node_min_size) & (length(right_index)>data$node_min_size)){
      # Getting out of the while
      break
    } else {

      # Adding one to the counter
      valid_count = valid_count + 1

      # Stop trying to search for a valid cutpoint
      if(valid_count > 2) {
        valid_terminal_node = FALSE
        return(tree)
      }
    }
  }

  # For convinience we are going to avoid terminal nodes less than 2
  if( (length(left_index)<2) || (length(right_index) < 2)) {
    stop("Error of invalid terminal node")
  }


  # Getting the node_index var
  node_index_var <- c_node$pred_vars


  # Calculating loglikelihood for the new changed nodes and the old ones
  c_loglike_left <- nodeLogLike(curr_part_res = curr_part_res,
                                index_node = tree[[c_node$left]]$train_index,
                                j_ = node_index_var,
                                data = data)


  c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                  index_node = tree[[c_node$right]]$train_index,
                                  j_ =  node_index_var,
                                  data = data)

  # Calculating a new ancestors left and right
  old_p_var <- tree[[c_node$left]]$node_var

  # Storing new left and right ancestors
  new_left_ancestors <- tree[[c_node$left]]$ancestors
  new_left_ancestors[length(new_left_ancestors)] <- p_var

  new_right_ancestors <- tree[[c_node$right]]$ancestors
  new_right_ancestors[length(new_right_ancestors)] <- p_var


  new_c_loglike_left <-  nodeLogLike(curr_part_res = curr_part_res,
                                     index_node = left_index,
                                     j = node_index_var,
                                     data = data)

  new_c_loglike_right <-  nodeLogLike(curr_part_res = curr_part_res,
                                      index_node = right_index,
                                      j =  node_index_var,
                                      data = data)


  # Calculating the acceptance probability
  acceptance <- exp(new_c_loglike_left+new_c_loglike_right-c_loglike_left-c_loglike_right)

  # Getting the training the left and the right index for the the grown node
  if(stats::runif(n = 1,min = 0,max = 1)<acceptance){

    # Updating the left and the right node
    # === Left =====
    tree[[c_node$left]]$node_var <- p_var
    tree[[c_node$left]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$left]]$train_index <- left_index
    tree[[c_node$left]]$test_index <- left_test_index
    tree[[c_node$left]]$ancestors <- new_left_ancestors

    #==== Right ====
    tree[[c_node$right]]$node_var <- p_var
    tree[[c_node$right]]$node_cutpoint_index <- sample_cutpoint
    tree[[c_node$right]]$train_index <- right_index
    tree[[c_node$right]]$test_index <- right_test_index
    tree[[c_node$right]]$ancestors <- new_right_ancestors

  } else {
    # Do nothing
  }

  # Return the new tree
  return(tree)

}


# ============
# Update Betas
# ============
updateBetas <- function(tree,
                        curr_part_res,
                        data,
                        trees_fit,
                        tree_number){


  # Getting the terminals
  t_nodes_names <- get_terminals(tree)


  # Getting the current prediction for that tree
  y_hat_train <- matrix(0,nrow = nrow(data$x_train),ncol = length(data$basis_subindex))
  y_hat_test <- matrix(0,nrow = nrow(data$x_test),ncol = length(data$basis_subindex))

  for(i in 1:length(t_nodes_names)){


    # Select the current terminal node
    cu_t <- tree[[t_nodes_names[i]]]

    # The lines above are summarised here
    node_index_var <- cu_t$pred_vars

    # Selecting the actually parameters subsetting
    basis_dim <- NCOL(data$P)
    basis_dim_interaction <- NCOL(data$P_interaction)
    n_leaf <- length(cu_t$train_index)
    diag_leaf <- diag(nrow = n_leaf)
    diag_basis <- diag(nrow = basis_dim)


    #  Calculating the quantities need to the posterior of \beta
    # == Starting to iterate over those coefficients ==========#
    for(jj in 1:length(node_index_var)){


      leaf_basis_subindex <- unlist(data$basis_subindex[node_index_var[jj]]) # Recall to the unique() here too

      # RES_LEAF also need to updated here from the new_curr_part_res
      old_betas <- matrix(tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex],nrow = 1)

      res_leaf <- matrix(curr_part_res[cu_t$train_index], ncol=1) - (trees_fit[tree_number,cu_t$train_index] - tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],old_betas))

      # Getting the index for the vector of betas
      b_ <- crossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],res_leaf)

      data_tau_beta_diag <- rep(data$tau_beta[node_index_var], NCOL(data$B_train[[node_index_var[jj]]])) # Don't really use this
      if(node_index_var[jj]<=length(data$dummy_x$continuousVars)){
        U_ <- data$P*data$tau_beta[node_index_var[jj]]
      } else {
        U_ <- data$P_interaction*data$tau_beta[node_index_var[jj]]
      }

      Q_ <- (crossprod(data$B_train[[node_index_var[jj]]]) + data$tau^(-1)*U_)
      Q_inv_ <- chol2inv(chol(Q_))

      # Storing the old betas
      # See that I also creating a vector with the new betas
      new_betas <- mvnfast::rmvn(n = 1,mu = Q_inv_%*%b_,sigma = (data$tau^(-1))*Q_inv_)
      tree[[t_nodes_names[i]]]$betas_vec[leaf_basis_subindex] <- new_betas
      new_betas <- matrix(new_betas,nrow = 1)
      # Updating the residuals
      new_partial_pred <- tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],new_betas)
      # Need to update the trees fit!
      trees_fit[tree_number,cu_t$train_index] <- trees_fit[tree_number,cu_t$train_index] - tcrossprod(data$B_train[[node_index_var[jj]]][cu_t$train_index,,drop=FALSE],old_betas) + new_partial_pred

      y_hat_train[cu_t$train_index,node_index_var[jj]] <- new_partial_pred
      y_hat_test[cu_t$test_index,node_index_var[jj]] <- tcrossprod(data$B_test[[node_index_var[jj]]][cu_t$test_index,,drop=FALSE],new_betas)
    }

  }

  # Returning the tree
  return(list(tree = tree,
              y_hat_train = y_hat_train,
              y_hat_test = y_hat_test))

}


# =================
# Update \tau_betas
# =================
update_tau_betas_j <- function(forest,
                               data){

  # Setting some default hyperparameters
  a_tau_beta <- data$a_tau_beta_j
  d_tau_beta <- data$d_tau_beta_j

  tau_b_shape <- 0.0
  tau_b_rate <- 0.0


  if(data$interaction_term){
    tau_b_shape <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_b_rate <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
    tau_beta_vec_aux <- numeric(NCOL(data$x_train)+NCOL(data$interaction_list))
  } else{
    tau_b_shape <- numeric(NCOL(data$x_train))
    tau_b_rate <- numeric(NCOL(data$x_train))
    tau_beta_vec_aux_proposal <- tau_beta_vec_aux <- numeric(NCOL(data$x_train))
  }

  # Iterating over all trees
  for(i in 1:length(forest)){

    # Getting terminal nodes
    t_nodes_names <- get_terminals(forest[[i]])
    n_t_nodes <- length(t_nodes_names)

    # Iterating over the terminal nodes
    for(j in 1:length(t_nodes_names)){

      cu_t <- forest[[i]][[t_nodes_names[j]]]


      # All the information from var_ now is summarised inside the element from ht enode pred_vars
      var_ <- cu_t$pred_vars


      # Getting ht leaf basis
      for(kk in 1:length(var_)){
        leaf_basis_subindex <- unlist(data$basis_subindex[var_[kk]]) # Recall to the unique() function here
        p_ <- length(leaf_basis_subindex)
        betas_mat_ <- matrix(cu_t$betas_vec[leaf_basis_subindex],nrow = p_)
        tau_b_shape[var_[kk]] <- tau_b_shape[var_[kk]] + p_

        if(var_[kk] <= NCOL(data$x_train)){
          tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P,betas_mat_)))
        } else {
          tau_b_rate[var_[kk]] <- tau_b_rate[var_[kk]] + c(crossprod(betas_mat_,crossprod(data$P_interaction,betas_mat_)))

        }
      }
      # }

    }


  }

  if(data$interaction_term){

    if(data$linero_sampler){

      # Use Linero sampler to solve this matter
      for(j in 1:(NCOL(data$x_train)+NCOL(data$interaction_list)) ){

        if(length(d_tau_beta)>1){
          tau_beta_vec_aux_proposal <- rgamma(n = 1,
                                        shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                        rate = 0.5*tau_b_rate[j] + d_tau_beta[j])

          # Getting the values from the proposal
          acceptance_tau_beta_ <- exp(extraDistr::dhcauchy(x = tau_beta_vec_aux_proposal^(-1/2),sigma = data$tau_mu^(-1/2),log = TRUE)+(-3/2)*log(tau_beta_vec_aux_proposal) - extraDistr::dhcauchy(x = data$tau_beta[j]^(-1/2),sigma = data$tau_mu^(-1/2),log = TRUE)-(-3/2)*log(data$tau_beta[j]))

          if(stats::runif(n = 1) < acceptance_tau_beta_){
            tau_beta_vec_aux[j] <- tau_beta_vec_aux_proposal
          } else {
            tau_beta_vec_aux[j] <- data$tau_beta[j]
          }

        } else {
          tau_beta_vec_aux_proposal <- rgamma(n = 1,
                                        shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                        rate = 0.5*tau_b_rate[j] + d_tau_beta)

          # Getting the values from the proposal
          acceptance_tau_beta_ <- exp(extraDistr::dhcauchy(x = tau_beta_vec_aux_proposal^(-1/2),sigma = data$tau_mu^(-1/2),log = TRUE)+(-3/2)*log(tau_beta_vec_aux_proposal) - extraDistr::dhcauchy(x = data$tau_beta[j]^(-1/2),sigma = data$tau_mu^(-1/2),log = TRUE)-(-3/2)*log(data$tau_beta[j]))

          if(stats::runif(n = 1) < acceptance_tau_beta_){
            tau_beta_vec_aux[j] <- tau_beta_vec_aux_proposal
          } else {
            tau_beta_vec_aux[j] <- data$tau_beta[j]
          }

        }
      }


    } else { # This else is regarding "Linero sampler" (the chunk above wouldn't use it)

        for(j in 1:(NCOL(data$x_train)+NCOL(data$interaction_list)) ){

          if(length(d_tau_beta)>1){
              tau_beta_vec_aux[j] <- rgamma(n = 1,
                                            shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                            rate = 0.5*tau_b_rate[j] + d_tau_beta[j])
          } else {
            tau_beta_vec_aux[j] <- rgamma(n = 1,
                                          shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                          rate = 0.5*tau_b_rate[j] + d_tau_beta)
          }
        }

    }
  } else {

    # Adding the "Linero sampler option
    if(data$linero_sampler){

      for(j in 1:NCOL(data$x_train)){
        tau_beta_vec_aux_proposal[j] <- rgamma(n = 1,
                                               shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                               rate = 0.5*tau_b_rate[j] + d_tau_beta)

        # Getting the values from the proposal
        acceptance_tau_beta_ <- exp(extraDistr::dhcauchy(x = tau_beta_vec_aux_proposal^(-1/2),sigma = data$tau_mu^(-1/2),log = TRUE)+(-3/2)*log(tau_beta_vec_aux_proposal) - extraDistr::dhcauchy(x = data$tau_beta[j]^(-1/2),sigma = data$tau_mu^(-1/2), log = TRUE)-(-3/2)*log(data$tau_beta[j]))

        if(stats::runif(n = 1) < acceptance_tau_beta_){
          tau_beta_vec_aux[j] <- tau_beta_vec_aux_proposal
        } else {
          tau_beta_vec_aux[j] <- data$tau_beta[j]
        }
      }
    } else {
      tau_beta_vec_aux[j] <- rgamma(n = 1,
                                    shape = 0.5*tau_b_shape[j] + a_tau_beta,
                                    rate = 0.5*tau_b_rate[j] + d_tau_beta)
    }

  }

  return(tau_beta_vec_aux)

}


# ===================
# Updating the \delta
# ===================

# A function to get predictions
getPredictions <- function(tree,
                           data){

  # Creating the vector to hold the values of the prediction
  if(data$interaction_term){
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = NCOL(data$x_train)+NCOL(data$interaction_list))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = NCOL(data$x_test)+NCOL(data$interaction_list))
  } else {
    y_hat <- matrix(0, nrow = nrow(data$x_train), ncol = ncol(data$x_train))
    y_hat_test <- matrix(0,nrow(data$x_test), ncol = ncol(data$x_test))
  }

  # Getting terminal nodes
  t_nodes <- get_terminals(tree = tree)
  n_t_nodes <- length(t_nodes)

  for(i in 1:n_t_nodes){


    # Getting the current terminal node
    cu_t <- tree[[t_nodes[[i]]]]
    leaf_train_index <- cu_t$train_index
    leaf_test_index <- cu_t$test_index

    # Getting the variables used in the model
    node_index_var <- cu_t$pred_vars

    leaf_ancestors <- node_index_var # here isnt really the ancestors, but the variables that are being used

    leaf_basis_subindex <- data$basis_subindex[leaf_ancestors]

    # This test doesn't make sense anymore
    # # Test unit
    # if(length(leaf_ancestors)!=length(leaf_basis_subindex)){
    #   stop("Error on the getPredictions function")
    # }

    # Only add the marginal effects if the variables are within that terminal node
    if(length(leaf_basis_subindex)!=0){
      for(k in 1:length(leaf_basis_subindex)){

        y_hat[leaf_train_index,leaf_ancestors[k]] <- y_hat[leaf_train_index,leaf_ancestors[k]] + data$D_train[leaf_train_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]
        y_hat_test[leaf_test_index,leaf_ancestors[k]] <- y_hat_test[leaf_test_index,leaf_ancestors[k]] + data$D_test[leaf_test_index,leaf_basis_subindex[[k]], drop = FALSE]%*%tree[[t_nodes[i]]]$betas_vec[leaf_basis_subindex[[k]]]

      }
    }

  }

  # Returning both training and test set predictions
  return(list(y_train_hat = y_hat,
              y_hat_test = y_hat_test))

}

# Updating tau
update_tau <- function(y_train_hat,
                       data){

  # Sampling a tau value
  n_ <- nrow(data$x_train)
  tau_sample <- stats::rgamma(n = 1,shape = 0.5*n_+data$a_tau,rate = 0.5*crossprod((data$y_train-y_train_hat))+data$d_tau)

  return(tau_sample)

}



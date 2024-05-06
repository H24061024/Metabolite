rescale_function = function(x,f){
  xdim = dim(x[[1]])
  rescale_matrix = matrix(NA,ncol = length(x[[1]]),nrow = length(x))
  for (i in 1:length(x)) {
    rescale_matrix[i,] = as.matrix(x[[i]],nrow = 1)
  }
  rescale_matrix = f(rescale_matrix)
  for (i in 1:length(x)) {
    x[[i]] = array(rescale_matrix[i,],dim=xdim)
  }
  return(x)
}

rescale_01 = function(x){
  for (i in 1:length(x[1,])) {
    max_x = max(x[,i])
    min_x = min(x[,i])
    if(max_x==min_x){
      x[,i] = x[,i]-min_x
    }else{
      x[,i] = (x[,i]-min_x)/(max_x-min_x)
    }

  }

  return(x)
}

rescale_log2 = function(x){
  return(log2(x))
}




rescale_standard = function(x){
  for (i in 1:length(x[1,])) {
    mean_x = mean(x[,i])
    sd_x = sd(x[,i])
    if(sd_x==0){
      x[,i] = x[,i]-mean_x
    }else{
      x[,i] = (x[,i]-mean_x)/(sd_x)
    }

  }

  return(x)
}

t_adjten = function (x, z, y, testx = NULL, testz = NULL, is.centered = FALSE)
{
  n = length(x)
  dimen = dim(x[[1]])
  nvars = prod(dimen)
  vec_x = matrix(0, nvars, n)
  for (i in 1:n) {
    vec_x[, i] = matrix(x[[i]], nvars, 1)
  }
  vec_x = t(vec_x)
  if (!is.null(testx)) {
    tesize = length(testx)
    vec_testx = matrix(0, nvars, tesize)
    for (i in 1:tesize) {
      vec_testx[, i] = matrix(testx[[i]], nvars, 1)
    }
    vec_testx = t(vec_testx)
  }
  z = as.matrix(z)
  cz = z
  n = dim(z)[1]
  q = dim(z)[2]
  cvecx = vec_x
  nclass <- as.integer(length(unique(y)))
  if (is.centered == FALSE) {
    for (i in 1:nclass) {
      if (q > 1) {
        cz[y == i, ] = sweep(z[y == i, ], 2, colMeans(z[y ==
                                                          i, ]))
      }
      else {
        cz[y == i] = z[y == i] - mean(z[y == i])
      }
      cvecx[y == i, ] = sweep(vec_x[y == i, ], 2, colMeans(vec_x[y == i, ]))
    }
  }
  c = solve(t(cz) %*% cz) %*% t(cz)
  c = c %*% cvecx
  vec_xres = vec_x - z %*% c
  xres = array(list(), n)
  for (i in 1:n) {
    xres[[i]] = array(vec_xres[i, ], dim = dimen)
  }
  if (!is.null(testx)) {
    vec_testxres = vec_testx - testz %*% c
    testxres = array(list(), tesize)
    for (i in 1:tesize) {
      testxres[[i]] = array(vec_testxres[i, ], dim = dimen)
    }
  }
  else {
    testxres = NULL
  }
  muz = matrix(0, nrow = q, ncol = (nclass - 1))
  for (i in 2:nclass) {
    if (q > 1) {
      muz[, (i - 1)] = apply(z[y == i, ], 2, mean) - apply(z[y ==
                                                               1, ], 2, mean)
    }
    else {
      muz[i - 1] = mean(z[y == i]) - mean(z[y == 1])
    }
  }
  gamma = solve(cov(z)) %*% muz
  outlist = list(gamma = gamma, xres = xres, testxres = testxres,alpha = c,vec_xres=vec_xres)
  outlist
}
predict.catch.prob = function (object, newx, z = NULL, ztest = NULL, gamma = NULL)
{
  if (is.null(gamma)) {
    pred = predict.tsda(object, newx)
  }
  else {
    thetatm <- object$beta
    theta = array(list(), length(thetatm))
    mu <- object$mu
    prior <- object$prior
    nclass <- length(prior)
    dimen <- dim(newx[[1]])
    nvars <- prod(dimen)
    nlambda <- length(theta)
    gamma = as.matrix(gamma)
    q = dim(gamma)[1]
    z = as.matrix(z)
    ztest = as.matrix(ztest)
    for (i in 1:nlambda) {
      theta[[i]] = matrix(0, nrow = dim(thetatm[[i]])[1] +
                            q, ncol = nclass - 1)
      for (j in 1:nclass - 1) {
        theta[[i]][1:nvars, j] = matrix(thetatm[[i]][,
                                                     j], ncol = 1)
        for (qq in 1:q) {
          theta[[i]][nvars + qq, j] = gamma[qq, j]
        }
      }
    }
    mubar = matrix(list(), nclass - 1, 1)
    for (i in 1:(nclass - 1)) {
      mubar[[i]] = (mu[[i + 1]] + mu[[1]])/2
    }
    n <- length(newx)
    nn <- length(object$x)
    x.train <- object$x
    vecx.train = matrix(0, ncol = nn, nrow = nvars + q)
    vecnewx = matrix(0, ncol = n, nrow = nvars + q)
    for (i in 1:nn) {
      vecx.train[1:nvars, i] <- matrix(x.train[[i]], ncol = 1)
      for (qq in 1:q) {
        vecx.train[nvars + qq, i] = z[i, qq]
      }
    }
    vecx.train = t(vecx.train)
    for (i in 1:(length(newx))) {
      vecnewx[1:nvars, i] <- matrix(newx[[i]], ncol = 1)
      for (qq in 1:q) {
        vecnewx[nvars + qq, i] = ztest[i, qq]
      }
    }
    vecnewx = t(vecnewx)
    y.train <- object$y
    pred <- matrix(0, n, nlambda)
    pred[1] <- which.max(prior)
    for (i in 1:nlambda) {
      nz <- sum(theta[[i]][, 1] != 0)
      if (nz == 0) {
        pred[, i] <- which.max(prior)
      }
      else {
        xfit <- vecx.train %*% theta[[i]][, 1:(min(nclass -
                                                     1, nz)), drop = FALSE]
        xfit.sd <- matrix(0, nclass, ncol(xfit))
        for (j in 1:nclass) {
          xfit.sd[j, ] <- apply(xfit[y.train == j, ,
                                     drop = FALSE], 2, sd)
        }
        xfit.sd <- apply(xfit.sd, 2, min)
        if (min(xfit.sd) < 1e-04) {
          print("xfit.sd less then 1e-04,can't output probability")
          pred = NULL
        }
        else {
          l <- lda(xfit, y.train)
          pred <- predict(l, vecnewx %*% theta[[i]][,
                                                    1:(min(nclass - 1, nz))])$posterior
        }
      }
    }
  }
  pred
}

preprocessing = function(x,y,z=NULL,testx=NULL,testy=NULL,testz=NULL,
                         baseline.factor=NULL,rescale_x=NULL,rescale_z=NULL,
                         tensor_namelist=NULL,cov_matrix_namelist=NULL){
  x_train_as_test = is.null(testx)
  y_train_as_test = is.null(testy)
  z_train_as_test = is.null(testz)
  train_as_test=FALSE
  if(is.null(z)){
    if(x_train_as_test&y_train_as_test){
      train_as_test = TRUE
    }else if(!((!x_train_as_test)&(!y_train_as_test))){
      stop("testx and testy must all need or non.")
    }
  }else{
    if(x_train_as_test&y_train_as_test&z_train_as_test){
      train_as_test = TRUE
    }else if(!((!x_train_as_test)&(!y_train_as_test)&(!z_train_as_test))){
      stop("testx,testy and testz must all need or non.")
    }
  }
  if(is.numeric(y)){
    orig_y = sort(unique(y))
  }else{
    orig_y = unique(y)
  }
  if(!is.null(baseline.factor)){
    if(baseline.factor %in% orig_y){
      orig_y = orig_y[orig_y!=baseline.factor]
      orig_y = c(baseline.factor,orig_y)
    }else{
      warning("baseline.factor not in y,ignore baseline.factor.")
    }
  }


  trans_y = c(1:length(orig_y))
  new_y = rep(NA,length(y))
  for (i in 1:length(orig_y)) {
    new_y[y==orig_y[i]] = trans_y[i]
  }
  if(!is.null(testy)){
    new_test_y = rep(NA,length(testy))
    for (i in 1:length(orig_y)) {
      new_test_y[testy==orig_y[i]] = trans_y[i]
    }
  }else if(train_as_test){
    testy = y
    new_test_y = new_y
  }
  if(class(x)=="data.frame"){
    y_n = length(new_y)
    actual_x_length = length(colnames(x))
    matrix_n = ceiling(length(x[1,])^(1/2))
    model_x_length = matrix_n*matrix_n
    x_dim = c(matrix_n,matrix_n)
    while(length(x[1,]) < (matrix_n*matrix_n)){x = cbind(x,0)}
    tensor.array = array(list(),y_n)
    for (i in 1:y_n){
      tensor.array[[i]] =  array(as.numeric(x[i,]),dim=x_dim)
    }
    if(is.null(tensor_namelist)){
      tensor_namelist = colnames(x)
    }else if (length(tensor_namelist) != colnames(x)){
      tensor_namelist = colnames(x)
      warning("length of tensor_namelist not equal to length of x ,ignore tensor_namelist.")
    }
    name_array = NULL
    if(train_as_test){
      test.tensor.array = tensor.array
    }else{
      test_y_n = length(new_test_y)
      if(class(testx) != "data.frame"){
        stop("x is dataframe but testx isn't a dataframe.")
      }else{
        while(length(testx[1,]) < (matrix_n*matrix_n)){
          testx = cbind(testx,0)
        }
        test.tensor.array = array(list(),test_y_n)
        for (i in 1:test_y_n){
          test.tensor.array[[i]] =  array(as.numeric(testx[i,]),dim=x_dim)
        }
      }
    }

  }else{
    y_n = length(new_y)
    actual_x_length = length(x[[1]])
    model_x_length = length(x[[1]])
    x_dim = dim(x[[1]])
    tensor.array = x
    if(train_as_test){
      test.tensor.array = tensor.array
    }else{
      test.tensor.array = testx
    }
    if(is.null(tensor_namelist)){
      tensor_namelist = paste0("tensor-",1:actual_x_length)
    }else if(length(tensor_namelist) != actual_x_length){
      tensor_namelist = paste0("tensor-",1:actual_x_length)
      warning("length of tensor_namelist not equal to length of x ,ignore tensor_namelist.")
    }
    name_array = array(tensor_namelist,dim=x_dim)
  }
  if(is.null(rescale_x)){
  }else if(class(rescale_x)=="character"){
    if(rescale_x=="standard"){

      tensor.array = rescale_function(tensor.array,rescale_standard)
      if(train_as_test){
        test.tensor.array = tensor.array
      }else{
        test.tensor.array = rescale_function(test.tensor.array,rescale_standard)
      }

    }else if(rescale_x=="01"){

      tensor.array = rescale_function(tensor.array,rescale_01)
      if(train_as_test){
        test.tensor.array = tensor.array
      }else{
        test.tensor.array = rescale_function(test.tensor.array,rescale_01)
      }
    }else if(rescale_x=="log2"){

      tensor.array = rescale_function(tensor.array,rescale_log2)
      if(train_as_test){
        test.tensor.array = tensor.array
      }else{
        test.tensor.array = rescale_function(test.tensor.array,rescale_log2)
      }
    }else{stop("rescale_x must be NULL \"01\" \"standard\" or a rescale function.")}


  }else if(class(rescale_x)=="function"){

    tensor.array = rescale_function(tensor.array,rescale_x)
    if(train_as_test){
      test.tensor.array = tensor.array
    }else{
      test.tensor.array = rescale_function(test.tensor.array,rescale_x)
    }

  }else{stop("rescale_x must be NULL , \"01\" , \"standard\" or a rescale function.")}
  if((!is.null(z))&is.null(cov_matrix_namelist)){
    if(is.null(colnames(z))){
      cov_matrix_namelist = paste0("cov-",1:length(z[1,]))
    }else{cov_matrix_namelist = colnames(z)}
  }else if(length(cov_matrix_namelist) != length(z[1,])){
    warning("length of cov_matrix_namelist can't fit with matrix z,ignore cov_matrix_namelist.")
    if(is.null(colnames(z))){
      cov_matrix_namelist = paste0("cov-",1:length(z[1,]))
    }else{cov_matrix_namelist = colnames(z)}
  }
  if(!is.null(z)){
    z = as.matrix(z)
  }
  if(!is.null(testz)){
    testz = as.matrix(testz)
  }
  if(train_as_test){
    testz = z
  }
  if(is.null(rescale_z)){
  }else if(class(rescale_z)=="character"){
    if(rescale_z=="standard"){

      z = rescale_standard(z)
      if(train_as_test){
        testz = z
      }else{
        testz = rescale_standard(testz)
      }

    }else if(rescale_z=="01"){

      z = rescale_01(z)
      if(train_as_test){
        testz = z
      }else{
        testz = rescale_01(testz)
      }

    }else{stop("rescale_z must be NULL , \"01\" , \"standard\" or a rescale function.")}
  }else if(class(rescale_z)=="function"){

    z = rescale_z(z)
    if(train_as_test){
      testz = z
    }else{
      testz = rescale_z(testz)
    }

  }else{stop("rescale_z must be NULL , \"01\" , \"standard\" or a rescale function.")}




  output  = list()
  output$actual_x_length = actual_x_length
  output$x = x
  output$tensor.array = tensor.array
  output$testx = testx
  output$test.tensor.array = test.tensor.array

  output$trans_y = trans_y
  output$orig_y = orig_y
  output$y = y
  output$new_y = new_y
  output$testy = testy

  output$z = z
  output$testz = testz

  output$name_array = name_array
  output$tensor_namelist = tensor_namelist
  output$cov_matrix_namelist = cov_matrix_namelist
  return(output)
}


catch_model = function(ppo,lambda=NULL,max_p=ppo$actual_x_length,min_p=2,positive.factor=NULL,cv.seed=NULL,annotation_colors = NA,alpha_show_name=F,alpha_group=NA,min_abs_alpha=0){
  if(max_p<min_p){
    stop("max_p should't smaller then min_p")
  }
  if(!("catch" %in% (.packages()))){
    library("catch")
  }
  if(!("caret" %in% (.packages()))){
    library("caret")
  }
  if(!("pheatmap" %in% (.packages()))){
    library("pheatmap")
  }
  if(!("MASS" %in% (.packages()))){
    library("MASS")
  }
  in_x = c(rep(TRUE,ppo$actual_x_length),rep(FALSE,(length(ppo$tensor.array[[1]])-ppo$actual_x_length)))
  y_n = length(ppo$new_y)
  test_y_n = length(ppo$testy)
  if(!is.null(positive.factor)){
    if(!(positive.factor %in% ppo$orig_y)){
      positive.factor = NULL
      warning("positive.factor not in y,ignore baseline.factor.")
    }
  }

  set.seed(cv.seed)
  cv_model = cv.catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y)
  estimate_mcve = cv_model$cvm
  cv_model2 = catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y,
                    testx=ppo$tensor.array,testz=ppo$z)
  estimate_significant_n = cv_model2$df
  estimate_lambda = cv_model2$lambda
  estimate_y = cv_model2$pred
  estimate_precision = c()
  for(i in 1:length(estimate_y[1,])){
    pre = sum(estimate_y[,i]==ppo$new_y)/length(estimate_y[,i])
    estimate_precision = c(estimate_precision,pre)
  }
  if(sum(estimate_significant_n<max_p)==0){
    max_p = min(estimate_significant_n)
    warning("max_p is smaller then all of cross validation estimate,ignore it information.")
  }
  if(sum(estimate_significant_n>min_p)==0){
    min_p = max(estimate_significant_n)
    warning("min_p is larger then all of cross validation estimate,ignore it information.")
  }
  if(!is.null(lambda)){
    best.lambda = lambda
  }else{

    bool1 = estimate_significant_n<=max_p
    bool2 = estimate_significant_n>=min_p
    bool3 = bool1&bool2


    if(sum(bool3)>0){
      best.lambda = estimate_lambda[(estimate_mcve == min(estimate_mcve[bool3]))&bool3]
      best.lambda = best.lambda[length(best.lambda)]
    }else{
      up_find = FALSE
      down_find = FALSE
      while(!up_find){
        max_p = max_p+1
        if (max_p %in% estimate_significant_n) {
          up_lambda = estimate_lambda[estimate_significant_n==max_p]
          up_find = TRUE
        }
      }
      while(!down_find){
        min_p = min_p-1
        if (min_p %in% estimate_significant_n) {
          down_lambda = estimate_lambda[estimate_significant_n==min_p]
          down_find = TRUE
        }
      }
      best.lambda = (up_lambda+down_lambda)/2
    }

  }
  if(!is.null(ppo$z)){
    adj_model = t_adjten(x = ppo$tensor.array,z=ppo$z,y=ppo$new_y,testx=ppo$test.tensor.array,testz=ppo$testz)
    alpha = as.matrix(adj_model$alpha)
    colnames(alpha) = ppo$tensor_namelist
    row.names(alpha) = ppo$cov_matrix_namelist
    gamma = adj_model$gamma

    testxres = adj_model$testxres
    adj = adj_model$xres
    if(class(ppo$x)=="data.frame"){
      matrix_n = ceiling(length(ppo$x[1,])^(1/2))
      adj.matrix = data.frame(array(dim=c(length(ppo$testy),matrix_n*matrix_n)))
      colnames(adj.matrix) = ppo$tensor_namelist
      for (i in 1:test_y_n) {
        adj.matrix[i,] = array(testxres[[i]])
      }
      adj.matrix = adj.matrix[,1:ppo$actual_x_length]
    }else{adj.matrix = testxres}
  }
  model = catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y ,
                lambda=best.lambda,testx = ppo$test.tensor.array,testz = ppo$testz)

  model_pred = model$pred[,1]
  model_pred_prob = predict.catch.prob(model,adj_model$testxres,ppo$z,ppo$testz,adj_model$gamma)
  model_pred2 = rep(NA,length(model_pred))
  for (i in 1:length(ppo$orig_y)) {
    model_pred2[model_pred==ppo$trans_y[i]] = ppo$orig_y[i]
  }
  if(is.null(positive.factor)){
    pred_table = confusionMatrix(factor(model_pred2),factor(ppo$testy),positive=levels(factor(ppo$testy))[length(levels(factor(ppo$testy)))])
  }else{
    pred_table = confusionMatrix(factor(model_pred2),factor(ppo$testy),positive=positive.factor)
  }

  if(length(unique(ppo$y))==2){
    pred_table = list(pred_table,pred_table$byClass[7])
  }
  beta = model$beta$`1`
  if(is.null(beta)){
    beta = model$beta[[1]]
  }
  rownames(beta) = ppo$tensor_namelist
  significant = rep(FALSE,length(beta[,1]))
  for (i in 1:(length(ppo$orig_y)-1)) {
    significant = significant|(beta[,i] !=0)
  }
  select_n = paste0("select significant beta(=",sum(significant),")")
  par(mar=c(5, 5, 5, 5) )
  max_y = max(estimate_significant_n)
  min_y = min(estimate_significant_n)
  range_n = max_y-min_y

  plot(estimate_lambda,estimate_significant_n,main="The numbers of significant beta\nand precision\nfor each lambda",
       type="l",ylab="",ylim=c(min_y,max_y+(max_y-min_y)*0.6))
  mtext("The number of significant beta", side=2, line=3, cex.lab=1)


  lines(c(best.lambda,best.lambda),c(min_y-range_n,max_y+range_n),col="red",lty=2)

  range_lambda = max(estimate_lambda)-min(estimate_lambda)
  lines(c(min(estimate_lambda)-range_lambda,best.lambda),
        c(sum(significant),sum(significant)),lty=2)

  par(new=T)
  max_y = max(estimate_precision)
  min_y = min(estimate_precision)
  if(is.null(lambda)){
    select_precision = estimate_precision[estimate_lambda==best.lambda]
  }else{
    pltmodel = catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y ,
                     lambda=best.lambda,testx = ppo$tensor.array,testz = ppo$z)
    select_precision = sum(pltmodel$pred==ppo$new_y)/length(pltmodel$pred)
  }



  plot(estimate_lambda,estimate_precision,type="l",ann=F,bty="n",xaxt="n",yaxt="n",col="blueviolet",
       ylim=c(min_y,max_y+(max_y-min_y)*0.6))

  lines(c(best.lambda,max(estimate_lambda)+range_lambda),
        c(select_precision,select_precision),lty=2,col="blueviolet")

  axis(4,col.axis="blueviolet",at=round(seq(from=min_y,to=1,length=4),2))

  mtext("precision", side=4, line=3, cex.lab=1,col="blueviolet")
  legend("top",
         lty = c(1,1,2,2,2),
         cex = 0.8 ,
         col = c("black","blueviolet","red","black","blueviolet"),
         ncol=2,
         legend = c("number of significant beta","precision",
                    paste0("select lambda(=",round(best.lambda,3),")"),select_n,paste0("select precision(p=",round(select_precision,2),")"))
  )
  p1 <- recordPlot()
  select_n = paste0("select significant beta(=",sum(significant),")")
  par(mar=c(5, 5, 5, 5) )
  max_y = max(estimate_significant_n)
  min_y = min(estimate_significant_n)
  range_n = max_y-min_y
  plot(estimate_lambda,estimate_significant_n,main="The numbers of significant beta\nand mean of cross validation errors\nfor each lambda",
       type="l",ylab="",ylim=c(min_y,max_y+(max_y-min_y)*0.6))
  mtext("The number of significant beta", side=2, line=3, cex.lab=1)
  lines(c(best.lambda,best.lambda),c(min_y-range_n,max_y+range_n),col="red",lty=2)
  range_lambda = max(estimate_lambda)-min(estimate_lambda)
  lines(c(min(estimate_lambda)-range_lambda,best.lambda),
        c(sum(significant),sum(significant)),lty=2)
  par(new=T)
  max_y = max(estimate_mcve)
  min_y = min(estimate_mcve)


  plot(estimate_lambda,estimate_mcve,type="l",ann=F,bty="n",xaxt="n",yaxt="n",col="blueviolet",
       ylim=c(min_y,max_y+(max_y-min_y)*0.6))


  if(is.null(lambda)){
    select_mcve = estimate_mcve[estimate_lambda==best.lambda]
    lines(c(best.lambda,max(estimate_lambda)+range_lambda),
          c(select_mcve,select_mcve),lty=2,col="blueviolet")
  }else{
    cv_model = cv.catch(x=ppo$tensor.array , z=ppo$z , y=ppo$new_y,lambda = lambda)
    select_mcve = cv_model$cvm
    lines(c(best.lambda,max(estimate_lambda)+range_lambda),
          c(select_mcve,select_mcve),lty=2,col="blueviolet")
  }
  axis(4,col.axis="blueviolet")
  mtext("mean_cv_error", side=4, line=3, cex.lab=1,col="blueviolet")
  legend("top",
         lty = c(1,1,2,2,2),
         cex = 0.8 ,
         col = c("black","blueviolet","red","black","blueviolet"),
         ncol=2,
         legend = c("number of significant beta","mean_cv_error",
                    paste0("select lambda(=",round(best.lambda,3),")"),select_n,paste0("select mean_cv_error(=",round(select_mcve,2),")")) # 顏色所對應的名稱
  )

  p2 <- recordPlot()
  output=list()
  output$baseline_class = ppo$orig_y[1]
  output$class_list = ppo$orig_y
  output$lambda = best.lambda
  output$name_array = ppo$name_array
  if(!is.null(ppo$z)){
    output$adj.array = adj.matrix
  }
  if(sum(significant)!=0){
    output$estimate_value = beta[significant&in_x,]
  }
  output$pred = model_pred2
  output$pred_table = pred_table
  output$precision_plot = p1
  output$mcve_plot = p2
  if(sum(significant)>1){
    heatmap_matrix = array(dim=c(test_y_n,length(testxres[[1]])))
    for (i in 1:test_y_n) {
      heatmap_matrix[i,] = testxres[[i]][1:length(testxres[[1]])]
    }
    colnames(heatmap_matrix) = ppo$tensor_namelist[1:length(ppo$test.tensor.array[[1]])]
    row.names(heatmap_matrix) = 1:length(heatmap_matrix[,1])
    heatmap_beta = beta[significant&in_x,1]
    output$adj_matrix=heatmap_matrix
    heatmap_matrix = heatmap_matrix[,significant&in_x]
    heatmap_matrix = heatmap_matrix[,order(heatmap_beta,decreasing=TRUE)]
    prob = matrix(ncol=length(ppo$trans_y),nrow = length(ppo$test.tensor.array))
    prob[,1] = 0
    phi = matrix(ncol=length(ppo$trans_y),nrow = length(ppo$z[1,]))
    mu = matrix(ncol=length(ppo$trans_y),nrow = length(matrix(model$mu[1,1][[1]],ncol=1)))
    for (i in 1:length(ppo$trans_y)) {
      if(!is.null(ppo$z)){
        phi[,i] = colMeans(x=ppo$z[ppo$new_y==i,])
        mu[,i] = matrix(model$mu[i,1][[1]],ncol=1)
      }else{mu[,i] = matrix(model$mu[i,1][[1]],ncol=1)  }

    }
    f_alpha = c()

    for (i in 2:length(ppo$trans_y)) {
      if(!is.null(ppo$testz)){
        f_alpha = c(f_alpha,log(model$prior[i]/model$prior[1]) - ((phi[,i]+phi[,1]) %*% gamma[,(i-1)])/2 - sum((mu[,i]+mu[,1])/2*beta[,i-1]))
      }else{
        f_alpha = c(f_alpha,log(model$prior[i]/model$prior[1])-sum((mu[,i]+mu[,1])/2*beta[,i-1]))
        gamma = matrix(0,nrow=i-1)
      }
    }

    for (i in 2:length(ppo$trans_y)) {
      for (j in 1:length(ppo$test.tensor.array)) {
        if(!is.null(ppo$z)){
          prob[j,i] = f_alpha[i-1]+sum(ppo$testz[j,]*gamma[,i-1])+sum(matrix(testxres[[j]],ncol=1)*beta[,i-1])
        }else{prob[j,i] = f_alpha[i-1]+sum(matrix(ppo$test.tensor.array[[j]],ncol=1)*beta[,i-1])}

      }
    }
    if(length(ppo$trans_y)==2){
      p0 = ppo$testy==ppo$orig_y[1]
      p1 = ppo$testy==ppo$orig_y[2]
      mt0 = heatmap_matrix[p0,,drop=F]
      mt0 = mt0[order(model_pred_prob[p0,2]),,drop=F]
      mt1 = heatmap_matrix[p1,,drop=F]
      mt1 = mt1[order(model_pred_prob[p1,2]),,drop=F]
      heatmap_matrix = rbind(mt0,mt1)
    }else if (length(ppo$trans_y)==3){
      pp0 = model_pred2==ppo$orig_y[1]
      pp1 = model_pred2==ppo$orig_y[2]
      pp2 = model_pred2==ppo$orig_y[3]
      p0 = ppo$testy==ppo$orig_y[1]
      p1 = ppo$testy==ppo$orig_y[2]
      p2 = ppo$testy==ppo$orig_y[3]

      mt0_0 = heatmap_matrix[p0&pp0,,drop=F]
      mt0_0 = mt0_0[order(model_pred_prob[p0&pp0,1],decreasing=T),,drop=F]
      mt0_12 = heatmap_matrix[p0&(!pp0),,drop=F]
      mt0_12 = mt0_12[order(model_pred_prob[p0&(!pp0),1],decreasing=T),,drop=F]

      mt1_0 = heatmap_matrix[p1&pp0,,drop=F]
      mt1_0 = mt1_0[order(model_pred_prob[p1&pp0,2]),,drop=F]
      mt1_1 = heatmap_matrix[p1&pp1,,drop=F]
      mt1_1 = mt1_1[order(model_pred_prob[p1&pp1,2]),,drop=F]
      mt1_2 = heatmap_matrix[p1&pp2,,drop=F]
      mt1_2 = mt1_2[order(model_pred_prob[p1&pp2,2]),,drop=F]

      mt2_01 = heatmap_matrix[p2&(!pp2),,drop=F]
      mt2_01 = mt2_01[order(model_pred_prob[p2&(!pp2),3]),,drop=F]
      mt2_2 = heatmap_matrix[p2&pp2,,drop=F]
      mt2_2 = mt2_2[order(model_pred_prob[p2&pp2,3]),,drop=F]

      heatmap_matrix = rbind(mt0_0,mt0_12,mt1_0,mt1_1,mt1_2,mt2_01,mt2_2)
    }

    output$prob = prob
    output$prob2 = model_pred_prob
    heatmap_y = data.frame(y_reference = factor(ppo$testy))
    output$heatmap_matrix = heatmap_matrix
    output$heatmap_y = heatmap_y
    output$tensor_heatmap = tryCatch({
      pheatmap(t(heatmap_matrix),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=heatmap_y,annotation_colors=annotation_colors,show_colnames=F)
    },error = function(err) {
      tryCatch({
        pheatmap(t(heatmap_matrix),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=heatmap_y,show_colnames=F)
        print("annotation_colors define error ,please ignore")
      },error = function(err) {past0("Error in plot pheatmap:","err")})
    })

    for (i in 1:length(heatmap_matrix[1,])) {
      if(sd(heatmap_matrix[,i])==0){
        heatmap_matrix[,i] = rep(0,length(heatmap_matrix[,i]))
      }else{
        heatmap_matrix[,i] = (heatmap_matrix[,i]-mean(heatmap_matrix[,i]))/sd(heatmap_matrix[,i])
      }

    }
    output$heatmap_matrix2 = heatmap_matrix
    output$tensor_heatmap2 = tryCatch({
      pheatmap(t(heatmap_matrix),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=heatmap_y,annotation_colors=annotation_colors,show_colnames=F)
    },error = function(err) {
      tryCatch({
        pheatmap(t(heatmap_matrix),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=heatmap_y,show_colnames=F)
        print("annotation_colors define error ,please ignore")
      },error = function(err) {past0("Error in plot pheatmap:","err")})
    })


    if(!is.null(ppo$z)){
      heatmap_matrix2 = alpha[,significant&in_x]
      heatmap_matrix2 = heatmap_matrix2[,order(heatmap_beta,decreasing=TRUE)]
      if(!is.null(alpha_group)){
        mt = heatmap_matrix2[alpha_group[,1]==unique(alpha_group[,1])[1],,drop=F]
        bool = rowSums(abs(mt)>min_abs_alpha)>0
        if((sum(bool)!=0)&(sum(!bool)!=0)){
          mtT = mt[bool,,drop=F]
          mtF = mt[!bool,,drop=F]
          mt = rbind(mtF,mtT)
        }
        for (i in 2:length(unique(alpha_group[,1]))) {
          mt0 = heatmap_matrix2[alpha_group[,1]==unique(alpha_group[,1])[i],,drop=F]
          bool = rowSums(abs(mt0)>min_abs_alpha)>0
          if((sum(bool)!=0)&(sum(!bool)!=0)){
            mtT = mt0[bool,,drop=F]
            mtF = mt0[!bool,,drop=F]
            mt0 = rbind(mtF,mtT)
          }
          mt = rbind(mt,mt0)
        }
        heatmap_matrix2 = mt
      }
      output$alpha_matrix = t(heatmap_matrix2)
      output$alpha_heatmap = pheatmap(t(heatmap_matrix2),cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=alpha_group,show_colnames = alpha_show_name)
    }
  }else{warning("Less then 2 significant tensor,no estimate_value、tensor_heatmap、alpha_heatmap")}
  return(output)
}
fit_and_catch_model <- function(x, y, z, positive_factor, pathway_MSI) {
  fit <- preprocessing(x = x, y = y, z = z,
                       testx = x, testy = y, testz = z,
                       rescale_x = NULL, rescale_z = NULL, baseline.factor = "NonMetastasis")
  model <- catch_model(ppo = fit, max_p = 15, min_p = 5, cv.seed = 1,
                       alpha_group = pathway_MSI, min_abs_alpha = 0, positive.factor = positive_factor)
  return(model)
}
make_supplementary1 <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney){
  GI_GYN_acc <- vector("list", length = 5)
  GI_GYN_acc <- lapply(list(GI_brain$pred_table[[1]],GI_lung$pred_table[[1]],GI_liver$pred_table[[1]],GI_bone$pred_table[[1]],GI_kidney$pred_table[[1]]), function(tbl) tbl$overall[[1]])
  GI_GYN_sen <- vector("list", length = 5)
  GI_GYN_sen <- lapply(list(GI_brain$pred_table[[1]],GI_lung$pred_table[[1]],GI_liver$pred_table[[1]],GI_bone$pred_table[[1]],GI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[1]])
  GI_GYN_spe <- vector("list", length = 5)
  GI_GYN_spe <- lapply(list(GI_brain$pred_table[[1]],GI_lung$pred_table[[1]],GI_liver$pred_table[[1]],GI_bone$pred_table[[1]],GI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[2]])
  GI_GYN_pre <- vector("list", length = 5)
  GI_GYN_pre <- lapply(list(GI_brain$pred_table[[1]],GI_lung$pred_table[[1]],GI_liver$pred_table[[1]],GI_bone$pred_table[[1]],GI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[5]])
  GI_GYN_f1 <- vector("list", length = 5)
  GI_GYN_f1 <- lapply(list(GI_brain$pred_table[[1]],GI_lung$pred_table[[1]],GI_liver$pred_table[[1]],GI_bone$pred_table[[1]],GI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[7]])
  GI_GYN_out_df=cbind(c("Brain","Lung","Liver","Bone","Kidney"),rep(c("GI_GYN"),5),GI_GYN_acc,GI_GYN_sen,GI_GYN_spe,GI_GYN_pre,GI_GYN_f1)
  colnames(GI_GYN_out_df)=c("Metastasis","Cancers","Accuracy","Sensitivity","Specificity","Precision","F1-score")
  Lung_acc <- vector("list", length = 5)
  Lung_acc <- lapply(list(Lung_brain$pred_table[[1]],Lung_lung$pred_table[[1]],Lung_liver$pred_table[[1]],Lung_bone$pred_table[[1]],Lung_kidney$pred_table[[1]]), function(tbl) tbl$overall[[1]])
  Lung_sen <- vector("list", length = 5)
  Lung_sen <- lapply(list(Lung_brain$pred_table[[1]],Lung_lung$pred_table[[1]],Lung_liver$pred_table[[1]],Lung_bone$pred_table[[1]],Lung_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[1]])
  Lung_spe <- vector("list", length = 5)
  Lung_spe <- lapply(list(Lung_brain$pred_table[[1]],Lung_lung$pred_table[[1]],Lung_liver$pred_table[[1]],Lung_bone$pred_table[[1]],Lung_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[2]])
  Lung_pre <- vector("list", length = 5)
  Lung_pre <- lapply(list(Lung_brain$pred_table[[1]],Lung_lung$pred_table[[1]],Lung_liver$pred_table[[1]],Lung_bone$pred_table[[1]],Lung_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[5]])
  Lung_f1 <- vector("list", length = 5)
  Lung_f1 <- lapply(list(Lung_brain$pred_table[[1]],Lung_lung$pred_table[[1]],Lung_liver$pred_table[[1]],Lung_bone$pred_table[[1]],Lung_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[7]])
  Lung_out_df=cbind(c("Brain","Lung","Liver","Bone","Kidney"),rep(c("Lung"),5),Lung_acc,Lung_sen,Lung_spe,Lung_pre,Lung_f1)
  colnames(Lung_out_df)=c("Metastasis","Cancers","Accuracy","Sensitivity","Specificity","Precision","F1-score")
  MSI_acc <- vector("list", length = 5)
  MSI_acc <- lapply(list(MSI_brain$pred_table[[1]],MSI_lung$pred_table[[1]],MSI_liver$pred_table[[1]],MSI_bone$pred_table[[1]],MSI_kidney$pred_table[[1]]), function(tbl) tbl$overall[[1]])
  MSI_sen <- vector("list", length = 5)
  MSI_sen <- lapply(list(MSI_brain$pred_table[[1]],MSI_lung$pred_table[[1]],MSI_liver$pred_table[[1]],MSI_bone$pred_table[[1]],MSI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[1]])
  MSI_spe <- vector("list", length = 5)
  MSI_spe <- lapply(list(MSI_brain$pred_table[[1]],MSI_lung$pred_table[[1]],MSI_liver$pred_table[[1]],MSI_bone$pred_table[[1]],MSI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[2]])
  MSI_pre <- vector("list", length = 5)
  MSI_pre <- lapply(list(MSI_brain$pred_table[[1]],MSI_lung$pred_table[[1]],MSI_liver$pred_table[[1]],MSI_bone$pred_table[[1]],MSI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[5]])
  MSI_f1 <- vector("list", length = 5)
  MSI_f1 <- lapply(list(MSI_brain$pred_table[[1]],MSI_lung$pred_table[[1]],MSI_liver$pred_table[[1]],MSI_bone$pred_table[[1]],MSI_kidney$pred_table[[1]]), function(tbl) tbl$byClass[[7]])
  MSI_out_df=cbind(c("Brain","Lung","Liver","Bone","Kidney"),rep(c("MSI"),5),MSI_acc,MSI_sen,MSI_spe,MSI_pre,MSI_f1)
  colnames(MSI_out_df)=c("Metastasis","Cancers","Accuracy","Sensitivity","Specificity","Precision","F1-score")
  table_1=rbind(GI_GYN_out_df,Lung_out_df,MSI_out_df)
  table_1[] <- lapply(table_1, function(x) ifelse(is.numeric(x), round(x, 3), x))
  table_1=data.frame(table_1)
  sort_order <- c("Brain", "Lung", "Liver", "Bone", "Kidney")
  table_1$Metastasis <- factor(table_1$Metastasis,levels = sort_order)
  table_1 <- table_1[order(table_1$Metastasis), ]
  table_1=data.frame(table_1)
  table_1$Cancers=unlist(table_1$Cancers)
  table_1$Accuracy=unlist(table_1$Accuracy)
  table_1$Sensitivity=unlist(table_1$Sensitivity)
  table_1$Specificity=unlist(table_1$Specificity)
  table_1$Precision=unlist(table_1$Precision)
  table_1$F1.score=unlist(table_1$F1.score)
  new_folder_path <- "supplementary"
  suppressWarnings({dir.create(new_folder_path)})
  write.csv(table_1,"supplementary/table_1.csv")
}
create_df <- function(data, metastasis_site, cancer_type) {
  df <- data.frame(data$estimate_value)
  rownames <- rownames(df)
  df <- cbind(rep(metastasis_site, length(data$estimate_value)),
              rep(cancer_type, length(data$estimate_value)),
              rownames,
              df)
  colnames(df) <- c("Metastasis.Sites", "Cancers", "Features", "Coefficient")
  return(df)
}
make_supplementary2 <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney){
  GI_GYN_brain_beta_df <- create_df(GI_brain, "Brain", "GI/GYN")
  lung_brain_beta_df <- create_df(Lung_brain, "Brain", "Lung")
  MSI_brain_beta_df <- create_df(MSI_brain, "Brain", "MSI")
  GI_GYN_lung_beta_df <- create_df(GI_lung, "Lung", "GI/GYN")
  lung_lung_beta_df <- create_df(Lung_lung, "Lung", "Lung")
  MSI_lung_beta_df <- create_df(MSI_lung, "Lung", "MSI")
  GI_GYN_liver_beta_df <- create_df(GI_liver, "Liver", "GI/GYN")
  lung_liver_beta_df <- create_df(Lung_liver, "Liver", "Lung")
  MSI_liver_beta_df <- create_df(MSI_liver, "Liver", "MSI")
  GI_GYN_bone_beta_df <- create_df(GI_bone, "Bone", "GI/GYN")
  lung_bone_beta_df <- create_df(Lung_bone, "Bone", "Lung")
  MSI_bone_beta_df <- create_df(MSI_bone, "Bone", "MSI")
  GI_GYN_kidney_beta_df <- create_df(GI_kidney, "Kidney", "GI/GYN")
  lung_kidney_beta_df <- create_df(Lung_kidney, "Kidney", "Lung")
  MSI_kidney_beta_df <- create_df(MSI_kidney, "Kidney", "MSI")
  df_table2_brain <- rbind(GI_GYN_brain_beta_df, lung_brain_beta_df, MSI_brain_beta_df)
  df_table2_lung <- rbind(GI_GYN_lung_beta_df, lung_lung_beta_df, MSI_lung_beta_df)
  df_table2_liver <- rbind(GI_GYN_liver_beta_df, lung_liver_beta_df, MSI_liver_beta_df)
  df_table2_bone <- rbind(GI_GYN_bone_beta_df, lung_bone_beta_df, MSI_bone_beta_df)
  df_table2_kidney <- rbind(GI_GYN_kidney_beta_df, lung_kidney_beta_df, MSI_kidney_beta_df)
  df_table2_brain <- df_table2_brain[order(abs(df_table2_brain$Coefficient), decreasing = TRUE), ]
  df_table2_lung <- df_table2_lung[order(abs(df_table2_lung$Coefficient), decreasing = TRUE), ]
  df_table2_liver <- df_table2_liver[order(abs(df_table2_liver$Coefficient), decreasing = TRUE), ]
  df_table2_bone <- df_table2_bone[order(abs(df_table2_bone$Coefficient), decreasing = TRUE), ]
  df_table2_kidney <- df_table2_kidney[order(abs(df_table2_kidney$Coefficient), decreasing = TRUE), ]
  df_table2=rbind(df_table2_brain, df_table2_lung,
                  df_table2_liver, df_table2_bone,
                  df_table2_kidney)
  df_table2_metastasis_4 <- rbind(df_table2_brain[1:4, ], df_table2_lung[1:4, ],
                                  df_table2_liver[1:4, ], df_table2_bone[1:4, ],
                                  df_table2_kidney[1:4, ])
  write.csv(df_table2,"supplementary/Figure_2and3_Upset_Plot_Supplementary_Data_4.csv",row.names = F)
  write.csv(df_table2_metastasis_4,"supplementary/Figure_1_Sankey_Supplementary_Data_4.csv",row.names = F)
  return(invisible(NULL))
}
make_supplementary4 <- function(pathway_MSI=pathway_MSI,pathway_GI=pathway_GI,pathway_LUNG=pathway_LUNG){
  amino_14 <- rownames(pathway_MSI)[pathway_MSI$pathway == "Metabolism of amino acids and derivatives"]
  carbohydrates_14 <- rownames(pathway_MSI)[pathway_MSI$pathway == "Metabolism of carbohydrates"]
  lipids_14 <- rownames(pathway_MSI)[pathway_MSI$pathway == "Metabolism of lipids"]
  nucleotides_14 <- rownames(pathway_MSI)[pathway_MSI$pathway == "Metabolism of nucleotides"]
  TCA_14 <- rownames(pathway_MSI)[pathway_MSI$pathway == "The citric acid (TCA) cycle and respiratory electron transport"]

  max_L14 <- max(length(amino_14), length(carbohydrates_14), length(lipids_14), length(nucleotides_14), length(TCA_14))

  amino_14 <- c(amino_14, rep("", max_L14 - length(amino_14)))
  carbohydrates_14 <- c(carbohydrates_14, rep("", max_L14 - length(carbohydrates_14)))
  lipids_14 <- c(lipids_14, rep("", max_L14 - length(lipids_14)))
  nucleotides_14 <- c(nucleotides_14, rep("", max_L14 - length(nucleotides_14)))
  TCA_14 <- c(TCA_14, rep("", max_L14 - length(TCA_14)))

  supplementary_4_14 <- data.frame("amino acids and derivatives" = amino_14,
                                   "carbohydrates" = carbohydrates_14,
                                   "lipids" = lipids_14,
                                   "nucleotides" = nucleotides_14,
                                   "The citric acid (TCA)" = TCA_14)
  vitamins_34 <- rownames(pathway_GI)[pathway_GI$pathway == "Metabolism of vitamins and cofactors"]
  amino_34 <- rownames(pathway_GI)[pathway_GI$pathway == "Metabolism of amino acids and derivatives"]
  carbohydrates_34 <- rownames(pathway_GI)[pathway_GI$pathway == "Metabolism of carbohydrates"]
  lipids_34 <- rownames(pathway_GI)[pathway_GI$pathway == "Metabolism of lipids"]
  Biological_34 <- rownames(pathway_GI)[pathway_GI$pathway == "Biological oxidations"]
  nucleotides_34 <- rownames(pathway_GI)[pathway_GI$pathway == "Metabolism of nucleotides"]
  TCA_34 <- rownames(pathway_GI)[pathway_GI$pathway == "The citric acid (TCA) cycle and respiratory electron transport"]

  max_L34 <- max(length(amino_34), length(carbohydrates_34), length(lipids_34), length(nucleotides_34),
                 length(TCA_34), length(Biological_34), length(vitamins_34))

  amino_34 <- c(amino_34, rep("", max_L34 - length(amino_34)))
  carbohydrates_34 <- c(carbohydrates_34, rep("", max_L34 - length(carbohydrates_34)))
  lipids_34 <- c(lipids_34, rep("", max_L34 - length(lipids_34)))
  nucleotides_34 <- c(nucleotides_34, rep("", max_L34 - length(nucleotides_34)))
  TCA_34 <- c(TCA_34, rep("", max_L34 - length(TCA_34)))
  Biological_34 <- c(Biological_34, rep("", max_L34 - length(Biological_34)))
  vitamins_34 <- c(vitamins_34, rep("", max_L34 - length(vitamins_34)))

  supplementary_4_34 <- data.frame("amino acids and derivatives" = amino_34,
                                   "carbohydrates" = carbohydrates_34,
                                   "lipids" = lipids_34,
                                   "nucleotides" = nucleotides_34,
                                   "The citric acid (TCA)" = TCA_34,
                                   "vitamins and cofactors" = vitamins_34,
                                   "Biological oxidations" = Biological_34)
  vitamins_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "Metabolism of vitamins and cofactors"]
  amino_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "Metabolism of amino acids and derivatives"]
  carbohydrates_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "Metabolism of carbohydrates"]
  lipids_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "Metabolism of lipids"]
  Biological_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "Biological oxidations"]
  nucleotides_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "Metabolism of nucleotides"]
  TCA_42 <- rownames(pathway_LUNG)[pathway_LUNG$pathway == "The citric acid (TCA) cycle and respiratory electron transport"]

  max_L42 <- max(length(amino_42), length(carbohydrates_42), length(lipids_42), length(nucleotides_42),
                 length(TCA_42), length(Biological_42), length(vitamins_42))

  amino_42 <- c(amino_42, rep("", max_L42 - length(amino_42)))
  carbohydrates_42 <- c(carbohydrates_42, rep("", max_L42 - length(carbohydrates_42)))
  lipids_42 <- c(lipids_42, rep("", max_L42 - length(lipids_42)))
  nucleotides_42 <- c(nucleotides_42, rep("", max_L42 - length(nucleotides_42)))
  TCA_42 <- c(TCA_42, rep("", max_L42 - length(TCA_42)))
  Biological_42 <- c(Biological_42, rep("", max_L42 - length(Biological_42)))
  vitamins_42 <- c(vitamins_42, rep("", max_L42 - length(vitamins_42)))
  library(openxlsx)
  supplementary_4_42 <- data.frame("amino acids and derivatives" = amino_42,
                                   "carbohydrates" = carbohydrates_42,
                                   "lipids" = lipids_42,
                                   "nucleotides" = nucleotides_42,
                                   "The citric acid (TCA)" = TCA_42,
                                   "vitamins and cofactors" = vitamins_42,
                                   "Biological oxidations" = Biological_42)
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "14 Gene")
  writeData(wb, sheet = "14 Gene", supplementary_4_14)
  addWorksheet(wb, sheetName = "34 Gene")
  writeData(wb, sheet = "34 Gene", supplementary_4_34)
  addWorksheet(wb, sheetName = "42 Gene")
  writeData(wb, sheet = "42 Gene", supplementary_4_42)
  saveWorkbook(wb, "supplementary/supplementary_4.xlsx",overwrite = TRUE)
  return(invisible(NULL))
}
make_supplementary7 <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,pathway_df=pathway_df){
  pathway_full <- data.frame(pathway = pathway_df$pathway, row.names = pathway_df$gene)
  colname_brain <- union(union(colnames(t(GI_brain$alpha_matrix)), colnames(t(Lung_brain$alpha_matrix))), colnames(t(MSI_brain$alpha_matrix)))

  GI_brain_s6 <- cbind(rep("GI and GYN", length(rownames(t(GI_brain$alpha_matrix)))), rownames(t(GI_brain$alpha_matrix)), t(GI_brain$alpha_matrix))
  colnames(GI_brain_s6) <- c("Cancers", "Gene", colnames(t(GI_brain$alpha_matrix)))

  Lung_brain_s6 <- cbind(rep("Lung", length(rownames(t(Lung_brain$alpha_matrix)))), rownames(t(Lung_brain$alpha_matrix)), t(Lung_brain$alpha_matrix))
  colnames(Lung_brain_s6) <- c("Cancers", "Gene", colnames(t(Lung_brain$alpha_matrix)))

  MSI_brain_s6 <- cbind(rep("MSI", length(rownames(t(MSI_brain$alpha_matrix)))), rownames(t(MSI_brain$alpha_matrix)), t(MSI_brain$alpha_matrix))
  colnames(MSI_brain_s6) <- c("Cancers", "Gene", colnames(t(MSI_brain$alpha_matrix)))

  colname_brain <- c("Cancers", "Gene", colname_brain)

  GI_brain_s6_cname <- colnames(GI_brain_s6)
  Lung_brain_s6_cname <- colnames(Lung_brain_s6)
  MSI_brain_s6_cname <- colnames(MSI_brain_s6)

  for (i in 1:length(colname_brain)) {
    if (!(colname_brain[i] %in% colnames(GI_brain_s6))) {
      GI_brain_s6_cname <- c(GI_brain_s6_cname, colname_brain[i])
      GI_brain_s6 <- cbind(GI_brain_s6, rep(0, length(GI_brain_s6[,1])))
    }
  }

  colnames(GI_brain_s6) <- GI_brain_s6_cname

  for (i in 1:length(colname_brain)) {
    if (!(colname_brain[i] %in% colnames(Lung_brain_s6))) {
      Lung_brain_s6_cname <- c(Lung_brain_s6_cname, colname_brain[i])
      Lung_brain_s6 <- cbind(Lung_brain_s6, rep(0, length(Lung_brain_s6[,1])))
    }
  }

  colnames(Lung_brain_s6) <- Lung_brain_s6_cname

  for (i in 1:length(colname_brain)) {
    if (!(colname_brain[i] %in% colnames(MSI_brain_s6))) {
      MSI_brain_s6_cname <- c(MSI_brain_s6_cname, colname_brain[i])
      MSI_brain_s6 <- cbind(MSI_brain_s6, rep(0, length(MSI_brain_s6[,1])))
    }
  }

  colnames(MSI_brain_s6) <- MSI_brain_s6_cname

  GI_brain_s6 <- GI_brain_s6[, colname_brain]
  Lung_brain_s6 <- Lung_brain_s6[, colname_brain]
  MSI_brain_s6 <- MSI_brain_s6[, colname_brain]

  brain_s6 <- rbind(GI_brain_s6, Lung_brain_s6, MSI_brain_s6)

  Pathway_col <- sapply(brain_s6[,2], function(x) pathway_full[match(x, rownames(pathway_full)), 1])
  brain_s6 <- cbind(brain_s6, Pathway_col)
  colnames(brain_s6) <- c(colname_brain, "Pathway")
  colname_lung <- union(union(colnames(t(GI_lung$alpha_matrix)), colnames(t(Lung_lung$alpha_matrix))), colnames(t(MSI_lung$alpha_matrix)))

  GI_lung_s6 <- cbind(rep("GI and GYN", length(rownames(t(GI_lung$alpha_matrix)))), rownames(t(GI_lung$alpha_matrix)), t(GI_lung$alpha_matrix))
  colnames(GI_lung_s6) <- c("Cancers", "Gene", colnames(t(GI_lung$alpha_matrix)))

  Lung_lung_s6 <- cbind(rep("Lung", length(rownames(t(Lung_lung$alpha_matrix)))), rownames(t(Lung_lung$alpha_matrix)), t(Lung_lung$alpha_matrix))
  colnames(Lung_lung_s6) <- c("Cancers", "Gene", colnames(t(Lung_lung$alpha_matrix)))

  MSI_lung_s6 <- cbind(rep("MSI", length(rownames(t(MSI_lung$alpha_matrix)))), rownames(t(MSI_lung$alpha_matrix)), t(MSI_lung$alpha_matrix))
  colnames(MSI_lung_s6) <- c("Cancers", "Gene", colnames(t(MSI_lung$alpha_matrix)))

  colname_lung <- c("Cancers", "Gene", colname_lung)

  GI_lung_s6_cname <- colnames(GI_lung_s6)
  Lung_lung_s6_cname <- colnames(Lung_lung_s6)
  MSI_lung_s6_cname <- colnames(MSI_lung_s6)

  for (i in 1:length(colname_lung)) {
    if (!(colname_lung[i] %in% colnames(GI_lung_s6))) {
      GI_lung_s6_cname <- c(GI_lung_s6_cname, colname_lung[i])
      GI_lung_s6 <- cbind(GI_lung_s6, rep(0, length(GI_lung_s6[,1])))
    }
  }

  colnames(GI_lung_s6) <- GI_lung_s6_cname

  for (i in 1:length(colname_lung)) {
    if (!(colname_lung[i] %in% colnames(Lung_lung_s6))) {
      Lung_lung_s6_cname <- c(Lung_lung_s6_cname, colname_lung[i])
      Lung_lung_s6 <- cbind(Lung_lung_s6, rep(0, length(Lung_lung_s6[,1])))
    }
  }

  colnames(Lung_lung_s6) <- Lung_lung_s6_cname

  for (i in 1:length(colname_lung)) {
    if (!(colname_lung[i] %in% colnames(MSI_lung_s6))) {
      MSI_lung_s6_cname <- c(MSI_lung_s6_cname, colname_lung[i])
      MSI_lung_s6 <- cbind(MSI_lung_s6, rep(0, length(MSI_lung_s6[,1])))
    }
  }

  colnames(MSI_lung_s6) <- MSI_lung_s6_cname

  GI_lung_s6 <- GI_lung_s6[, colname_lung]
  Lung_lung_s6 <- Lung_lung_s6[, colname_lung]
  MSI_lung_s6 <- MSI_lung_s6[, colname_lung]

  lung_s6 <- rbind(GI_lung_s6, Lung_lung_s6, MSI_lung_s6)

  Pathway_col <- sapply(lung_s6[,2], function(x) pathway_full[match(x, rownames(pathway_full)), 1])
  lung_s6 <- cbind(lung_s6, Pathway_col)
  colnames(lung_s6) <- c(colname_lung, "Pathway")
  pathway_full <- data.frame(pathway = pathway_df$pathway, row.names = pathway_df$gene)
  colname_liver <- union(union(colnames(t(GI_liver$alpha_matrix)), colnames(t(Lung_liver$alpha_matrix))), colnames(t(MSI_liver$alpha_matrix)))

  GI_liver_s6 <- cbind(rep("GI and GYN", length(rownames(t(GI_liver$alpha_matrix)))), rownames(t(GI_liver$alpha_matrix)), t(GI_liver$alpha_matrix))
  colnames(GI_liver_s6) <- c("Cancers", "Gene", colnames(t(GI_liver$alpha_matrix)))

  Lung_liver_s6 <- cbind(rep("Lung", length(rownames(t(Lung_liver$alpha_matrix)))), rownames(t(Lung_liver$alpha_matrix)), t(Lung_liver$alpha_matrix))
  colnames(Lung_liver_s6) <- c("Cancers", "Gene", colnames(t(Lung_liver$alpha_matrix)))

  MSI_liver_s6 <- cbind(rep("MSI", length(rownames(t(MSI_liver$alpha_matrix)))), rownames(t(MSI_liver$alpha_matrix)), t(MSI_liver$alpha_matrix))
  colnames(MSI_liver_s6) <- c("Cancers", "Gene", colnames(t(MSI_liver$alpha_matrix)))

  colname_liver <- c("Cancers", "Gene", colname_liver)

  GI_liver_s6_cname <- colnames(GI_liver_s6)
  Lung_liver_s6_cname <- colnames(Lung_liver_s6)
  MSI_liver_s6_cname <- colnames(MSI_liver_s6)

  for (i in 1:length(colname_liver)) {
    if (!(colname_liver[i] %in% colnames(GI_liver_s6))) {
      GI_liver_s6_cname <- c(GI_liver_s6_cname, colname_liver[i])
      GI_liver_s6 <- cbind(GI_liver_s6, rep(0, length(GI_liver_s6[,1])))
    }
  }

  colnames(GI_liver_s6) <- GI_liver_s6_cname

  for (i in 1:length(colname_liver)) {
    if (!(colname_liver[i] %in% colnames(Lung_liver_s6))) {
      Lung_liver_s6_cname <- c(Lung_liver_s6_cname, colname_liver[i])
      Lung_liver_s6 <- cbind(Lung_liver_s6, rep(0, length(Lung_liver_s6[,1])))
    }
  }

  colnames(Lung_liver_s6) <- Lung_liver_s6_cname

  for (i in 1:length(colname_liver)) {
    if (!(colname_liver[i] %in% colnames(MSI_liver_s6))) {
      MSI_liver_s6_cname <- c(MSI_liver_s6_cname, colname_liver[i])
      MSI_liver_s6 <- cbind(MSI_liver_s6, rep(0, length(MSI_liver_s6[,1])))
    }
  }

  colnames(MSI_liver_s6) <- MSI_liver_s6_cname

  GI_liver_s6 <- GI_liver_s6[, colname_liver]
  Lung_liver_s6 <- Lung_liver_s6[, colname_liver]
  MSI_liver_s6 <- MSI_liver_s6[, colname_liver]

  liver_s6 <- rbind(GI_liver_s6, Lung_liver_s6, MSI_liver_s6)

  Pathway_col <- sapply(liver_s6[,2], function(x) pathway_full[match(x, rownames(pathway_full)), 1])
  liver_s6 <- cbind(liver_s6, Pathway_col)
  colnames(liver_s6) <- c(colname_liver, "Pathway")
  pathway_full <- data.frame(pathway = pathway_df$pathway, row.names = pathway_df$gene)
  colname_bone <- union(union(colnames(t(GI_bone$alpha_matrix)), colnames(t(Lung_bone$alpha_matrix))), colnames(t(MSI_bone$alpha_matrix)))

  GI_bone_s6 <- cbind(rep("GI and GYN", length(rownames(t(GI_bone$alpha_matrix)))), rownames(t(GI_bone$alpha_matrix)), t(GI_bone$alpha_matrix))
  colnames(GI_bone_s6) <- c("Cancers", "Gene", colnames(t(GI_bone$alpha_matrix)))

  Lung_bone_s6 <- cbind(rep("Lung", length(rownames(t(Lung_bone$alpha_matrix)))), rownames(t(Lung_bone$alpha_matrix)), t(Lung_bone$alpha_matrix))
  colnames(Lung_bone_s6) <- c("Cancers", "Gene", colnames(t(Lung_bone$alpha_matrix)))

  MSI_bone_s6 <- cbind(rep("MSI", length(rownames(t(MSI_bone$alpha_matrix)))), rownames(t(MSI_bone$alpha_matrix)), t(MSI_bone$alpha_matrix))
  colnames(MSI_bone_s6) <- c("Cancers", "Gene", colnames(t(MSI_bone$alpha_matrix)))

  colname_bone <- c("Cancers", "Gene", colname_bone)

  GI_bone_s6_cname <- colnames(GI_bone_s6)
  Lung_bone_s6_cname <- colnames(Lung_bone_s6)
  MSI_bone_s6_cname <- colnames(MSI_bone_s6)

  for (i in 1:length(colname_bone)) {
    if (!(colname_bone[i] %in% colnames(GI_bone_s6))) {
      GI_bone_s6_cname <- c(GI_bone_s6_cname, colname_bone[i])
      GI_bone_s6 <- cbind(GI_bone_s6, rep(0, length(GI_bone_s6[,1])))
    }
  }

  colnames(GI_bone_s6) <- GI_bone_s6_cname

  for (i in 1:length(colname_bone)) {
    if (!(colname_bone[i] %in% colnames(Lung_bone_s6))) {
      Lung_bone_s6_cname <- c(Lung_bone_s6_cname, colname_bone[i])
      Lung_bone_s6 <- cbind(Lung_bone_s6, rep(0, length(Lung_bone_s6[,1])))
    }
  }

  colnames(Lung_bone_s6) <- Lung_bone_s6_cname

  for (i in 1:length(colname_bone)) {
    if (!(colname_bone[i] %in% colnames(MSI_bone_s6))) {
      MSI_bone_s6_cname <- c(MSI_bone_s6_cname, colname_bone[i])
      MSI_bone_s6 <- cbind(MSI_bone_s6, rep(0, length(MSI_bone_s6[,1])))
    }
  }

  colnames(MSI_bone_s6) <- MSI_bone_s6_cname

  GI_bone_s6 <- GI_bone_s6[, colname_bone]
  Lung_bone_s6 <- Lung_bone_s6[, colname_bone]
  MSI_bone_s6 <- MSI_bone_s6[, colname_bone]

  bone_s6 <- rbind(GI_bone_s6, Lung_bone_s6, MSI_bone_s6)

  Pathway_col <- sapply(bone_s6[,2], function(x) pathway_full[match(x, rownames(pathway_full)), 1])
  bone_s6 <- cbind(bone_s6, Pathway_col)
  colnames(bone_s6) <- c(colname_bone, "Pathway")
  pathway_full <- data.frame(pathway = pathway_df$pathway, row.names = pathway_df$gene)
  colname_kidney <- union(union(colnames(t(GI_kidney$alpha_matrix)), colnames(t(Lung_kidney$alpha_matrix))), colnames(t(MSI_kidney$alpha_matrix)))

  GI_kidney_s6 <- cbind(rep("GI and GYN", length(rownames(t(GI_kidney$alpha_matrix)))), rownames(t(GI_kidney$alpha_matrix)), t(GI_kidney$alpha_matrix))
  colnames(GI_kidney_s6) <- c("Cancers", "Gene", colnames(t(GI_kidney$alpha_matrix)))

  Lung_kidney_s6 <- cbind(rep("Lung", length(rownames(t(Lung_kidney$alpha_matrix)))), rownames(t(Lung_kidney$alpha_matrix)), t(Lung_kidney$alpha_matrix))
  colnames(Lung_kidney_s6) <- c("Cancers", "Gene", colnames(t(Lung_kidney$alpha_matrix)))

  MSI_kidney_s6 <- cbind(rep("MSI", length(rownames(t(MSI_kidney$alpha_matrix)))), rownames(t(MSI_kidney$alpha_matrix)), t(MSI_kidney$alpha_matrix))
  colnames(MSI_kidney_s6) <- c("Cancers", "Gene", colnames(t(MSI_kidney$alpha_matrix)))

  colname_kidney <- c("Cancers", "Gene", colname_kidney)

  GI_kidney_s6_cname <- colnames(GI_kidney_s6)
  Lung_kidney_s6_cname <- colnames(Lung_kidney_s6)
  MSI_kidney_s6_cname <- colnames(MSI_kidney_s6)

  for (i in 1:length(colname_kidney)) {
    if (!(colname_kidney[i] %in% colnames(GI_kidney_s6))) {
      GI_kidney_s6_cname <- c(GI_kidney_s6_cname, colname_kidney[i])
      GI_kidney_s6 <- cbind(GI_kidney_s6, rep(0, length(GI_kidney_s6[,1])))
    }
  }

  colnames(GI_kidney_s6) <- GI_kidney_s6_cname

  for (i in 1:length(colname_kidney)) {
    if (!(colname_kidney[i] %in% colnames(Lung_kidney_s6))) {
      Lung_kidney_s6_cname <- c(Lung_kidney_s6_cname, colname_kidney[i])
      Lung_kidney_s6 <- cbind(Lung_kidney_s6, rep(0, length(Lung_kidney_s6[,1])))
    }
  }

  colnames(Lung_kidney_s6) <- Lung_kidney_s6_cname

  for (i in 1:length(colname_kidney)) {
    if (!(colname_kidney[i] %in% colnames(MSI_kidney_s6))) {
      MSI_kidney_s6_cname <- c(MSI_kidney_s6_cname, colname_kidney[i])
      MSI_kidney_s6 <- cbind(MSI_kidney_s6, rep(0, length(MSI_kidney_s6[,1])))
    }
  }

  colnames(MSI_kidney_s6) <- MSI_kidney_s6_cname

  GI_kidney_s6 <- GI_kidney_s6[, colname_kidney]
  Lung_kidney_s6 <- Lung_kidney_s6[, colname_kidney]
  MSI_kidney_s6 <- MSI_kidney_s6[, colname_kidney]

  kidney_s6 <- rbind(GI_kidney_s6, Lung_kidney_s6, MSI_kidney_s6)

  Pathway_col <- sapply(kidney_s6[,2], function(x) pathway_full[match(x, rownames(pathway_full)), 1])
  kidney_s6 <- cbind(kidney_s6, Pathway_col)
  colnames(kidney_s6) <- c(colname_kidney, "Pathway")
  s6_xlsx_df=list(brain=brain_s6,lung=lung_s6,liver=liver_s6,bone=bone_s6,kidney=kidney_s6)
  library(openxlsx)

  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Brain Metastasis")
  writeData(wb, sheet = "Brain Metastasis", s6_xlsx_df$brain)
  addWorksheet(wb, sheetName = "Lung Metastasis")
  writeData(wb, sheet = "Lung Metastasis", s6_xlsx_df$lung)
  addWorksheet(wb, sheetName = "Liver Metastasis")
  writeData(wb, sheet = "Liver Metastasis", s6_xlsx_df$liver)
  addWorksheet(wb, sheetName = "Bone Metastasis")
  writeData(wb, sheet = "Bone Metastasis", s6_xlsx_df$bone)
  addWorksheet(wb, sheetName = "Kidney Metastasis")
  writeData(wb, sheet = "Kidney Metastasis", s6_xlsx_df$kidney)
  #saveWorkbook(wb, "data/add/supplementary_6.xlsx",overwrite = TRUE)
  saveWorkbook(wb, "supplementary/Figure_5_Supplementary_Data_7_Metabolic_Pathway_(Metabolite)_versus_Metabolic_Metapathway (Gene).xlsx",overwrite = TRUE)
  return(invisible(NULL))
}
simplify_data <- function(data, ID,lineage, organ) {
  colname_df <- c("depMapID", "Lineage", "Metastasis", colnames(data))
  simplified_data <- cbind(ID, rep(lineage, length(ID)), rep(organ, length(ID)), data)
  colnames(simplified_data) <- colname_df
  return(simplified_data)
}
make_supplementary5 <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI){
  datasets <- list(
    GI_brain = list(ID = X_GI$X, lineage="GI and GYN",organ = "Brain", data = GI_brain$adj_matrix),
    Lung_brain = list(ID = X_Lung$X, lineage="Lung", organ = "Brain", data = Lung_brain$adj_matrix),
    MSI_brain = list(ID = X_MSI$X, lineage="MSI", organ = "Brain", data = MSI_brain$adj_matrix),
    GI_lung = list(ID = X_GI$X, lineage="GI and GYN", organ = "Lung", data = GI_lung$adj_matrix),
    Lung_lung = list(ID = X_Lung$X, lineage="Lung", organ = "Lung", data = Lung_lung$adj_matrix),
    MSI_lung = list(ID = X_MSI$X, lineage="MSI", organ = "Lung", data = MSI_lung$adj_matrix),
    GI_liver = list(ID = X_GI$X, lineage="GI and GYN", organ = "Liver", data = GI_liver$adj_matrix),
    Lung_liver = list(ID = X_Lung$X, lineage="Lung", organ = "Liver", data = Lung_liver$adj_matrix),
    MSI_liver = list(ID = X_MSI$X, lineage="MSI", organ = "Liver", data = MSI_liver$adj_matrix),
    GI_bone = list(ID = X_GI$X, lineage="GI and GYN", organ = "Bone", data = GI_bone$adj_matrix),
    Lung_bone = list(ID = X_Lung$X, lineage="Lung", organ = "Bone", data = Lung_bone$adj_matrix),
    MSI_bone = list(ID = X_MSI$X, lineage="MSI", organ = "Bone", data = MSI_bone$adj_matrix),
    GI_kidney = list(ID = X_GI$X, lineage="GI and GYN", organ = "Kidney", data = GI_kidney$adj_matrix),
    Lung_kidney = list(ID = X_Lung$X, lineage="Lung", organ = "Kidney", data = Lung_kidney$adj_matrix),
    MSI_kidney = list(ID = X_MSI$X, lineage="MSI", organ = "Kidney", data = MSI_kidney$adj_matrix)
  )

  simplified_data <- lapply(datasets, function(dataset) {
    ID <- dataset$ID
    organ <- dataset$organ
    lineage <- dataset$lineage
    data <- dataset$data
    simplify_data(data, ID, lineage,organ)
  })

  S5_df <- do.call(rbind, simplified_data)
  write.csv(S5_df,"supplementary/Figure_4_Vertical_Dumbbell_Supplementary_Data_5.csv",row.names = F)
  return(invisible(NULL))
}
make_table_3 <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney){
  GI_brain_t3 <- cbind("Brain", "GI and GYN", rownames(data.frame(GI_brain$estimate_value)))
  Lung_brain_t3 <- cbind("Brain", "Lung", rownames(data.frame(Lung_brain$estimate_value)))
  MSI_brain_t3 <- cbind("Brain", "MSI", rownames(data.frame(MSI_brain$estimate_value)))
  GI_lung_t3 <- cbind("Lung", "GI and GYN", rownames(data.frame(GI_lung$estimate_value)))
  Lung_lung_t3 <- cbind("Lung", "Lung", rownames(data.frame(Lung_lung$estimate_value)))
  MSI_lung_t3 <- cbind("Lung", "MSI", rownames(data.frame(MSI_lung$estimate_value)))
  GI_liver_t3 <- cbind("Liver", "GI and GYN", rownames(data.frame(GI_liver$estimate_value)))
  Lung_liver_t3 <- cbind("Liver", "Lung", rownames(data.frame(Lung_liver$estimate_value)))
  MSI_liver_t3 <- cbind("Liver", "MSI", rownames(data.frame(MSI_liver$estimate_value)))
  GI_bone_t3 <- cbind("Bone", "GI and GYN", rownames(data.frame(GI_bone$estimate_value)))
  Lung_bone_t3 <- cbind("Bone", "Lung", rownames(data.frame(Lung_bone$estimate_value)))
  MSI_bone_t3 <- cbind("Bone", "MSI", rownames(data.frame(MSI_bone$estimate_value)))
  GI_kidney_t3 <- cbind("Kidney", "GI and GYN", rownames(data.frame(GI_kidney$estimate_value)))
  Lung_kidney_t3 <- cbind("Kidney", "Lung", rownames(data.frame(Lung_kidney$estimate_value)))
  MSI_kidney_t3 <- cbind("Kidney", "MSI", rownames(data.frame(MSI_kidney$estimate_value)))
  table3_df <- rbind(
    GI_brain_t3, Lung_brain_t3, MSI_brain_t3,
    GI_lung_t3, Lung_lung_t3, MSI_lung_t3,
    GI_liver_t3, Lung_liver_t3, MSI_liver_t3,
    GI_bone_t3, Lung_bone_t3, MSI_bone_t3,
    GI_kidney_t3, Lung_kidney_t3, MSI_kidney_t3
  )
  pathway_t3=c()
  Metabolite=c()
  metabolite_225=read.xlsx("Pathway/metabolitesandMetabolicGenes2.xlsx",sheet = 2)
  metabolite_225_ori_name=metabolite_225$Metabolite
  metabolite_225_name=colnames(GI_brain$adj_matrix)
  for(i in 1:length(table3_df[,1])){
    if(!(table3_df[i,3]%in%metabolite_225_name)){
      pathway_t3=c(pathway_t3,"error")
    }
    for(j in 1:length(metabolite_225_name)){
      if(table3_df[i,3]==metabolite_225_name[j]){
        if(metabolite_225$Pathway[j]=="尼古丁代謝"){
          pathway_t3=c(pathway_t3,"Metabolism of Nicotine")
        }
        else{
          pathway_t3=c(pathway_t3,metabolite_225$Pathway[j])
        }

        Metabolite=c(Metabolite,metabolite_225_ori_name[j])
      }
    }
  }
  table3_df=cbind(table3_df,pathway_t3)
  table3_df=data.frame(table3_df)
  colnames(table3_df)=c("Metastasis","Cancers","Metabolite","Pathway")
  table3_df$Metabolite=Metabolite
  write.csv(table3_df,"supplementary/Table_3.csv",row.names = F)
  return(invisible(NULL))
}
Before_adjusted_matrix <- function(X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,metabolite_name_trans=metabolite_name_trans){
  Brain_before=rbind(data.frame(depMapID=X_GI$X,Lineage=rep(c("GI_GYN"),length(X_GI$X)),Metastasis=rep(c("Brain"),length(X_GI$X)),X_GI[,-1]),data.frame(depMapID=X_Lung$X,Lineage=rep(c("Lung"),length(X_Lung$X)),Metastasis=rep(c("Brain"),length(X_Lung$X)),X_Lung[,-1]),data.frame(depMapID=X_MSI$X,Lineage=rep(c("MSI"),length(X_MSI$X)),Metastasis=rep(c("Brain"),length(X_MSI$X)),X_MSI[,-1]))
  Lung_before=rbind(data.frame(depMapID=X_GI$X,Lineage=rep(c("GI_GYN"),length(X_GI$X)),Metastasis=rep(c("Lung"),length(X_GI$X)),X_GI[,-1]),data.frame(depMapID=X_Lung$X,Lineage=rep(c("Lung"),length(X_Lung$X)),Metastasis=rep(c("Lung"),length(X_Lung$X)),X_Lung[,-1]),data.frame(depMapID=X_MSI$X,Lineage=rep(c("MSI"),length(X_MSI$X)),Metastasis=rep(c("Lung"),length(X_MSI$X)),X_MSI[,-1]))
  Liver_before=rbind(data.frame(depMapID=X_GI$X,Lineage=rep(c("GI_GYN"),length(X_GI$X)),Metastasis=rep(c("Liver"),length(X_GI$X)),X_GI[,-1]),data.frame(depMapID=X_Lung$X,Lineage=rep(c("Lung"),length(X_Lung$X)),Metastasis=rep(c("Liver"),length(X_Lung$X)),X_Lung[,-1]),data.frame(depMapID=X_MSI$X,Lineage=rep(c("MSI"),length(X_MSI$X)),Metastasis=rep(c("Liver"),length(X_MSI$X)),X_MSI[,-1]))
  Bone_before=rbind(data.frame(depMapID=X_GI$X,Lineage=rep(c("GI_GYN"),length(X_GI$X)),Metastasis=rep(c("Bone"),length(X_GI$X)),X_GI[,-1]),data.frame(depMapID=X_Lung$X,Lineage=rep(c("Lung"),length(X_Lung$X)),Metastasis=rep(c("Bone"),length(X_Lung$X)),X_Lung[,-1]),data.frame(depMapID=X_MSI$X,Lineage=rep(c("MSI"),length(X_MSI$X)),Metastasis=rep(c("Bone"),length(X_MSI$X)),X_MSI[,-1]))
  Kidney_before=rbind(data.frame(depMapID=X_GI$X,Lineage=rep(c("GI_GYN"),length(X_GI$X)),Metastasis=rep(c("Kidney"),length(X_GI$X)),X_GI[,-1]),data.frame(depMapID=X_Lung$X,Lineage=rep(c("Lung"),length(X_Lung$X)),Metastasis=rep(c("Kidney"),length(X_Lung$X)),X_Lung[,-1]),data.frame(depMapID=X_MSI$X,Lineage=rep(c("MSI"),length(X_MSI$X)),Metastasis=rep(c("Kidney"),length(X_MSI$X)),X_MSI[,-1]))
  before_adj_matrix=rbind(Brain_before,Lung_before,Liver_before,Bone_before,Kidney_before)
  colnames(before_adj_matrix)=c("depMapID","Lineage","Metastasis",metabolite_name_trans)
  return(before_adj_matrix)
}
Sankey <- function(df_table2_metastasis_4 = df_table2_metastasis_4){
  library(ggsankey)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  df=df_table2_metastasis_4
  df$Cancers <- ifelse(df$Cancers=="Lung", "LUNG", df$Cancers)
  df <- df %>%
    make_long(Cancers,Features,`Metastasis.Sites`)
  node_counts <- table(df$node)
  sorted_nodes <- names(node_counts)[order(-node_counts)]
  sorted_nodes_metabolite=sorted_nodes[sorted_nodes%in%df_table2_metastasis_4$Features]
  sorted_nodes <- c("MSI","LUNG","GI/GYN","Kidney","Bone","Liver","Lung","Brain",rev(sorted_nodes_metabolite))
  df$node <- factor(df$node, levels = sorted_nodes)

  pl <- ggplot(df, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill = factor(node)
                       , label = node)
  )
  pl <- pl +geom_sankey(flow.alpha = 0.5
                        , node.color = "black"
                        ,show.legend = FALSE)
  pl <- pl +geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
  pl <- pl +  theme_bw()
  pl <- pl + theme(legend.position = "none")
  pl <- pl +  theme(axis.title = element_blank()
                    , axis.text.y = element_blank()
                    , axis.ticks = element_blank()
                    , panel.grid = element_blank())
  library(RColorBrewer)
  colors_1 <- brewer.pal(8, "Set2")
  colors_2 <- brewer.pal(12, "Set3")
  colors_3 <- brewer.pal(9, "Set1")
  colors_4 <- brewer.pal(11, "Paired")
  colors <- c(colors_3,colors_2,colors_1,colors_4)
  pl <- pl + scale_fill_manual(values = colors)
  pl <- pl + labs(title = "Sankey diagram using ggplot")
  pl <- pl + labs(subtitle = "using  David Sjoberg's ggsankey package")
  pl <- pl + labs(caption = "@techanswers88")
  pl <- pl + labs(fill = 'Nodes')
  sankey_with_text=pl
  pl2 <- ggplot(df, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill = factor(node)
                       , label = node)
  )
  pl2 <- pl2 +geom_sankey(flow.alpha = 0.5
                        , node.color = "black"
                        ,show.legend = FALSE)
  pl2 <- pl2 +  theme_bw()
  pl2 <- pl2 + theme(legend.position = "none")
  pl2 <- pl2 +  theme(axis.title = element_blank()
                    , axis.text.y = element_blank()
                    , axis.ticks = element_blank()
                    , panel.grid = element_blank())
  library(RColorBrewer)
  pl2 <- pl2 + scale_fill_manual(values = colors)
  sankey_without_text=pl2
  ggsave("supplementary/sankey_with_text.png", plot = sankey_with_text, width = 6, height = 4, dpi = 300)
  ggsave("supplementary/sankey_without_text.png", plot = sankey_without_text, width = 6, height = 4, dpi = 300)
  return(invisible(NULL))
}
make_heatmap <- function(df_table2=df_table2){
  figure_2_brain=df_table2[df_table2$`Metastasis.Sites`=="Brain",]
  metabolite=unique(figure_2_brain$Features)
  cancer=unique(figure_2_brain$Cancers)
  heatmap_matrix=matrix(0,nrow = length(metabolite),ncol=length(cancer))
  colnames(heatmap_matrix)=cancer
  rownames(heatmap_matrix)=metabolite
  for(i in 1:length(figure_2_brain[,1])){
    for(j in 1:length(metabolite)){
      for(k in 1:length(cancer)){
        if((figure_2_brain[i,3]==metabolite[j])&&(figure_2_brain[i,2]==cancer[k])){
          heatmap_matrix[j,k]=figure_2_brain[i,4]
        }
      }
    }
  }
  cols = colorRampPalette(c("blue","white", "red"))(100)
  max_value <- max(abs(heatmap_matrix))
  num_colors <- 100
  color_range <- seq(-max_value, max_value, length.out = num_colors+1)
  my_col <- colorRampPalette(c("blue", "white", "red"))(num_colors+1)
  my_col[51] <- "white"
  myBreaks <- c(seq(min(heatmap_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(heatmap_matrix)/num_colors, max(heatmap_matrix), length.out=floor(num_colors/2)))

  heatmap_brain=pheatmap(heatmap_matrix,cluster_cols = T,color = my_col,breaks = myBreaks)
  ggsave("supplementary/heatmap_brain.png", plot = heatmap_brain, width = 6, height = 4, dpi = 300)
  figure_2_Lung=df_table2[df_table2$`Metastasis.Sites`=="Lung",]
  metabolite=unique(figure_2_Lung$Features)
  cancer=unique(figure_2_Lung$Cancers)
  heatmap_matrix=matrix(0,nrow = length(metabolite),ncol=length(cancer))
  colnames(heatmap_matrix)=cancer
  rownames(heatmap_matrix)=metabolite
  for(i in 1:length(figure_2_Lung[,1])){
    for(j in 1:length(metabolite)){
      for(k in 1:length(cancer)){
        if((figure_2_Lung[i,3]==metabolite[j])&&(figure_2_Lung[i,2]==cancer[k])){
          heatmap_matrix[j,k]=figure_2_Lung[i,4]
        }
      }
    }
  }
  max_value <- max(abs(heatmap_matrix))
  num_colors <- 100
  color_range <- seq(-max_value, max_value, length.out = num_colors+1)
  my_col <- colorRampPalette(c("blue", "white", "red"))(num_colors+1)
  my_col[51] <- "white"
  myBreaks <- c(seq(min(heatmap_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(heatmap_matrix)/num_colors, max(heatmap_matrix), length.out=floor(num_colors/2)))
  heatmap_lung=pheatmap(heatmap_matrix,cluster_cols = T,color = my_col,breaks = myBreaks)
  ggsave("supplementary/heatmap_lung.png", plot = heatmap_lung, width = 6, height = 4, dpi = 300)
  figure_2_Liver=df_table2[df_table2$`Metastasis.Sites`=="Liver",]
  metabolite=unique(figure_2_Liver$Features)
  cancer=unique(figure_2_Liver$Cancers)
  heatmap_matrix=matrix(0,nrow = length(metabolite),ncol=length(cancer))
  colnames(heatmap_matrix)=cancer
  rownames(heatmap_matrix)=metabolite
  for(i in 1:length(figure_2_Liver[,1])){
    for(j in 1:length(metabolite)){
      for(k in 1:length(cancer)){
        if((figure_2_Liver[i,3]==metabolite[j])&&(figure_2_Liver[i,2]==cancer[k])){
          heatmap_matrix[j,k]=figure_2_Liver[i,4]
        }
      }
    }
  }
  max_value <- max(abs(heatmap_matrix))
  num_colors <- 100
  color_range <- seq(-max_value, max_value, length.out = num_colors+1)
  my_col <- colorRampPalette(c("blue", "white", "red"))(num_colors+1)
  my_col[51] <- "white"
  myBreaks <- c(seq(min(heatmap_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(heatmap_matrix)/num_colors, max(heatmap_matrix), length.out=floor(num_colors/2)))
  heatmap_liver=pheatmap(heatmap_matrix,cluster_cols = T,color = my_col,breaks = myBreaks)
  ggsave("supplementary/heatmap_liver.png", plot = heatmap_liver, width = 6, height = 4, dpi = 300)
  figure_2_Bone=df_table2[df_table2$`Metastasis.Sites`=="Bone",]
  metabolite=unique(figure_2_Bone$Features)
  cancer=unique(figure_2_Bone$Cancers)
  heatmap_matrix=matrix(0,nrow = length(metabolite),ncol=length(cancer))
  colnames(heatmap_matrix)=cancer
  rownames(heatmap_matrix)=metabolite
  for(i in 1:length(figure_2_Bone[,1])){
    for(j in 1:length(metabolite)){
      for(k in 1:length(cancer)){
        if((figure_2_Bone[i,3]==metabolite[j])&&(figure_2_Bone[i,2]==cancer[k])){
          heatmap_matrix[j,k]=figure_2_Bone[i,4]
        }
      }
    }
  }
  max_value <- max(abs(heatmap_matrix))
  num_colors <- 100
  color_range <- seq(-max_value, max_value, length.out = num_colors+1)
  my_col <- colorRampPalette(c("blue", "white", "red"))(num_colors+1)
  my_col[51] <- "white"
  myBreaks <- c(seq(min(heatmap_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(heatmap_matrix)/num_colors, max(heatmap_matrix), length.out=floor(num_colors/2)))
  heatmap_bone=pheatmap(heatmap_matrix,cluster_cols = T,color = my_col,breaks = myBreaks)
  ggsave("supplementary/heatmap_bone.png", plot = heatmap_bone, width = 6, height = 4, dpi = 300)
  figure_2_Kidney=df_table2[df_table2$`Metastasis.Sites`=="Kidney",]
  metabolite=unique(figure_2_Kidney$Features)
  cancer=unique(figure_2_Kidney$Cancers)
  heatmap_matrix=matrix(0,nrow = length(metabolite),ncol=length(cancer))
  colnames(heatmap_matrix)=cancer
  rownames(heatmap_matrix)=metabolite
  for(i in 1:length(figure_2_Kidney[,1])){
    for(j in 1:length(metabolite)){
      for(k in 1:length(cancer)){
        if((figure_2_Kidney[i,3]==metabolite[j])&&(figure_2_Kidney[i,2]==cancer[k])){
          heatmap_matrix[j,k]=figure_2_Kidney[i,4]
        }
      }
    }
  }
  max_value <- max(abs(heatmap_matrix))
  num_colors <- 100
  color_range <- seq(-max_value, max_value, length.out = num_colors+1)
  my_col <- colorRampPalette(c("blue", "white", "red"))(num_colors+1)
  my_col[51] <- "white"
  myBreaks <- c(seq(min(heatmap_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(heatmap_matrix)/num_colors, max(heatmap_matrix), length.out=floor(num_colors/2)))
  heatmap_kidney=pheatmap(heatmap_matrix,cluster_cols = T,color = my_col,breaks = myBreaks)
  ggsave("supplementary/heatmap_kidney.png", plot = heatmap_kidney, width = 6, height = 4, dpi = 300)
  return(invisible(NULL))
}
UpSet_plot <- function(df_table2){
  library(UpSetR)
  library(grid)
  library(gridSVG)
  library(ggplotify)
  figure_2_brain=df_table2[df_table2$`Metastasis.Sites`=="Brain",]
  figure_2_Lung=df_table2[df_table2$`Metastasis.Sites`=="Lung",]
  figure_2_Liver=df_table2[df_table2$`Metastasis.Sites`=="Liver",]
  figure_2_Bone=df_table2[df_table2$`Metastasis.Sites`=="Bone",]
  figure_2_Kidney=df_table2[df_table2$`Metastasis.Sites`=="Kidney",]
  MSI_brain_feature=figure_2_brain[figure_2_brain$Cancers=="MSI",3]
  GI_GYN_brain_feature=figure_2_brain[figure_2_brain$Cancers=="GI/GYN",3]
  Lung_brain_feature=figure_2_brain[figure_2_brain$Cancers=="Lung",3]
  listInput<-list(MSI=MSI_brain_feature,GI_GYN=GI_GYN_brain_feature,Lung=Lung_brain_feature)
  upset_brain=upset(fromList(listInput), keep.order = T,sets = c("MSI","Lung","GI_GYN"),text.scale = 1.8)
  gg_pl <- as.ggplot(upset_brain)
  ggsave("supplementary/upset_brain.png", plot = gg_pl, width = 6, height = 4, dpi = 300)
  brain_importance_list=c(GI_GYN_brain_feature,Lung_brain_feature,MSI_brain_feature)
  brain_importance_list=unique(brain_importance_list)
  GI_GYN_brain_upset=c()
  Lung_brain_upset=c()
  MSI_brain_upset=c()
  for(i in 1:length(brain_importance_list)){
    if(brain_importance_list[i]%in%GI_GYN_brain_feature){
      GI_GYN_brain_upset=c(GI_GYN_brain_upset,1)
    }
    else{
      GI_GYN_brain_upset=c(GI_GYN_brain_upset,0)
    }
    if(brain_importance_list[i]%in%Lung_brain_feature){
      Lung_brain_upset=c(Lung_brain_upset,1)
    }
    else{
      Lung_brain_upset=c(Lung_brain_upset,0)
    }
    if(brain_importance_list[i]%in%MSI_brain_feature){
      MSI_brain_upset=c(MSI_brain_upset,1)
    }
    else{
      MSI_brain_upset=c(MSI_brain_upset,0)
    }
  }
  brain_importance_df=cbind(GI_GYN_brain_upset,Lung_brain_upset,MSI_brain_upset)
  colnames(brain_importance_df)=c("GI_GYN","Lung","MSI")
  rownames(brain_importance_df)=brain_importance_list
  MSI_Lung_feature=figure_2_Lung[figure_2_Lung$Cancers=="MSI",3]
  GI_GYN_Lung_feature=figure_2_Lung[figure_2_Lung$Cancers=="GI/GYN",3]
  Lung_Lung_feature=figure_2_Lung[figure_2_Lung$Cancers=="Lung",3]
  listInput<-list(MSI=MSI_Lung_feature,GI_GYN=GI_GYN_Lung_feature,Lung=Lung_Lung_feature)
  upset_lung=upset(fromList(listInput), keep.order = T,sets = c("MSI","Lung","GI_GYN"),text.scale = 1.8)
  gg_pl <- as.ggplot(upset_lung)
  ggsave("supplementary/upset_lung.png", plot = gg_pl, width = 6, height = 4, dpi = 300)
  Lung_importance_list=c(GI_GYN_Lung_feature,Lung_Lung_feature,MSI_Lung_feature)
  Lung_importance_list=unique(Lung_importance_list)
  GI_GYN_Lung_upset=c()
  Lung_Lung_upset=c()
  MSI_Lung_upset=c()
  for(i in 1:length(Lung_importance_list)){
    if(Lung_importance_list[i]%in%GI_GYN_Lung_feature){
      GI_GYN_Lung_upset=c(GI_GYN_Lung_upset,1)
    }
    else{
      GI_GYN_Lung_upset=c(GI_GYN_Lung_upset,0)
    }
    if(Lung_importance_list[i]%in%Lung_Lung_feature){
      Lung_Lung_upset=c(Lung_Lung_upset,1)
    }
    else{
      Lung_Lung_upset=c(Lung_Lung_upset,0)
    }
    if(Lung_importance_list[i]%in%MSI_Lung_feature){
      MSI_Lung_upset=c(MSI_Lung_upset,1)
    }
    else{
      MSI_Lung_upset=c(MSI_Lung_upset,0)
    }
  }
  Lung_importance_df=cbind(GI_GYN_Lung_upset,Lung_Lung_upset,MSI_Lung_upset)
  colnames(Lung_importance_df)=c("GI_GYN","Lung","MSI")
  rownames(Lung_importance_df)=Lung_importance_list
  MSI_Liver_feature=figure_2_Liver[figure_2_Liver$Cancers=="MSI",3]
  GI_GYN_Liver_feature=figure_2_Liver[figure_2_Liver$Cancers=="GI/GYN",3]
  Lung_Liver_feature=figure_2_Liver[figure_2_Liver$Cancers=="Lung",3]
  listInput<-list(MSI=MSI_Liver_feature,GI_GYN=GI_GYN_Liver_feature,Lung=Lung_Liver_feature)
  upset_liver=upset(fromList(listInput), keep.order = T,sets = c("MSI","Lung","GI_GYN"),text.scale = 1.8)
  gg_pl <- as.ggplot(upset_liver)
  ggsave("supplementary/upset_liver.png", plot = gg_pl, width = 6, height = 4, dpi = 300)
  Liver_importance_list=c(GI_GYN_Liver_feature,Lung_Liver_feature,MSI_Liver_feature)
  Liver_importance_list=unique(Liver_importance_list)
  GI_GYN_Liver_upset=c()
  Lung_Liver_upset=c()
  MSI_Liver_upset=c()
  for(i in 1:length(Liver_importance_list)){
    if(Liver_importance_list[i]%in%GI_GYN_Liver_feature){
      GI_GYN_Liver_upset=c(GI_GYN_Liver_upset,1)
    }
    else{
      GI_GYN_Liver_upset=c(GI_GYN_Liver_upset,0)
    }
    if(Liver_importance_list[i]%in%Lung_Liver_feature){
      Lung_Liver_upset=c(Lung_Liver_upset,1)
    }
    else{
      Lung_Liver_upset=c(Lung_Liver_upset,0)
    }
    if(Liver_importance_list[i]%in%MSI_Liver_feature){
      MSI_Liver_upset=c(MSI_Liver_upset,1)
    }
    else{
      MSI_Liver_upset=c(MSI_Liver_upset,0)
    }
  }
  Liver_importance_df=cbind(GI_GYN_Liver_upset,Lung_Liver_upset,MSI_Liver_upset)
  colnames(Liver_importance_df)=c("GI_GYN","Lung","MSI")
  rownames(Liver_importance_df)=Liver_importance_list
  MSI_Bone_feature=figure_2_Bone[figure_2_Bone$Cancers=="MSI",3]
  GI_GYN_Bone_feature=figure_2_Bone[figure_2_Bone$Cancers=="GI/GYN",3]
  Lung_Bone_feature=figure_2_Bone[figure_2_Bone$Cancers=="Lung",3]
  listInput<-list(MSI=MSI_Bone_feature,GI_GYN=GI_GYN_Bone_feature,Lung=Lung_Bone_feature)
  upset_bone=upset(fromList(listInput), keep.order = T,sets = c("MSI","Lung","GI_GYN"),text.scale = 1.8)
  gg_pl <- as.ggplot(upset_bone)
  ggsave("supplementary/upset_bone.png", plot = gg_pl, width = 6, height = 4, dpi = 300)
  Bone_importance_list=c(GI_GYN_Bone_feature,Lung_Bone_feature,MSI_Bone_feature)
  Bone_importance_list=unique(Bone_importance_list)
  GI_GYN_Bone_upset=c()
  Lung_Bone_upset=c()
  MSI_Bone_upset=c()
  for(i in 1:length(Bone_importance_list)){
    if(Bone_importance_list[i]%in%GI_GYN_Bone_feature){
      GI_GYN_Bone_upset=c(GI_GYN_Bone_upset,1)
    }
    else{
      GI_GYN_Bone_upset=c(GI_GYN_Bone_upset,0)
    }
    if(Bone_importance_list[i]%in%Lung_Bone_feature){
      Lung_Bone_upset=c(Lung_Bone_upset,1)
    }
    else{
      Lung_Bone_upset=c(Lung_Bone_upset,0)
    }
    if(Bone_importance_list[i]%in%MSI_Bone_feature){
      MSI_Bone_upset=c(MSI_Bone_upset,1)
    }
    else{
      MSI_Bone_upset=c(MSI_Bone_upset,0)
    }
  }
  Bone_importance_df=cbind(GI_GYN_Bone_upset,Lung_Bone_upset,MSI_Bone_upset)
  colnames(Bone_importance_df)=c("GI_GYN","Lung","MSI")
  rownames(Bone_importance_df)=Bone_importance_list
  MSI_Kidney_feature=figure_2_Kidney[figure_2_Kidney$Cancers=="MSI",3]
  GI_GYN_Kidney_feature=figure_2_Kidney[figure_2_Kidney$Cancers=="GI/GYN",3]
  Lung_Kidney_feature=figure_2_Kidney[figure_2_Kidney$Cancers=="Lung",3]
  listInput<-list(MSI=MSI_Kidney_feature,GI_GYN=GI_GYN_Kidney_feature,Lung=Lung_Kidney_feature)
  upset_kidney=upset(fromList(listInput), keep.order = T,sets = c("MSI","Lung","GI_GYN"),text.scale = 1.8)
  gg_pl <- as.ggplot(upset_kidney)
  ggsave("supplementary/upset_kidney.png", plot = gg_pl, width = 6, height = 4, dpi = 300)
  Kidney_importance_list=c(GI_GYN_Kidney_feature,Lung_Kidney_feature,MSI_Kidney_feature)
  Kidney_importance_list=unique(Kidney_importance_list)
  GI_GYN_Kidney_upset=c()
  Lung_Kidney_upset=c()
  MSI_Kidney_upset=c()
  for(i in 1:length(Kidney_importance_list)){
    if(Kidney_importance_list[i]%in%GI_GYN_Kidney_feature){
      GI_GYN_Kidney_upset=c(GI_GYN_Kidney_upset,1)
    }
    else{
      GI_GYN_Kidney_upset=c(GI_GYN_Kidney_upset,0)
    }
    if(Kidney_importance_list[i]%in%Lung_Kidney_feature){
      Lung_Kidney_upset=c(Lung_Kidney_upset,1)
    }
    else{
      Lung_Kidney_upset=c(Lung_Kidney_upset,0)
    }
    if(Kidney_importance_list[i]%in%MSI_Kidney_feature){
      MSI_Kidney_upset=c(MSI_Kidney_upset,1)
    }
    else{
      MSI_Kidney_upset=c(MSI_Kidney_upset,0)
    }
  }
  Kidney_importance_df=cbind(GI_GYN_Kidney_upset,Lung_Kidney_upset,MSI_Kidney_upset)
  colnames(Kidney_importance_df)=c("GI_GYN","Lung","MSI")
  rownames(Kidney_importance_df)=Kidney_importance_list
  return(invisible(NULL))
}
make_heatmap2 <- function(metabolite_225=metabolite_225,meta_225_col=meta_225_col,metabolite_name_trans=metabolite_name_trans,supplementary_7=supplementary_7){
  brain_s6 = read.xlsx(supplementary_7,sheet = "Brain Metastasis")
  lung_s6 = read.xlsx(supplementary_7,sheet = "Lung Metastasis")
  liver_s6 = read.xlsx(supplementary_7,sheet = "Liver Metastasis")
  bone_s6 = read.xlsx(supplementary_7,sheet = "Bone Metastasis")
  kidney_s6 = read.xlsx(supplementary_7,sheet = "Kidney Metastasis")
  figure_5_brain=data.frame(brain_s6)
  pathway_df=cbind(metabolite_225,meta_225_col$Pathway)
  pathway_df=data.frame(pathway_df)
  colnames(pathway_df)=c("metabolite","pathway")
  f5_meta=colnames(figure_5_brain)[-c(1,2)]
  f5_meta= f5_meta[f5_meta != "Pathway"]
  pathway_df <- pathway_df[order(pathway_df$pathway), ]
  brain_pathway=pathway_df[pathway_df$metabolite%in%f5_meta,]
  sort_brain=c("Cancers","Gene",brain_pathway$metabolite,"Pathway")
  figure_5_brain=figure_5_brain[,sort_brain]
  figure_5_brain=figure_5_brain[order(figure_5_brain$Pathway),]
  row_name_figure_5_brain=figure_5_brain$Gene
  matrix_figure_5_brain=figure_5_brain[,3:39]
  col_name_figure_5_brain=colnames(figure_5_brain[,3:39])
  colnames(matrix_figure_5_brain)=col_name_figure_5_brain
  t_matrix_figure_5_brain=t(matrix_figure_5_brain)
  colnames(t_matrix_figure_5_brain)=row_name_figure_5_brain
  rownames(t_matrix_figure_5_brain)=col_name_figure_5_brain
  write.csv(t_matrix_figure_5_brain,"f5_brain.csv",row.names = T)
  f5_brain_matrix=read.csv("f5_brain.csv",row.names = 1)
  file.remove('f5_brain.csv')
  brain_annotation=data.frame(brain_pathway$pathway,row.names = brain_pathway$metabolite)
  colnames(brain_annotation)=c("Metabolite_pathway")
  brain_annotation_col=data.frame(figure_5_brain$Pathway,row.names = colnames(f5_brain_matrix))
  colnames(brain_annotation_col)=c("Gene_pathway")
  num_colors <- 100
  cols = colorRampPalette(c("blue","white", "red"))(100)
  myBreaks <- c(seq(min(f5_brain_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(f5_brain_matrix)/num_colors, max(f5_brain_matrix), length.out=floor(num_colors/2)))
  annotation_colors <- list(
    Metabolite_pathway = setNames(brewer.pal(length(unique(brain_annotation$Metabolite_pathway)), "Set1"),
                                  unique(brain_annotation$Metabolite_pathway)),
    Gene_pathway = setNames(brewer.pal(length(unique(brain_annotation_col$Gene_pathway)), "Set2"),
                            unique(brain_annotation_col$Gene_pathway))
  )
  levels(brain_annotation$Metabolite_pathway) <- unique(brain_annotation$Metabolite_pathway)
  levels(brain_annotation_col$Gene_pathway) <- unique(brain_annotation_col$Gene_pathway)
  temp_name=c()
  for(i in rownames(f5_brain_matrix)){
    for(j in 1:length(metabolite_name_trans)){
      if(i==meta_225_col$Metabolite_trans[j]){
        temp_name=c(temp_name,metabolite_name_trans[j])
      }
    }
  }
  rownames(f5_brain_matrix)=temp_name
  brain_heatmap=pheatmap(f5_brain_matrix,cluster_rows = F,cluster_cols = F,color = cols,annotation_row = brain_annotation,annotation_col = brain_annotation_col,breaks = myBreaks,annotation_colors = annotation_colors)
  ggsave("supplementary/brain_heatmap.png", plot = brain_heatmap, width = 20, height = 12, dpi = 300)
  figure_5_Lung=data.frame(lung_s6)
  pathway_df=cbind(metabolite_225,meta_225_col$Pathway)
  pathway_df=data.frame(pathway_df)
  colnames(pathway_df)=c("metabolite","pathway")
  f5_meta=colnames(figure_5_Lung)[-c(1,2)]
  f5_meta= f5_meta[f5_meta != "Pathway"]
  pathway_df <- pathway_df[order(pathway_df$pathway), ]
  Lung_pathway=pathway_df[pathway_df$metabolite%in%f5_meta,]
  sort_Lung=c("Cancers","Gene",Lung_pathway$metabolite,"Pathway")
  figure_5_Lung=figure_5_Lung[,sort_Lung]
  figure_5_Lung=figure_5_Lung[order(figure_5_Lung$Pathway),]
  row_name_figure_5_Lung=figure_5_Lung$Gene
  matrix_figure_5_Lung=figure_5_Lung[,3:24]
  col_name_figure_5_Lung=colnames(figure_5_Lung[,3:24])
  colnames(matrix_figure_5_Lung)=col_name_figure_5_Lung
  t_matrix_figure_5_Lung=t(matrix_figure_5_Lung)
  colnames(t_matrix_figure_5_Lung)=row_name_figure_5_Lung
  rownames(t_matrix_figure_5_Lung)=col_name_figure_5_Lung
  write.csv(t_matrix_figure_5_Lung,"f5_Lung.csv",row.names = T)
  f5_Lung_matrix=read.csv("f5_Lung.csv",row.names = 1)
  file.remove('f5_Lung.csv')
  Lung_annotation=data.frame(Lung_pathway$pathway,row.names = Lung_pathway$metabolite)
  colnames(Lung_annotation)=c("Metabolite_pathway")
  Lung_annotation_col=data.frame(figure_5_Lung$Pathway,row.names = colnames(f5_Lung_matrix))
  colnames(Lung_annotation_col)=c("Gene_pathway")
  myBreaks <- c(seq(min(f5_Lung_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(f5_Lung_matrix)/num_colors, max(f5_Lung_matrix), length.out=floor(num_colors/2)))
  annotation_colors <- list(
    Metabolite_pathway = setNames(brewer.pal(length(unique(Lung_annotation$Metabolite_pathway)), "Set1"),
                                  unique(Lung_annotation$Metabolite_pathway)),
    Gene_pathway = setNames(brewer.pal(length(unique(Lung_annotation_col$Gene_pathway)), "Set2"),
                            unique(Lung_annotation_col$Gene_pathway))
  )
  levels(Lung_annotation$Metabolite_pathway) <- unique(Lung_annotation$Metabolite_pathway)
  levels(Lung_annotation_col$Gene_pathway) <- unique(Lung_annotation_col$Gene_pathway)
  temp_name=c()
  for(i in rownames(f5_Lung_matrix)){
    for(j in 1:length(metabolite_name_trans)){
      if(i==meta_225_col$Metabolite_trans[j]){
        temp_name=c(temp_name,metabolite_name_trans[j])
      }
    }
  }
  rownames(f5_Lung_matrix)=temp_name
  lung_heatmap=pheatmap(f5_Lung_matrix,cluster_rows = F,cluster_cols = F,color = cols,annotation_row = Lung_annotation,annotation_col = Lung_annotation_col,breaks = myBreaks,annotation_colors = annotation_colors)
  ggsave("supplementary/lung_heatmap.png", plot = lung_heatmap, width = 20, height = 12, dpi = 300)
  figure_5_Liver=data.frame(liver_s6)
  pathway_df=cbind(metabolite_225,meta_225_col$Pathway)
  pathway_df=data.frame(pathway_df)
  colnames(pathway_df)=c("metabolite","pathway")
  f5_meta=colnames(figure_5_Liver)[-c(1,2)]
  f5_meta= f5_meta[f5_meta != "Pathway"]
  pathway_df <- pathway_df[order(pathway_df$pathway), ]
  Liver_pathway=pathway_df[pathway_df$metabolite%in%f5_meta,]
  sort_Liver=c("Cancers","Gene",Liver_pathway$metabolite,"Pathway")
  figure_5_Liver=figure_5_Liver[,sort_Liver]
  figure_5_Liver=figure_5_Liver[order(figure_5_Liver$Pathway),]
  row_name_figure_5_Liver=figure_5_Liver$Gene
  matrix_figure_5_Liver=figure_5_Liver[,3:29]
  col_name_figure_5_Liver=colnames(figure_5_Liver[,3:29])
  colnames(matrix_figure_5_Liver)=col_name_figure_5_Liver
  t_matrix_figure_5_Liver=t(matrix_figure_5_Liver)
  colnames(t_matrix_figure_5_Liver)=row_name_figure_5_Liver
  rownames(t_matrix_figure_5_Liver)=col_name_figure_5_Liver
  write.csv(t_matrix_figure_5_Liver,"f5_Liver.csv",row.names = T)
  f5_Liver_matrix=read.csv("f5_Liver.csv",row.names = 1)
  file.remove('f5_Liver.csv')
  Liver_annotation=data.frame(Liver_pathway$pathway,row.names = Liver_pathway$metabolite)
  colnames(Liver_annotation)=c("Metabolite_pathway")
  Liver_annotation_col=data.frame(figure_5_Liver$Pathway,row.names = colnames(f5_Liver_matrix))
  colnames(Liver_annotation_col)=c("Gene_pathway")
  myBreaks <- c(seq(min(f5_Liver_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(f5_Liver_matrix)/num_colors, max(f5_Liver_matrix), length.out=floor(num_colors/2)))
  annotation_colors <- list(
    Metabolite_pathway = setNames(brewer.pal(length(unique(Liver_annotation$Metabolite_pathway)), "Set1"),
                                  unique(Liver_annotation$Metabolite_pathway)),
    Gene_pathway = setNames(brewer.pal(length(unique(Liver_annotation_col$Gene_pathway)), "Set2"),
                            unique(Liver_annotation_col$Gene_pathway))
  )
  levels(Liver_annotation$Metabolite_pathway) <- unique(Liver_annotation$Metabolite_pathway)
  levels(Liver_annotation_col$Gene_pathway) <- unique(Liver_annotation_col$Gene_pathway)
  temp_name=c()
  for(i in rownames(f5_Liver_matrix)){
    for(j in 1:length(metabolite_name_trans)){
      if(i==meta_225_col$Metabolite_trans[j]){
        temp_name=c(temp_name,metabolite_name_trans[j])
      }
    }
  }
  rownames(f5_Liver_matrix)=temp_name
  liver_heatmap=pheatmap(f5_Liver_matrix,cluster_rows = F,cluster_cols = F,color = cols,annotation_row = Liver_annotation,annotation_col = Liver_annotation_col,breaks = myBreaks,annotation_colors = annotation_colors)
  ggsave("supplementary/liver_heatmap.png", plot = liver_heatmap, width = 20, height = 12, dpi = 300)
  figure_5_Bone=data.frame(bone_s6)
  pathway_df=cbind(metabolite_225,meta_225_col$Pathway)
  pathway_df=data.frame(pathway_df)
  colnames(pathway_df)=c("metabolite","pathway")
  f5_meta=colnames(figure_5_Bone)[-c(1,2)]
  f5_meta= f5_meta[f5_meta != "Pathway"]
  pathway_df <- pathway_df[order(pathway_df$pathway), ]
  Bone_pathway=pathway_df[pathway_df$metabolite%in%f5_meta,]
  sort_Bone=c("Cancers","Gene",Bone_pathway$metabolite,"Pathway")
  figure_5_Bone=figure_5_Bone[,sort_Bone]
  figure_5_Bone=figure_5_Bone[order(figure_5_Bone$Pathway),]
  row_name_figure_5_Bone=figure_5_Bone$Gene
  matrix_figure_5_Bone=figure_5_Bone[,3:27]
  col_name_figure_5_Bone=colnames(figure_5_Bone[,3:27])
  colnames(matrix_figure_5_Bone)=col_name_figure_5_Bone
  t_matrix_figure_5_Bone=t(matrix_figure_5_Bone)
  colnames(t_matrix_figure_5_Bone)=row_name_figure_5_Bone
  rownames(t_matrix_figure_5_Bone)=col_name_figure_5_Bone
  write.csv(t_matrix_figure_5_Bone,"f5_Bone.csv",row.names = T)
  f5_Bone_matrix=read.csv("f5_Bone.csv",row.names = 1)
  file.remove('f5_Bone.csv')
  Bone_annotation=data.frame(Bone_pathway$pathway,row.names = Bone_pathway$metabolite)
  colnames(Bone_annotation)=c("Metabolite_pathway")
  Bone_annotation_col=data.frame(figure_5_Bone$Pathway,row.names = colnames(f5_Bone_matrix))
  colnames(Bone_annotation_col)=c("Gene_pathway")
  myBreaks <- c(seq(min(f5_Bone_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(f5_Bone_matrix)/num_colors, max(f5_Bone_matrix), length.out=floor(num_colors/2)))
  annotation_colors <- list(
    Metabolite_pathway = setNames(brewer.pal(length(unique(Bone_annotation$Metabolite_pathway)), "Set1"),
                                  unique(Bone_annotation$Metabolite_pathway)),
    Gene_pathway = setNames(brewer.pal(length(unique(Bone_annotation_col$Gene_pathway)), "Set2"),
                            unique(Bone_annotation_col$Gene_pathway))
  )
  levels(Bone_annotation$Metabolite_pathway) <- unique(Bone_annotation$Metabolite_pathway)
  levels(Bone_annotation_col$Gene_pathway) <- unique(Bone_annotation_col$Gene_pathway)
  temp_name=c()
  for(i in rownames(f5_Bone_matrix)){
    for(j in 1:length(metabolite_name_trans)){
      if(i==meta_225_col$Metabolite_trans[j]){
        temp_name=c(temp_name,metabolite_name_trans[j])
      }
    }
  }
  rownames(f5_Bone_matrix)=temp_name
  bone_heatmap=pheatmap(f5_Bone_matrix,cluster_rows = F,cluster_cols = F,color = cols,annotation_row = Bone_annotation,annotation_col = Bone_annotation_col,breaks = myBreaks,annotation_colors = annotation_colors)
  ggsave("supplementary/bone_heatmap.png", plot = bone_heatmap, width = 20, height = 12, dpi = 300)
  figure_5_Kidney=data.frame(kidney_s6)
  pathway_df=cbind(metabolite_225,meta_225_col$Pathway)
  pathway_df=data.frame(pathway_df)
  colnames(pathway_df)=c("metabolite","pathway")
  f5_meta=colnames(figure_5_Kidney)[-c(1,2)]
  f5_meta= f5_meta[f5_meta != "Pathway"]
  pathway_df <- pathway_df[order(pathway_df$pathway), ]
  Kidney_pathway=pathway_df[pathway_df$metabolite%in%f5_meta,]
  sort_Kidney=c("Cancers","Gene",Kidney_pathway$metabolite,"Pathway")
  figure_5_Kidney=figure_5_Kidney[,sort_Kidney]
  figure_5_Kidney=figure_5_Kidney[order(figure_5_Kidney$Pathway),]
  row_name_figure_5_Kidney=figure_5_Kidney$Gene
  matrix_figure_5_Kidney=figure_5_Kidney[,3:27]
  col_name_figure_5_Kidney=colnames(figure_5_Kidney[,3:27])
  colnames(matrix_figure_5_Kidney)=col_name_figure_5_Kidney
  t_matrix_figure_5_Kidney=t(matrix_figure_5_Kidney)
  colnames(t_matrix_figure_5_Kidney)=row_name_figure_5_Kidney
  rownames(t_matrix_figure_5_Kidney)=col_name_figure_5_Kidney
  write.csv(t_matrix_figure_5_Kidney,"f5_Kidney.csv",row.names = T)
  f5_Kidney_matrix=read.csv("f5_Kidney.csv",row.names = 1)
  file.remove('f5_Kidney.csv')
  Kidney_annotation=data.frame(Kidney_pathway$pathway,row.names = Kidney_pathway$metabolite)
  colnames(Kidney_annotation)=c("Metabolite_pathway")
  Kidney_annotation_col=data.frame(figure_5_Kidney$Pathway,row.names = colnames(f5_Kidney_matrix))
  colnames(Kidney_annotation_col)=c("Gene_pathway")
  myBreaks <- c(seq(min(f5_Kidney_matrix), 0, length.out=ceiling(num_colors/2) + 1),
                seq(max(f5_Kidney_matrix)/num_colors, max(f5_Kidney_matrix), length.out=floor(num_colors/2)))
  annotation_colors <- list(
    Metabolite_pathway = setNames(brewer.pal(length(unique(Kidney_annotation$Metabolite_pathway)), "Set1"),
                                  unique(Kidney_annotation$Metabolite_pathway)),
    Gene_pathway = setNames(brewer.pal(length(unique(Kidney_annotation_col$Gene_pathway)), "Set2"),
                            unique(Kidney_annotation_col$Gene_pathway))
  )
  levels(Kidney_annotation$Metabolite_pathway) <- unique(Kidney_annotation$Metabolite_pathway)
  levels(Kidney_annotation_col$Gene_pathway) <- unique(Kidney_annotation_col$Gene_pathway)
  temp_name=c()
  for(i in rownames(f5_Kidney_matrix)){
    for(j in 1:length(metabolite_name_trans)){
      if(i==meta_225_col$Metabolite_trans[j]){
        temp_name=c(temp_name,metabolite_name_trans[j])
      }
    }
  }
  rownames(f5_Kidney_matrix)=temp_name
  kidney_heatmap=pheatmap(f5_Kidney_matrix,cluster_rows = F,cluster_cols = F,color = cols,annotation_row = Kidney_annotation,annotation_col = Kidney_annotation_col,breaks = myBreaks,annotation_colors = annotation_colors)
  ggsave("supplementary/kidney_heatmap.png", plot = kidney_heatmap, width = 20, height = 12, dpi = 300)
  return(invisible(NULL))
}
dumbbell_overall <- function(before_adj_matrix=before_adj_matrix,S5_df=S5_df){
  library(ggplot2)
  before_adj_matrix=before_adj_matrix
  after_adj_matrix=S5_df
  before_adj_matrix_numeric=before_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=after_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=data.frame(after_adj_matrix_numeric)
  after_adj_matrix_numeric <- as.data.frame(sapply(after_adj_matrix_numeric, function(x) as.numeric(as.character(x))))
  mean_for_before=colMeans(before_adj_matrix_numeric)
  mean_for_after=colMeans(after_adj_matrix_numeric)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_df=cbind(colnames(before_adj_matrix_numeric),dumbbell_df)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(before_adj_matrix_numeric)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)
  dumbbell_overall_with_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), linewidth = 1, linetype = "dashed") +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270)+
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(2, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  dumbbell_overall_without_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), linewidth = 1, linetype = "dashed") +
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(2, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  ggsave("supplementary/dumbbell_overall_with_text.png", plot = dumbbell_overall_with_text, width = 12, height = 10, dpi = 300)
  ggsave("supplementary/dumbbell_overall_without_text.png", plot = dumbbell_overall_without_text, width = 12, height = 10, dpi = 300)
  return(invisible(NULL))
}
dumbbell_with_sankey_significant <- function(before_adj_matrix=before_adj_matrix,S5_df=S5_df,df_table2_metastasis_4 = df_table2_metastasis_4){
  before_adj_matrix=before_adj_matrix
  after_adj_matrix=S5_df
  before_adj_matrix_numeric=before_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=after_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=data.frame(after_adj_matrix_numeric)
  after_adj_matrix_numeric <- as.data.frame(sapply(after_adj_matrix_numeric, function(x) as.numeric(as.character(x))))
  mean_for_before=colMeans(before_adj_matrix_numeric)
  mean_for_after=colMeans(after_adj_matrix_numeric)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_color=rep("Non-Significant",length(mean_for_after))
  dumbbell_df=cbind(colnames(before_adj_matrix_numeric),dumbbell_df,dumbbell_color)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after","color")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(before_adj_matrix_numeric)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  dumbbell_df[dumbbell_df$Metabolite%in%df_table2_metastasis_4$Features,]$color="Significant"
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_overall_with_text = ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = color, shape = "Before"), size = 2) +
    geom_point(aes(x = order, y = Mean_for_after, color = color, shape = "After"), size = 2) +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270) +
    labs(x = NULL, y = "Value", color = "Significant") +
    scale_color_manual(values = c("Non-Significant" = "grey", "Significant" = "red")) +
    scale_shape_manual(name = "Adjustment", labels = c("After", "Before"), values = c("Before" = 1, "After" = 3)) +
    theme_bw() +
    ylim(2, 10) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(1, 1, 3, 1), "lines"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.5, "lines")
    ) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  ggsave("supplementary/dumbbell_overall_with_text.png", plot = dumbbell_overall_with_text, width = 12, height = 10, dpi = 300)
  before_adj_matrix=before_adj_matrix
  after_adj_matrix=S5_df
  before_adj_matrix_numeric=before_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=after_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=data.frame(after_adj_matrix_numeric)
  after_adj_matrix_numeric <- as.data.frame(sapply(after_adj_matrix_numeric, function(x) as.numeric(as.character(x))))
  mean_for_before=colMeans(before_adj_matrix_numeric)
  mean_for_after=colMeans(after_adj_matrix_numeric)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_color=rep("Non-Significant",length(mean_for_after))
  dumbbell_df=cbind(colnames(before_adj_matrix_numeric),dumbbell_df,dumbbell_color)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after","color")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(before_adj_matrix_numeric)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  dumbbell_df[dumbbell_df$Metabolite%in%df_table2_metastasis_4$Features,]$color="Significant"
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_overall_without_text = ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = color, shape = "Before"), size = 2) +
    geom_point(aes(x = order, y = Mean_for_after, color = color, shape = "After"), size = 2) +
    labs(x = NULL, y = "Value", color = "Significant") +
    scale_color_manual(values = c("Non-Significant" = "grey", "Significant" = "red")) +
    scale_shape_manual(name = "Adjustment", labels = c("After", "Before"), values = c("Before" = 1, "After" = 3)) +
    theme_bw() +
    ylim(2, 10) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(1, 1, 3, 1), "lines"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.5, "lines")
    ) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  write.csv(dumbbell_df,"supplementary/dumbbell_sorting.csv")
  ggsave("supplementary/dumbbell_overall_without_text.png", plot = dumbbell_overall_without_text, width = 12, height = 10, dpi = 300)
  return(invisible(NULL))
}
PC_plot <- function(before_adj_matrix=before_adj_matrix,S5_df=S5_df){
  library(ggplot2)
  library(dplyr)
  library(factoextra)
  before_adj_matrix_numeric=before_adj_matrix[,-c(1,2,3)]
  set.seed(3)
  pca <- prcomp(before_adj_matrix_numeric, scale. = TRUE)

  pca_data <- as.data.frame(pca$x[, 1:2])

  kmeans_cluster <- kmeans(pca_data, centers = 3)

  pca_data$cluster <- as.factor(kmeans_cluster$cluster)
  pca_data_before=pca_data
  pca_before=ggplot(pca_data_before, aes(PC1, PC2, color = cluster)) +
    geom_point(size=0.5) +
    labs(title = "PCA Clustering") +
    theme_bw()+
    theme(legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size = 14))
  ggsave("supplementary/pca_before.png", plot = pca_before, width = 8, height = 6, dpi = 300)
  set.seed(3)
  after_adj_matrix=S5_df
  after_adj_matrix_numeric=after_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=data.frame(after_adj_matrix_numeric)
  after_adj_matrix_numeric <- as.data.frame(sapply(after_adj_matrix_numeric, function(x) as.numeric(as.character(x))))
  pca <- prcomp(after_adj_matrix_numeric, scale. = TRUE)

  pca_data <- as.data.frame(pca$x[, 1:2])

  kmeans_cluster <- kmeans(pca_data, centers = 3)

  pca_data$cluster <- as.factor(kmeans_cluster$cluster)
  pca_data_after=pca_data
  levels(pca_data_after$cluster) <- c("MSI", "GI and GYN", "Lung")
  after_adj_matrix=S5_df
  after_adj_matrix=data.frame(after_adj_matrix)
  pca_data_after$Metastasis=after_adj_matrix$Metastasis
  pca_after_MSI=ggplot(pca_data_after[pca_data_after$cluster=="MSI",], aes(PC1, PC2, color = cluster,shape = Metastasis)) +
    geom_point(size=2.5) +
    labs(title = "PCA Clustering(MSI)") +
    theme_bw()+
    theme(legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size = 14))
  ggsave("supplementary/pca_after_MSI.png", plot = pca_after_MSI, width = 8, height = 6, dpi = 300)
  pca_after_Lung=ggplot(pca_data_after[pca_data_after$cluster=="Lung",], aes(PC1, PC2, color = cluster,shape = Metastasis)) +
    geom_point(size=2.5) +
    labs(title = "PCA Clustering(Lung)") +
    theme_bw()+
    theme(legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size = 14))
  ggsave("supplementary/pca_after_Lung.png", plot = pca_after_Lung, width = 8, height = 6, dpi = 300)
  pca_after_GI_GYN=ggplot(pca_data_after[pca_data_after$cluster=="GI and GYN",], aes(PC1, PC2, color = cluster,shape = Metastasis)) +
    geom_point(size=2.5) +
    labs(title = "PCA Clustering(GI and GYN)") +
    theme_bw()+
    theme(legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size = 14))
  ggsave("supplementary/pca_after_GI_GYN.png", plot = pca_after_GI_GYN, width = 8, height = 6, dpi = 300)
  pca_after=ggplot(pca_data_after, aes(PC1, PC2, color = cluster,shape = Metastasis)) +
    geom_point(size=2.5) +
    labs(title = "PCA Clustering") +
    theme_bw()+
    theme(legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size = 14))
  ggsave("supplementary/pca_after.png", plot = pca_after, width = 8, height = 6, dpi = 300)
  return(invisible(NULL))
}
tSNE_plot <- function(S5_df=S5_df,before_adj_matrix=before_adj_matrix){
  library(Rtsne)
  library(ggplot2)


  set.seed(3)
  after_adj_matrix=S5_df
  after_adj_matrix=data.frame(after_adj_matrix)
  after_adj_matrix_numeric=after_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=data.frame(after_adj_matrix_numeric)
  after_adj_matrix_numeric <- as.data.frame(sapply(after_adj_matrix_numeric, function(x) as.numeric(as.character(x))))
  tsne_data <- Rtsne(after_adj_matrix_numeric, perplexity = 30, dims = 2, check_duplicates=F)

  tsne_coords <- tsne_data$Y

  kmeans_cluster <- kmeans(tsne_coords, centers = 3)

  tsne_df <- data.frame(x = tsne_coords[, 1], y = tsne_coords[, 2])
  tsne_df$cluster <- as.factor(kmeans_cluster$cluster)
  tsne_df_after=tsne_df
  levels(tsne_df_after$cluster) <- c("GI and GYN", "MSI", "Lung")
  tsne_df_after$Metastasis=after_adj_matrix$Metastasis
  tsne_after = ggplot(tsne_df_after, aes(x, y, color = cluster, shape = Metastasis)) +
    geom_point(size = 1.5, alpha=0.9) +
    labs(title = "Adjusted Clustering") +
    guides(color = guide_legend(override.aes = list(size = 6)))+
    guides(shape = guide_legend(override.aes = list(size = 6)))+
    scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))+
    theme_bw() +
    theme(
      legend.key.size = unit(0.6, "cm"),
      legend.title=element_text(size=20),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 16),
      panel.grid = element_blank(),
      plot.title = element_text(size = 24)
    )+
    scale_shape_manual(values = c("Brain" = 1, "Lung" = 2, "Liver" = 3, "Bone" = 4, "Kidney" = 6))
  ggsave("supplementary/tsne_after.png", plot = tsne_after, width = 8, height = 6, dpi = 300)
  library(Rtsne)
  library(ggplot2)

  set.seed(3)
  before_adj_matrix_numeric=before_adj_matrix[,-c(1,2,3)]
  tsne_data <- Rtsne(before_adj_matrix_numeric, perplexity = 5, dims = 2,check_duplicates =F)

  tsne_coords <- tsne_data$Y

  kmeans_cluster <- kmeans(tsne_coords, centers = 3)

  tsne_df <- data.frame(x = tsne_coords[, 1], y = tsne_coords[, 2])
  tsne_df$cluster <- as.factor(kmeans_cluster$cluster)
  tsne_df_before=tsne_df
  set.seed(3)

  tsne_data <- Rtsne(after_adj_matrix_numeric, perplexity = 50, dims = 2,check_duplicates=F)

  tsne_coords <- tsne_data$Y

  kmeans_cluster <- kmeans(tsne_coords, centers = 3)

  tsne_df <- data.frame(x = tsne_coords[, 1], y = tsne_coords[, 2])
  tsne_df$cluster <- as.factor(kmeans_cluster$cluster)
  tsne_df_after=tsne_df
  levels(tsne_df_after$cluster) <- c("GI and GYN", "MSI", "Lung")
  tsne_df_after$Metastasis=after_adj_matrix$Metastasis
  tsne_after_old=ggplot(tsne_df_after, aes(x, y, color = Metastasis,shape = cluster)) +
    scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))+
    guides(color = guide_legend(override.aes = list(size = 6)))+
    guides(shape = guide_legend(override.aes = list(size = 6)))+
    geom_point(size=1) +
    labs(title = "t-SNE Clustering") +
    theme_bw()+
    theme(
      legend.key.size = unit(0.6, "cm"),
      legend.title=element_text(size=20),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 16),
      panel.grid = element_blank(),
      plot.title = element_text(size = 24)
    )
  ggsave("supplementary/tsne_after_ver2.png", plot = tsne_after_old, width = 8, height = 6, dpi = 300)
  return(invisible(NULL))
}
dumbbell_each <- function(df_table2=df_table2,S5_df=S5_df,before_adj_matrix=before_adj_matrix,GI_brain=GI_brain){
  colnames(before_adj_matrix)=c("depMapID","Lineage","Metastasis",colnames(GI_brain$adj_matrix))
  figure_2_brain=df_table2[df_table2$`Metastasis.Sites`=="Brain",]
  figure_2_Lung=df_table2[df_table2$`Metastasis.Sites`=="Lung",]
  figure_2_Liver=df_table2[df_table2$`Metastasis.Sites`=="Liver",]
  figure_2_Bone=df_table2[df_table2$`Metastasis.Sites`=="Bone",]
  figure_2_Kidney=df_table2[df_table2$`Metastasis.Sites`=="Kidney",]
  MSI_brain_feature=figure_2_brain[figure_2_brain$Cancers=="MSI",3]
  GI_GYN_brain_feature=figure_2_brain[figure_2_brain$Cancers=="GI/GYN",3]
  Lung_brain_feature=figure_2_brain[figure_2_brain$Cancers=="Lung",3]
  MSI_Lung_feature=figure_2_Lung[figure_2_Lung$Cancers=="MSI",3]
  GI_GYN_Lung_feature=figure_2_Lung[figure_2_Lung$Cancers=="GI/GYN",3]
  Lung_Lung_feature=figure_2_Lung[figure_2_Lung$Cancers=="Lung",3]
  MSI_Liver_feature=figure_2_Liver[figure_2_Liver$Cancers=="MSI",3]
  GI_GYN_Liver_feature=figure_2_Liver[figure_2_Liver$Cancers=="GI/GYN",3]
  Lung_Liver_feature=figure_2_Liver[figure_2_Liver$Cancers=="Lung",3]
  MSI_Bone_feature=figure_2_Bone[figure_2_Bone$Cancers=="MSI",3]
  GI_GYN_Bone_feature=figure_2_Bone[figure_2_Bone$Cancers=="GI/GYN",3]
  Lung_Bone_feature=figure_2_Bone[figure_2_Bone$Cancers=="Lung",3]
  MSI_Kidney_feature=figure_2_Kidney[figure_2_Kidney$Cancers=="MSI",3]
  GI_GYN_Kidney_feature=figure_2_Kidney[figure_2_Kidney$Cancers=="GI/GYN",3]
  Lung_Kidney_feature=figure_2_Kidney[figure_2_Kidney$Cancers=="Lung",3]
  after_adj_matrix=S5_df
  after_adj_matrix=data.frame(after_adj_matrix)
  after_adj_matrix_numeric=after_adj_matrix[,-c(1,2,3)]
  after_adj_matrix_numeric=data.frame(after_adj_matrix_numeric)
  after_adj_matrix_numeric <- as.data.frame(sapply(after_adj_matrix_numeric, function(x) as.numeric(as.character(x))))
  before_adj_matrix_numeric=before_adj_matrix[,-c(1,2,3)]




  GI_GYN_brain_feature_t=t(GI_GYN_brain_feature)
  colnames(GI_GYN_brain_feature_t)=GI_GYN_brain_feature

  GI_GYN_brain_feature_name=GI_GYN_brain_feature
  Lung_brain_feature_t=t(Lung_brain_feature)
  colnames(Lung_brain_feature_t)=Lung_brain_feature

  Lung_brain_feature_name=Lung_brain_feature
  MSI_brain_feature_t=t(MSI_brain_feature)
  colnames(MSI_brain_feature_t)=MSI_brain_feature

  MSI_brain_feature_name=MSI_brain_feature
  brain_before_adj_matrix=before_adj_matrix_numeric[,unique(c(GI_GYN_brain_feature_name,Lung_brain_feature_name,MSI_brain_feature_name))]
  colnames(after_adj_matrix_numeric)=colnames(before_adj_matrix_numeric)
  brain_after_adj_matrix=after_adj_matrix_numeric[,unique(c(GI_GYN_brain_feature_name,Lung_brain_feature_name,MSI_brain_feature_name))]


  mean_for_before=colMeans(brain_before_adj_matrix)
  mean_for_after=colMeans(brain_after_adj_matrix)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_df=cbind(colnames(brain_before_adj_matrix),dumbbell_df)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(brain_before_adj_matrix)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_brain_with_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270)+
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(2, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  dumbbell_brain_without_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(2, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  GI_GYN_Lung_feature_t=t(GI_GYN_Lung_feature)
  colnames(GI_GYN_Lung_feature_t)=GI_GYN_Lung_feature

  GI_GYN_Lung_feature_name=GI_GYN_Lung_feature
  Lung_Lung_feature_t=t(Lung_Lung_feature)
  colnames(Lung_Lung_feature_t)=Lung_Lung_feature

  Lung_Lung_feature_name=Lung_Lung_feature
  MSI_Lung_feature_t=t(MSI_Lung_feature)
  colnames(MSI_Lung_feature_t)=MSI_Lung_feature

  MSI_Lung_feature_name=MSI_Lung_feature
  Lung_before_adj_matrix=before_adj_matrix_numeric[,unique(c(GI_GYN_Lung_feature_name,Lung_Lung_feature_name,MSI_Lung_feature_name))]
  Lung_after_adj_matrix=after_adj_matrix_numeric[,unique(c(GI_GYN_Lung_feature_name,Lung_Lung_feature_name,MSI_Lung_feature_name))]

  mean_for_before=colMeans(Lung_before_adj_matrix)
  mean_for_after=colMeans(Lung_after_adj_matrix)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_df=cbind(colnames(Lung_before_adj_matrix),dumbbell_df)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(Lung_before_adj_matrix)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_lung_with_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270)+
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(3, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  dumbbell_lung_without_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(3, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  GI_GYN_Liver_feature_t=t(GI_GYN_Liver_feature)
  colnames(GI_GYN_Liver_feature_t)=GI_GYN_Liver_feature

  GI_GYN_Liver_feature_name=GI_GYN_Liver_feature
  Lung_Liver_feature_t=t(Lung_Liver_feature)
  colnames(Lung_Liver_feature_t)=Lung_Liver_feature

  Lung_Liver_feature_name=Lung_Liver_feature
  MSI_Liver_feature_t=t(MSI_Liver_feature)
  colnames(MSI_Liver_feature_t)=MSI_Liver_feature

  MSI_Liver_feature_name=MSI_Liver_feature
  Liver_before_adj_matrix=before_adj_matrix_numeric[,unique(c(GI_GYN_Liver_feature_name,Lung_Liver_feature_name,MSI_Liver_feature_name))]
  Liver_after_adj_matrix=after_adj_matrix_numeric[,unique(c(GI_GYN_Liver_feature_name,Lung_Liver_feature_name,MSI_Liver_feature_name))]

  mean_for_before=colMeans(Liver_before_adj_matrix)
  mean_for_after=colMeans(Liver_after_adj_matrix)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_df=cbind(colnames(Liver_before_adj_matrix),dumbbell_df)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(Liver_before_adj_matrix)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_liver_with_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270)+
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(2, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  dumbbell_liver_without_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(2, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  GI_GYN_Bone_feature_t=t(GI_GYN_Bone_feature)
  colnames(GI_GYN_Bone_feature_t)=GI_GYN_Bone_feature

  GI_GYN_Bone_feature_name=GI_GYN_Bone_feature
  Lung_Bone_feature_t=t(Lung_Bone_feature)
  colnames(Lung_Bone_feature_t)=Lung_Bone_feature

  Lung_Bone_feature_name=Lung_Bone_feature
  MSI_Bone_feature_t=t(MSI_Bone_feature)
  colnames(MSI_Bone_feature_t)=MSI_Bone_feature

  MSI_Bone_feature_name=MSI_Bone_feature
  Bone_before_adj_matrix=before_adj_matrix_numeric[,unique(c(GI_GYN_Bone_feature_name,Lung_Bone_feature_name,MSI_Bone_feature_name))]
  Bone_after_adj_matrix=after_adj_matrix_numeric[,unique(c(GI_GYN_Bone_feature_name,Lung_Bone_feature_name,MSI_Bone_feature_name))]

  mean_for_before=colMeans(Bone_before_adj_matrix)
  mean_for_after=colMeans(Bone_after_adj_matrix)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_df=cbind(colnames(Bone_before_adj_matrix),dumbbell_df)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(Bone_before_adj_matrix)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_bone_with_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270)+
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(3, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  dumbbell_bone_without_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(3, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  GI_GYN_Kidney_feature_t=t(GI_GYN_Kidney_feature)
  colnames(GI_GYN_Kidney_feature_t)=GI_GYN_Kidney_feature

  GI_GYN_Kidney_feature_name=GI_GYN_Kidney_feature
  Lung_Kidney_feature_t=t(Lung_Kidney_feature)
  colnames(Lung_Kidney_feature_t)=Lung_Kidney_feature

  Lung_Kidney_feature_name=Lung_Kidney_feature
  MSI_Kidney_feature_t=t(MSI_Kidney_feature)
  colnames(MSI_Kidney_feature_t)=MSI_Kidney_feature

  MSI_Kidney_feature_name=MSI_Kidney_feature
  Kidney_before_adj_matrix=before_adj_matrix_numeric[,unique(c(GI_GYN_Kidney_feature_name,Lung_Kidney_feature_name,MSI_Kidney_feature_name))]
  Kidney_after_adj_matrix=after_adj_matrix_numeric[,unique(c(GI_GYN_Kidney_feature_name,Lung_Kidney_feature_name,MSI_Kidney_feature_name))]
  mean_for_before=colMeans(Kidney_before_adj_matrix)
  mean_for_after=colMeans(Kidney_after_adj_matrix)
  dumbbell_df=cbind(mean_for_before,mean_for_after)
  dumbbell_df=cbind(colnames(Kidney_before_adj_matrix),dumbbell_df)
  colnames(dumbbell_df)=c("Metabolite","Mean_for_before","Mean_for_after")
  dumbbell_df=data.frame(dumbbell_df)
  dumbbell_df$Mean_for_before=as.numeric(dumbbell_df$Mean_for_before)
  dumbbell_df$Mean_for_after=as.numeric(dumbbell_df$Mean_for_after)
  dumbbell_df$Metabolite=as.factor(dumbbell_df$Metabolite)
  dumbbell_df2=data.frame(Mean_for_after=mean_for_after,Mean_for_before=mean_for_before,Metabolite=as.factor(colnames(Kidney_before_adj_matrix)))
  dumbbell_df <- dumbbell_df[order(dumbbell_df$Mean_for_after, decreasing = TRUE), ]
  dumbbell_df$order=factor(seq(1:length(dumbbell_df$Metabolite)))
  label_start=c()
  for(i in 1:length(dumbbell_df$Mean_for_after)){
    if(dumbbell_df$Mean_for_before[i]<=dumbbell_df$Mean_for_after[i]){
      label_start=c(label_start,5)
    }
    else{
      label_start=c(label_start,10)
    }
  }
  dumbbell_df=cbind(dumbbell_df,label_start)

  dumbbell_kidney_with_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    geom_text(aes(x = c(order), y = label_start, label = Metabolite), hjust = -0.2, vjust = 0.3, angle = 270)+
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(3, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  dumbbell_kidney_without_text=ggplot(dumbbell_df) +
    geom_point(aes(x = order, y = Mean_for_before, color = "Before Adjustment"), size = 3) +
    geom_point(aes(x = order, y = Mean_for_after, color = "After Adjustment"), size = 3) +
    geom_segment(aes(x = order, xend = order, y = Mean_for_before, yend = Mean_for_after), size = 1, linetype = "dashed") +
    labs(x = NULL, y = "Value", color = "Adjustment") +
    scale_color_manual(values = c("Before Adjustment" = "red", "After Adjustment" = "blue")) +
    theme_bw() +
    ylim(3, 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(1, 1, 3, 1), "lines"),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank()) +
    coord_cartesian(xlim = c(min(as.numeric(dumbbell_df$order)), max(as.numeric(dumbbell_df$order))))
  ggsave("supplementary/dumbbell_brain_with_text.png", plot = dumbbell_brain_with_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_brain_without_text.png", plot = dumbbell_brain_without_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_lung_with_text.png", plot = dumbbell_lung_with_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_lung_without_text.png", plot = dumbbell_lung_without_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_liver_with_text.png", plot = dumbbell_liver_with_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_liver_without_text.png", plot = dumbbell_liver_without_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_bone_with_text.png", plot = dumbbell_bone_with_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_bone_without_text.png", plot = dumbbell_bone_without_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_kidney_with_text.png", plot = dumbbell_kidney_with_text, width = 8, height = 6, dpi = 300)
  ggsave("supplementary/dumbbell_kidney_without_text.png", plot = dumbbell_kidney_without_text, width = 8, height = 6, dpi = 300)
  return(invisible(NULL))
}
Generate_importance <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney){
  listInput <- list(brain = rownames(data.frame(GI_brain$estimate_value, drop = FALSE)), lung = rownames(data.frame(GI_lung$estimate_value, drop = FALSE)), liver = rownames(data.frame(GI_liver$estimate_value, drop = FALSE)), bone = rownames(data.frame(GI_bone$estimate_value, drop = FALSE)), kidney = rownames(data.frame(GI_kidney$estimate_value, drop = FALSE)))
  GI_GYN_listInput <- listInput
  upset(fromList(listInput), keep.order = TRUE, sets = c("kidney", "bone", "liver", "lung", "brain"))

  GI_GYN_importance_list <- unique(c(rownames(data.frame(GI_brain$estimate_value, drop = FALSE)), rownames(data.frame(GI_lung$estimate_value, drop = FALSE)), rownames(data.frame(GI_liver$estimate_value, drop = FALSE)), rownames(data.frame(GI_bone$estimate_value, drop = FALSE)), rownames(data.frame(GI_kidney$estimate_value, drop = FALSE))))

  GI_GYN_importance_df <- data.frame(matrix(0, nrow = length(GI_GYN_importance_list), ncol = 5))
  colnames(GI_GYN_importance_df) <- c("Brain", "Lung", "Liver", "Bone", "Kidney")
  rownames(GI_GYN_importance_df) <- GI_GYN_importance_list

  for (i in 1:length(GI_GYN_importance_list)) {
    if (GI_GYN_importance_list[i] %in% rownames(data.frame(GI_brain$estimate_value, drop = FALSE))) {
      GI_GYN_importance_df[i, "Brain"] <- 1
    }
    if (GI_GYN_importance_list[i] %in% rownames(data.frame(GI_lung$estimate_value, drop = FALSE))) {
      GI_GYN_importance_df[i, "Lung"] <- 1
    }
    if (GI_GYN_importance_list[i] %in% rownames(data.frame(GI_liver$estimate_value, drop = FALSE))) {
      GI_GYN_importance_df[i, "Liver"] <- 1
    }
    if (GI_GYN_importance_list[i] %in% rownames(data.frame(GI_bone$estimate_value, drop = FALSE))) {
      GI_GYN_importance_df[i, "Bone"] <- 1
    }
    if (GI_GYN_importance_list[i] %in% rownames(data.frame(GI_kidney$estimate_value, drop = FALSE))) {
      GI_GYN_importance_df[i, "Kidney"] <- 1
    }
  }
  write.csv(GI_GYN_importance_df, "supplementary/GI_GYN_importance.csv")
  listInput <- list(brain = rownames(data.frame(Lung_brain$estimate_value, drop = FALSE)), lung = rownames(data.frame(Lung_lung$estimate_value, drop = FALSE)), liver = rownames(data.frame(Lung_liver$estimate_value, drop = FALSE)), bone = rownames(data.frame(Lung_bone$estimate_value, drop = FALSE)), kidney = rownames(data.frame(Lung_kidney$estimate_value, drop = FALSE)))
  Lung_listInput <- listInput
  upset(fromList(listInput), keep.order = TRUE, sets = c("kidney", "bone", "liver", "lung", "brain"))

  Lung_importance_list <- unique(c(rownames(data.frame(Lung_brain$estimate_value, drop = FALSE)), rownames(data.frame(Lung_lung$estimate_value, drop = FALSE)), rownames(data.frame(Lung_liver$estimate_value, drop = FALSE)), rownames(data.frame(Lung_bone$estimate_value, drop = FALSE)), rownames(data.frame(Lung_kidney$estimate_value, drop = FALSE))))

  Lung_importance_df <- data.frame(matrix(0, nrow = length(Lung_importance_list), ncol = 5))
  colnames(Lung_importance_df) <- c("Brain", "Lung", "Liver", "Bone", "Kidney")
  rownames(Lung_importance_df) <- Lung_importance_list

  for (i in 1:length(Lung_importance_list)) {
    if (Lung_importance_list[i] %in% rownames(data.frame(Lung_brain$estimate_value, drop = FALSE))) {
      Lung_importance_df[i, "Brain"] <- 1
    }
    if (Lung_importance_list[i] %in% rownames(data.frame(Lung_lung$estimate_value, drop = FALSE))) {
      Lung_importance_df[i, "Lung"] <- 1
    }
    if (Lung_importance_list[i] %in% rownames(data.frame(Lung_liver$estimate_value, drop = FALSE))) {
      Lung_importance_df[i, "Liver"] <- 1
    }
    if (Lung_importance_list[i] %in% rownames(data.frame(Lung_bone$estimate_value, drop = FALSE))) {
      Lung_importance_df[i, "Bone"] <- 1
    }
    if (Lung_importance_list[i] %in% rownames(data.frame(Lung_kidney$estimate_value, drop = FALSE))) {
      Lung_importance_df[i, "Kidney"] <- 1
    }
  }
  write.csv(Lung_importance_df, "supplementary/Lung_importance.csv")
  listInput <- list(brain = rownames(data.frame(MSI_brain$estimate_value, drop = FALSE)), lung = rownames(data.frame(MSI_lung$estimate_value, drop = FALSE)), liver = rownames(data.frame(MSI_liver$estimate_value, drop = FALSE)), bone = rownames(data.frame(MSI_bone$estimate_value, drop = FALSE)), kidney = rownames(data.frame(MSI_kidney$estimate_value, drop = FALSE)))
  MSI_listInput <- listInput
  upset(fromList(listInput), keep.order = TRUE, sets = c("kidney", "bone", "liver", "lung", "brain"))

  MSI_importance_list <- unique(c(rownames(data.frame(MSI_brain$estimate_value, drop = FALSE)), rownames(data.frame(MSI_lung$estimate_value, drop = FALSE)), rownames(data.frame(MSI_liver$estimate_value, drop = FALSE)), rownames(data.frame(MSI_bone$estimate_value, drop = FALSE)), rownames(data.frame(MSI_kidney$estimate_value, drop = FALSE))))

  MSI_importance_df <- data.frame(matrix(0, nrow = length(MSI_importance_list), ncol = 5))
  colnames(MSI_importance_df) <- c("Brain", "Lung", "Liver", "Bone", "Kidney")
  rownames(MSI_importance_df) <- MSI_importance_list

  for (i in 1:length(MSI_importance_list)) {
    if (MSI_importance_list[i] %in% rownames(data.frame(MSI_brain$estimate_value, drop = FALSE))) {
      MSI_importance_df[i, "Brain"] <- 1
    }
    if (MSI_importance_list[i] %in% rownames(data.frame(MSI_lung$estimate_value, drop = FALSE))) {
      MSI_importance_df[i, "Lung"] <- 1
    }
    if (MSI_importance_list[i] %in% rownames(data.frame(MSI_liver$estimate_value, drop = FALSE))) {
      MSI_importance_df[i, "Liver"] <- 1
    }
    if (MSI_importance_list[i] %in% rownames(data.frame(MSI_bone$estimate_value, drop = FALSE))) {
      MSI_importance_df[i, "Bone"] <- 1
    }
    if (MSI_importance_list[i] %in% rownames(data.frame(MSI_kidney$estimate_value, drop = FALSE))) {
      MSI_importance_df[i, "Kidney"] <- 1
    }
  }
  write.csv(MSI_importance_df, "supplementary/MSI_importance.csv")
  return(invisible(NULL))
}
create_plot <- function(target, group, title, subtitle,j) {
  p <- ggplot(data.frame(target), aes(x = group, y = target, fill = group)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(color = "black", size = 0.1, alpha = 0.8) +
    scale_x_discrete("") +
    scale_y_continuous(j) +
    labs(title = title, subtitle = subtitle) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  return(p)
}
create_plot_without_legend <- function(target, group, title, subtitle, j) {
  p <- ggplot(data.frame(target), aes(x = group, y = target, fill = group)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(color = "black", size = 0.1, alpha = 0.8) +
    scale_x_discrete("") +
    scale_y_continuous(j) +
    labs(title = title, subtitle = subtitle) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none")
  return(p)
}
box_plots_generate <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,metabolite_name_trans = Metabolite_pathway,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,U_GI=U_GI,U_Lung=U_Lung,U_MSI=U_MSI,Y_GI=Y_GI,Y_LUNG=Y_Lung,Y_MSI=Y_MSI){
  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$brain))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$brain)")
  n <- sum(Y_GI$brain == "Metastasis")
  GI_brain_box=list()
  j=1
  for (i in row.names(data.frame(GI_brain$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_GI$brain == "Metastasis"], target[Y_GI$brain == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$brain == "Metastasis"], target[Y_GI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_GI$brain), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- GI_brain$adj_matrix[, colnames(GI_brain$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_GI$brain == "Metastasis"], target3[Y_GI$brain == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$brain == "Metastasis"], target3[Y_GI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_GI$brain), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","GI_","brain_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    GI_brain_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$lung))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$lung)")
  n <- sum(Y_GI$lung == "Metastasis")
  GI_lung_box=list()
  j=1
  for (i in row.names(data.frame(GI_lung$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_GI$lung == "Metastasis"], target[Y_GI$lung == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$lung == "Metastasis"], target[Y_GI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_GI$lung), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- GI_lung$adj_matrix[, colnames(GI_lung$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_GI$lung == "Metastasis"], target3[Y_GI$lung == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$lung == "Metastasis"], target3[Y_GI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_GI$lung), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","GI_","lung_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    GI_lung_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$liver))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$liver)")
  n <- sum(Y_GI$liver == "Metastasis")
  GI_liver_box=list()
  j=1
  for (i in row.names(data.frame(GI_liver$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_GI$liver == "Metastasis"], target[Y_GI$liver == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$liver == "Metastasis"], target[Y_GI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_GI$liver), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- GI_liver$adj_matrix[, colnames(GI_liver$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_GI$liver == "Metastasis"], target3[Y_GI$liver == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$liver == "Metastasis"], target3[Y_GI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_GI$liver), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","GI_","liver_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    GI_liver_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$bone))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$bone)")
  n <- sum(Y_GI$bone == "Metastasis")
  GI_bone_box=list()
  j=1
  for (i in row.names(data.frame(GI_bone$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]


    statistic <- wilcox.test(target[Y_GI$bone == "Metastasis"], target[Y_GI$bone == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$bone == "Metastasis"], target[Y_GI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_GI$bone), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- GI_bone$adj_matrix[, colnames(GI_bone$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_GI$bone == "Metastasis"], target3[Y_GI$bone == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$bone == "Metastasis"], target3[Y_GI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_GI$bone), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","GI_","bone_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    GI_bone_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$kidney))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$kidney)")
  n <- sum(Y_GI$kidney == "Metastasis")
  GI_kidney_box=list()
  j=1
  for (i in row.names(data.frame(GI_kidney$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_GI$kidney == "Metastasis"], target[Y_GI$kidney == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$kidney == "Metastasis"], target[Y_GI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_GI$kidney), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- GI_kidney$adj_matrix[, colnames(GI_kidney$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_GI$kidney == "Metastasis"], target3[Y_GI$kidney == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$kidney == "Metastasis"], target3[Y_GI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_GI$kidney), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","GI_","kidney_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    GI_kidney_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$brain))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$brain)")
  n <- sum(Y_LUNG$brain == "Metastasis")
  Lung_brain_box=list()
  j=j+1
  for (i in row.names(data.frame(Lung_brain$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_LUNG$brain == "Metastasis"], target[Y_LUNG$brain == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$brain == "Metastasis"], target[Y_LUNG$brain == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_LUNG$brain), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- Lung_brain$adj_matrix[, colnames(Lung_brain$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_LUNG$brain == "Metastasis"], target3[Y_LUNG$brain == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$brain == "Metastasis"], target3[Y_LUNG$brain == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_LUNG$brain), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","Lung_","brain_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    Lung_brain_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$lung))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$lung)")
  n <- sum(Y_LUNG$lung == "Metastasis")
  Lung_lung_box=list()
  j=1
  for (i in row.names(data.frame(Lung_lung$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_LUNG$lung == "Metastasis"], target[Y_LUNG$lung == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$lung == "Metastasis"], target[Y_LUNG$lung == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_LUNG$lung), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- Lung_lung$adj_matrix[, colnames(Lung_lung$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_LUNG$lung == "Metastasis"], target3[Y_LUNG$lung == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$lung == "Metastasis"], target3[Y_LUNG$lung == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_LUNG$lung), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","Lung_","lung_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    Lung_lung_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$liver))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$liver)")
  n <- sum(Y_LUNG$liver == "Metastasis")
  Lung_liver_box=list()
  j=j+1
  for (i in row.names(data.frame(Lung_liver$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_LUNG$liver == "Metastasis"], target[Y_LUNG$liver == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$liver == "Metastasis"], target[Y_LUNG$liver == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_LUNG$liver), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- Lung_liver$adj_matrix[, colnames(Lung_liver$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_LUNG$liver == "Metastasis"], target3[Y_LUNG$liver == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$liver == "Metastasis"], target3[Y_LUNG$liver == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_LUNG$liver), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","Lung_","liver_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    Lung_liver_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$bone))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$bone)")
  n <- sum(Y_LUNG$bone == "Metastasis")
  Lung_bone_box=list()
  j=j+1
  for (i in row.names(data.frame(Lung_bone$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]


    statistic <- wilcox.test(target[Y_LUNG$bone == "Metastasis"], target[Y_LUNG$bone == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$bone == "Metastasis"], target[Y_LUNG$bone == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_LUNG$bone), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- Lung_bone$adj_matrix[, colnames(Lung_bone$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_LUNG$bone == "Metastasis"], target3[Y_LUNG$bone == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$bone == "Metastasis"], target3[Y_LUNG$bone == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_LUNG$bone), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","Lung_","bone_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    Lung_bone_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$kidney))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$kidney)")
  n <- sum(Y_LUNG$kidney == "Metastasis")
  Lung_kidney_box=list()
  j=j+1
  for (i in row.names(data.frame(Lung_kidney$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- wilcox.test(target[Y_LUNG$kidney == "Metastasis"], target[Y_LUNG$kidney == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$kidney == "Metastasis"], target[Y_LUNG$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_LUNG$kidney), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- Lung_kidney$adj_matrix[, colnames(Lung_kidney$adj_matrix) == i]
    statistic3 <- wilcox.test(target3[Y_LUNG$kidney == "Metastasis"], target3[Y_LUNG$kidney == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$kidney == "Metastasis"], target3[Y_LUNG$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_LUNG$kidney), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","Lung_","kidney_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    Lung_kidney_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$brain))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$brain)")
  n <- sum(Y_MSI$brain == "Metastasis")
  MSI_brain_box=list()
  j=1
  for (i in row.names(data.frame(MSI_brain$estimate_value))) {
    i=row.names(data.frame(MSI_brain$estimate_value))[1]
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$brain == "Metastasis"], target[Y_MSI$brain == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$brain == "Metastasis"], target[Y_MSI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_MSI$brain), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)

    target3 <- MSI_brain$adj_matrix[, colnames(MSI_brain$adj_matrix) == i]
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$brain == "Metastasis"], target3[Y_MSI$brain == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$brain == "Metastasis"], target3[Y_MSI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_MSI$brain), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","MSI_","brain_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    MSI_brain_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$lung))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$lung)")
  n <- sum(Y_MSI$lung == "Metastasis")
  MSI_lung_box=list()
  j=1
  for (i in row.names(data.frame(MSI_lung$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$lung == "Metastasis"], target[Y_MSI$lung == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$lung == "Metastasis"], target[Y_MSI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_MSI$lung), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)

    target3 <- MSI_lung$adj_matrix[, colnames(MSI_lung$adj_matrix) == i]
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$lung == "Metastasis"], target3[Y_MSI$lung == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$lung == "Metastasis"], target3[Y_MSI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_MSI$lung), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","MSI_","lung_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    MSI_lung_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$liver))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$liver)")
  n <- sum(Y_MSI$liver == "Metastasis")
  MSI_liver_box=list()
  j=1
  for (i in row.names(data.frame(MSI_liver$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$liver == "Metastasis"], target[Y_MSI$liver == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$liver == "Metastasis"], target[Y_MSI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_MSI$liver), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)

    target3 <- MSI_liver$adj_matrix[, colnames(MSI_liver$adj_matrix) == i]
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$liver == "Metastasis"], target3[Y_MSI$liver == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$liver == "Metastasis"], target3[Y_MSI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_MSI$liver), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","MSI_","liver_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    MSI_liver_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$bone))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$bone)")
  n <- sum(Y_MSI$bone == "Metastasis")
  MSI_bone_box=list()
  j=1
  for (i in row.names(data.frame(MSI_bone$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]


    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$bone == "Metastasis"], target[Y_MSI$bone == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$bone == "Metastasis"], target[Y_MSI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_MSI$bone), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- MSI_bone$adj_matrix[, colnames(MSI_bone$adj_matrix) == i]
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$bone == "Metastasis"], target3[Y_MSI$bone == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$bone == "Metastasis"], target3[Y_MSI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_MSI$bone), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1  + p3
    filename=paste0("supplementary/","MSI_","bone_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    MSI_bone_box[[j]]=ggp_all
    j=j+1
  }

  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$kidney))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$kidney)")
  n <- sum(Y_MSI$kidney == "Metastasis")
  MSI_kidney_box=list()
  j=1
  for (i in row.names(data.frame(MSI_kidney$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]

    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$kidney == "Metastasis"], target[Y_MSI$kidney == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$kidney == "Metastasis"], target[Y_MSI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    p1 <- create_plot_without_legend(target, factor(Y_MSI$kidney), "(Non-adjusted)", paste0(expression("statistic"), statistic, expression(",p"), p.value),i)


    target3 <- MSI_kidney$adj_matrix[, colnames(MSI_kidney$adj_matrix) == i]
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$kidney == "Metastasis"], target3[Y_MSI$kidney == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$kidney == "Metastasis"], target3[Y_MSI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    p3 <- create_plot(target3, factor(Y_MSI$kidney), "(Catch-adjusted)", paste0(expression("statistic"), statistic3, expression(",p"), p.value3),i)

    ggp_all <- p1 +  p3
    filename=paste0("supplementary/","MSI_","kidney_",i,"_boxplot.png")
    ggsave(filename, plot = ggp_all, width = 8, height = 3, dpi = 300)
    MSI_kidney_box[[j]]=ggp_all
    j=j+1
  }

  return(invisible(NULL))
}
make_meta_vs_non_meta <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,U_GI=U_GI,U_Lung=U_Lung,U_MSI=U_MSI,Y_GI=Y_GI,Y_LUNG=Y_Lung,Y_MSI=Y_MSI){
  table_meta_nmeta=c()
  num_GI_brain=c(sum(Y_GI$brain=="Metastasis"),sum(Y_GI$brain=="NonMetastasis"),sum(Y_GI$brain=="Metastasis")/length(Y_GI$brain),sum(Y_GI$brain=="NonMetastasis")/length(Y_GI$brain))
  num_GI_lung=c(sum(Y_GI$lung=="Metastasis"),sum(Y_GI$lung=="NonMetastasis"),sum(Y_GI$lung=="Metastasis")/length(Y_GI$lung),sum(Y_GI$lung=="NonMetastasis")/length(Y_GI$lung))
  num_GI_liver=c(sum(Y_GI$liver=="Metastasis"),sum(Y_GI$liver=="NonMetastasis"),sum(Y_GI$liver=="Metastasis")/length(Y_GI$liver),sum(Y_GI$liver=="NonMetastasis")/length(Y_GI$liver))
  num_GI_bone=c(sum(Y_GI$bone=="Metastasis"),sum(Y_GI$bone=="NonMetastasis"),sum(Y_GI$bone=="Metastasis")/length(Y_GI$bone),sum(Y_GI$bone=="NonMetastasis")/length(Y_GI$bone))
  num_GI_kidney=c(sum(Y_GI$kidney=="Metastasis"),sum(Y_GI$kidney=="NonMetastasis"),sum(Y_GI$kidney=="Metastasis")/length(Y_GI$kidney),sum(Y_GI$kidney=="NonMetastasis")/length(Y_GI$kidney))
  num_LUNG_brain=c(sum(Y_LUNG$brain=="Metastasis"),sum(Y_LUNG$brain=="NonMetastasis"),sum(Y_LUNG$brain=="Metastasis")/length(Y_LUNG$brain),sum(Y_LUNG$brain=="NonMetastasis")/length(Y_LUNG$brain))
  num_LUNG_lung=c(sum(Y_LUNG$lung=="Metastasis"),sum(Y_LUNG$lung=="NonMetastasis"),sum(Y_LUNG$lung=="Metastasis")/length(Y_LUNG$lung),sum(Y_LUNG$lung=="NonMetastasis")/length(Y_LUNG$lung))
  num_LUNG_liver=c(sum(Y_LUNG$liver=="Metastasis"),sum(Y_LUNG$liver=="NonMetastasis"),sum(Y_LUNG$liver=="Metastasis")/length(Y_LUNG$liver),sum(Y_LUNG$liver=="NonMetastasis")/length(Y_LUNG$liver))
  num_LUNG_bone=c(sum(Y_LUNG$bone=="Metastasis"),sum(Y_LUNG$bone=="NonMetastasis"),sum(Y_LUNG$bone=="Metastasis")/length(Y_LUNG$bone),sum(Y_LUNG$bone=="NonMetastasis")/length(Y_LUNG$bone))
  num_LUNG_kidney=c(sum(Y_LUNG$kidney=="Metastasis"),sum(Y_LUNG$kidney=="NonMetastasis"),sum(Y_LUNG$kidney=="Metastasis")/length(Y_LUNG$kidney),sum(Y_LUNG$kidney=="NonMetastasis")/length(Y_LUNG$kidney))
  num_MSI_brain=c(sum(Y_MSI$brain=="Metastasis"),sum(Y_MSI$brain=="NonMetastasis"),sum(Y_MSI$brain=="Metastasis")/length(Y_MSI$brain),sum(Y_MSI$brain=="NonMetastasis")/length(Y_MSI$brain))
  num_MSI_lung=c(sum(Y_MSI$lung=="Metastasis"),sum(Y_MSI$lung=="NonMetastasis"),sum(Y_MSI$lung=="Metastasis")/length(Y_MSI$lung),sum(Y_MSI$lung=="NonMetastasis")/length(Y_MSI$lung))
  num_MSI_liver=c(sum(Y_MSI$liver=="Metastasis"),sum(Y_MSI$liver=="NonMetastasis"),sum(Y_MSI$liver=="Metastasis")/length(Y_MSI$liver),sum(Y_MSI$liver=="NonMetastasis")/length(Y_MSI$liver))
  num_MSI_bone=c(sum(Y_MSI$bone=="Metastasis"),sum(Y_MSI$bone=="NonMetastasis"),sum(Y_MSI$bone=="Metastasis")/length(Y_MSI$bone),sum(Y_MSI$bone=="NonMetastasis")/length(Y_MSI$bone))
  num_MSI_kidney=c(sum(Y_MSI$kidney=="Metastasis"),sum(Y_MSI$kidney=="NonMetastasis"),sum(Y_MSI$kidney=="Metastasis")/length(Y_MSI$kidney),sum(Y_MSI$kidney=="NonMetastasis")/length(Y_MSI$kidney))
  table_meta_nmeta=rbind(num_GI_brain,num_GI_lung,num_GI_liver,num_GI_bone,num_GI_kidney,num_LUNG_brain,num_LUNG_lung,num_LUNG_liver,num_LUNG_bone,num_LUNG_kidney,num_MSI_brain,num_MSI_lung,num_MSI_liver,num_MSI_bone,num_MSI_kidney)
  table_meta_nmeta=round(table_meta_nmeta,2)
  colnames(table_meta_nmeta)=c("Metastasis","Nonmetastasis","Metastasis(%)","Nonmetastasis(%)")
  rownames(table_meta_nmeta)=c("GI_brain","GI_lung","GI_liver","GI_bone","GI_kidney","Lung_brain","Lung_lung","Lung_liver","Lung_bone","Lung_kidney","MSI_brain","MSI_lung","MSI_liver","MSI_bone","MSI_kidney")
  write.csv(table_meta_nmeta,"supplementary/table_meta_or_nmeta.csv",row.names = T)
  return(invisible(NULL))
}
create_raincloudplot_without_legend <- function(value, group, title, subtitle, j) {
  group_new=c()
  for(k in 1:length(group)){
    if(group[k]=="Metastasis"){
      group_new=c(group_new,"(+)")
    }else{
      group_new=c(group_new,"(-)")
    }
  }
  group=group_new
  group = factor(group_new, levels = c("(+)", "(-)"))
  rain_cloud_plot <- ggplot(data.frame(value), aes(x = group, y = value, fill = group, colour = group)) +
    geom_flat_violin(position = position_nudge(x = 0.1), adjust = 2, trim = FALSE) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.1), size = 0.25) +
    geom_boxplot(aes(x = group, y = value), position = position_nudge(x = 0), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    xlab('') + ylab('') +
    theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
    scale_fill_manual(values = c("(+)" = "#FB8072", "(-)" = "#80B1D3")) +
    scale_color_manual(values = c("(+)" = "#FB8072", "(-)" = "#80B1D3")) +
    labs(title = title, subtitle = subtitle) +
    theme(axis.text.x = element_text(size = 15, hjust = 0.2),
          plot.title = element_text(hjust = 0,size=18),
          plot.subtitle = element_text(hjust = 0, face = "italic",size=18),
          legend.position = "none",
          axis.title.y = element_text(size = 17),
          text = element_text(family = "Arial"),
          axis.ticks = element_blank())
  return(rain_cloud_plot)
}
create_raincloudplot_with_legend <- function(value, group, title, subtitle, j) {
  group_new=c()
  for(k in 1:length(group)){
    if(group[k]=="Metastasis"){
      group_new=c(group_new,"Metastasis (+)")
    }else{
      group_new=c(group_new,"              (-)")
    }
  }
  group=group_new
  group = factor(group_new, levels = c("Metastasis (+)", "              (-)"))
  rain_cloud_plot <- ggplot(data.frame(value), aes(x = group, y = value, fill = group, colour = group)) +
    geom_flat_violin(position = position_nudge(x = 0.1), adjust = 2, trim = FALSE) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.1), size = 0.25) +
    geom_boxplot(aes(x = group, y = value), position = position_nudge(x = 0), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    xlab('') + ylab(j) +
    theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
    scale_fill_manual(values = c("Metastasis (+)" = "#FB8072", "              (-)" = "#80B1D3")) +
    scale_color_manual(values = c("Metastasis (+)" = "#FB8072", "              (-)" = "#80B1D3")) +
    labs(title = title, subtitle = subtitle) +
    theme(axis.text.x = element_text(size = 15, hjust = 0.85),
          plot.title = element_text(hjust = 0,size=18),
          plot.subtitle = element_text(hjust = 0, face = "italic",size=18),
          legend.position = "none",
          axis.title.y = element_text(size = 17),
          text = element_text(family = "Arial"),
          axis.ticks = element_blank())
  return(rain_cloud_plot)
}
raincloud_plot_generate <- function(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,U_GI=U_GI,U_Lung=U_Lung,U_MSI=U_MSI,Y_GI=Y_GI,Y_LUNG=Y_Lung,Y_MSI=Y_MSI,metabolite_name_trans=metabolite_name_trans){
  library(PupillometryR)
  library(cowplot)
  library(ggdist)
  library(ggplot2)
  library(patchwork)
  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$brain))
  colnames(dfx) <- c(metabolite_name_trans, "factor(Y_GI$brain)")
  n <- sum(Y_GI$brain == "Metastasis")
  GI_brain_rain=list()
  j=1
  for (i in row.names(data.frame(GI_brain$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_GI$brain)`
    target3 <- GI_brain$adj_matrix[, colnames(GI_brain$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_GI$brain == "Metastasis"], target[Y_GI$brain == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$brain == "Metastasis"], target[Y_GI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_GI$brain == "Metastasis"], target3[Y_GI$brain == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$brain == "Metastasis"], target3[Y_GI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot <- rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","GI_","brain_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    GI_brain_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$lung))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$lung)")
  n <- sum(Y_GI$lung == "Metastasis")
  GI_lung_rain=list()
  j=1
  for (i in row.names(data.frame(GI_lung$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_GI$lung)`
    target3 <- GI_lung$adj_matrix[, colnames(GI_lung$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_GI$lung == "Metastasis"], target[Y_GI$lung == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$lung == "Metastasis"], target[Y_GI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_GI$lung == "Metastasis"], target3[Y_GI$lung == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$lung == "Metastasis"], target3[Y_GI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","GI_","lung_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    GI_lung_rain[[j]]=rain_cloud_plot
    j=j+1

  }
  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$liver))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$liver)")
  n <- sum(Y_GI$liver == "Metastasis")
  GI_liver_rain=list()
  j=1
  for (i in row.names(data.frame(GI_liver$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_GI$liver)`
    target3 <- GI_liver$adj_matrix[, colnames(GI_liver$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_GI$liver == "Metastasis"], target[Y_GI$liver == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$liver == "Metastasis"], target[Y_GI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_GI$liver == "Metastasis"], target3[Y_GI$liver == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$liver == "Metastasis"], target3[Y_GI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","GI_","liver_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    GI_liver_rain[[j]]=rain_cloud_plot
    j=j+1

  }
  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$bone))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$bone)")
  n <- sum(Y_GI$bone == "Metastasis")
  GI_bone_rain=list()
  j=1
  for (i in row.names(data.frame(GI_bone$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_GI$bone)`
    target3 <- GI_bone$adj_matrix[, colnames(GI_bone$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_GI$bone == "Metastasis"], target[Y_GI$bone == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$bone == "Metastasis"], target[Y_GI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_GI$bone == "Metastasis"], target3[Y_GI$bone == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$bone == "Metastasis"], target3[Y_GI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","GI_","bone_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    GI_bone_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_GI[, -1])
  dfx <- cbind(dfx, factor(Y_GI$kidney))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_GI$kidney)")
  n <- sum(Y_GI$kidney == "Metastasis")
  GI_kidney_rain=list()
  j=1
  for (i in row.names(data.frame(GI_kidney$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_GI$kidney)`
    target3 <- GI_kidney$adj_matrix[, colnames(GI_kidney$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_GI$kidney == "Metastasis"], target[Y_GI$kidney == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_GI$kidney == "Metastasis"], target[Y_GI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_GI$kidney == "Metastasis"], target3[Y_GI$kidney == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_GI$kidney == "Metastasis"], target3[Y_GI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","GI_","kidney_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    GI_kidney_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$brain))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$brain)")
  n <- sum(Y_LUNG$brain == "Metastasis")
  Lung_brain_rain=list()
  j=1
  for (i in row.names(data.frame(Lung_brain$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_LUNG$brain)`
    target3 <- Lung_brain$adj_matrix[, colnames(Lung_brain$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_LUNG$brain == "Metastasis"], target[Y_LUNG$brain == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$brain == "Metastasis"], target[Y_LUNG$brain == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_LUNG$brain == "Metastasis"], target3[Y_LUNG$brain == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$brain == "Metastasis"], target3[Y_LUNG$brain == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","Lung_","brain_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    Lung_brain_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$lung))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$lung)")
  n <- sum(Y_LUNG$lung == "Metastasis")
  Lung_lung_rain=list()
  j=1
  for (i in row.names(data.frame(Lung_lung$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_LUNG$lung)`
    target3 <- Lung_lung$adj_matrix[, colnames(Lung_lung$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_LUNG$lung == "Metastasis"], target[Y_LUNG$lung == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$lung == "Metastasis"], target[Y_LUNG$lung == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_LUNG$lung == "Metastasis"], target3[Y_LUNG$lung == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$lung == "Metastasis"], target3[Y_LUNG$lung == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","Lung_","lung_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    Lung_lung_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$liver))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$liver)")
  n <- sum(Y_LUNG$liver == "Metastasis")
  Lung_liver_rain=list()
  j=1
  for (i in row.names(data.frame(Lung_liver$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_LUNG$liver)`
    target3 <- Lung_liver$adj_matrix[, colnames(Lung_liver$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_LUNG$liver == "Metastasis"], target[Y_LUNG$liver == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$liver == "Metastasis"], target[Y_LUNG$liver == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_LUNG$liver == "Metastasis"], target3[Y_LUNG$liver == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$liver == "Metastasis"], target3[Y_LUNG$liver == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","Lung_","liver_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    Lung_liver_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$bone))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$bone)")
  n <- sum(Y_LUNG$bone == "Metastasis")
  Lung_bone_rain=list()
  j=1
  for (i in row.names(data.frame(Lung_bone$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_LUNG$bone)`
    target3 <- Lung_bone$adj_matrix[, colnames(Lung_bone$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_LUNG$bone == "Metastasis"], target[Y_LUNG$bone == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$bone == "Metastasis"], target[Y_LUNG$bone == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_LUNG$bone == "Metastasis"], target3[Y_LUNG$bone == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$bone == "Metastasis"], target3[Y_LUNG$bone == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","Lung_","bone_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    Lung_bone_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_Lung[, -1])
  dfx <- cbind(dfx, factor(Y_LUNG$kidney))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_LUNG$kidney)")
  n <- sum(Y_LUNG$kidney == "Metastasis")
  Lung_kidney_rain=list()
  j=1
  for (i in row.names(data.frame(Lung_kidney$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_LUNG$kidney)`
    target3 <- Lung_kidney$adj_matrix[, colnames(Lung_kidney$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- wilcox.test(target[Y_LUNG$kidney == "Metastasis"], target[Y_LUNG$kidney == "NonMetastasis"])$statistic
    p.value <- round(wilcox.test(target[Y_LUNG$kidney == "Metastasis"], target[Y_LUNG$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- wilcox.test(target3[Y_LUNG$kidney == "Metastasis"], target3[Y_LUNG$kidney == "NonMetastasis"])$statistic
    p.value3 <- round(wilcox.test(target3[Y_LUNG$kidney == "Metastasis"], target3[Y_LUNG$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","Lung_","kidney_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    Lung_kidney_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$brain))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$brain)")
  n <- sum(Y_MSI$brain == "Metastasis")
  MSI_brain_rain=list()
  j=1
  for (i in row.names(data.frame(MSI_brain$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_MSI$brain)`
    target3 <- MSI_brain$adj_matrix[, colnames(MSI_brain$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$brain == "Metastasis"], target[Y_MSI$brain == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$brain == "Metastasis"], target[Y_MSI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$brain == "Metastasis"], target3[Y_MSI$brain == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$brain == "Metastasis"], target3[Y_MSI$brain == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","MSI_","brain_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    MSI_brain_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$lung))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$lung)")
  n <- sum(Y_MSI$lung == "Metastasis")
  MSI_lung_rain=list()
  j=1
  for (i in row.names(data.frame(MSI_lung$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_MSI$lung)`
    target3 <- MSI_lung$adj_matrix[, colnames(MSI_lung$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$lung == "Metastasis"], target[Y_MSI$lung == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$lung == "Metastasis"], target[Y_MSI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$lung == "Metastasis"], target3[Y_MSI$lung == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$lung == "Metastasis"], target3[Y_MSI$lung == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","MSI_","lung_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    MSI_lung_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$liver))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$liver)")
  n <- sum(Y_MSI$liver == "Metastasis")
  MSI_liver_rain=list()
  j=1
  for (i in row.names(data.frame(MSI_liver$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_MSI$liver)`
    target3 <- MSI_liver$adj_matrix[, colnames(MSI_liver$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$liver == "Metastasis"], target[Y_MSI$liver == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$liver == "Metastasis"], target[Y_MSI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$liver == "Metastasis"], target3[Y_MSI$liver == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$liver == "Metastasis"], target3[Y_MSI$liver == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","MSI_","liver_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    MSI_liver_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$bone))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$bone)")
  n <- sum(Y_MSI$bone == "Metastasis")
  MSI_bone_rain=list()
  j=1
  for (i in row.names(data.frame(MSI_bone$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_MSI$bone)`
    target3 <- MSI_bone$adj_matrix[, colnames(MSI_bone$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$bone == "Metastasis"], target[Y_MSI$bone == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$bone == "Metastasis"], target[Y_MSI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$bone == "Metastasis"], target3[Y_MSI$bone == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$bone == "Metastasis"], target3[Y_MSI$bone == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","MSI_","bone_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    MSI_bone_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  dfx <- data.frame(X_MSI[, -1])
  dfx <- cbind(dfx, factor(Y_MSI$kidney))
  colnames(dfx)=c(metabolite_name_trans,"factor(Y_MSI$kidney)")
  n <- sum(Y_MSI$kidney == "Metastasis")
  MSI_kidney_rain=list()
  j=1
  for (i in row.names(data.frame(MSI_kidney$estimate_value))) {
    set.seed(1)
    target <- dfx[, colnames(dfx) == i]
    dfy <- dfx$`factor(Y_MSI$kidney)`
    target3 <- MSI_kidney$adj_matrix[, colnames(MSI_kidney$adj_matrix) == i]
    before_df_rain <- data.frame(value = target, group = dfy)
    after_df_rain <- data.frame(value = target3, group = dfy)
    statistic <- suppressWarnings(wilcox.test(target[Y_MSI$kidney == "Metastasis"], target[Y_MSI$kidney == "NonMetastasis"])$statistic)
    p.value <- round(wilcox.test(target[Y_MSI$kidney == "Metastasis"], target[Y_MSI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value == 0) {
      p.value <- "<0.001"
    }else{
      p.value <- paste0(" = ",p.value)
    }
    rain_cloud_plot_before = create_raincloudplot_with_legend(before_df_rain$value,before_df_rain$group,"Non-adjusted",paste0(expression("p"), p.value),i)
    statistic3 <- suppressWarnings(wilcox.test(target3[Y_MSI$kidney == "Metastasis"], target3[Y_MSI$kidney == "NonMetastasis"])$statistic)
    p.value3 <- round(wilcox.test(target3[Y_MSI$kidney == "Metastasis"], target3[Y_MSI$kidney == "NonMetastasis"])$p.value, 3)
    if (p.value3 == 0) {
      p.value3 <- "<0.001"
    }else{
      p.value3 <- paste0(" = ",p.value3)
    }
    rain_cloud_plot_after = create_raincloudplot_without_legend(after_df_rain$value,after_df_rain$group,"Adjusted",paste0(expression("p"), p.value3),i)
    rain_cloud_plot = rain_cloud_plot_before+rain_cloud_plot_after
    i <- gsub("/", "_", i)
    i <- gsub(":", "_", i)
    i <- gsub(" ", "_", i)
    filename=paste0("supplementary/","MSI_","kidney_",i,"_raincloud.png")
    ggsave(filename, plot = rain_cloud_plot, width = 8, height = 4, dpi = 300)
    MSI_kidney_rain[[j]]=rain_cloud_plot
    j=j+1
  }
  return(invisible(NULL))
}
Generate_Three_cancer_type_data <- function(X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,Sample_ID=Sample_ID){
  library(openxlsx)
  GI_info_df=Sample_ID[Sample_ID$DepMap_ID%in%X_GI$X,]
  colnames_GI=colnames(GI_info_df)
  Lung_info_df=Sample_ID[Sample_ID$DepMap_ID%in%X_Lung$X,]
  colnames_Lung=colnames(Lung_info_df)
  MSI_info_df=Sample_ID[Sample_ID$DepMap_ID%in%X_MSI$X,]
  colnames_MSI=colnames(MSI_info_df)
  df_GI_Y=Y_GI[match(GI_info_df$DepMap_ID, Y_GI$DepMap_ID),]
  GI_info_df=cbind(GI_info_df,df_GI_Y$brain,df_GI_Y$lung,df_GI_Y$liver,df_GI_Y$bone,df_GI_Y$kidney)
  df_Lung_Y=Y_LUNG[match(Lung_info_df$DepMap_ID, Y_LUNG$DepMap_ID),]
  Lung_info_df=cbind(Lung_info_df,df_Lung_Y$brain,df_Lung_Y$lung,df_Lung_Y$liver,df_Lung_Y$bone,df_Lung_Y$kidney)
  df_MSI_Y=Y_MSI[match(MSI_info_df$DepMap_ID, Y_MSI$DepMap_ID),]
  MSI_info_df=cbind(MSI_info_df,df_MSI_Y$brain,df_MSI_Y$lung,df_MSI_Y$liver,df_MSI_Y$bone,df_MSI_Y$kidney)
  colnames(GI_info_df)=c(colnames_GI,"Brain","Lung","Liver","Bone","Kidney")
  colnames(Lung_info_df)=c(colnames_Lung,"Brain","Lung","Liver","Bone","Kidney")
  colnames(MSI_info_df)=c(colnames_MSI,"Brain","Lung","Liver","Bone","Kidney")
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "GI_GYN")
  writeData(wb, sheet = "GI_GYN", GI_info_df)
  addWorksheet(wb, sheetName = "Lung")
  writeData(wb, sheet = "Lung", Lung_info_df)
  addWorksheet(wb, sheetName = "MSI")
  writeData(wb, sheet = "MSI", MSI_info_df)
  saveWorkbook(wb, "Supplementary_1_three_cancer_type.xlsx",overwrite = TRUE)
}


























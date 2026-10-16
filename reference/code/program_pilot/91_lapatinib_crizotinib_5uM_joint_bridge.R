#!/usr/bin/env Rscript
# v82 joint molecular bridge restricted to lapatinib/crizotinib, both at 5 uM.
# The v82 decoder/objective below are preserved. No pair or other-assay labels
# exist in this fit; their unused coefficients remain zero. No PLS is used.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(Matrix);library(glmnet);library(digest);library(ggplot2)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
root <- '/work/results/program_pilot/91_lapatinib_crizotinib_5uM_joint_bridge'
for(folder in c('data','models','figures')) dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
P <- sprintf('P%02d',1:40)
AU <- 'ACH-000248'
selected_drugs <- c('Lapatinib ditosylate','crizotinib')
native_doses <- 5
baseline_rank <- 8L
molecular_ridge <- 1
molecular_weight <- 10
observed_function_weight <- 1
within_function_weight <- .8
# All decoder and optimization settings are the fixed v82 settings.
decoder_rank <- 4L
latent_pairs <- which(upper.tri(matrix(0,decoder_rank,decoder_rank),diag=TRUE),arr.ind=TRUE)
nonlinear_features <- 80L+nrow(latent_pairs)
interaction_features <- 40L+nrow(latent_pairs)
decoder_penalty <- .003
interaction_penalty <- .01
maximum_iterations <- 500L
families <- c('Frozen molecular bridge','Linked molecular bridge')
assays <- c('secondary','primary','nci_FF','nci_FG')
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
status <- function(x) {cat(format(Sys.time()),x,'\n');flush.console();writeLines(x,file.path(root,'STATUS.txt'))}
pearson <- function(x,y) if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
macro_weights <- function(group) {t <- table(group);1/as.numeric(t[group])/length(t)}
getmatrix <- function(dat,ids,cols) {i<-match(ids,dat$depmap_id);stopifnot(!anyNA(i));as.matrix(dat[i,cols])}
key <- function(cell,drug) paste(cell,drug,sep='|')
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
archive <- readRDS(archive_path)
PREP <- archive$prep
PREP$m <- as.data.table(PREP$m)[drug %in% selected_drugs & dose==5 & depmap_id!=AU]
# Match actual PRISM curves to 5 uM: log-dose interpolation, no extrapolation,
# first within screen/compound, then replicate-weighted across screens (v78/v82).
prism <- '/work/cache/tahoe_prism_full_chain/prism'
treatment <- fread(file.path(prism,'secondary-screen-replicate-treatment-info.csv'))[
  tolower(name) %in% c('lapatinib','crizotinib') & perturbation_type=='experimental_treatment']
treatment[,drug:=ifelse(tolower(name)=='lapatinib',selected_drugs[1],selected_drugs[2])]
raw_path <- file.path(prism,'secondary-screen-logfold-change.csv')
header <- names(fread(raw_path,nrows=0))
raw <- fread(raw_path,select=c(header[1],treatment$column_name))
setnames(raw,1L,'depmap_id')
raw <- raw[depmap_id %in% PREP$raw$depmap_id & depmap_id!=AU]
z <- melt(raw,id.vars='depmap_id',variable.name='column_name',value.name='y',variable.factor=FALSE)
z <- merge(z,treatment[,.(column_name,drug,dose,screen_id,broad_id)],by='column_name')
curve <- z[is.finite(y),.(y=mean(y),replicates=.N),by=.(depmap_id,drug,screen_id,broad_id,dose)]
native <- curve[,{o<-order(dose);.(dose=5,y=approx(log(dose[o]),y[o],xout=log(5),rule=1)$y,
  weight=approx(log(dose[o]),replicates[o],xout=log(5),rule=1)$y)},by=.(depmap_id,drug,screen_id,broad_id)]
SINGLE <- native[is.finite(y)&is.finite(weight),.(y=weighted.mean(y,weight)),by=.(depmap_id,drug,dose)]
SINGLE[,`:=`(assay='secondary',source='prism_secondary',group=drug)]
setorder(SINGLE,drug,depmap_id)
SINGLE[,evaluation_id:=.I]
old <- fread('/work/results/program_pilot/82_joint_program_rules_and_pair_interactions/data/recipient_folds.tsv')
recipient_fold <- setNames(old$fold,old$depmap_id)
missing_cells <- setdiff(SINGLE$depmap_id,names(recipient_fold))
set.seed(77012)
if(length(missing_cells)) recipient_fold <- c(recipient_fold,setNames(sample(rep(1:3,length.out=length(missing_cells))),missing_cells))
SINGLE[,fold:=unname(recipient_fold[depmap_id])]
stopifnot(!anyNA(SINGLE$fold),uniqueN(SINGLE$drug)==2L,!anyDuplicated(SINGLE[,.(depmap_id,drug)]),all(SINGLE$dose==5))
pair_available <- as.data.table(archive$pair)[drugA %in% selected_drugs & drugB %in% selected_drugs & doseA==5 & doseB==5]
stopifnot(nrow(pair_available)==0L)
write_table(SINGLE,'individual_training_outcomes')
write_table(data.table(depmap_id=names(recipient_fold),fold=recipient_fold),'recipient_folds')
write_table(SINGLE[,.(cells=.N,observed_below_50pct=sum(y < -1)),by=.(drug,dose)],'coverage')
write_table(pair_available,'available_5uM_combinations')
status(paste('Inputs ready:',nrow(SINGLE),'PRISM outcomes;',nrow(PREP$m),'RNA profiles; zero matching combinations.'))
rm(raw,z,curve,archive)
invisible(gc())

# Same v82 preparation with one molecular anchor per drug and no pair rows.
prepare <- function(excluded) {
  molecular <- PREP$m[!depmap_id %in% excluded]
  ref_ids <- sort(unique(molecular$depmap_id))
  baseline <- getmatrix(PREP$raw,ref_ids,P)
  baseline_mean <- colMeans(baseline)
  baseline_sd <- apply(baseline,2,sd)
  centered <- molecular[,lapply(.SD,function(v)v-mean(v)),by=drug,.SDcols=P]
  response_sd <- vapply(centered[,..P],sd,numeric(1))
  stopifnot(all(baseline_sd>0),all(response_sd>0))
  training_baseline <- sweep(sweep(baseline,2,baseline_mean),2,baseline_sd,'/')
  baseline_projection <- svd(training_baseline,nu=0,nv=baseline_rank)$v
  baseline_projection_scale <- sqrt(colMeans((training_baseline%*%baseline_projection)^2))
  ss <- SINGLE[!depmap_id %in% excluded]
  query <- unique(rbind(ss[,.(depmap_id,drug)],molecular[,.(depmap_id,drug)]))
  query[,dose:=5]
  nquery <- nrow(query)
  query_key <- key(query$depmap_id,query$drug)
  s <- sweep(sweep(sweep(getmatrix(PREP$raw,query$depmap_id,P),2,baseline_mean),2,baseline_sd,'/')%*%
    baseline_projection,2,baseline_projection_scale,'/')
  molecular_response <- sweep(as.matrix(molecular[,..P]),2,response_sd,'/')
  molecular_index <- match(key(molecular$depmap_id,molecular$drug),query_key)
  group_rows <- lapply(selected_drugs,function(dr)which(query$drug==dr))
  group_molecular_rows <- lapply(selected_drugs,function(dr)which(molecular$drug==dr))
  group_baseline_mean <- matrix(0,2,baseline_rank)
  group_response_mean <- matrix(0,2,40)
  coefficient <- array(0,c(baseline_rank,40,2))
  initial_response <- matrix(0,nquery,40)
  for(g in 1:2) {
    ii <- group_rows[[g]]
    mi <- group_molecular_rows[[g]]
    stopifnot(length(mi)>=baseline_rank+2)
    group_baseline_mean[g,] <- colMeans(s[molecular_index[mi],,drop=FALSE])
    group_response_mean[g,] <- colMeans(molecular_response[mi,,drop=FALSE])
    s[ii,] <- sweep(s[ii,,drop=FALSE],2,group_baseline_mean[g,])
    fit <- glmnet(s[molecular_index[mi],,drop=FALSE],sweep(molecular_response[mi,,drop=FALSE],2,group_response_mean[g,]),
      family='mgaussian',alpha=0,lambda=molecular_ridge,intercept=FALSE,standardize=FALSE,thresh=1e-10)
    coefficient[,,g] <- vapply(coef(fit),function(v)as.numeric(v)[-1],numeric(baseline_rank))
    initial_response[ii,] <- sweep(s[ii,,drop=FALSE]%*%coefficient[,,g],2,group_response_mean[g,],'+')
  }
  u <- eigen(crossprod(molecular_response),symmetric=TRUE)$vectors[,1:decoder_rank,drop=FALSE]
  u <- sweep(u,2,sqrt(colMeans((molecular_response%*%u)^2)),'/')
  obs <- merge(ss,molecular[,.(depmap_id,drug,molecular_row=.I)],by=c('depmap_id','drug'))
  geometry <- function(y,weight,group) {
    g <- match(group,unique(group));gw <- as.numeric(rowsum(weight,g,reorder=FALSE))
    ym <- as.numeric(rowsum(y*weight,g,reorder=FALSE))/gw
    list(y=y,weight=weight,group=g,group_weight=gw,within_scale=max(sum(weight*(y-ym[g])^2),1e-8),absolute_scale=max(sum(weight*y^2),1e-8))
  }
  empty_map <- sparseMatrix(i=integer(),j=integer(),dims=c(0L,nquery*4L))
  list(s=s,initial_response=initial_response,initial_coefficient=coefficient,
    group_response_mean=group_response_mean,group_molecular_rows=group_molecular_rows,
    query=query,nquery=nquery,group_rows=group_rows,ngroup=2L,u=u,
    molecular_response=molecular_response,molecular_index=molecular_index,
    observed_response=molecular_response[obs$molecular_row,,drop=FALSE],
    ms=sparseMatrix(i=seq_len(nrow(ss)),j=match(key(ss$depmap_id,ss$drug),query_key),x=1,dims=c(nrow(ss),nquery*4L)),
    ma=empty_map,mb=empty_map,mp=sparseMatrix(i=integer(),j=integer(),dims=c(0L,0L)),
    pair_points=data.table(a=integer(),b=integer(),assay=integer()),ns=nrow(ss),np=0L,metadata=ss,
    function_geometry=geometry(ss$y,macro_weights(ss$drug),ss$drug),
    observed_geometry=geometry(obs$y,macro_weights(obs$drug),obs$drug),observed_cells=unique(obs$depmap_id),
    context=list(response_sd=response_sd,baseline_mean=baseline_mean,baseline_sd=baseline_sd,
      baseline_projection=baseline_projection,baseline_projection_scale=baseline_projection_scale,
      group_baseline_mean=group_baseline_mean,excluded_cells=excluded,reference_cells=ref_ids),
    decoder_rms=pmax(sqrt(colMeans(initial_response^2)),1e-6),
    nonlinear_rms=pmax(sqrt(colMeans(phi(initial_response,u)^2)),1e-6),interaction_rms=rep(1,interaction_features))
}
phi <- function(x, u) {
  latent <- tanh(x %*% u)
  cbind(tanh(x), x^2 / (1 + x^2),
    latent[, latent_pairs[, 1], drop = FALSE] * latent[, latent_pairs[, 2], drop = FALSE])
}
decoder <- function(x, par, assay = NULL) {
  z <- x %*% par$u
  h <- phi(x, par$u)
  common <- as.numeric(x %*% par$w + h %*% par$a)
  if (is.null(assay)) {
    output <- matrix(common, nrow(x), length(assays))
    output[, -1L] <- output[, -1L] + x %*% par$deviation
  } else {
    output <- common
    for (j in 2:length(assays)) {
      ii <- which(assay == j)
      if (length(ii)) output[ii] <- output[ii] + as.numeric(x[ii, , drop = FALSE] %*% par$deviation[, j - 1L])
    }
  }
  list(output = output, z = z, h = h, x = x)
}
decoder_back <- function(node, gradient, par, assay = NULL) {
  if (is.null(assay)) {
    g <- rowSums(gradient)
    gd <- gradient[, -1L, drop = FALSE]
  } else {
    g <- as.numeric(gradient)
    gd <- vapply(2:length(assays), function(j) g * as.numeric(assay == j), numeric(length(g)))
  }
  dz <- matrix(0, nrow(node$x), decoder_rank)
  latent <- tanh(node$z)
  for (j in seq_len(nrow(latent_pairs))) {
    a <- latent_pairs[j, 1]
    b <- latent_pairs[j, 2]
    dz[, a] <- dz[, a] + g * par$a[80L + j] * latent[, b] * (1 - latent[, a]^2)
    dz[, b] <- dz[, b] + g * par$a[80L + j] * latent[, a] * (1 - latent[, b]^2)
  }
  dx <- g * (matrix(par$w, nrow(node$x), 40L, byrow = TRUE) +
    matrix(par$a[1:40], nrow(node$x), 40L, byrow = TRUE) * (1 - tanh(node$x)^2) +
    matrix(par$a[40L + 1:40], nrow(node$x), 40L, byrow = TRUE) * (2 * node$x / (1 + node$x^2)^2))
  list(x = dx + dz %*% t(par$u) + gd %*% t(par$deviation),
    w = as.numeric(crossprod(node$x, g)), a = as.numeric(crossprod(node$h, g)),
    u = crossprod(node$x, dz), deviation = crossprod(node$x, gd))
}
function_loss <- function(predicted, geometry) {
  error <- predicted - geometry$y
  mean_error <- as.numeric(rowsum(error * geometry$weight, geometry$group, reorder = FALSE)) / geometry$group_weight
  centered <- error - mean_error[geometry$group]
  value <- within_function_weight * sum(geometry$weight * centered^2) / geometry$within_scale +
    (1 - within_function_weight) * sum(geometry$weight * error^2) / geometry$absolute_scale
  gradient <- 2 * geometry$weight * (within_function_weight * centered / geometry$within_scale +
    (1 - within_function_weight) * error / geometry$absolute_scale)
  list(value = value, gradient = gradient)
}

# Both inputs are retained. Products are symmetric and vanish at a zero input;
# no drug, pair identity, baseline, or independent mean predictor enters J.
pair_features <- function(a, b, u) {
  ha <- tanh(a %*% u)
  hb <- tanh(b %*% u)
  h <- matrix(0, nrow(a), nrow(latent_pairs))
  for (j in seq_len(nrow(latent_pairs))) {
    p <- latent_pairs[j, 1]
    q <- latent_pairs[j, 2]
    h[, j] <- if (p == q) ha[, p] * hb[, p] else
      (ha[, p] * hb[, q] + ha[, q] * hb[, p]) / sqrt(2)
  }
  cbind(tanh(a) * tanh(b), h)
}
pair_back <- function(a, b, gradient, par) {
  g <- as.numeric(gradient)
  ha <- tanh(a %*% par$u)
  hb <- tanh(b %*% par$u)
  ga <- gb <- matrix(0, nrow(a), decoder_rank)
  for (j in seq_len(nrow(latent_pairs))) {
    p <- latent_pairs[j, 1]
    q <- latent_pairs[j, 2]
    c <- g * par$interaction[40L + j]
    if (p == q) {
      ga[, p] <- ga[, p] + c * hb[, p]
      gb[, p] <- gb[, p] + c * ha[, p]
    } else {
      ga[, p] <- ga[, p] + c * hb[, q] / sqrt(2)
      ga[, q] <- ga[, q] + c * hb[, p] / sqrt(2)
      gb[, p] <- gb[, p] + c * ha[, q] / sqrt(2)
      gb[, q] <- gb[, q] + c * ha[, p] / sqrt(2)
    }
  }
  a_gradient <- g * sweep((1 - tanh(a)^2) * tanh(b), 2L, par$interaction[1:40], '*') +
    (ga * (1 - ha^2)) %*% t(par$u)
  b_gradient <- g * sweep((1 - tanh(b)^2) * tanh(a), 2L, par$interaction[1:40], '*') +
    (gb * (1 - hb^2)) %*% t(par$u)
  list(a = a_gradient, b = b_gradient,
    interaction = as.numeric(crossprod(pair_features(a, b, par$u), g)))
}
decoder_vector <- function(par) c(par$w, par$a, par$deviation, par$interaction)
decoder_unpack <- function(theta, u) {
  list(w = theta[1:40], a = theta[40 + seq_len(nonlinear_features)], u = u,
    deviation = matrix(theta[40L + nonlinear_features + seq_len(120L)], 40L, 3L),
    interaction = tail(theta, interaction_features))
}
ndecoder <- 40L + nonlinear_features + 120L + interaction_features

# Warm start uses one shared regularized decoder. Measured and inferred RNA
# supervision have equal total weights, with within-treatment contrasts weighted
# more than absolute effects. No context-specific calibration heads are fitted.
initialize_decoder <- function(p, family) {
  par <- list(w = rep(0, 40L), a = rep(0, nonlinear_features), u = p$u,
    deviation = matrix(0, 40L, 3L), interaction = rep(0, interaction_features))
  basis <- cbind(p$initial_response, phi(p$initial_response, par$u))
  native_design <- lapply(seq_along(assays), function(a)
    cbind(basis, do.call(cbind, lapply(2:length(assays), function(j)
      p$initial_response * as.numeric(j == a)))))
  native_design <- cbind(do.call(rbind, native_design), matrix(0, p$nquery * length(assays), interaction_features))
  pair_a <- p$initial_response[p$pair_points$a, , drop = FALSE]
  pair_b <- p$initial_response[p$pair_points$b, , drop = FALSE]
  pair_design <- cbind(matrix(0, nrow(pair_a), ndecoder - interaction_features),
    pair_features(pair_a, pair_b, par$u))
  xs <- as.matrix(p$ms %*% native_design)
  xa <- as.matrix(p$ma %*% native_design)
  xb <- as.matrix(p$mb %*% native_design)
  xj <- as.matrix(p$mp %*% pair_design)
  x <- rbind(xs, xa + xb + xj, xj)
  observed_x <- cbind(p$observed_response, phi(p$observed_response, par$u),
    matrix(0, nrow(p$observed_response), 120L + interaction_features))
  views <- function(x, geo) {
    xm <- rowsum(x * geo$weight, geo$group, reorder = FALSE) / geo$group_weight
    ym <- as.numeric(rowsum(geo$y * geo$weight, geo$group, reorder = FALSE)) / geo$group_weight
    list(x = rbind(x - xm[geo$group, , drop = FALSE], x), y = c(geo$y - ym[geo$group], geo$y),
      weight = c(geo$weight * within_function_weight / geo$within_scale,
        geo$weight * (1 - within_function_weight) / geo$absolute_scale))
  }
  observed <- views(observed_x, p$observed_geometry)
  inferred <- views(x, p$function_geometry)
  fitting <- list(x = rbind(inferred$x, observed$x), y = c(inferred$y, observed$y),
    weight = c(inferred$weight, observed$weight * observed_function_weight))
  fit <- glmnet(fitting$x, fitting$y, weights = fitting$weight, alpha = 0,
    lambda = .003, intercept = FALSE, standardize = TRUE, thresh = 1e-7, maxit = 100000L)
  stopifnot(fit$jerr == 0L)
  coefficient <- as.numeric(coef(fit))[-1L]
  par$w <- coefficient[seq_len(40L)]
  par$a <- coefficient[40L + seq_len(nonlinear_features)]
  par$deviation <- matrix(coefficient[40L + nonlinear_features + seq_len(120L)], 40L, 3L)
  par$interaction <- tail(coefficient, interaction_features)
  par
}

fit_model <- function(p, family, job_tag, check_gradient = FALSE) {
  linked <- family == families[2]
  par0 <- initialize_decoder(p, family)
  initial <- decoder_vector(par0)
  na <- baseline_rank * 40L * p$ngroup
  if (linked) initial <- c(initial, as.numeric(p$initial_coefficient))
  evaluations <- 0L
  cache_theta <- cache_result <- NULL
  objective <- function(theta) {
    if (!is.null(cache_theta) && identical(theta, cache_theta)) return(cache_result)
    par <- decoder_unpack(theta[seq_len(ndecoder)], p$u)
    if (linked) {
      coefficient <- array(theta[ndecoder + seq_len(na)], c(baseline_rank, 40L, p$ngroup))
      response <- matrix(0, p$nquery, 40L)
      for (g in seq_len(p$ngroup)) {
        ii <- p$group_rows[[g]]
        response[ii, ] <- sweep(p$s[ii, , drop = FALSE] %*% coefficient[, , g],
          2L, p$group_response_mean[g, ], '+')
      }
    } else response <- p$initial_response
    n <- decoder(response, par)
    pair_a <- response[p$pair_points$a, , drop = FALSE]
    pair_b <- response[p$pair_points$b, , drop = FALSE]
    pair_h <- pair_features(pair_a, pair_b, par$u)
    single_prediction <- as.numeric(p$ms %*% as.vector(n$output))
    a_prediction <- as.numeric(p$ma %*% as.vector(n$output))
    b_prediction <- as.numeric(p$mb %*% as.vector(n$output))
    interaction_prediction <- as.numeric(p$mp %*% (pair_h %*% par$interaction))
    predicted <- c(single_prediction, a_prediction + b_prediction + interaction_prediction, interaction_prediction)
    fl <- function_loss(predicted, p$function_geometry)
    obs_node <- decoder(p$observed_response, par, rep(1L, nrow(p$observed_response)))
    ol <- function_loss(obs_node$output, p$observed_geometry)
    ns <- p$ns
    np <- p$np
    gs <- fl$gradient[seq_len(ns)]
    gp <- fl$gradient[ns + seq_len(np)]
    gi <- fl$gradient[ns + np + seq_len(np)]
    native_gradient <- as.numeric(crossprod(p$ms, gs) + crossprod(p$ma, gp) + crossprod(p$mb, gp))
    native_gradient <- matrix(native_gradient, p$nquery, length(assays))
    pair_gradient <- as.numeric(crossprod(p$mp, gp + gi))
    ng <- decoder_back(n, native_gradient, par)
    pg <- pair_back(pair_a, pair_b, pair_gradient, par)
    og <- decoder_back(obs_node, ol$gradient * observed_function_weight, par, rep(1L, nrow(p$observed_response)))
    grad_decoder <- list(w = ng$w + og$w, a = ng$a + og$a,
      deviation = ng$deviation + og$deviation, interaction = pg$interaction)
    penalty_value <- decoder_penalty * (sum((par$w * p$decoder_rms)^2) + sum((par$a * p$nonlinear_rms)^2) +
      8 * sum((par$deviation * p$decoder_rms)^2)) +
      interaction_penalty * sum((par$interaction * p$interaction_rms)^2)
    grad_decoder$w <- grad_decoder$w + 2 * decoder_penalty * par$w * p$decoder_rms^2
    grad_decoder$a <- grad_decoder$a + 2 * decoder_penalty * par$a * p$nonlinear_rms^2
    grad_decoder$deviation <- grad_decoder$deviation + 16 * decoder_penalty * par$deviation * p$decoder_rms^2
    grad_decoder$interaction <- grad_decoder$interaction +
      2 * interaction_penalty * par$interaction * p$interaction_rms^2
    value <- fl$value + observed_function_weight * ol$value + penalty_value
    gradient <- decoder_vector(grad_decoder)
    if (linked) {
      grad_response <- ng$x
      for (which_input in c('a', 'b')) {
        index <- p$pair_points[[which_input]]
        contribution <- rowsum(pg[[which_input]], index, reorder = FALSE)
        ii <- as.integer(rownames(contribution))
        grad_response[ii, ] <- grad_response[ii, , drop = FALSE] + contribution
      }
      grad_coefficient <- coefficient * (2 * molecular_weight * molecular_ridge / (p$ngroup * 40L))
      for (g in seq_len(p$ngroup)) {
        ii <- p$group_rows[[g]]
        mi <- p$group_molecular_rows[[g]]
        qi <- p$molecular_index[mi]
        molecular_error <- response[qi, , drop = FALSE] - p$molecular_response[mi, , drop = FALSE]
        value <- value + molecular_weight / p$ngroup *
          (mean(molecular_error^2) + molecular_ridge * sum(coefficient[, , g]^2) / 40L)
        grad_response[qi, ] <- grad_response[qi, , drop = FALSE] +
          2 * molecular_weight * molecular_error / (p$ngroup * length(molecular_error))
        grad_coefficient[, , g] <- grad_coefficient[, , g] +
          crossprod(p$s[ii, , drop = FALSE], grad_response[ii, , drop = FALSE])
      }
      gradient <- c(gradient, grad_coefficient)
    }
    stopifnot(is.finite(value), all(is.finite(gradient)))
    evaluations <<- evaluations + 1L
    cache_theta <<- theta
    cache_result <<- list(value = value, gradient = gradient)
    cache_result
  }
  if (check_gradient && linked) {
    probe <- initial
    probe[ndecoder + 1L] <- .013
    positions <- unique(c(1L, 41L, 49L, 100L, ndecoder - c(0L, 2L, 10L),
      ndecoder + c(1L, 40L, 159L, 2301L), length(probe) - c(0L, 39L, 120L)))
    positions <- positions[positions <= length(probe)]
    analytical <- objective(probe)$gradient[positions]
    numerical <- vapply(positions, function(j) {
      step <- 1e-5 * max(1, abs(probe[j]))
      plus <- minus <- probe
      plus[j] <- plus[j] + step
      minus[j] <- minus[j] - step
      (objective(plus)$value - objective(minus)$value) / (2 * step)
    }, numeric(1))
    relative <- abs(numerical - analytical) / pmax(1, abs(numerical), abs(analytical))
    write_table(data.table(position = positions, analytical, numerical, relative), 'analytical_gradient_check')
    stopifnot(max(relative) < 1e-5)
    if (nrow(p$pair_points) >= 20L) {
    a <- p$initial_response[p$pair_points$a[1:20], , drop = FALSE]
    b <- p$initial_response[p$pair_points$b[1:20], , drop = FALSE]
    stopifnot(max(abs(pair_features(a, b, p$u) - pair_features(b, a, p$u))) < 1e-12,
      max(abs(pair_features(a, a * 0, p$u))) == 0,
      max(abs(pair_features(b * 0, b, p$u))) == 0)
    }
  }
  status(paste(job_tag, '| optimizing', family))
  # Program-coordinate units differ substantially. Parameter scaling conditions
  # the same objective; it does not change the representation or supervision.
  decoder_scale <- c(1 / p$decoder_rms, 1 / p$nonlinear_rms, rep(1 / p$decoder_rms, 3L), 1 / p$interaction_rms)
  parameter_scale <- if (linked) c(decoder_scale, rep(sqrt(p$ngroup / molecular_weight), na)) else decoder_scale
  opt <- optim(initial, function(theta) objective(theta)$value, function(theta) objective(theta)$gradient,
    method = 'L-BFGS-B', control = list(maxit = maximum_iterations, factr = 1e8,
      pgtol = 1e-4, lmm = 7, parscale = parameter_scale))
  if (opt$convergence == 1L && max(abs(objective(opt$par)$gradient * parameter_scale)) > 1e-3) {
    status(paste(job_tag, '| continuing optimization of', family))
    opt <- optim(opt$par, function(theta) objective(theta)$value, function(theta) objective(theta)$gradient,
      method = 'L-BFGS-B', control = list(maxit = maximum_iterations, factr = 1e8,
        pgtol = 1e-4, lmm = 7, parscale = parameter_scale))
  }
  par <- decoder_unpack(opt$par[seq_len(ndecoder)], p$u)
  model <- list(family = family, decoder = par, context = p$context,
    molecular_coefficient = if (linked) array(opt$par[ndecoder + seq_len(na)],
      c(baseline_rank, 40L, p$ngroup)) else p$initial_coefficient,
    initial_coefficient = p$initial_coefficient, group_response_mean = p$group_response_mean,
    selected_drugs = selected_drugs,
    functional_features = P, functional_baseline_inputs = FALSE,
    training_cells = unique(p$metadata$depmap_id),
    training_molecular_cells = p$context$reference_cells,
    optimization = list(convergence = opt$convergence, counts = opt$counts,
      maximum_gradient = max(abs(objective(opt$par)$gradient)),
      maximum_scaled_gradient = max(abs(objective(opt$par)$gradient * parameter_scale)), value = opt$value))
  status(paste(job_tag, '| finished', family, '| solver status', opt$convergence))
  model
}

# Inference remains baseline -> programs -> shared decoder, as in v82.
predict_displacement <- function(model,query) {
  s <- sweep(sweep(sweep(getmatrix(PREP$raw,query$depmap_id,P),2,model$context$baseline_mean),2,
    model$context$baseline_sd,'/')%*%model$context$baseline_projection,2,model$context$baseline_projection_scale,'/')
  response <- matrix(0,nrow(query),40)
  for(g in 1:2) {
    ii <- which(query$drug==selected_drugs[g])
    response[ii,] <- sweep(sweep(s[ii,,drop=FALSE],2,model$context$group_baseline_mean[g,])%*%
      model$molecular_coefficient[,,g],2,model$group_response_mean[g,],'+')
  }
  response
}
predictions <- measured_predictions <- audits <- list()
for(k in 1:4) {
  job <- if(k==4) 'final' else paste0('recipient_',k)
  held <- if(k==4) character() else names(recipient_fold)[recipient_fold==k]
  status(paste(job,'preparing'))
  p <- prepare(union(AU,held))
  test <- if(k==4) SINGLE else SINGLE[depmap_id %in% held]
  measured <- if(k==4) PREP$m else PREP$m[depmap_id %in% held]
  for(family in families) {
    model <- fit_model(p,family,job,check_gradient=(k==1))
    stopifnot(!any(model$training_cells %in% held),!any(model$training_molecular_cells %in% held),
      max(abs(model$decoder$interaction))<1e-10,max(abs(model$decoder$deviation))<1e-10)
    model$held_cells <- held
    model$native_doses <- 5
    model$combination_training_rows <- 0L
    qsave(model,file.path(root,'models',paste0(job,'_',gsub(' ','_',family),'.qs')))
    x <- predict_displacement(model,test)
    out <- test[,.(depmap_id,drug,dose,fold,observed_log2=y)]
    out[,`:=`(family=family,evaluation=if(k==4)'fitted' else 'held_cell',predicted_log2=decoder(x,model$decoder,rep(1L,.N))$output)]
    predictions[[paste(job,family)]] <- out
    mx <- sweep(as.matrix(measured[,..P]),2,model$context$response_sd,'/')
    mt <- merge(measured[,.(depmap_id,drug,dose,molecular_row=.I)],SINGLE[,.(depmap_id,drug,observed_log2=y)],by=c('depmap_id','drug'))
    mt[,`:=`(family=family,evaluation=if(k==4)'fitted' else 'held_cell',
      predicted_log2=decoder(mx[molecular_row,,drop=FALSE],model$decoder,rep(1L,.N))$output)]
    measured_predictions[[paste(job,family)]] <- mt
    audits[[paste(job,family)]] <- data.table(job,family,convergence=model$optimization$convergence,
      max_scaled_gradient=model$optimization$maximum_scaled_gradient,function_evaluations=model$optimization$counts[1])
    write_table(rbindlist(predictions),'single_predictions_frozen')
    write_table(rbindlist(measured_predictions),'measured_displacement_predictions_frozen')
    write_table(rbindlist(audits),'fit_audit')
  }
}
pred <- rbindlist(predictions)
metrics <- pred[,.(cells=.N,pearson_log2=pearson(observed_log2,predicted_log2),
  observed_sensitive=sum(observed_log2< -1),predicted_sensitive=sum(predicted_log2< -1)),by=.(drug,dose,family,evaluation)]
write_table(metrics,'single_correlations')
mp <- rbindlist(measured_predictions)
write_table(mp[,.(cells=.N,pearson_log2=pearson(observed_log2,predicted_log2)),by=.(drug,dose,family,evaluation)],'measured_displacement_correlations')
print(metrics[evaluation=='held_cell'])
held <- pred[evaluation=='held_cell']
held[,Drug:=ifelse(drug==selected_drugs[1],'Lapatinib','Crizotinib')]
lab <- metrics[evaluation=='held_cell']
lab[,`:=`(Drug=ifelse(drug==selected_drugs[1],'Lapatinib','Crizotinib'),label=sprintf('r = %.3f; n = %d',pearson_log2,cells))]
plot <- ggplot(held,aes(observed_log2,predicted_log2))+
  geom_point(aes(color=factor(fold)),size=.9,alpha=.65)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
  geom_label(data=lab,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.1,size=3,label.size=0)+
  facet_wrap(~family+Drug,ncol=2,scales='free')+
  theme_bw(base_size=10)+theme(aspect.ratio=1,legend.position='bottom')+
  labs(title='Lapatinib and crizotinib at 5 µM: v82 joint molecular bridge',
    x='Observed PRISM log2 abundance / vehicle',y='Held-cell predicted log2 abundance / vehicle',color='Held-cell fold',
    caption='Whole cell lines excluded from molecular and functional fitting. Fixed v82 hyperparameters.\nOnly these two drugs at 5 µM; no PLS, baseline functional shortcut, or combination training.')
ggsave(file.path(root,'figures/held_cell_predictions.pdf'),plot,width=8,height=8,device=cairo_pdf)
ggsave(file.path(root,'figures/held_cell_predictions.png'),plot,width=8,height=8,dpi=130)
writeLines(c('v82 fixed decoder, loss, penalties and baseline/program architecture; two drugs, one dose.',
  'Zero supported combination rows at 5+5 uM. No combination or HGF prediction is supported by this fit.',
  'Baseline rank 8; fixed molecular eigenvectors rank 4 in nonlinear decoder. No PLS.',
  'Mean molecular responses fixed; only baseline-conditioned slopes are jointly adjusted.',
  'Functional loss: 80% within-treatment variation, 20% absolute endpoint. No mean-only candidate or drug intercept.',
  'PRISM interpolated in log-dose to 5 uM; no extrapolation. Same original recipient folds, extended for new cells.',
  capture.output(sessionInfo())),file.path(root,'RUN_DETAILS.txt'))
status('Complete')

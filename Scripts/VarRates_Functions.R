
require(ape)
require(phytools)

###############################################################################

# Generate a control file for fabric/VarRates analysis that tracks the ancestral trait values
# for all nodes of a tree

BTControl <- function(phy, burnin, iterations, beta.prior, model=c("Fabric","VarRates","BM","Delta","Kappa","OU"), 
  traits=NULL, sigma=NULL, stones=NULL, file.name="BayesTraits_Control.txt", fabric.limits){
  
  # establish the control file and basic parameters
  writeLines("7", file.name) # independent contrasts
  control <- file(file.name, "a") # Open connection to append
  writeLines("2", control) # MCMC
  # if model is fabric or VR, specify that
  if(model %in% c("Fabric","VarRates","Delta","Kappa","OU")){writeLines(model, control)} # model of interest (e.g. 'Fabric', 'VarRates')
  # if we are using the fabric or a derivative, write model specifics
  if(model=="Fabric"){writeLines(paste("Prior", beta.prior), control)}
  if(model=="Fabric" && fabric.limits=="noV"){writeLines(paste("RJThreshold Node -9999"), control)}
  if(model=="Fabric" && fabric.limits=="noB"){writeLines(paste("RJThreshold LandscapeBL -9999"), control)}
  if(model=="Fabric" && fabric.limits=="noBnoV"){writeLines(paste("RJThreshold Node -9999"), control);writeLines(paste("RJThreshold LandscapeBL -9999"), control)}


  # if we have traits and want to specify priors based on that information
  if(model %in% c("Fabric","VarRates")){
    if(!is.null(traits)){
      # fit a constant rate BM model
      bm.fit <- geiger::fitContinuous(phy, traits, model="BM")
      # write the Alpha prior to file
      writeLines(paste("Prior", "Alpha-1", "normal", round(bm.fit$opt$z0,6), 0.75), control)
      if(sigma=="estimate"){
        writeLines(paste("Prior", "Sigma-1", "gamma", 2, round(bm.fit$opt$sigsq,6)), control) # visualize here: https://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html a=2 b=ML fit
      }
    }
  }

  if(length(sigma)==2){writeLines(paste("Prior", paste0("Sigma-",j), "uniform", sigma[1], sigma[2]), control)}

  # add 'tags' to keep track of ancestral trait values
  for(k in (Ntip(phy)+1):(Nnode(phy)+Ntip(phy))){
    tip.nos <- getDescendants(phy, k)
    tip.nos <- tip.nos[tip.nos<=Ntip(phy)]
    tip.nos <- phy$tip.label[tip.nos]
    writeLines(paste("AddTag", paste0("Tag0",k), paste(tip.nos, collapse=" ")), control)
    writeLines(paste("AddMRCA", paste0("n",k), paste0("Tag0",k)), control)
  }

  # add stepping stone sampler if you want to estimate the marginal likelihood
  if(!is.null(stones)){writeLines(paste("Stones", stones[1], stones[2]), control)}
  
  # specify the burnin, iterations, and run commands
  writeLines(paste("Burnin", burnin), control)
  writeLines(paste("Iterations", iterations), control)
  writeLines("Run", control)
  close(control)                        #Close connection
}
# NOTE: iterations and burnin must be in quotes to avoid being translated to scientific notation


###############################################################################


# Make a shell script to run lots of BayesTraits runs

BT.shell.old <- function(in.files, tree.name, dir.path){
  curr.dir <- getwd()
  
  if(is.null(dir.path)){
    for (k in 1:length(in.files)){
      curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4", 
                         paste0(getwd(),"/",tree.name), 
                         paste0(getwd(),"/",in.files[[k]]),
                         "<", 
                         paste0(getwd(),"/fabric_Control.txt"))
      write.table(curr.call, file=paste0(curr.dir,"/","fabric_Shell.sh"), append=T, row.names=F, col.names=F, quote=F)
    }
  }
  if(!is.null(dir.path)){
    for (k in 1:length(in.files)){
      curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4", 
                         tree.name, 
                         in.files[[k]],
                         "<", 
                         paste0(getwd(),"/Control.txt"))
      write.table(curr.call, file=paste0(curr.dir,"/","BTShell.sh"), append=T, row.names=F, col.names=F, quote=F)
    }    
  }
}

# Might be smart to flip the for-loop around so that it does it by model and not by trait, but whatever
# the drawback is that this only accepts a single tree file so I'll make a separate version below
BT.shell <- function(in.files, ctl.files, tree.name){
  curr.dir <- getwd()
  
  tfs <- sapply(in.files, function(x) strex::str_after_last(x, "/"))
  tnames <- sapply(tfs, function(x) strex::str_before_first(x, "_TRAIT")); names(tnames) <- NULL
  
  mds <- sapply(ctl.files, function(x) strex::str_after_last(x,"/"))
  mnames <- unique(sapply(mds, function(x) strex::str_before_first(x, "_Control"))); names(mnames) <- NULL
  
  for (k in 1:length(mnames)){
    model.in <- ctl.files[grep(mnames[[k]], ctl.files)]
    for (j in 1:length(model.in)){
      curr.trait <- tnames[which(sapply(tnames, function(x) grepl(x, model.in[[j]]))=="TRUE")]
      curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4",
                         paste0(getwd(),"/",tree.name),
                         in.files[grep(curr.trait, in.files)],
                         "<",
                         model.in[[j]])
      write.table(curr.call, file=paste0(curr.dir,"/SHELL_",mnames[[k]],".sh"), append=T, row.names=F, col.names=F, quote=F)
    }
  }
  
#  for (k in 1:length(tnames)){
#    ctl.in <- ctl.files[grep(tnames[[k]], ctl.files)]
#    for (j in 1:length(ctl.in)){
#      curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4",
#                         paste0(getwd(),"/",tree.name),
#                         in.files[[k]],
#                         "<",
#                         ctl.in[[j]])
#      write.table(curr.call, file=paste0(curr.dir,"/SHELL_",tnames[[k]],".sh"), append=T, row.names=F, col.names=F, quote=F)
#    }
#  }
}

# Run the fabric_Shell.sh file in parallel with GNU Parallel
# $ parallel -j [#parallel] --bar :::: [shell_file.sh]

# This version is very similar to the above, but matches a specific tree to the data
BT.shell2 <- function(in.files, ctl.files, tr.files){
  curr.dir <- getwd()
  
  tfs <- sapply(in.files, function(x) strex::str_after_last(x, "/"))
  tnames <- sapply(tfs, function(x) strex::str_before_first(x, "_TRAIT")); names(tnames) <- NULL
  
  mds <- sapply(ctl.files, function(x) strex::str_after_last(x,"/"))
  mnames <- unique(sapply(mds, function(x) strex::str_before_first(x, "_Control"))); names(mnames) <- NULL
  

  for (k in 1:length(mnames)){
    model.in <- ctl.files[grep(mnames[[k]], ctl.files)]
    for (j in 1:length(model.in)){
      curr.trait <- tnames[which(sapply(tnames, function(x) grepl(x, model.in[[j]]))=="TRUE")]
      curr.tree  <- tr.files[grep(curr.trait, tr.files)]
      curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4",
                         curr.tree,
                         in.files[grep(curr.trait, in.files)],
                         "<",
                         model.in[[j]])
      write.table(curr.call, file=paste0(curr.dir,"/SHELL_",mnames[[k]],".sh"), append=T, row.names=F, col.names=F, quote=F)
    }
  }
}

# this version is the same as the original BT.shell() function but works with traits that are numbered
# e.g. "peak1" "peak2" ... "peak10". This otherwise causes a problem with grep'ing the trait names
BT.shell.numbers <- function(in.files, ctl.files, tree.name){
  curr.dir <- getwd()
  
  tfs <- sapply(in.files, function(x) strex::str_after_last(x, "/"))
  tnames <- sapply(tfs, function(x) strex::str_before_first(x, "_TRAIT")); names(tnames) <- NULL
  
  mds <- sapply(ctl.files, function(x) strex::str_after_last(x,"/"))
  mnames <- unique(sapply(mds, function(x) strex::str_before_first(x, "_Control"))); names(mnames) <- NULL
  
  for (k in 1:length(mnames)){
    model.in <- ctl.files[grep(mnames[[k]], ctl.files)]
    for (j in 1:length(model.in)){
      #curr.trait <- tnames[which(sapply(tnames, function(x) grepl(x, model.in[[j]]))=="TRUE")]
      curr.trait <- strex::str_before_first(strex::str_after_last(strex::str_after_last(model.in[[j]],"/"),"_"),".txt")
      curr.call <- paste("/Applications/BayesTraitsV4.1.2/BayesTraitsV4",
                         paste0(getwd(),"/",tree.name),
                         in.files[grep(paste0(curr.trait,"_TRAIT"), in.files)],
                         "<",
                         model.in[[j]])
      write.table(curr.call, file=paste0(curr.dir,"/SHELL_",mnames[[k]],".sh"), append=T, row.names=F, col.names=F, quote=F)
    }
  }
}

# And a shell script to run the VarRates PostProcessor
PP.shell <- function(in.files, out.path){
  curr.dir <- getwd()
  
  tfs <- sapply(in.files, function(x) strex::str_after_last(x, "/"))
  tnames <- sapply(tfs, function(x) strex::str_before_first(x, "_TRAIT")); names(tnames) <- NULL
  
  for(k in 1:length(tnames)){
    curr.trait <- tnames[[k]]
    curr.call <- paste("/Applications/BayesTraitsV4/PPPostProcess/PPPostProcess",
                       in.files[grep(curr.trait,in.files)],
                       ">",
                       paste0(out.path,"/",curr.trait,"_Results.txt"))
    write.table(curr.call, file=paste0(curr.dir,"/SHELL_PP.sh"), append=T, row.names=F, col.names=F, quote=F)
  }
}


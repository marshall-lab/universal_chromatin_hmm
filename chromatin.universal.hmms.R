# Chromatin state processing script (part 1: HMM fitting)
# Copyright © 2014-19, Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

# NB: The viterbi function will only work when R is run with the --max-ppsize=500000 command-line option (no rstudio for you!)

### version
my.version = "0.9.10 (2019-11-25)"
cat(paste("Universal Chromatin HMM script
version ",my.version,"

** Please run with --max-ppsize=500000 to avoid protection stack overflow errors

",sep=""))

### CLI options
if (!dir.exists("~/.config/chromatin.universal.rscripts")) {
  cat("Creating config dir at ~/.config/chromatin.universal.rscripts/hmm.config\n")
  dir.create("~/.config/chromatin.universal.rscripts")
}

op.args = list(
  "name" = "test",
  "genes.file" = "/mnt/genomes/dmel_release/DmR6/DmR6.genes.gff",
  "tfs" = "/mnt/genomes/dmel_release/dros.transcription.factors.txt",
  "states.matrix" = "/mnt/genomes/state.means.csv",
  "nStates" = 25,
  "rhmm.iter" = 2000,
  "rhmm.iter.init" = 10,
  "gff.path" = getwd(),
  "expected.state.clusters" = 8,
  "mc.cores" = 8,
  "chr.model" = "all",
  "fit.hmm" = T,
  "fit.viterbi" = T,
  "gene.expression" = T,
  "load.mcore.seeds" = "",
  "save.defaults" = F
)

op.args.defaults.to.save = c("genes.file","tfs","states.matrix","nStates","rhmm.iter","rhmm.iter.init","expected.state.clusters","mc.cores","chr.model")

if (file.exists("~/.config/chromatin.universal.rscripts/hmm.defaults")) {
  cat("Loading saved defaults ...\n")
  load("~/.config/chromatin.universal.rscripts/hmm.defaults")
  for (n in names(op.args.saved)) {
    op.args[[n]] = op.args.saved[[n]]
  }
}

op.notes = list(
  "name" = "Sample name (no spaces)",
  "genes.file" = "Genes file in GFF format",
  "tfs" = "List of all transcription factors",
  "states.matrix" = "Prediction matrix for chromatin state calling",
  "nStates" = "Number of states to fit HMM",
  "rhmm.iter" = "Number of random starts for HMM fitting",
  "rhmm.iter.init" = "Number of iterations to run per random start",
  "gff.path" = "Path to binding tracks in bedgraph or GFF format",
  "expected.state.clusters" = "Number of clusters to split states tree on",
  "mc.cores" = "Number of cores to use for HMM fitting",
  "chr.model" = "Chromosomes to fit HMM model to (separate by commas, no spaces)",
  "fit.hmm" = "Whether to fit HMM",
  "fit.viterbi" = "Whether to fit Viterbi path",
  "gene.expression" = "Whether to profile gene expression (needs RNA pol binding profile)",
  "load.mcore.seeds" = "If specified, previously saved multicore seeds will be loaded from this file",
  "save.defaults" = "Save current config as default\n  (excludes gff.path, fit.hmm, fit.viterbi, gene.expression, load.mcore.seeds)"
)

### save random seed for future reproducibility
test = runif(1)
my.seed  = .Random.seed
write.table(my.seed,".randomseed")

### Read CLI options
input.args = commandArgs(trailingOnly = TRUE)

read.ops = function (x) {
  for (op in x) {
    op = gsub("^--","",op)
    y = unlist(strsplit(op,"="))

    if (y[1] == "help") {
      cat("Options:\n")
      for (n in names(op.args)) {
        cat(paste("  ",op.notes[[n]],":\n",sep=""))
        cat(paste("  --",n,"=",op.args[[n]],"\n\n",sep=""))
      }
      q()
    }

    if (!is.null(op.args[[ y[1] ]])) {
      op.args[[ y[1] ]] <<- y[2]
    } else {
      cat("Error: Option",y[1],"not recognised ...\n")
      q()
    }
  }
}

write.ops = function () {
  oldw = getOption("warn")
  options(warn = -1)
  out.df = data.frame(option="version",value=my.version)
  for (n in names(op.args)) {
    v <<- op.args[[n]]
    df.line = data.frame(
      option=n,
      value=v
    )
    out.df = rbind(out.df, df.line)
  }
  write.table(out.df,"input.args.txt",row.names=F)
  options(warn = oldw)
}

read.gff = function (x,name="score") {
  fn.ext = file_ext(x)
  
  if (grepl("gff",ignore.case=T,fn.ext)) {
	temp.data <- read.table(x,row.names=NULL)
  	if (ncol(temp.data) > 5) {
  	  # GFF
  	  trim.data = temp.data[,c(1,4,5,6)]
  	} else {
  		cat("Error: file does not appear to be in GFF format\n\n")
  		quit("no",1)
  	}
  } else if (grepl("bed",ignore.case=T,fn.ext)) {
  	temp.data = read.table(x,row.names=NULL,skip=1)
  	if (ncol(temp.data) == 4) {
  		# bedgraph
  		trim.data = temp.data
  	} else {
  		cat("Error: file does not appear to be in bedGraph format\n\n")
  		quit("no",1)
  	}
  } else {
  	cat("Error: input file does not appear to be in bedGraph or GFF format ...\n\n")
  	quit("no",1)
  }
  
  names(trim.data) = c("chr","start","end",name)
  trim.data$chr = gsub("^chr","",trim.data$chr,perl=T)
  
  return(trim.data)
}

build.dataframes = function () {
  # find bedgraphs in working dir
  bedgraphs <<- Sys.glob(paste(base.path,"*.bed*",sep="/"))
  if (length(bedgraphs)<1) {
    bedgraphs <<- Sys.glob(paste(base.path,"*.gff",sep="/"))
  }
  polii.name = bedgraphs[grep("II",bedgraphs)]
  
  # put PolII first if exists
  if (length(polii.name)) {
    bedgraphs[c(1,grep("II",bedgraphs))] <<- bedgraphs[c(grep("II",bedgraphs),1)]
  }
  
  # profiled proteins 
  prot.names <<- regmatches(bedgraphs,regexpr("(?<=/)(?!.*/).*?(?=(-|_))",bedgraphs,perl=T))
  
  cat("Building dataframes:\n")
  file = bedgraphs[1]
  prot = regmatches(file,regexpr("(?<=/)(?!.*/).*?(?=[-\\.])",file,perl=T))
  cat(paste("  reading",prot,"...\n"))
  data.all = read.gff(bedgraphs[1],tolower(prot))
  
  for (i in 2:length(bedgraphs)) {
    file = bedgraphs[i]
    prot = regmatches(file,regexpr("(?<=/)(?!.*/).*?(?=[-\\.])",file,perl=T))
    cat(paste("  reading",prot,"...\n"))
    tempin = read.gff(file,tolower(prot))
    data.all = merge(data.all,tempin,by=c("chr","start","end"))
  }
  
  # order by chr and fragment
  data.all = data.all[order(data.all$chr,data.all$start),]
  return(data.all)
}

order.clusters = function () {
  # establish order from clusering
  names(data.na.ex) = str_remove(names(data.na.ex),paste("\\.{0,1}",name,sep=""))
  corr = cor(data.na.ex[4:length(names(data.na.ex))])
  corr.hc = hclust(dist(corr))
  ggcorrplot(corr, hc.order=T, lab_size=2, outline.col="white")
  ggsave(paste("bedgraph_correlation_plot",name,"pdf",sep="."))
  
  names(data.na.ex) = str_remove(names(data.na.ex),"_.*")
  corr2 = cor(data.na.ex[4:length(names(data.na.ex))])
  ggcorrplot(corr2, hc.order=T, lab_size=2, outline.col="white")
  ggsave(paste("bedgraph_correlation_plot_minimal_names",name,"pdf",sep="."))
  
  return(corr.hc)
}

fit.hmm = function (data = data.na) {
  cat(paste("  Fitting",nStates,"states ...\n"))
  model = HMMFit(
      data.matrix(data[,4:ncol(data)]),
      nStates=nStates,
      control=list(
                   verbose=2,
                   nInit=rhmm.n.iter,
                   nIterInit=rhmm.n.iter.init
                   )
  )
  return(model)
}

fit.hmm.mc = function (data = data.na, cores = mc.cores, chrs = NULL) {
  # multi.core HMM to speed up model fitting
  
  # Allow restriction of training dataset to individual chromosome or chromosomes
  # Reducing the training set will decrease time to fit model (but use with some caution)
  model.input = data.frame()
  if (chrs == "all") {
    model.input = data
  } else {
    for (c in chrs) {
      model.input = rbind(model.input,data[data$chr == c,])
    }
  }
  
  cat(paste("  Fitting",nStates,"states ...\n"))
  
  # We need different random seeds for each thread, but want this operation to be reproducible
  seeds = sample(0:2147483647,mc.cores,replace=F)
  if (op.args[["load.mcore.seeds"]] != "") {
    if (file.exists(paste("../",op.args[["load.mcore.seeds"]],sep=""))) {
      # load previously saved random seeds
      cat("  Loading seeds file ...\n")
      seeds = read.table(paste("../",op.args[["load.mcore.seeds"]],sep=""))[[1]]
    } else {
      cat("  Error, cannot read seeds file.\n")
    }
  }
  write.table(seeds,"multicore.random.seeds",sep="\n",quote=F,row.names=F,col.names=F)
  
  # Fit multi-threaded
  mc.model = mclapply(seeds, function (x) {
      set.seed(x)
      HMMFit(
      data.matrix(model.input[,4:ncol(model.input)]),
      nStates=nStates,
      control=list(
                   verbose=2,
                   nInit=as.integer(rhmm.n.iter/length(seeds)),
                   nIterInit=rhmm.n.iter.init
                   )
      )},
      mc.cores=cores
  )
  
  hmm.best = vector()
  for (i in 1:length(mc.model)) {
    if (is.finite(mc.model[[i]]$LLH)) {
      hmm.best=c(hmm.best,mc.model[[i]]$LLH)
    }
  }
  
  cat("Model fitted.\n")
  print(hmm.best)
  print( grep(max(hmm.best),hmm.best) )
  
  return( mc.model[[ grep(max(hmm.best),hmm.best) ]] )
}

viterbi.path = function (data = data.na.ex, hmm.model = chrom.states[[name]]) {
  chr.names = unique(data$chr)
  vit.gff = data.frame()
  chr.num = 0
  for (cn in chr.names) {
    chr.num = chr.num + 1
    data.chr = na.exclude(data[data$chr == cn,])
    if (!nrow(data.chr)) {next}
    cat(paste(cn," (",chr.num,"/",length(chr.names),") ",sep=""))
    vit.states = viterbi(hmm.model,data.matrix(data.chr[,4:ncol(data.chr)]))
    vit.gff = rbind(vit.gff, data.frame(chr=data.chr$chr,start=data.chr$start,end=data.chr$end,state=vit.states$states))
  }
  return(vit.gff)
}

get.prot.order = function (prot.names=prot.names,corr.hc=corr.hc) {
  # order can go either way
  if (any(grepl("II", prot.names[corr.hc$order]) == T)) {
    # Pol II sample present
    cat("Found Pol II sample and reordering (Pol II comes first) ...\n")
    if (grep("II",prot.names[corr.hc$order]) > length(prot.names)/2) {
      prot.order = rev(corr.hc$order)
    } else {
      prot.order = corr.hc$order
    }
    # make sure polii always comes first
    prot.order[c(1,grep("^1$",prot.order))] = prot.order[c(grep("^1$",prot.order),1)]
  } else {
    prot.order = corr.hc$order
  }
  
  #cat("Prot names order:\n")
  #print(prot.names[prot.order])
  return(prot.order)
}

genomic.coverage = function (vit.gff = vit.gff) {
  vit.genomic = data.frame(stringsAsFactors=F)
  for (i in 1:nStates) {
    vit.genomic = rbind(vit.genomic,
                        data.frame(
                                   state=i,
                                   coverage=sum(subset(vit.gff,state==i)$end - subset(vit.gff,state==i)$start)
                                   )
                        )
  }
  write.table(vit.genomic,"genomic.coverage.txt",quote=F,row.names=F,sep="\t")
  return(vit.genomic)
}

cluster.transMat = function () {
  hclust(dist(chrom.states[[name]]$HMM$transMat),method='ward.D2')
}

plot.transMat = function () {
  pdf(paste("transMat.heatmap",name,nStates,".pdf",sep="."))
  named_tmat = chrom.states[[name]]$HMM$transMat
  rownames(named_tmat)=c(1:nStates) 
  colnames(named_tmat)=c(1:nStates)
  print(Heatmap(
      named_tmat,
      cluster_columns=hc,
      cluster_rows=hc,
      col=circlize::colorRamp2(breaks=c(0,0.5,1),colors=c(rgb(0.2,0.3,0.6),rgb(1,1,1),rgb(0.9,0.1,0.1))),
      column_dend_height=unit(20,"mm"),
      row_dend_width=unit(20,"mm"),
      name="P(trans)",
    ))
  dev.off()
}

plot.chrom.phylo = function (hc = hc) {
  pdf(paste("hclust.transMat",name,nStates,"pdf",sep="."))
  plot( as.phylo(hc), type="unrooted", lab4ut="axial")
  dev.off()
}

Mode = function(x) {
  ux = unique(x)
  tx = tabulate(match(x, ux))
  mx = ux[which.max(tx)]
  cx = tx[which.max(tx)]
  return(c(mx,cx,length(x)))
}

gene.exp.modal = function (input.df, genes.file=op.args[["genes.file"]], pad=0, subset.genes=T) {  
  cat("Reading genes data file ...\n")
  genes = read.table(genes.file, comment.char="#", sep="\t", quote="", stringsAsFactors=F)
  names(genes) =  c('chr','source','type','start','end','score','strand','c','details')
    
  # only subset if there is a type termed "gene"
  if ( (subset.genes) & (any(genes$type == 'gene'))) {
    genes = subset(genes, type=='gene')
  }
  
  total = length(genes$chr)
  
  genes$name = unlist(lapply(genes$details,function(x) str_extract(x,"(?<=Name=).*?(?=;)")))
  genes$FBgnID = unlist(lapply(genes$details,function(x) str_extract(x,"(?<=ID=).*?(?=;)")))
  genes = genes[,c('chr','start','end','strand','name','FBgnID')]
  
  avg.exp = data.frame()
  avg = vector(length=1)
  avg.exp = avg.exp[0,]
  
  count = 0
    
  # unroll chromosomes for speed:
  for (chromo in unique(genes$chr)) {
    input.chr = subset(input.df, chr==chromo)
        
    if (is.null(input.chr)) { next }
    
    genes.chr = subset(genes, chr==chromo)
    
    for (i in c(1:length(genes.chr$name))) {
      expr = data.frame(input.chr[ (input.chr$start <= genes.chr[i,"end"]) 
                                      & (input.chr$end >= genes.chr[i,"start"]), ] )
      if (length(expr[,1]) == 0) {next}
      
      # trim to gene boundaries ...
      expr$start[1] = genes.chr[i,"start"]-pad
      expr$end[length(expr[,1])] = genes.chr[i,"end"]+pad
      
      # gene length (not required here)
      gene.len = genes.chr[i,"end"]-genes.chr[i,"start"]
      
      # roll through each row
      avg = vector()
      for (j in c(1:nrow(expr))) {
        avg = append(avg,rep(expr[j,4],(expr$end[j]-expr$start[j])))
      }
      
      if (is.null(avg)) { cat("Error!",chromo,genes.chr[i,"start"],genes.chr[i,"end"],avg,"\n");next }
      
      m = Mode(avg)
      modeavg = m[1]
      mprop = m[2]/m[3]
      
      new.df = data.frame(name=as.character(genes.chr[i,"name"]),FBgnID=genes.chr[i,"FBgnID"], state=modeavg, prop=mprop, gene.size=gene.len)
      
      avg.exp = rbind(avg.exp,new.df);
      count = count+1
      
      if (count %% 50 == 0) {cat(paste(count,"/",total,"genes averaged ...\r"))}
    }
  }
  
  cat("\nAll done.\n\n");
  return(list("exp" = avg.exp, "count" = total));
}

gene.exp = function (input.df, buffer=0, iter=50000, debug=F, genes.file=op.args[["genes.file"]]) {
  cat("Reading genes data file ...\n")
  genes = read.table(genes.file, comment.char="#", sep="\t", quote="")
  names(genes) =  c('chr','source','type','start','end','score','strand','c','details')

  # only subset if there is a type termed "gene"
  if (any(genes$type == 'gene')) {
    genes = subset(genes, type=='gene')
  }

  genes$name = sapply(genes$details, FUN = function (x) {regmatches(x,gregexpr("(?<=Name=).*?(?=;)", x, perl=T))} )
  genes = genes[,c('chr','start','end','strand','name')]

  avg.exp = data.frame(input.df[1,c(4:(length(names(input.df))))])
  avg = vector(length=(length(names(input.df)) - 4))
  avg.exp = avg.exp[0,]

  ### FDR calcs ###
  # Method based off perl scripts by Tony Southall (TDS) as published in
  # Southall et al. (2013). Dev Cell, 26(1), 101–12. doi:10.1016/j.devcel.2013.05.020
  #
  # Significant modifications to the original methodology include:
  # * taking a linear regression of log data rather than trial-and-error curve fitting of non-log data
  # * using a linear regression for the final intercept value rather than using the average intercept value for all conditions
  # -- both of these should increase the accuracy of the final FDR value.

  input.len = length(input.df[,1])
  frag.samp = c(1,2,3,4,6,8,10,12,15)
  thres.samp = c(0.1,0.2,0.3,0.4,0.5,0.65,0.8,1.0,1.5,2.0)

  rand = list()
  for (thres in thres.samp) {
    cat(paste("  Calculating FDR for threshold",thres,"\n",sep=" "))

    # init vars
    for (f in frag.samp) {
      # e.g: frag 1, thres 0.2: rand[[thres.1.0.2]]
      rand[[paste("thres.",f,".",thres,sep="")]] = 0;
    }

    for (i in 1:iter) {
      if (i %% 200 == 0) {cat(paste("  iter",i,"\r"))}
      # get random sample for different fragment lengths

      rand.samp = list()
      for (f in frag.samp) {
        # Using the fourth column as we're only calculating FDR for one sample ...
        rand.samp[[paste("rand.",f,sep="")]] = mean(input.df[runif(f,1,input.len),4])
      }

      # count number of times exp > thres
      for (f in frag.samp) {
        if (rand.samp[[paste("rand.",f,sep="")]] > thres) {rand[[paste("thres.",f,".",thres,sep="")]] = rand[[paste("thres.",f,".",thres,sep="")]] + 1}
      }
    }
  }

  rand.fdr = list()
  for (thres in thres.samp) {
    for (f in frag.samp) {
      rand.fdr[[paste("thres.",f,".",thres,sep="")]] = rand[[paste("thres.",f,".",thres,sep="")]]/iter
    }
  }

  cat("Fitting curves ...\n")

  # curve fit: fdr vs thresholds
  var.thres = list()
  for (thres in thres.samp) {
    for (f in frag.samp) {
      var.thres[[paste("frags.",f,sep="")]] = append(var.thres[[paste("frags.",f,sep="")]], rand.fdr[[paste("thres.",f,".",thres,sep="")]])
    }
  }

  inf.log.lm = function (v) {
    non.inf = log(v) != -Inf
    ret = lm(log(v)[non.inf] ~ thres.samp[non.inf])
    return(ret)
  }

  # The relationship is exponential, so we need log data for a linear regression
  # (in R, linear regression is: y = lm$coefficients[[2]]x + lm$coefficients[[1]] ... )
  var.lm = list()
  for (f in frag.samp) {
    var.lm[[paste("frags.",f,sep="")]] = inf.log.lm(var.thres[[paste("frags.",f,sep="")]])
  }

  # ... and now we do a linear regression on the slopes and intercepts of our previous regressions
  #
  # (This is the clever bit, and it actually seems to work.  The correlation of slope to fragment size is linear ...
  # By doing this on the slope and intercept, we can now predict the FDR for any number of fragments with any expression.)
  slope = vector()
  for (f in frag.samp) {
    slope = append(slope, var.lm[[paste("frags.",f,sep="")]]$coefficients[[2]])
  }

  # slope regression predicts the average slope
  slope.lm = lm(slope ~ frag.samp)

  # TDS used an average intercept value for the intercept, however ...
  inter = vector()
  for (f in frag.samp) {
    inter = append(inter, var.lm[[paste("frags.",f,sep="")]]$coefficients[[1]])
  }

  # ... there's actually quite a bit of variation of the intercept with real data,
  # so we're going to do a lin regression on the intercept instead.
  #
  # (I'm not convinced it's a linear relationship.  But it's close to linear,
  # and will certainly perform better than taking the mean intercept ...)
  #
  # If you're interested, set the debug flag to TRUE and take a look at the plots generated below ...
  inter.lm = lm(inter ~ frag.samp)

  if (debug == T) {
    # plots for debugging/checking
    plot.debug = function (y,x,l,name="Debug plot") {
      plot(y ~ x)
      abline(l)

      lsum = summary(l)
      r2 = lsum$r.squared

      legend("topleft",legend=r2,bty='n')
      title(name)

      dev.copy(png,paste(name,".png",sep=""));dev.off()
    }

    plot.debug(slope,frag.samp,slope.lm,name="correlation of slope")
    plot.debug(inter,frag.samp,inter.lm,name="correlation of intercepts")
  }

  # ok, so putting that all together ...
  fdr = function (frags, expr) {
    inter.test = inter.lm$coefficients[[2]] * frags + inter.lm$coefficients[[1]]
    slope.test = slope.lm$coefficients[[2]] * frags + slope.lm$coefficients[[1]]

    fdr.out = exp(slope.test * expr + inter.test)
    return(fdr.out)
  }

  ### Gene expression values ###
  cat("Calculating gene values ...\n")

  count = 0

  # unroll chromosomes for speed:
  for (chromo in unique(genes$chr)) {
    input.chr = subset(input.df, chr==chromo)
    genes.chr = subset(genes, chr==chromo)
    for (i in 1:length(genes.chr$name)) {
      # Roll through each gene

      # Note: the original script calculated expression values for a table of proteins,
      # here we just limit ourselves to the FDR for the first column past "chr", "start" and "end"
      # so if you're thinking of looking at chromatin, for example, put PolII as column 4 in your table!

      gene.start = genes.chr[i,"start"] - buffer
      gene.end = genes.chr[i,"end"] + buffer

      gene.start = ifelse(gene.start < 1, 1, gene.start)

      # Create data frames for all gatc fragments covering current gene
      exp = data.frame(input.chr[ (input.chr$start <= gene.end)
                                 & (input.chr$end >= gene.start)
                                ,] )

      gatc.num = length(exp[,1])

      # skip if no gatc fragments cover gene :(
      if (gatc.num == 0) {next}

      # trim to gene boundaries ...
      exp$start[1] = gene.start
      exp$end[length(exp[,1])] = gene.end

      # gene length covered by gatc fragments
      len = sum(exp$end-exp$start)

      # calculate weighted score for each column (representing different proteins)
      for (j in 4:length(names(input.chr))) {
        avg[j] = (sum((exp$end-exp$start)*exp[j]))/len
      }

      # make data.frame of averages (to be appended to avg.exp)
      df = cbind(avg[1])
      for (k in 2:length(avg)) {
        df = cbind(df,avg[k])
      }
      df = cbind(df,gatc.num)

      # only fdr for first column for now ...
      gene.fdr = fdr(gatc.num,avg[4])
      df = cbind(df, gene.fdr)

      # append current gene to list
      avg.exp = rbind(avg.exp,data.frame(name=as.character(genes.chr[i,"name"]), df))
      count = count+1
      if (count %% 50 == 0) {cat(paste(count,"genes averaged ...\r"))}
    }
  }

  avg.exp = avg.exp[,c(1,5:(length(names(avg.exp))))]
  names(avg.exp) = c("name",names(input.df)[c(4:(length(names(input.df))))],"gatc.num","FDR")

  cat("\nAll done.\n\n")
  return(avg.exp)
}


generate.state.template = function () {
  ### single state template
  detected.clusters.out=""
  clustercut = cutree(hc,expected.state.clusters)
  
  for (k in 1:expected.state.clusters) {
    cluster.states = hc$order[(1:nStates)[clustercut[hc$order]==k]]
    detected.clusters.out = paste(paste(detected.clusters.out,collapse=""), "#", k,":",paste(cluster.states,collapse=","),"\n")
  }
  write.table(detected.clusters.out,"detected.state.clusters.txt",row.names=F,col.names=F,quote=F)
  
  # default state names and colours
  chrom.cols = list()
  chrom.cols[["TrxG"]]   = "#dc3912"
  chrom.cols[["TrxGR"]]  = "#860000"
  chrom.cols[["Yellow"]] = "#FF9900"
  chrom.cols[["YellowR"]]= "#9a5c02"
  chrom.cols[["PcGM"]]   = "#3366cc"
  chrom.cols[["PcGR"]]    = "#2050aa"
  chrom.cols[["HP1"]]    = "#4e9a06"
  chrom.cols[["Black"]]  = "#666666"
  
  state.means =  chrom.states[[name]]$HMM$distribution$mean

  test = data.frame(matrix(unlist(state.means), nrow=length(state.means), byrow=T))
  names(test) = prot.names
  
  cluster.means = data.frame()
  for (i in unique(clustercut)) {
    temp = as.data.frame(t(colMeans(test[grep(i,clustercut),])))
    cluster.means = rbind(cluster.means,temp)
  }
  
  # compatibilty tweaks
  names(cluster.means) = gsub(pattern="RPII18",replacement="polII",x=names(cluster.means),perl=T,ignore.case=T)
  names(cluster.means) = gsub(pattern="HP1$",replacement="HP1a",x=names(cluster.means),perl=T,ignore.case=T)
  names(cluster.means) = tolower(names(cluster.means))
  
  pred.chrom.groups = read.table(op.args[["states.matrix"]])
  row.names(pred.chrom.groups) = gsub(pattern="HP1$",replacement="HP1a",x=row.names(pred.chrom.groups),perl=T,ignore.case=T)
  row.names(pred.chrom.groups) = tolower(row.names(pred.chrom.groups))
  
  # make matrices compatible
  test.cm = cluster.means[ ,names(cluster.means) %in% row.names(pred.chrom.groups)]
  test.pcg = pred.chrom.groups[row.names(pred.chrom.groups) %in% names(test.cm), ]
  
  # predict states based on best correlation to established means
  cor.pred = cor(test.pcg, t(test.cm))
  pred.states = rownames(cor.pred)[apply(cor.pred, 2, which.max)]
  
  chrom.clusters = list()
  
  for (i in unique(clustercut)) {
    chrom.clusters[[ pred.states[i] ]] = c(chrom.clusters[[ pred.states[i] ]],hc$order[hc$order %in% grep(i,clustercut)])
  }
  
  names.list = ""
  states.list = ""
  cols.list = ""
  state.count = 0
  for (i in names(chrom.cols)) {
    if (!length(chrom.clusters[[i]]>0)) {
      # nothing to see here ... comment out this block
      names.list = paste(names.list,"#",sep="")
      states.list = paste(states.list, "#",sep="")
      cols.list = paste(cols.list, "#",sep="")
    } else {
      state.count=state.count+1
      for (j in chrom.clusters[[i]]) {
        hb.annot[j] <<- i
      }
      hb.cols[["states"]] <<- c(hb.cols[["states"]], c(i = chrom.cols[[i]]))
      names(hb.cols[["states"]])[state.count] <<- i
    }
      
  
    names.list = paste(names.list,paste('\t\"',i,'\",',"\n",sep=""),sep="")
    states.list = paste(states.list,paste('\tc(',paste(chrom.clusters[[i]],collapse=","),"),  #",i,"\n",sep=""),sep="")
    cols.list = paste(cols.list,paste('\t\"',chrom.cols[[i]],'\",  #',i,"\n",sep=""),sep="")
  }	

  single.state.template <- paste("
## Enter the target state and sample name
nStates =",nStates,"
name =\"",name,"\"

### assign names and colours to states, and simplify
# (if not simplifying, then just identify all states)

# state names
simp.names[[name]] = c(
",names.list,"
)

# state colours (generally leave as is)
simp.cols[[name]] = c(
",cols.list,"
)

# list of simplified states, with each element being a vector of the states in the simplified state
simp.states[[name]] = list(
",states.list,"
)
",sep="")
  
  write.table(single.state.template,paste("single.state",nStates,"r",sep="."),quote=F,row.names=F,col.names=F)
}

plot.universal.heatmap = function () {
  ### Universal heatmap
  state.mean.hmap <- data.frame(name=prot.names)
  for (state in 1:nStates) {
    state.mean.hmap[,paste("s",state)] <- chrom.states[[name]]$HMM$distribution$mean[[state]]
  }
  
  named.smh = state.mean.hmap[,2:ncol(state.mean.hmap)]
  rownames(named.smh) = prot.names
  annot_df = data.frame(
              genomic.cov=vit.genomic$coverage,
              num.genes=state.table.all.states$num.genes,
              tfs.pc=state.table.all.states$pc.tfs.in.state
              )
  
  col = list(
    genomic.cov = circlize::colorRamp2(c(0, max(na.omit(annot_df$genomic.cov))), 
                                         c("white", rgb(0.102,0.102,0.349))),
    num.genes = circlize::colorRamp2(c(0, max(na.omit(annot_df$num.genes))), 
                                         c("white", "purple")),
    tfs.pc		  = circlize::colorRamp2(c(0, max(na.omit(annot_df$tfs))), 
                                         c("white", "darkgreen"))
    )
  
  if (op.args[["gene.expression"]]) {
    annot_df$exp.pc = state.table.all.states$pc.of.state
    col[["exp.pc"]] = circlize::colorRamp2(c(0, max(na.omit(annot_df$exp.pc))), 
                                         c("white", "orange"))
  }
  
  ha = HeatmapAnnotation(df = annot_df, col = col, show_annotation_name = F)
  
  hb = HeatmapAnnotation(states = hb.annot,
    col = hb.cols
  )
  
  hc = hclust(
        d = dist(chrom.states[[name]]$HMM$transMat),
        method = "ward.D2"
        )
  
  data.min = min(named.smh)
  data.max = max(named.smh)
  
  data.mid = 0
  if (data.min>0) {
    data.mid = (data.min+data.max)/2
  }
  
  pdf.output.name = paste("Universal_heatmap",name,nStates,sep=".")
  pdf(paste(pdf.output.name,"pdf",sep="."),
      title=pdf.output.name,
      height=0.5*length(prot.names)+1,
      width=0.5*nStates+1
      )
  draw(Heatmap(
      as.matrix(named.smh),
      cluster_columns=hc,
      cluster_rows=F,
      row_order=prot.order,
      col=circlize::colorRamp2(c(data.min,data.mid,data.max),c(rgb(0.2,0.3,0.6),rgb(1,1,1),rgb(0.9,0.1,0.1))),
      column_dend_height=unit(20,"mm"),
      top_annotation = ha,
      bottom_annotation = hb,
      name="log enrichment",
      row_names_side="left",
      column_split=expected.state.clusters,
      column_gap=unit(3,"mm")
    ))
  dev.off()
}

output.state.table <- function (states,state.names, use.exp=T, exp.genes, state.order=c(1:length(state.names)), fdr=0.01) {
    all.tfs <- read.table(op.args[["tfs"]],quote="")
    out.df <- data.frame()

    for (i in state.order) {
        gene.state.names <- subset(states, state %in% i, name)[[1]]
        
        all.tfs.in.state <- gene.state.names[gene.state.names %in% all.tfs[[1]] ]
        tfs.present.in.state <- length(gene.state.names[gene.state.names %in% all.tfs[[1]] ])
        
        num.genes.in.state <- length(gene.state.names)
        
        pc.tot.tfs <- length(all.tfs.in.state)/length(all.tfs[[1]]) * 100
        pc.tfs.in.state <- length(all.tfs.in.state)/length(gene.state.names) * 100
        
        if (use.exp) {
            genes.expressed <- gene.state.names[gene.state.names %in% exp.genes$name[ exp.genes$FDR < fdr]]
            tfs.expressed.in.state <- length(genes.expressed[ genes.expressed %in% all.tfs[[1]] ])
            
            expression.value <- mean(exp.genes[exp.genes$name %in% genes.expressed,2])
            expression.state <- mean(exp.genes[exp.genes$name %in% gene.state.names,2])
            
            pc.exp.st <- length(genes.expressed)/length(gene.state.names) * 100
            pc.exp.tot <- length(genes.expressed)/length(exp.genes[,1]) * 100
        } else {
            genes.expressed = NA
            tfs.expressed.in.state = NA
            
            expression.value = NA
            expression.state = NA
            
            pc.exp.st = NA
            pc.exp.tot = NA
            
            num.tfs = NA
            num.genes = NA
        }
        
        num.tfs <- length(all.tfs[,1])
        num.genes <- length(exp.genes[,1])
        
        fisher.p <- fisher.test(matrix(c(
                        tfs.present.in.state,
                        num.tfs - tfs.present.in.state,
                        num.genes.in.state - tfs.present.in.state,
                        num.genes - num.tfs - num.genes.in.state + tfs.present.in.state
                    ),2,2),alternative="greater")["p.value"]
        
        df.line <- data.frame(
            state = i,
            state.name = state.names[i],
            num.genes = length(gene.state.names),
            num.exp = length(genes.expressed),
            pc.of.state = pc.exp.st,
            pc.of.total = pc.exp.tot,
            mean.exp.exp = expression.value,
            mean.exp.all = expression.state,
            tfs.in.state = tfs.present.in.state,
            pc.tfs.in.state = pc.tfs.in.state,
            pc.all.tfs = pc.tot.tfs,
            fisher.exact.test = fisher.p,
            tfs.expressed = tfs.expressed.in.state
        )

        out.df <- rbind(out.df, df.line)
    }

    return(out.df)
}

generate.state.table = function () {
  state.table.all.states = output.state.table(states=all.modal[["exp"]],state.names=c(1:nStates),exp.genes=exp.genes,state.order=c(1:nStates),fdr=0.01, use.exp=op.args[["gene.expression"]])
  write.table(state.table.all.states,"state.table.all.states.txt", row.names=F, col.names=T, quote=F)
  capture.output( print(state.table.all.states, print.gap=3), file="state.table.all.states.capture.txt")
  
  write.table(state.table.all.states[hc$order[(1:nStates)],],"state.table.all.states.predicted.txt", row.names=F, col.names=T, quote=F)
  capture.output( print(state.table.all.states[hc$order,], print.gap=3), file="state.table.all.states.predicted.capture.txt")
  
  return(state.table.all.states)
}

save.final.data = function () {
  write.table(all.modal[["exp"]],paste(curr.path,"/","single.state.dat",sep=""))

  for (i in c(1:nStates)) {
    out.names = all.modal[["exp"]]$name[ all.modal[["exp"]]$state == i ]
    write.table(out.names,paste(curr.path,"/","state.",i,".genes.txt",sep=""),row.names=F,col.names=F,quote=F)
  }
  
  save.image("all.states.RData")
}

### Process CLI arguments
read.ops(input.args)
if (is.na(op.args[["name"]])) {
  cat("Error: name must be set (use --name=[name])\n")
  q()
}
if (is.na(op.args[["gff.path"]])) {
  cat("Error: path to chromatin GFFs must be set (use --gff.path=[name])\n")
  q()
}

# save defaults if requested
if (op.args[["save.defaults"]]) {
  op.args.saved = op.args[op.args.defaults.to.save]
  cat("Saving defaults ...\n")
  save(list=c("op.args.saved"),file="~/.config/chromatin.universal.rscripts/hmm.defaults")
  cat("Defaults saved to ~/.config/chromatin.universal.rscripts/hmm.defaults\n\n")
  q()
}

### Load libs
library(tools)
library(parallel)
library(RHmm,quietly=T)
library(ggplot2,quietly=T)
library(ggcorrplot,quietly=T)
library(stringr,quietly=T)
library(ComplexHeatmap,quietly=T)

### read input data
name = op.args[["name"]] # must conform to R variable naming conventions -- no spaces, etc
nStates = as.integer(op.args[["nStates"]]) # number of HMM states to fit
rhmm.n.iter = as.integer(op.args[["rhmm.iter"]])
rhmm.n.iter.init = as.integer(op.args[["rhmm.iter.init"]])
base.path = op.args[["gff.path"]]
expected.state.clusters = as.integer(op.args[["expected.state.clusters"]])
mc.cores = as.integer(op.args[["mc.cores"]])
chr.model = strsplit(op.args[["chr.model"]],",")[[1]]

# Variables to save for HMM modelling
save.vars=c("chrom.states","nStates","rhmm.n.iter","rhmm.n.iter.init","name","chr.model","base.path","data.all","data.na","data.na.ex","prot.names","bedgraphs")

### Other globals
prot.names= vector()
bedgraphs = vector()
hmm.fitted.this.session=F

##########################
### Pipeline begins here
#

if (op.args[["fit.hmm"]]) {
  # Make working dir and change to it
  dir.name = paste("analysis",
                   name,
                   format(Sys.time(),"%Y-%m-%d"),
                   "iter",rhmm.n.iter,
                   "init",rhmm.n.iter.init,
                   "st",nStates,
                   sep=".")
  dir.create(dir.name)
  setwd(dir.name)
  
  # Write options file
  write.ops()
  data.all = build.dataframes()
    data.na = data.all
    data.na[data.na == 0] = NA
    #data.na[data.na == "NA"] = NA
    data.na.ex = na.exclude(data.na)
  
  #save.image("test.RData",compress=T)
  #q()
  
  cat("\nFitting HMMs ...\n")
  chrom.states = list()
  chrom.states[[name]] = fit.hmm.mc(data = data.na.ex, cores = mc.cores, chrs = chr.model)
  
  cat("\nSaving data to hmm.RData ...\n")
  oldw = getOption("warn")
  options(warn = -1)
  save(list=save.vars,file="hmm.RData",compress=T)
  options(warn = oldw)
  hmm.fitted.this.session = T

} else if (!file.exists("hmm.RData")) {
  # Not fitting HMM but no current version HMM exists (legacy sanity check)
  # Try to create file from older .RData format ...
  
  if (file.exists(".RData")) {
    cat("  Warning: HMM fit not requested but hmm.RData data not found.\n  Trying to generate hmm.RData from .RData ...\n")
    load(".RData")
    cat("  Saving hmm.RData ...\n")
    save(list=save.vars,file="hmm.RData",compress=T,precheck=T)
    
    if (exists("vit.gff")) {
      cat("\n  Found viterbi path.\n  Saving vit.RData ...\n")
      save(list=c("vit.gff"),file="vit.RData",compress=T)
    } else {
      cat("  Warning: viterbi path not found.\n")
    }
    
    cat("\n  If files saved correctly, please re-run.\n\n")
    q()
  } else {
    cat("  Error: HMM fit not requested but hmm.RData data not found and .RData not found.  Please run with --fit.hmm=T.\n\n")
    q()
  }
  
}

if (is.null(op.args[["fit.viterbi"]]) || as.logical(op.args[["fit.viterbi"]])) {
  if (!hmm.fitted.this.session) {
    cat("Loading previously saved data from hmm.RData ...\n")
    load("hmm.RData")
  }
  
  cat("Fitting Viterbi path ...\n")
  vit.gff = viterbi.path(data = data.na.ex, hmm.model = chrom.states[[name]])
  
  cat("\nSaving data to vit.RData ...\n")
  save(list=c("vit.gff"),file="vit.RData",compress=T)
  cat("HMM fitted.\n")
  
} else {
  cat("Loading model from hmm.RData ...\n")
  load("hmm.RData")
  
  if (!file.exists("vit.RData")) {
    cat("Error: viterbi path not fitted, but no existing data found.  Please run with --fit.viterbi=T.  Exiting ...\n")
    q()
  } else {
    cat("Loading viterbi path from vit.RData ...\n")
    load("vit.RData")
  }
}


#################
### Post-HMMs
#

corr.hc = order.clusters()
prot.order = get.prot.order(prot.names=prot.names,corr.hc=corr.hc)

# Genomic coverage
vit.genomic = genomic.coverage(vit.gff)

# Plot state means for all states
curr.path = "state.mean.plots"
dir.create(curr.path)

# State transitions plotting, clustering and trees
hc = cluster.transMat()
plot.transMat()

# Generate template file for analysis
hb.annot = vector()
hb.cols = list()
generate.state.template()

### genes in chromatin states
cat("\nCalculating chromatin state coverage ...\n")
all.modal = gene.exp.modal(vit.gff, genes.file=op.args[["genes.file"]])

if (op.args[["gene.expression"]]) {
  cat("\nCalculating gene expression ...\n")
  exp.genes = gene.exp(data.na.ex, genes.file=op.args[["genes.file"]])
  cat("Saving gene.exp.RData ...\n")
  save(list=c("exp.genes"),file="gene.exp.RData",compress=T)
} else if (!file.exists("gene.exp.RData")) {
  exp.genes = all.modal$exp
} else {
  cat("Loading saved gene expression data from gene.exp.RData ...\n")
  load("gene.exp.RData")
}

# state tables and predictions
state.table.all.states = generate.state.table()

# Universal heatmap
plot.universal.heatmap()

### save data
save.final.data()

cat("All done.\n")

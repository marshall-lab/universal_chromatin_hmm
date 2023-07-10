# Universal chromatin state analysis script (the single-state bit)
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

# This script should be run after and inherit all variables from the chromatin.universal.hmms script.


library(ape)
library(ggplot2)
library(stringr)

### version
my.version = "0.9.10 (2019-11-25)"
cat(paste("\nUniversal Chromatin HMM analysis script
version ",my.version,"\n\n",sep=""))

### Read CLI options
input.args = commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

read.ops = function (x) {
  for (op in x) {
    op = gsub("^--","",op)
    y = unlist(strsplit(op,"="))

    if (y[1] == "help") {
      cat("Run this script in the analysis directory created by the chromatin.universal.hmms script\n\nOptions:\n")
      for (n in names(single.op.args)) {
        cat(paste("  ",n,"=",single.op.args[[n]],"\n",sep=""))
      }
      cat("\n")
      q()
    }

    if (!is.null(single.op.args[[ y[1] ]])) {
	  print(y[2])
	  if ( grepl("^[[:digit:]]*$", y[2])) {
		single.op.args[[ y[1] ]] <<- as.integer(y[2])
	  } else {
		single.op.args[[ y[1] ]] <<- y[2]
	  }
    } else {
      cat("Error: Option",y[1],"not recognised ...\n")
      q()
    }
  }
}

write.ops = function () {
  out.df = data.frame()
  for (n in names(single.op.args)) {
	v <<- single.op.args[[n]]
	df.line = data.frame(
	  option=n,
	  value=v
	)
	out.df = rbind(out.df, df.line)
  }
  write.table(out.df,"input.args.single.txt",row.names=F)
}

convert.hex.rgb = function (x) {
  hex.r = substr(x,2,3)
  hex.g = substr(x,4,5)
  hex.b = substr(x,6,7)
  dec.r = strtoi(hex.r,16L)
  dec.g = strtoi(hex.g,16L)
  dec.b = strtoi(hex.b,16L)
  out = paste(dec.r,",",dec.g,",",dec.b,sep="")
  return(out)
}

create.wd = function () {
  # Date -- for working directory
  adate = format(Sys.time(), "%Y-%m-%d_%H:%M%P")
  wd = paste(name,adate,sep=".")
  dir.create(wd)
  setwd(wd)
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

gene.pie.simplified = function (genes.df,cols=NULL,names=NULL,simplified.states=NULL,state.order=c(1:length(simplified.states)),num.genes=17718) {
    state.len = vector()
    props = vector()
    j = 0
    out.cols=vector()
    n.states = length(simplified.states)
    
    for (i in state.order) {
        cat("state ",i,":",simplified.states[[i]],"\n")
        j = j+1
        cat(i,"\n")
        out.cols[j] = cols[i]
        stlen = length(genes.df$name[genes.df$state %in% simplified.states[[i]] ])
        state.len[j]=stlen
        fnum=format(round((stlen/num.genes*100),2),nsmall=2)
        props[j]=paste(names[i],"\n",fnum,"%")
    }
    
    rem = num.genes-length(genes.df$name)
    state.len[length(state.len)+1] = rem
    props[length(props)+1] = paste("N/A\n",format(round((rem/num.genes*100),2),nsmall=2),"%")
    out.cols[length(out.cols)+1] = "white"

    pie(state.len,col=out.cols,labels=props,border="white")
    return(state.len)
}

condense.gff = function (input.df) {
	total = nrow(input.df)
	merged = list()
	
	chr.last = ""
	state.start = NA
	state.end = NA
	state.current = NA
	rcount = 1
		
	for (i in c(1:total)) {
	  if (i%%1000 == 0) {cat(paste("Processing row",i,"of",total,"             \r"))}
	  
	  chr = .subset2(input.df,1)[i]
	  start = .subset2(input.df,2)[i]
	  end = .subset2(input.df,3)[i]
	  state = .subset2(input.df,4)[i]
	  
	  if (chr == chr.last) {
		if (state == state.current && start == state.end) {
		  state.end = end
		} else {
		  # save state block
		  merged[[rcount]] = data.frame(chr=chr, start=state.start, end=state.end, state=state.current)
		  rcount = rcount+1
		  
		  state.start = start
		  state.end = end
		  state.current = state
		}
	  } else {
		# new chromosome
		state.start = start
		state.end = end
		state.current = state
	  }
	  
	  chr.last = chr
	}
	merged[[rcount]] = data.frame(chr=chr, start=state.start, end=state.end, state=state.current)
	
	merged.df = do.call('rbind',merged)
	return(merged.df)
}

make.states.bed = function (datf, simp) {
  # expect "states gff" input: chr start end state
  datf$zero = 0
  datf$plus = "+"
  datf$cols = sapply(datf$state, function(x) convert.hex.rgb(simp.cols[[name]][x]) )
  
  states.bed = data.frame(datf$chr,datf$start,datf$end,datf$state,datf$zero,datf$plus,datf$start,datf$end,datf$cols)
  return(states.bed)
}

reverse.states = function () {
  #reverse states lookup
  states.reverse = vector()
  for (i in 1:length(simp.states[[name]])) {
    for (j in unlist(simp.states[[name]][i])) {
      states.reverse[j] = i
    }
  }
  return(states.reverse)
}

create.simp.vit = function () {
  simp.vit.gff = vit.gff
  simp.vit.gff$state = states.reverse[simp.vit.gff$state]
  write.table(data.frame(simp.vit.gff$chr,".",".",simp.vit.gff$start,simp.vit.gff$end,simp.vit.gff$state),paste("simplified.states.gff",sep=""), col.names=F, row.names=F, quote=F)
  return(simp.vit.gff)
}

write.states.data = function () {
  write.table(simp.vit.genes,paste("single.state.simp.dat"))
  write.table(simp.vit.genes$name,"all.analysed.genes.txt",row.names=F,col.names=F,quote=F)
  write.table(chrom.states.genes,paste("single.state.dat"))
}

write.state.table = function (state.table,filename) {
  write.table(state.table,file=paste(filename,"txt",sep="."), row.names=F, col.names=T, quote=F)
  capture.output( print(state.table, print.gap=3), file=paste(filename,"capture","txt",sep="."))
}

plot.pie.charts = function () {
  ###########################################
  ### Pie charts
  # Simplified states
  pdf(paste("simplified.state.pie.pdf",sep=""))
  gene.pie.simplified(
                      data.frame(
                                name=chrom.states.genes$name,
                                state=chrom.states.genes$state
                                ),
                      cols=simp.cols[[name]],
                      names=simp.names[[name]],
                      simplified.states=simp.states[[name]],
                      state.order=c(1:length(simp.states[[name]])),
                      num.genes=gene.count
                    )
  #dev.copy(png,paste("simplified.state.pie.png",sep=""),width=600,height=600);dev.off()
  dev.off()
  
  # All states pie
  full.cols = vector()
  for (i in 1:length(simp.states[[name]])) {
    for (j in unlist(simp.states[[name]][i])) {
      full.cols[j] = simp.cols[[name]][i]
    }
  }
  
  pdf(paste("all.states.pie.pdf",sep=""))
  gene.pie.simplified(
                      chrom.states.genes,
                      names=c(1:target.state),
                      simplified.states=c(1:target.state),
                      state.order=unlist(simp.states[[name]]),
                      cols=full.cols,
                      num.genes=gene.count
                      )
  dev.off()
  
  pdf(paste("simp.vit.pie.pdf",sep=""))
  gene.pie.simplified(
                      data.frame(
                                name=simp.vit.genes$name,
                                state=simp.vit.genes$state
                                ),
                      cols=simp.cols[[name]],
                      names=simp.names[[name]],
                      simplified.states=c(1:length(simp.names[[name]])),
                      state.order=c(1:length(simp.states[[name]])),
                      num.genes=gene.count
                    )
  #dev.copy(png,paste("simp.vit.pie.png",sep=""),width=600,height=600);dev.off()
  dev.off()
}

write.gene.lists = function () {
  ### Gene lists
  # simplified states
  for ( i in c(1:length(simp.states[[name]])) ) {
    out.names = chrom.states.genes$name[ chrom.states.genes$state %in% simp.states[[name]][[i]] ]
    write.table(out.names,paste(simp.names[[name]][i],".genes.txt",sep=""),row.names=F,col.names=F,quote=F)
  }
  
  # all states
  for (i in c(1:target.state)) {
    out.names = chrom.states.genes$name[ chrom.states.genes$state == i ]
    write.table(out.names,paste("state",i,"genes.txt",sep="."),row.names=F,col.names=F,quote=F)
  }
}

plot.barplots = function () {
  ### barplots
  ### reverse lookup lists for colours
  all.states = list()
  all.cols = list()
  all.name = list()
  real.cols = list()
  real.names = list()
  for (i in 1:length(simp.states[[name]])) {
    all.states[[name]] = append(all.states[[name]], simp.states[[name]][[i]])
    for (j in 1: length(simp.states[[name]][[i]])) {
      all.cols[[name]] = append(all.cols[[name]], simp.cols[[name]][i])
      all.name[[name]] = append(all.name[[name]], simp.names[[name]][i])
    }
  }
  
  # proportion of genes expressed
  pdf(paste("barpolot.prop.of",name,"genes.expr.by.state.pdf",sep="."),height=5,width=0.5*target.state)
  barplot(state.table.all.states$pc.of.state,col=all.cols[[name]],main=paste("Proportion of",name,"genes expressed by state"),names.arg=state.table.all.states$state)
  dev.off()
  
  pdf(paste("barpolot.prop.of",name,"simp.vit.genes.expr.by.state.pdf",sep="."),height=5,width=5)
  par(las=1)
  par(mar=c(5,8,4,2))
  barplot(rev(state.table.simp.vit$pc.of.state),col=rev(simp.cols[[name]]),main=paste("Proportion of",name,"genes expressed by state"),names.arg=rev(state.table.simp.vit$state.name), horiz=T)
  dev.off()
  
  # Enrichment of TFs
  pdf(paste("barpolot.enrichment.of.tfs.log-p-value",name,"simp.vit.pdf",sep="."),height=5,width=0.5*target.state)
  par(las=1)
  par(mar=c(5,8,4,2))
  barplot(rev(-log(state.table.simp.vit$p.value)),col=rev(simp.cols[[name]]),main=paste("Enrichment of TFs in state"),names.arg=rev(state.table.simp.vit$state.name), horiz=T,xlab="-log(P value) Fischer's exact test"); abline(v=-log(0.05/length(simp.cols[[name]])));
  dev.off()
  
  # Percentage of TFs
  pdf(paste("barpolot.percentage.of.tfs",name,"simp.vit.pdf",sep="."),height=5,width=0.5*target.state)
  par(las=1)
  par(mar=c(5,8,4,2))
  barplot(rev(state.table.simp.vit$pc.tfs.in.state),col=rev(simp.cols[[name]]),main=paste("Proportion of TFs in state"),names.arg=rev(state.table.simp.vit$state.name), horiz=T,xlab="Percentage of TFs")
  dev.off()
}

plot.boxplots = function () {
  ### protein binding boxplots -- unordered
  for (i in c(4:(3+length(prot.names)))) {
    prot.name = names(data.na.ex)[i]
    binding = list()
    for (j in c(1:target.state)) {
    binding[[j]] = data.na.ex[ vit.gff$states == j, i]
    }
    pdf(paste("prot.binding.boxplot",prot.name,"pdf",sep="."))
    title = prot.name
    boxplot(binding,ylim=c(-5,5), main=title,names=c(1:target.state))
    abline(0,0)
    dev.off()
  }
  
  ### protein binding boxplots -- ordered
  for (i in c(4:(3+length(prot.names)))) {
    prot.name = names(data.na.ex)[i]
    binding = list()
    k = 1;
    for (j in all.states[[name]]) {
    binding[[k]] = data.na.ex[ vit.gff$states == j, i]
    k = k+1
    }
    pdf(paste("prot.binding.boxplot",prot.name,"ordered.pdf",sep="."))
    title = prot.name
    boxplot(binding,ylim=c(-5,5), main=title,names=all.states[[name]])
    abline(0,0)
    dev.off()
  }
}

draw.prop.hists = function () {
  # histograms of proportional coverage
  for ( i in c( 1:length(simp.names[[name]]) ) ) {
    simp.vit.state = subset(simp.vit.genes, state == i)
    
    if (length(simp.vit.state$prop) < 2) { next }
    pdf(paste("simp.vit.state.coverage.",simp.names[[name]][i],".pdf",sep=""))
    hist(simp.vit.state$prop, col=simp.cols[[name]][i], main=simp.names[[name]][i], xlim=c(0,1), xlab="Proportion of gene body covered by modal state")
    dev.off()
  }
  
  pdf(paste("simp.vit.state.coverage.all.pdf",sep=""))
  hist(simp.vit.genes$prop, col="#666666", main="All states", xlim=c(0,1), xlab="Proportion of gene body covered by modal state")
  dev.off()
}

write.simp.vit.genelists = function () {
  for (i in 1:length(simp.states[[name]])) {
    out.names = simp.vit.genes$name[ simp.vit.genes$state == i ]
    write.table(out.names,paste("simp.vit.genes.state",simp.names[[name]][i],"genes.txt",sep="."),row.names=F,col.names=F,quote=F)
  }
}

plot.gene.state.boxplots = function () {
  ### Gene state boxplots
  for (i in c(1:target.state)) {
    data.cut = data.na.ex[ vit.gff$state==i, c(4:ncol(data.na.ex))]
    pdf(paste("detailed.boxplot.state",i,"pdf",sep="."))
    title = paste("State ",i)
    boxplot(data.cut, ylim=c(min(data.cut),max(data.cut)),main=title,names=prot.names)
    abline(0,0)
    dev.off()
  }
}
  
write.coloured.bedtracks = function () {
  ### reverse lookup lists for colours
  all.states = list()
  all.cols = list()
  all.name = list()
  real.cols = list()
  real.names = list()
  for (i in 1:length(simp.states[[name]])) {
    all.states[[name]] = append(all.states[[name]], simp.states[[name]][[i]])
    for (j in 1: length(simp.states[[name]][[i]])) {
      all.cols[[name]] = append(all.cols[[name]], simp.cols[[name]][i])
      all.name[[name]] = append(all.name[[name]], simp.names[[name]][i])
    }
  }
  
  for (i in 1:length(all.states[[name]])) {
    real.cols[[name]] = append(real.cols[[name]], all.cols[[name]][ all.states[[name]] == i ])
    real.names[[name]] = append(real.names[[name]], all.name[[name]][ all.states[[name]] == i ] )
  }
  
  write.table(real.cols[[name]],"state.mean.hmap.sc.cols.dat")
  
  ### coloured IGV tracks
  write('track name="" description="" visibility=2 itemRgb="On"', paste(name,"igv.cols.bed",sep="."))
  write.table(make.states.bed(simp.vit.gff, simp.cols), paste(name,"igv.cols.bed",sep="."), append=T, quote=F, row.names=F,col.names=F)
  
  write('track name="" description="" visibility=2 itemRgb="On"', paste(name,"igv.cols.cond.bed",sep="."))
  write.table(make.states.bed(condense.gff(simp.vit.gff), simp.cols), paste(name,"igv.cols.cond.bed",sep="."), append=T, quote=F, row.names=F,col.names=F)
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
  
  if (	op.args[["gene.expression"]]) {
    annot_df$exp.pc = state.table.all.states$pc.of.state
    col[["exp.pc"]] = circlize::colorRamp2(c(0, max(na.omit(annot_df$exp.pc))), 
                                         c("white", "orange"))
  }
  
  ha = HeatmapAnnotation(df = annot_df, col = col, show_annotation_name = F)
  
  hc = hclust(
        d = dist(chrom.states[[name]]$HMM$transMat),
        method = "ward.D2"
        )
  
  pdf(paste("Universal_heatmap",name,nStates,"pdf",sep="."),height=0.5*length(prot.names)+1,width=0.5*nStates+1)
  print(Heatmap(
      as.matrix(named.smh),
      cluster_columns=hc,
      cluster_rows=F,
      row_order=prot.order,
      col=circlize::colorRamp2(c(-5,0,5),c(rgb(0.2,0.3,0.6),rgb(1,1,1),rgb(0.9,0.1,0.1))),
      column_dend_height=unit(20,"mm"),
      top_annotation = ha,
      bottom_annotation = hb,
      name="log enrichment",
      row_names_side="left",
      column_split=expected.state.clusters,
      column_gap=unit(3,"mm")
      #rect_gp = gpar(col = "white", lwd = 1),
      #column_split=cutree(hc,9),
    ))
  dev.off()
}

single.op.args = list(
  "input" = Sys.glob("single*.r"),
  "data" = "all.states.RData",
  "genes.file" = "/mnt/genomes/dmel_release/DmR6/DmR6.genes.gff",
  "tfs" = "/home/owen/dros.transcription.factors.tfs-update.txt",
  "pad" = 0,
  "process.exons" = F,
  "analyse.goterms" = T
)

read.ops(input.args)

if (is.na(single.op.args[["input"]])) {
  cat("Error: input file must be set (use --input=[single.state.r file])\n\n")
  q()
}

write.ops()

simp.names = list()
simp.cols = list()
simp.states = list()
load(single.op.args[["data"]])

data.file = tryCatch({
    load(single.op.args[["data"]])
  },  warning = function (e) {
    cat("Error: unable to read .RData file (use --data=[.RData file] to set)\n\n")
    q()
  }
)


###########################################
#### analysis pipeline starts here


# input file is read from the commandline option:
source(single.op.args[["input"]])

create.wd()

# read TFs data
all.tfs = read.table(single.op.args[["tfs"]],quote="")
target.state = nStates

states.reverse = reverse.states()
simp.vit.gff = create.simp.vit()
simp.vit.genes = gene.exp.modal(simp.vit.gff,
                                genes.file=single.op.args[["genes.file"]],
                                pad=single.op.args[["pad"]]
                                )$exp

chrom.states.genes = all.modal[["exp"]]
gene.count = all.modal[["count"]]

### save data
write.states.data()

# exon scatter plots
if (single.op.args[["process.exons"]]) {
  process.exons()
}

#cat("Saving analysis.debug.RData ...\n")
#save.image("analysis.debug.RData")

plot.pie.charts()
write.gene.lists()

write.table(exp.genes, paste(name,"gene.expression.dat",sep="."))
write.table(exp.genes$name[exp.genes$FDR < 0.01], paste(name,"expressed.genes.txt",sep="."),row.names=F,col.names=F,quote=F)

all.states.simp.genes = chrom.states.genes
all.states.simp.genes$state = states.reverse[all.states.simp.genes$state]

# expressed genes in states
state.table.all.states = output.state.table(states=chrom.states.genes,state.names=simp.names[[name]][states.reverse],exp.genes=exp.genes,state.order=unlist(simp.states[[name]]),fdr=0.01)
write.state.table(state.table=state.table.all.states,filename="state.table.all.states")

# expressed genes in simp.vit states
state.table.simp.vit = output.state.table(states=simp.vit.genes,state.names=simp.names[[name]],exp.genes=exp.genes,fdr=0.01)
write.state.table(state.table=state.table.simp.vit,filename="state.table.simp.vit.states")

draw.prop.hists()
write.simp.vit.genelists()
plot.gene.state.boxplots()
plot.barplots()
write.coloured.bedtracks()

save.image("single.Rdata")

if (single.op.args[["analyse.goterms"]]) {
  cat("Running goterm analysis ...")
  system("perl ~/Dropbox/perl/flymine.goterms.pl *.genes.txt")
}

setwd("../")
cat("Copying single.state.r file to working directory ...")
file.copy(single.op.args[["input"]],wd)

cat("All done.\n")

library(MCRestimate)
library(randomForest)
library(pamr)
library(e1071)

MyLoad <- function(x){
 A <- get(load(x))
 return(A)
}


the.expression.set <- MyLoad(DataSet)


if(RF){
if (is.null(mtry.range)){
  if(is.null(ntree.range)){
    list.of.parameter <- c(parameter.for.preprocessing)
  }else{
    list.of.parameter <- c(list(ntree=ntree.range),parameter.for.preprocessing)
  }
}else{
 if(is.null(ntree.range)){
    list.of.parameter <- c(list(mtry=mtry.range),parameter.for.preprocessing)
  }else{
    list.of.parameter <- c(list(mtry=mtry.range,ntree=ntree.range),parameter.for.preprocessing)
  }
}
 r.forest <- MCRestimate(the.expression.set,
                         class.colum,
                         classification.fun="RF.wrap",
                         poss.parameters=list.of.parameter,
			 thePreprocessingMethods=thePreprocessingfunktionsRF,
                         cross.outer=cross.outer,
                         cross.inner=cross.inner,
                         cross.repeat=cross.repeat,
                         reference.class=ref.class,
                         plot.label=plot.label,
                         block.column=block.column,
			 stratify=Strat,
			 rand=SEED)
 save(r.forest, file=paste("backRF",SEED,".RData",sep=""))
    }


if(GPLS)
   {r.gpls <- MCRestimate(the.expression.set,
                          class.colum,
                          classification.fun="GPLS.wrap",
                          poss.parameter=parameter.for.preprocessing,
			  thePreprocessingMethods=thePreprocessingfunktionsGPLS,
                          cross.outer=cross.outer,
                          cross.repeat=cross.repeat,
                          cross.inner=cross.inner,
                          reference.class=ref.class,
                          plot.label=plot.label,
                          block.column=block.column,
			  stratify=Strat,
			 rand=SEED)
    save(r.gpls, file=paste("backGPLS",seed,".RData",sep=""))
  }



if(PAM){
list.of.parameter <- c(list(threshold=thresholds),parameter.for.preprocessing)

r.pam <- MCRestimate(the.expression.set,
                     class.colum,
                     classification.fun="PAM.wrap",
                     poss.parameter=list.of.parameter,
		     thePreprocessingMethods=thePreprocessingfunktionsPAM,
                     cross.outer=cross.outer,
                     cross.repeat=cross.repeat,
                     cross.inner=cross.inner,
                     reference.class=ref.class,
                     plot.label=plot.label,
                     block.column=block.column,
		     stratify=Strat,
                     rand=SEED)
save(r.pam,file=paste("backPAM",SEED,".RData",sep=""))
}


if(PLR){
list.of.parameter <- c(list(kappa=kappa.range),parameter.for.preprocessing)
r.logReg <- MCRestimate(the.expression.set,
                          class.colum,
                          classification.fun="PLR.wrap",
                          poss.parameter=list.of.parameter,
			  thePreprocessingMethods=thePreprocessingfunktionsPLR,
                          cross.outer=cross.outer,
                          cross.repeat=cross.repeat,
                          cross.inner=cross.inner,
                          reference.class=ref.class,
                          plot.label=plot.label,
                        block.column=block.column,
			stratify=Strat,
			 rand=SEED)
save(r.logReg,file=paste("backlogReg",SEED,".RData",sep=""))
}

if(SVM){
 list.of.parameter <- c(list (gamma=gamma.range,cost=cost.range),parameter.for.preprocessing)
 r.svm <- MCRestimate(the.expression.set,
                      class.colum,
                      classification.fun="SVM.wrap",
                      poss.parameter=list.of.parameter,
		      thePreprocessingMethods=thePreprocessingfunktionsSVM,
                      cross.outer=cross.outer,
                      cross.repeat=cross.repeat,
                      cross.inner=cross.inner,
                      reference.class=ref.class,
                      plot.label=plot.label,
                      block.column=block.column,
		      stratify=Strat,
			 rand=SEED)
 save(r.svm, file=paste("backSVM",SEED,".RData",sep=""))
}

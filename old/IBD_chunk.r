#### prep data (don't forget "as.matrix")
setwd('/Volumes/SeaDrive/Users/Gabe/Desktop/Books/knights/input.v4')

 table=read.delim('map-subset-imputed.txt',sep='\t',comment='',head=T,row.names=1,check.names=F);
#attach(table)

 Mergetaxa=as.matrix(read.delim('merged-taxa-subset.txt',sep='\t',comment='',head=T,row.names=1,check.names=F));
 taxa=as.matrix(read.delim('taxa-subset/otutable-min500-subset_L7.txt',sep='\t',comment='',head=T,row.names=1,check.names=F));


tax2tab = match(rownames(table),colnames(taxa))
taxa = taxa[,tax2tab] # put table in taxa's order

# weed out OTUs that have only 0's
taxThres = 0.000001;
XT = taxa[rowMeans(taxa) > taxThres,] # 1,000,000 times faster!

########## t test fusion version: setup
Groups = table$Disease;
num = dim(XT)[1]
# whole = sample(seq(1,dim(table)[1]),dim(table)[1]); 
# part = 3/4; revP = 1-part;
# half = whole[seq(1,length(whole)*part)]; 
#XT = XT[,half]; Groups = Groups[half];
CDix = Groups == "CD"; UCix = Groups == "UC";

# extremely fast hybrid apply+for versions
library(polycor)
baseline.tt <- apply(XT,1,function(xx) t.test(xx[CDix], xx[UCix])$p.value)
#baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Disease)) # at end

pt = proc.time();
mat.tt = matrix(rep(0,num*num),num); #mat.ps = mat.tt;
for (i in 1:num) for (j in i:num) {
 mat.tt[i,j] = t.test((XT[i,CDix] + XT[j,CDix])*.5, (XT[i,UCix] + XT[j,UCix])*.5)$p.value;
 mat.tt[j,i] = mat.tt[i,j];
 #mat.ps[i,j] = polyserial(XT[i,] + XT[j,], Disease);
 #mat.ps[j,i] = mat.ps[i,j];
}
proc.time() - pt

# save(whole, half, XT, baseline.tt, mat.tt, file="halfsies.RData")

##### gradient descent version (average)
doBlocking = TRUE; thresh = 0.9; maxThresh = 0.5; path = NULL;
ptm <- proc.time()
iter = 1; found = 0; mdim=num; numTests = num*(num-1)/2 + num;
blocked = rep(FALSE,num);
while (TRUE) { #(min(mat.tt) < Inf) { # (dim(mat.tt)[1] > 10) {
#print(paste("Iteration:",iter));

# here's the one with minimum p value:
pair = which.min(mat.tt);
px = ceiling(pair/mdim); py = pair - mdim * (px - 1);
low = min(px,py); high = max(px, py);
if (mat.tt[pair] == Inf) break;
tt.both = mat.tt[low,high];
#ps.both = mat.ps[low,high];
if (low == high) {
  mat.tt[low,] = Inf;
  mat.tt[,low] = Inf;
  if (tt.both > maxThresh) { print("Finalizing all blocks..."); break; }
  if (doBlocking) { blocked[low] = TRUE;
  print(paste("finalizing block:",low)); }
  next;
}
if (tt.both > baseline.tt[low]*thresh || tt.both > baseline.tt[high]*thresh) {
  mat.tt[low,high] = Inf;
  mat.tt[high,low] = Inf;
  #mdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
  #print(paste('skipping',low,high,"as max difference:",mdiff));
  next; 
}

found = found + 1; mdim = mdim - 1;
print(paste("found:",found,"at:",low,high,"NewP, low, high:",tt.both,baseline.tt[low], baseline.tt[high]));
path = c(path,low,high);
XT[low,] = (XT[low,] + XT[high,])*.5;
rownames(XT)[low] = paste(rownames(XT)[low],rownames(XT)[high],sep='\n')

# remove the trailing extra row (and column, for the matrices based on XT)
XT = XT[-high,];
#mat.ps = mat.ps[-high,-high];
mat.tt = mat.tt[-high,-high];
#baseline.ps = baseline.ps[-high];
baseline.tt = baseline.tt[-high];
blocked = blocked[-high];


# recalculate the new row and col's tt and ps
baseline.tt[low] = tt.both;
for (j in 1:mdim) {
  #if (sum(mat.tt[j,]==Inf) < mdim - 2) 
  if (!blocked[j]) 
  { mat.tt[low,j] = t.test((XT[low,CDix] + XT[j,CDix])*.5,(XT[low,UCix] + XT[j,UCix])*.5)$p.value; }
  else mat.tt[low,j] = Inf;
}
mat.tt[,low] = mat.tt[low,];

numTests = numTests + mdim - sum(blocked);
iter = iter + 1; 
}

baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Groups))
names(baseline.tt) = rownames(XT);
proc.time() - ptm

###### recreate a combined microbe dataset from given Path
#XT.old = XT;
half2 = whole[seq(length(whole)*part + 1,length(whole) + 1)]; 
Groups2 = table$Disease[half2];
CDix2 = Groups2 == "CD"; UCix2 = Groups2 == "UC";
YT = taxa[rowMeans(taxa) > taxThres,half2]

endp = length(path) / 2; 
for (i in 1:endp) {
 high = path[i*2]; low = path[i*2 - 1];
 YT[low,] = (YT[low,] + YT[high,])*.5;
 rownames(YT)[low] = paste(rownames(YT)[low],rownames(YT)[high],sep='\n');
 YT = YT[-high,];
}

######## random forests
library(randomForest)
RF = randomForest(t(XT),Groups,importance=TRUE,ntree=2000); RF
#RFb = randomForest(t(XT[order(baseline.tt)[1:sum(baseline.tt*numTests < 1)],]),Groups,importance=TRUE,ntree=1000); RFb
RFo = randomForest(t(taxa[,half]),Groups,importance=TRUE,ntree=2000); RFo



RF2 = randomForest(t(YT),Groups2,importance=TRUE,ntree=2000); RF2
RFo2 = randomForest(t(taxa[,half2]),Groups2,importance=TRUE,ntree=2000); RFo2

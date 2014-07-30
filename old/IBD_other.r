# old apply setup
#mat.tt = apply(XT,1,function(yy) apply(XT,1,function(xx) t.test((xx + yy)[CDix], (xx + yy)[UCix])$p.value))
#mat.ps = apply(XT,1,function(yy) apply(XT,1,function(xx) polyserial(xx + yy, Disease)))


######## polyserial version (too slow)
# get baseline values
num = dim(taxa)[1]
baseline = rep(0,num)
# build L1 comparison matrix
for (i in 1:num) baseline[i] = polyserial(as.numeric(taxa[i,]),Disease)
mat = matrix(rep(0,num*num),num)
#res = rep(0,num*(num-1)/2) # alternative version
ix = 1;
for (i in 1:num) {
  for (j in i+1:num) {
    if (j > num) break;
    mat[i,j] = polyserial(as.numeric(taxa[i,]+taxa[j,]),Disease)
  }
}

########## t test fusion version

num = dim(XT)[1]
baseline = rep(0,num)

# extremely fast hybrid apply+for versions
library(polycor)
baseline.tt <- apply(XT,1,function(xx) t.test(xx[CDix], xx[UCix])$p.value)
baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Disease))

baseline.tt = rep(0,num); baseline.ps = baseline.tt;
for (i in 1:num) {
 baseline.tt[i] = t.test(XT[i,CDix],XT[i,UCix])$p.value;
 baseline.ps[i] = polyserial(XT[i,],Disease);
}

#mat = matrix(rep(0,num*num),num)
#pt = proc.time();
#mat.tt = 0; mat.ps = 0;
#mat.tt = apply(XT,1,function(yy) apply(XT,1,function(xx) t.test((xx + yy)[CDix], (xx + yy)[UCix])$p.value))
#mat.ps = apply(XT,1,function(yy) apply(XT,1,function(xx) polyserial(xx + yy, Disease)))
#proc.time() - pt

pt = proc.time();
mat.tt = matrix(rep(0,num*num),num); mat.ps = mat.tt;
for (i in 1:num) for (j in i:num) {
 mat.tt[i,j] = t.test(XT[i,CDix] + XT[j,CDix], XT[i,UCix] + XT[j,UCix])$p.value;
 mat.tt[j,i] = mat.tt[i,j];
 mat.ps[i,j] = polyserial(XT[i,] + XT[j,], Disease);
 mat.ps[j,i] = mat.ps[i,j];
}
proc.time() - pt


# print 10 best correlations:
as.numeric(baseline.ps[order(abs(baseline.ps),decreasing=T)[1:10]])

# print 10 best p values:


#for (i in 1:num) baseline[i] = t.test(taxa[i,Disease=="CD"],taxa[i,Disease=="UC"])
CDix = Disease == "CD"; UCix = Disease == "UC";



######### test version -- pairwise diffs only
ptm <- proc.time()

iter = 1; found = 0; numTests = num*(num-1)/2 + num; mdim = num;
blocked = rep(FALSE,num);


## MAKE ANOTHER MATRIX
mat.tt.distpre = apply(mat.tt,1,function(xx) baseline.tt / xx)
mat.tt.dist = mat.tt.distpre;
for (j in 1:mdim) for (k in j:mdim) { 
 mat.tt.dist[j,k] = min(mat.tt.dist[j,k],mat.tt.dist[k,j]); 
 mat.tt.dist[k,j]=mat.tt.dist[j,k]; 
}

while (TRUE) { #(min(mat.tt) < Inf) { # (dim(mat.tt)[1] > 10) {
pair = which.max(mat.tt.dist);
px = ceiling(pair/mdim); py = pair - mdim * (px - 1);
low = min(px,py); high = max(px, py);
#if (mat.tt[pair] == Inf) break;
tt.both = mat.tt[low,high];
#ps.both = mat.ps[low,high];

print(paste("Iteration:",iter, "mat.tt.dist:",mat.tt.dist[px,py]));


if (low == high) {
  mat.tt[low,] = Inf; 
  mat.tt[,low] = Inf; 
  blocked[low] = TRUE;
  print(paste("match exact at:",low));
  break;
}
if (tt.both >= baseline.tt[low]*0.5 || tt.both >= baseline.tt[high]*0.5) {
  mat.tt[low,high] = Inf;
  mat.tt[high,low] = Inf;
  mdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
  print(paste('skipping',low,high,"as max difference:",mdiff));
  break; 
}
found = found + 1; mdim = mdim - 1;

 print(paste("found:",found,"at:",low,high,"Diff:",tt.both / min(baseline.tt[low],baseline.tt[high]), "New, l, h:",tt.both,baseline.tt[low], baseline.tt[high]));

XT[low,] = XT[low,] + XT[high,];
rownames(XT)[low] = paste(rownames(XT)[low],rownames(XT)[high],sep='\n')

# remove the trailing extra row (and column, for the matrices based on XT)
XT = XT[-high,];
#mat.ps = mat.ps[-high,-high];
mat.tt = mat.tt[-high,-high];
#mat.tt.distpre = mat.tt.distpre[-high,-high];
mat.tt.dist = mat.tt.dist[-high,-high];
#baseline.ps = baseline.ps[-high];
baseline.tt = baseline.tt[-high];


# recalculate the new row and col's tt # and ps
baseline.tt[low] = tt.both;
#baseline.ps[low] = ps.both;
#low_s = as.numeric(XT[low,])
##mat.tt[low,] = apply(XT,1,function(xx) t.test((low_s + xx)[CDix], (low_s + xx)[UCix])$p.value)
# new version
for (j in 1:mdim) {
  #if (sum(mat.tt[j,]==Inf) < mdim - 2) 
  if (!blocked[j]) 
  { mat.tt[low,j] = t.test(XT[low,CDix] + XT[j,CDix],XT[low,UCix] + XT[j,UCix])$p.value; }
  else mat.tt[low,j] = Inf;
}
mat.tt[,low] = mat.tt[low,] # row and col the same

#mat.ps[low,] = apply(XT,1,function(xx) polyserial(low_s + xx, Disease))
#mat.ps[,low] = mat.ps[low,]

for (j in 1:mdim) {
 if (!blocked[j]) 
 { mat.tt.dist[low,j] = min(baseline.tt[j]/mat.tt[low,j],baseline.tt[low]/mat.tt[j,low]); }
 else { print(paste("blocked",low,j)); mat.tt.dist[low,j] = 0; }
}
mat.tt.dist[,low] = mat.tt.dist[low,];

numTests = numTests + mdim - 1; # really -1?
iter = iter + 1; #mdim = num - found;

}
baseline.ps <- apply(XT,1,function(xx) polyserial(xx, table$Disease))
names(baseline.tt) = rownames(XT);
proc.time() - ptm

#### correlation version
ptm <- proc.time()
iter = 1; found = 0; numTests = num*(num-1)/2 + num;
while (TRUE) { #(min(mat.tt) < Inf) { # (dim(mat.tt)[1] > 10) {
#print(paste("Iteration:",iter));
iter = iter + 1;
# here's the one with minimum correlation value:
# Faster edition
pair = which.max(abs(mat.ps));
px = ceiling(pair/(num-found)); py = pair - (num-found) * (px - 1);
low = min(px,py); high = max(px, py);
if (mat.ps[pair] == 0) break;
tt.both = mat.tt[low,high];
ps.both = mat.ps[low,high];
if (low == high) {
  mat.ps[low,] = 0; 
  mat.ps[,low] = 0; 
  print(paste("match exact at:",low,high));
  next;
}
if (ps.both > 0) if (ps.both < baseline.ps[low]*1.1 || ps.both < baseline.ps[high]*1.1) {
  mat.ps[low,high] = 0;
  mat.ps[high,low] = 0;
  mdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
  #print(paste('skipping',low,high,"as max difference too low:",mdiff));
  next; 
}
if (ps.both < 0) if (ps.both > baseline.ps[low]*0.9 || ps.both > baseline.ps[high]*0.9) {
  mat.ps[low,high] = 0;
  mat.ps[high,low] = 0;
  mdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
  #print(paste('skipping',low,high,"as max difference too high:",mdiff));
  next;
}

found = found + 1;
print(paste("found:",found,"at:",low,high,"NewP, low, high:",tt.both,baseline.tt[low], baseline.tt[high]));

XT[low,] = XT[low,] + XT[high,];
rownames(XT)[low] = paste(rownames(XT)[low],rownames(XT)[high],sep='\n')

# remove the trailing extra row (and column, for the matrices based on XT)
XT = XT[-high,];
mat.ps = mat.ps[-high,-high];
mat.tt = mat.tt[-high,-high];
baseline.ps = baseline.ps[-high];
baseline.tt = baseline.tt[-high];

# transfer the new names over (not strictly necessary)
#colnames(XT) = rownames(XT); #columns are sample names
rownames(mat.ps) = rownames(XT);
colnames(mat.ps) = rownames(XT);
rownames(mat.tt) = rownames(XT);
colnames(mat.tt) = rownames(XT);
names(baseline.ps) = rownames(XT);
names(baseline.tt) = rownames(XT);

# recalculate the new row and col's tt and ps
baseline.tt[low] = tt.both;
baseline.ps[low] = ps.both;
low_s = as.numeric(XT[low,])
mat.tt[low,] = apply(XT,1,function(xx) t.test((low_s + xx)[CDix], (low_s + xx)[UCix])$p.value)
mat.tt[,low] = mat.tt[low,] # row and col the same

mat.ps[low,] = apply(XT,1,function(xx) polyserial(low_s + xx, Disease))
mat.ps[,low] = mat.ps[low,]

numTests = numTests + length(mat.tt[1,]) - 1;

}

proc.time() - ptm

################## nightmare edition (not mean, not add -- max fill!)
### prep the data (redundant?)
Groups = table$Disease;
num = dim(XT)[1]
# whole = sample(seq(1,dim(table)[1]),dim(table)[1]); 
# half = whole[seq(1,length(whole)/2)]; 
#XT = XT[,half]; Groups = Groups[half];
CDix = Groups == "CD"; UCix = Groups == "UC";

# extremely fast hybrid apply+for versions
library(polycor)
baseline.tt <- apply(XT,1,function(xx) t.test(xx[CDix], xx[UCix])$p.value)
#baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Disease)) # at end

# combine two rows in an additive way but only up to the max in either
maxFill = function(grp1a, grp1b,grp2a,grp2b) {
  mx = max(max(grp1a),max(grp1b),max(grp2a),max(grp2b));
  temp1 = grp1a + grp1b;
  temp2 = grp2a + grp2b;
  #temp*(mx/max(temp));
  t.test(sapply(temp1,function(x) min(x,mx)),sapply(temp2,function(x) min(x,mx)))$p.value;
}

maxFillT = function(grp1, grp2) {
  mx = max(max(grp1),max(grp2));
  temp = grp1 + grp2;
  #temp*(mx/max(temp));
  sapply(temp,function(x) min(x,mx));
}

pt = proc.time();
mat.tt = matrix(rep(0,num*num),num); #mat.ps = mat.tt;
for (i in 1:num) for (j in i:num) {
 mat.tt[i,j] = maxFill(XT[i,CDix],XT[j,CDix],XT[i,UCix],XT[j,UCix]);
 mat.tt[j,i] = mat.tt[i,j];
 #mat.ps[i,j] = polyserial(XT[i,] + XT[j,], Disease);
 #mat.ps[j,i] = mat.ps[i,j];
}
proc.time() - pt

#### actually run the loop -- gradient descent version
doBlocking = TRUE; thresh = 0.99; maxThresh = 0.5; path = NULL;
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
XT[low,] = maxFillT(XT[low,],XT[high,]);
rownames(XT)[low] = paste(rownames(XT)[low],rownames(XT)[high],sep='\n')

# remove the trailing extra row (and column, for the matrices based on XT)
XT = XT[-high,];
print(XT)
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
  { mat.tt[low,j] = maxFill(XT[low,CDix],XT[j,CDix],XT[low,UCix], XT[j,UCix]); }
  else mat.tt[low,j] = Inf;
}
mat.tt[,low] = mat.tt[low,];

numTests = numTests + mdim - sum(blocked);
iter = iter + 1; 
}

baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Groups))
names(baseline.tt) = rownames(XT);
proc.time() - ptm


###################### original additive version
########## t test fusion version: setup
Groups = table$Disease;
num = dim(XT)[1]
# whole = sample(seq(1,dim(table)[1]),dim(table)[1]); 
# part = 3/4; revP = 1-part;
# half = whole[seq(1,length(whole)*part)]; 
# XT = XT[,half]; Groups = Groups[half];
CDix = Groups == "CD"; UCix = Groups == "UC";

# extremely fast hybrid apply+for versions
library(polycor)
baseline.tt <- apply(XT,1,function(xx) t.test(xx[CDix], xx[UCix])$p.value)
#baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Disease)) # at end

pt = proc.time();
mat.tt = matrix(rep(0,num*num),num); #mat.ps = mat.tt;
for (i in 1:num) for (j in i:num) {
 mat.tt[i,j] = t.test((XT[i,CDix] + XT[j,CDix])*1, (XT[i,UCix] + XT[j,UCix])*1)$p.value;
 mat.tt[j,i] = mat.tt[i,j];
 #mat.ps[i,j] = polyserial(XT[i,] + XT[j,], Disease);
 #mat.ps[j,i] = mat.ps[i,j];
}
proc.time() - pt

# save(whole, half, XT, baseline.tt, mat.tt, file="halfsies.RData")

##### gradient descent (sum)
doBlocking = TRUE; thresh = 0.5; maxThresh = 0.25; path = NULL;
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
#corThres = min(cor(XT[low,CDix],XT[high,CDix]), cor(XT[low,UCix],XT[low,UCix]));

if (tt.both > baseline.tt[low]*thresh || tt.both > baseline.tt[high]*thresh) # || is.na(corThres) || corThres > -0.25) 
{
  mat.tt[low,high] = Inf;
  mat.tt[high,low] = Inf;
  #mdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
  #print(paste('skipping',low,high,"as max difference:",mdiff));
  next; 
}

found = found + 1; mdim = mdim - 1;
print(paste("found:",found,"at:",low,high,"NewP, low, high:",tt.both,baseline.tt[low], baseline.tt[high]));
path = c(path,low,high);
XT[low,] = (XT[low,] + XT[high,])*1;
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
  { mat.tt[low,j] = t.test((XT[low,CDix] + XT[j,CDix])*1.0,(XT[low,UCix] + XT[j,UCix])*1.0)$p.value; }
  else mat.tt[low,j] = Inf;
}
mat.tt[,low] = mat.tt[low,];

numTests = numTests + mdim - sum(blocked);
#print(XT)
iter = iter + 1; 
}

baseline.ps <- apply(XT,1,function(xx) polyserial(xx, Groups))
names(baseline.tt) = rownames(XT);
proc.time() - ptm


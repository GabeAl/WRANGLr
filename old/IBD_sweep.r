############## pairwise SWEEP version
######### test version -- pairwise diffs only
ptm <- proc.time()
thresh = 0.5; shearing = -1; minThresh = 0.1; 
PThresh = 0.1; Pshearing = 1; minPThresh = 0.05;
maxBatch = 20; doBlocking = T;

batch = 1; iter = 1; found = 0; blocked = rep(FALSE,num); 
numTests = num*(num-1)/2 + num; mdim = num;
lastFound = 0; Tdiff = thresh; minP = PThresh;


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
tt.both = mat.tt[low,high];
#ps.both = mat.ps[low,high];

#print(paste("Iteration:",iter, "mat.tt.dist:",mat.tt.dist[px,py]));

if (blocked[low] || blocked[high]) { mat.tt.dist[pair] = 0; next; }
if (low == high) {
  #mat.tt[low,] = Inf; 
  #mat.tt[,low] = Inf; 
  blocked[low] = TRUE;
  print(paste("match exact at:",low));
  break;
}
if (tt.both >= baseline.tt[low]*thresh || tt.both >= baseline.tt[high]*thresh) {
  #mat.tt[low,high] = Inf;
  #mat.tt[high,low] = Inf;
  mdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
  print(paste('skipping',low,high,"as max difference:",mdiff));
  if (batch < maxBatch || sum(blocked) > length(blocked) - 2) {
    print(paste("End of pass",batch,"reached. Resetting blocks. mdim:",mdim));
    batch = batch + 1;
    if (found == lastFound) break;
    blocked = rep(F,length(blocked));
    #for (i in 1:mdim) for (j in i:mdim) {
    #  mat.tt[i,j] = t.test((XT[i,CDix] + XT[j,CDix])*.5, (XT[i,UCix] + XT[j,UCix])*.5)$p.value;
    #  mat.tt[j,i] = mat.tt[i,j];
    #}
    mat.tt.distpre = apply(mat.tt,1,function(xx) baseline.tt / xx)
    mat.tt.dist = mat.tt.distpre;
    for (j in 1:mdim) for (k in j:mdim) { 
      mat.tt.dist[j,k] = min(mat.tt.dist[j,k],mat.tt.dist[k,j]); 
      mat.tt.dist[k,j]=mat.tt.dist[j,k]; 
    }
    thresh = max(minThresh,ifelse(shearing > 0, thresh * shearing, Tdiff));
    PThresh = max(minPThresh,ifelse(Pshearing > 0, PThresh * Pshearing, max(minP,PThresh*0.01)));
    lastFound = found;
    next;
  }
  break; 
}
if (tt.both > PThresh) { mat.tt.dist[pair] = 0; next; }


found = found + 1; mdim = mdim - 1;
if (doBlocking) {
  blocked[low] = T;
  #blocked[high] = T;
}

Tdiff = tt.both / min(baseline.tt[low],baseline.tt[high]);
minP = min(tt.both,minP);
 print(paste("found:",found,"at:",low,high,"Diff:",Tdiff, "New, l, h:",tt.both,baseline.tt[low], baseline.tt[high]));


XT[low,] = XT[low,] + XT[high,];
rownames(XT)[low] = paste(rownames(XT)[low],'\n',rownames(XT)[high],'\n ',sep='')

# remove the trailing extra row (and column, for the matrices based on XT)
XT = XT[-high,];
#mat.ps = mat.ps[-high,-high];
mat.tt = mat.tt[-high,-high];
#mat.tt.distpre = mat.tt.distpre[-high,-high];
mat.tt.dist = mat.tt.dist[-high,-high];
#baseline.ps = baseline.ps[-high];
baseline.tt = baseline.tt[-high];
blocked = blocked[-high];


# recalculate the new row and col's tt # and ps
baseline.tt[low] = tt.both;
#baseline.ps[low] = ps.both;
#low_s = as.numeric(XT[low,])
##mat.tt[low,] = apply(XT,1,function(xx) t.test((low_s + xx)[CDix], (low_s + xx)[UCix])$p.value)
# new version
for (j in 1:mdim) {
  #if (sum(mat.tt[j,]==Inf) < mdim - 2) 
  #if (!blocked[j]) 
  { mat.tt[low,j] = t.test(XT[low,CDix] + XT[j,CDix],XT[low,UCix] + XT[j,UCix])$p.value; }
  #else mat.tt[low,j] = Inf;
}
mat.tt[,low] = mat.tt[low,] # row and col the same

#mat.ps[low,] = apply(XT,1,function(xx) polyserial(low_s + xx, Disease))
#mat.ps[,low] = mat.ps[low,]

for (j in 1:mdim) {
 if (!blocked[j]) 
 { mat.tt.dist[low,j] = min(baseline.tt[j]/mat.tt[low,j],baseline.tt[low]/mat.tt[j,low]); }
 else { mat.tt.dist[low,j] = 0; }
}
mat.tt.dist[,low] = mat.tt.dist[low,];

numTests = numTests + mdim - 1; # really -1?
iter = iter + 1; #mdim = num - found;

}
baseline.ps <- apply(XT,1,function(xx) polyserial(xx, table$Disease))
names(baseline.tt) = rownames(XT);
proc.time() - ptm

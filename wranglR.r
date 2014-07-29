################ Within-Rank Analysis of Novel Groups by Likeness in R ###############
###################### 1987 Nei-Saitou-inspired special edition ######################

# Read map and taxonomy files in old QIIME format (no comments), OR provide directly
read_map_taxa = function(mapPath, taxaPath, GroupCol, taxThres)  {
	if (length(mapPath) > 1) { map = mapPath; }
	else if (!length(mapPath)) { map = 0; }
	else map = read.delim(mapPath,sep='\t',comment='',head=T,row.names=1,check.names=F);
	if (length(taxaPath) > 1) { taxa = taxaPath; } else taxa = 
		read.delim(taxaPath,sep='\t',comment='',head=T,row.names=1,check.names=F); 
	if (map) {
		# map = map[!is.na(match(rownames(map),colnames(taxa))),]; # asymmetry
		taxa = as.matrix(taxa[,match(rownames(map),colnames(taxa))]); 
	}
	if (length(GroupCol) > 1) { Groups = GroupCol; } else Groups = map[[GroupCol]];
	return(list(map = map, XT = taxa[rowMeans(taxa) > taxThres,], Groups = factor(Groups)));
}

# Setup function; this generates initial matrices upon which later algorithms operate
setup = function(mapPath, taxaPath, GroupCol, taxThres = 0, randSubset = 0, mat.tt = NULL, 
  A=T, hardCut = F, shift = 1) {
	data = read_map_taxa(mapPath, taxaPath, GroupCol, taxThres);
	  map = data$map; XT = data$XT; Groups = data$Groups;
	num = dim(XT)[1]; area = num * num; mdim = num; whole = seq(1,dim(XT)[2]);
	G1ix = Groups == levels(Groups)[1]; G2ix = Groups == levels(Groups)[2];
	if (randSubset) {  # Randomly subset a portion of samples to examine
		whole = sample(seq(1,dim(XT)[2]),dim(XT)[2]); 
		part = whole[seq(1,length(whole)*randSubset)]; 
		  XT = XT[,part]; Groups = Groups[part];
	}
	
	# Generate baseline t-test p-values between (binary) groups for each feature
	baseline.tt = apply(XT,1,function(xx) t.test(xx[G1ix], xx[G2ix])$p.value);
	
	# Generate matrix of pairwise additive t-test p values for each feature pair
	if (is.null(mat.tt)) {
		print("Generating baseline matrix...");
		mat.tt = matrix(rep(0,area),num); count = 0; timer = 0;
		for (i in 1:num) { 
			for (j in i:num) {
				mat.tt[i,j] = t.test(XT[i,G1ix] + XT[j,G1ix], XT[i,G2ix] + XT[j,G2ix])$p.value;
				mat.tt[j,i] = mat.tt[i,j];
				q = ifelse(j == i, 1, 2); 
				count = count + q; timer = timer + q;
			}
			if (timer >= area*0.1) 
				{ print(paste("Progress:",count/area*100,"%...")); timer = 0; }
		}
	}
	print("Baseline matrix generated.");
	
	# Calculate pairwise distances based on ratio of improvement
	mat.tt.dist = matrix(rep(0,area),num);
	for (i in 1:mdim) for (j in i:mdim) {
		if (A) { mat.tt.dist[i,j] = mat.tt[i,j] / (.5 * (baseline.tt[i] + baseline.tt[j])); }
		else mat.tt.dist[i,j] = 0.5*(mat.tt[i,j]/baseline.tt[j]+mat.tt[i,j]/baseline.tt[i]);
		if (mat.tt.dist[i,j] > 1) mat.tt.dist[i,j] = shift + 
			ifelse(hardCut, 0, log10(mat.tt.dist[i,j]));
		mat.tt.dist[j,i] = mat.tt.dist[i,j];
	}
	# stor = mat.tt.dist;
	
	# mat.tt.dist = apply(mat.tt,1,function(xx) xx / baseline.tt)
	# for (j in 1:mdim) for (k in j:mdim) { 
	 # mat.tt.dist[j,k] = 0.5*(mat.tt.dist[j,k]+mat.tt.dist[k,j]); 
	 # if (mat.tt.dist[j,k] > 1) mat.tt.dist[j,k] = 1+log10(mat.tt.dist[j,k]);
	 # mat.tt.dist[k,j]=mat.tt.dist[j,k]; 
	# }
	# print(paste("premat similarity:",sum(mat.tt.dist==stor)/length(stor)))
	
	# Calculate D, the Nei-Saitou neighbor-finder, based on improvement distances
	D = matrix(rep(0,area),num);
	r = rowSums(mat.tt.dist);
	for (i in 1:mdim) for (j in i:mdim) {
	  d = mat.tt.dist[i,j];
	  ri = 1 / (mdim - 2) * r[i] - d;
	  rj = 1 / (mdim - 2) * r[j] - d;
	  D[i,j] = d - ri - rj;
	  D[j,i] = D[i,j];
	}
	
	## faster version, only use if not reusing ri, rj in later dist calculation
	# for (i in 1:mdim) for (j in i:mdim) {
	 # D[i,j] = (mdim - 2) * mat.tt.dist[i,j] - r[i] - r[j];
	 # D[j,i] = D[i,j];
	# }
	
	for (i in 1:mdim) D[i,i] = Inf; # Set diagonal to infinity to avert self-selection
	return(list(G1ix = G1ix, G2ix = G2ix, XT = XT, baseline.tt = baseline.tt, mat.tt = mat.tt,
			mat.tt.dist = mat.tt.dist, D = D, r = r, randSubset = randSubset, whole = whole));
}

# The main function to be called from R or commandline; only binary groups are supported
wranglR = function(mapPath, taxaPath, GroupCol = "Group", taxThres = 0, randSubset = 0, dat=NA,
  A=T, hardCut = F, shift = 1) {
	if (!is.na(dat)) { dt = dat; }
	else dt = setup(mapPath, taxaPath, GroupCol, taxThres, randSubset, NULL, A, hardCut, shift);
	  G1ix = dt$G1ix; G2ix = dt$G2ix; XT = dt$XT; baseline.tt = dt$baseline.tt; D = dt$D;
	  mat.tt = dt$mat.tt; mat.tt.dist = dt$mat.tt.dist; randSubset = dt$randSubset; r = dt$r;
	num = dim(XT)[1]; mdim = num;   # Cache the initial number of species/entries
	  iter = 1; found = 0; 
	
	numTests = num*(num-1)/2 + num; # The number of t-tests performed thus far
	while (mdim > 2) {
		pair = which.min(D);
		px = ceiling(pair/mdim); py = pair - mdim * (px - 1);
		low = min(px,py); high = max(px, py); # to shrink matrix later
		tt.both = mat.tt[low,high]; 
		if (low == high) { print(paste("match exact at:",low)); break; } # Error condition
		
		found = found + 1; 
		XT[low,] = XT[low,] + XT[high,]; # Merge rows into lower ix row
		# distG = 0.5 * (mat.tt.dist[low,high] + (1 / (mdim - 2) * r[low] - 
			# mat.tt.dist[low,high]) - (1 / (mdim - 2) * r[high] - mat.tt.dist[low,high]));
		distG = mat.tt.dist[low,high];
		print(paste("found:",found,"at:",low,high,"Dist:",mat.tt.dist[low,high], 
			"New, l, h:",tt.both,baseline.tt[low], baseline.tt[high]));
		# Assign the new row a Newick name with distances according to Nei-Saitou NJ
		rownames(XT)[low] = paste('(',rownames(XT)[low],':',distG,',',
			rownames(XT)[high],':',distG,')',sep='');
		# Remove the trailing extra row (and column, for the matrices based on XT)
		XT = XT[-high,];
		mat.tt = mat.tt[-high,-high];
		mat.tt.dist = mat.tt.dist[-high,-high];
		baseline.tt = baseline.tt[-high];
		mdim = mdim - 1;
		
		# The baseline of this row becomes this round's combined score
		baseline.tt[low] = tt.both;
		for (j in 1:mdim) { 
		  mat.tt[low,j] = t.test(XT[low,G1ix] + XT[j,G1ix],XT[low,G2ix] + XT[j,G2ix])$p.value; 
		  if (A) { mat.tt.dist[low,j] = mat.tt[low,j] / (.5 * (baseline.tt[j] + baseline.tt[low])); }
		  else mat.tt.dist[low,j] = 0.5*(mat.tt[low,j]/baseline.tt[j]+mat.tt[low,j]/baseline.tt[low]);
		  if (mat.tt.dist[low,j] > 1) mat.tt.dist[low,j] = shift + 
		  	ifelse(hardCut, 0, log10(mat.tt.dist[low,j]));
		}
		mat.tt[,low] = mat.tt[low,] # row and col the same
		mat.tt.dist[,low] = mat.tt.dist[low,];
		# stor = mat.tt.dist;
		
		# mat.tt.dist = matrix(rep(0,mdim*mdim),mdim);
		# for (i in 1:mdim) for (j in i:mdim) {
			# if (A) { mat.tt.dist[i,j] = mat.tt[i,j] / (.5 * (baseline.tt[i] + baseline.tt[j])); }
			# else mat.tt.dist[i,j] = 0.5*(mat.tt[i,j]/baseline.tt[j]+mat.tt[i,j]/baseline.tt[i]);
			# if (mat.tt.dist[i,j] > 1) mat.tt.dist[i,j] = 1 + log10(mat.tt.dist[i,j]);
			# mat.tt.dist[j,i] = mat.tt.dist[i,j];
		# }
		# print(sum(mat.tt.dist==stor)/length(stor))
		
		# This bit repeats the D calculation code. Functionalize?
		D = matrix(rep(0,mdim*mdim),mdim);
		r = rowSums(mat.tt.dist);
		for (i in 1:mdim) for (j in i:mdim) {
			d = mat.tt.dist[i,j];
			ri = 1 / (mdim - 2) * r[i] - d;
			rj = 1 / (mdim - 2) * r[j] - d;
			D[i,j] = d - ri - rj;
			D[j,i] = D[i,j];
		}
		# for (i in 1:mdim) for (j in i:mdim) {
			# D[i,j] = (mdim - 2) * mat.tt.dist[i,j] - r[i] - r[j];
			# D[j,i] = D[i,j];
		# }
		for (i in 1:mdim) D[i,i] = Inf;  # Set diagonal to infinity
		
		numTests = numTests + mdim - 1;  # Increment the number of tests performed
		iter = iter + 1; 
	}
	names(baseline.tt) = rownames(XT);   # Assign new names back to tt object
	
	distG = mat.tt.dist[1,2]; # Set the final distance to the remaining
	tree_pre = paste('(',rownames(XT)[1],':',distG,',',rownames(XT)[2],':',distG,')',sep='');
	
	tree = paste(gsub(";","-",tree_pre),";",sep="");
	tree = gsub('[[]',"<",tree);
	tree = gsub('[]]',">",tree);
	tree = gsub(" ","_",tree);
	return(list(tree = tree, baseline.tt = baseline.tt, combinedData = XT));
}

# Re-combine any dataset of same row order from given path (w/o branch lengths)
fitTree = function(YT, path) {
	endp = length(path) / 2; 
	for (i in 1:endp) {
	 high = path[i*2]; low = path[i*2 - 1];
	 YT[low,] = YT[low,] + YT[high,];
	 rownames(YT)[low] = paste('(',rownames(YT)[low],',',
			rownames(XT)[high],')',sep='');
	 YT = YT[-high,];
	}
}

# Cuts the tree at a given branch depth and returns subtrees with paths
cutTree = function (XT, depth, useDist = F) {
	
}

function encode(images,names,destFolder)
% *****************************************************
% Encode newly observed images by
% 		displaying activated templates.
% *****************************************************

showPartBoundingBox = true;

% make sure all input images are gray-scale
for i = 1:length(images)
	if size(images{i},3) == 3
		images{i} = rgb2gray(images{i});
	end
end

storeExponentialModelName = ['storedExponentialModel' num2str(1)];   
load(storeExponentialModelName);
load partLocConfig templateSize numCluster numIter nTransform numOrient ...
	rotationRange templateTransform localOrNot ...
	doubleOrNot thresholdFactor saturation  locationShiftLimit ...
	orientShiftLimit S1softthres

subsampleM1 = 1; S2Thres = 0;
subsampleS2 = 2; % adjust this for speed (larger subsample step will make it faster)
locationPerturbationFraction = .5; % how large is the inhibition neighborhood
partSizeX = templateSize(1);
partSizeY = templateSize(2);

%% preparation
% load and transform the templates
S2Templates = cell(numCluster,1);
iter = numIter;
for cc = 1:numCluster
	load(sprintf('working/learnedmodel%d_iter%d.mat',cc,iter), 'numElement', 'selectedOrient', 'selectedx', 'selectedy', 'selectedlambda', 'selectedLogZ');
	S2Templates{cc} = struct( 'selectedRow', single(selectedx -1 - floor(templateSize(1)/2)),...
		'selectedCol', single(selectedy -1 - floor(templateSize(2)/2)), ...
		'selectedOri', single(selectedOrient), 'selectedScale', zeros(length(selectedx),1,'single'), ...
		'selectedLambda', single(selectedlambda), 'selectedLogZ', single(selectedLogZ) );
end
TransformedTemplate = cell(nTransform,numCluster);
selectedScale = zeros(1,length(selectedx),'single');
for cc = 1:numCluster
	for iT = 1:nTransform
		templateScaleInd = templateTransform{iT}(1);
		rowScale = templateTransform{iT}(2);
		colScale = templateTransform{iT}(3);
		rotation = templateTransform{iT}(4);
		[tmpSelectedRow tmpSelectedCol tmpSelectedOri tmpSelectedScale] = ...
			mexc_TemplateAffineTransform( templateScaleInd, rowScale,...
			colScale, rotation, S2Templates{cc}.selectedRow, S2Templates{cc}.selectedCol,...
			S2Templates{cc}.selectedOri, selectedScale, numOrient );
		TransformedTemplate{iT,cc}.selectedRow = tmpSelectedRow;
		TransformedTemplate{iT,cc}.selectedCol = tmpSelectedCol;
		TransformedTemplate{iT,cc}.selectedOri = tmpSelectedOri;
		TransformedTemplate{iT,cc}.selectedScale = tmpSelectedScale;
		TransformedTemplate{iT,cc}.selectedLambda = S2Templates{cc}.selectedLambda;
		TransformedTemplate{iT,cc}.selectedLogZ = S2Templates{cc}.selectedLogZ;
	end
end



%% begin encoding
% compute SUM1 maps for all images
allSUM1map = ApplyFilterfftSame(images, allFilter, localOrNot, localHalfx, localHalfy, doubleOrNot, thresholdFactor);
for ii = 1:numel(allSUM1map)
	allSUM1map{ii} = single(allSUM1map{ii}); % soft thresholding
	allSUM1map{ii}(:) = max(0,allSUM1map{ii}(:)-S1softthres);
end
buf_length = 0;
for img = 1:length(images)

	for k = 1:buf_length
		fprintf(1,'\b');
	end
	str = sprintf('dealing with: %s',names{img});
	fprintf(1,str);
	buf_length = length(str);

	% compute MAX1 map
	SUM1map = allSUM1map(img,:);
	[MAX1map ARGMAX1map M1RowShift M1ColShift M1OriShifted] = ...
		mexc_ComputeMAX1( numOrient, SUM1map, locationShiftLimit,... %single(locationShiftLimit)/(halfFilterSize*2+1),...
		orientShiftLimit, subsampleM1 );
	for ii = 1:numel(MAX1map)
		MAX1map{ii} = single( saturation*(2./(1.+exp(-2.*MAX1map{ii}/saturation))-1.) );
	end
	
	% compute SUM2/MAX2 and inhibit (explain away) using matching pursuit
	SUM2map = mexc_ComputeSUM2( numOrient, MAX1map, TransformedTemplate, subsampleS2 );
	% random perturbation (to break ties arbitrarily for MAX2)
	for ii = 1:numel(SUM2map)
		SUM2map{ii}(:) = SUM2map{ii}(:) + 1e-3 * ( rand(numel(SUM2map{ii}),1) - .5 );
	end
	% exclude the near-boundary region
	for ii = 1:numel(SUM2map)
		SUM2map{ii}(:,[1:floor(templateSize(2)/2) end-floor(templateSize(2)/2):end]) = min(-1001,S2Thres-1);
		SUM2map{ii}([1:floor(templateSize(1)/2) end-floor(templateSize(1)/2):end],:) = min(-1001,S2Thres-1);
	end
	% find the activations [row; col; templateInd; score]
	activations = mexc_ComputeMAX2MP( SUM2map, int32(locationPerturbationFraction*sqrt(partSizeX*partSizeY)/subsampleS2), S2Thres );
	activations(1:2,:) = activations(1:2,:) * subsampleS2;
	
	% initialize the mask of bounding boxes
	matchedBoundingBox = zeros([size(images{img}),3]);
	% track Gabor elements of activated templates
	gaborXX = [];
    gaborYY = [];
    gaborOO = [];
    gaborMM = [];
    % Gabor basis elements locations
    for t = 1:size(activations,2) % for each activated template (or part)
    	iCluster = ceil( ( activations(3,t) + 1 ) / nTransform );
    	iTransform = activations(3,t) + 1 - (iCluster-1) * nTransform;
    	for j = 1:numElement
		    gaborX = floor(activations(1,t) + TransformedTemplate{iTransform,iCluster}.selectedRow(j));
		    gaborY = floor(activations(2,t) + TransformedTemplate{iTransform,iCluster}.selectedCol(j));
		    gaborO = TransformedTemplate{iTransform,iCluster}.selectedOri(j);
		    if gaborX > 0 && gaborX <= size(MAX1map{1},1) && gaborY > 0 && gaborY <= size(MAX1map{1},2)
		        trace = ARGMAX1map{gaborO+1}(gaborX,gaborY) + 1;
		        dx = M1RowShift{gaborO+1}(trace);
		        dy = M1ColShift{gaborO+1}(trace);
		        shiftedo = M1OriShifted{gaborO+1}(trace);
		        gaborX = floor(.5 + gaborX + single(dx));
		        gaborY = floor(.5 + gaborY + single(dy));
		        gaborO = single(shiftedo);
		    end
		    gaborXX = [gaborXX;gaborX];
		    gaborYY = [gaborYY;gaborY];
		    gaborOO = [gaborOO;gaborO];
		    if gaborX > 0 && gaborX <= size(MAX1map{1},1) && gaborY > 0 && gaborY <= size(MAX1map{1},2)
		        val = SUM1map{gaborO+1}(gaborX,gaborY);
		    else
		        val = 0;
		    end
		    gaborMM = [gaborMM; max(0,sqrt(val)-.2)];
        end
        
        if showPartBoundingBox
		    margin = 2;
		    xx = repmat((1:partSizeX),1,margin*2);
		    yy = [];
		    for y = [1:margin partSizeY-margin+1:partSizeY]
		        yy = [yy,ones(1,partSizeX)*y];
		    end
		    yy = [yy,repmat((1:partSizeY),1,margin*2)];
		    for x = [1:margin partSizeX-margin+1:partSizeX]
		        xx = [xx,ones(1,partSizeY)*x];
		    end
		    inRow = single(xx-floor(partSizeX/2)); inCol = single(yy-floor(partSizeY/2));
		    tScale = 0; rScale = 1; cScale = 1; inO = zeros(numel(inRow),1,'single'); inS = zeros(numel(inRow),1,'single');
		    actualPartRotation = iTransform-1;
		    [outRow, outCol] = ...
		        mexc_TemplateAffineTransform(tScale,rScale,cScale,...
		            actualPartRotation,inRow,inCol,inO,inS,numOrient);

		    % directly overwrite the corresponding pixels
		    for p = 1:length(outRow)
		        x = floor(.5 + outRow(p) + activations(1,t)); y = floor(.5 + outCol(p) + activations(2,t));
		        if x > 0 && x <= size(matchedBoundingBox,1) && y > 0 && y <= size(matchedBoundingBox,2)
		            matchedBoundingBox(x,y,:) = [.5 .9 .6];
		        end
		    end
		end
    end
    
	tmpMatchedSym = displayMatchedTemplate(size(images{img}),gaborXX,...
    	gaborYY,gaborOO,zeros(length(gaborXX),1,'single'),gaborMM,allSymbol,numOrient);
    
    % overlay
    matchedSym = repmat(-single(tmpMatchedSym),[1 1 3]);
    matchedSym = 1 * (matchedSym-min(matchedSym(:)))/(max(matchedSym(:))-min(matchedSym(:)));
    if showPartBoundingBox
		for y = 1:size(matchedSym,2)
		    for x = 1:size(matchedSym,1)
		        if sum(abs(matchedBoundingBox(x,y,:))) > 0
		            matchedSym(x,y,:) = matchedBoundingBox(x,y,:);
		        end
		    end
		end
    end
	
	tmp = single(repmat(images{img},[1 1 3])) / 255;
	alpha = .6;
    tmp = tmp * alpha + matchedSym * (1-alpha);
    imwrite( tmp, sprintf('%s/overlayed_%s.png',destFolder,names{img}) );
    imwrite( matchedSym, sprintf('%s/sketch_%s.png',destFolder,names{img}) );
end
fprintf(1,'\n');


% *****************************************************
% display activated templates on selected images
% *****************************************************

storeExponentialModelName = ['storedExponentialModel' num2str(1)];   
load(storeExponentialModelName);
load partLocConfig

load activations
showPartBoundingBox = true;
partSizeX = templateSize(1);
partSizeY = templateSize(2);

%% preparation
% transform the templates
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



%% begin displaying
buf_length = 0;
for img = selected_img

	for k = 1:buf_length
		fprintf(1,'\b');
	end
	str = sprintf('%d',img);
	fprintf(1,str);
	buf_length = length(str);

	% load SUM1/MAX1 map
	SUM1MAX1mapName = ['working/SUM1MAX1map' 'image' num2str(img) 'scale' num2str(1)];
	load(SUM1MAX1mapName, 'SUM1map', 'MAX1map', 'ARGMAX1map', 'M1RowShift', 'M1ColShift',...
		'M1OriShifted', 'J');
	% find activations
	ind = find( activations(1,:) == img );
	
	% initialize the mask of bounding boxes and sketched image
	matchedBoundingBox = zeros([size(J{end}),3]); % the bounding box map is of the highest resolution
	matchedSym = zeros(size(J{end})); % the sketched image is of the highest resolution
	
    % Gabor basis elements locations
    for t = ind % for each activated template (or codeword, or part)
	
		% track Gabor elements of activated templates for this template only
		gaborXX = [];
		gaborYY = [];
		gaborOO = [];
		gaborMM = [];
	
    	iCluster = ceil( ( activations(5,t) + 1 ) / nTransform );
    	iTransform = activations(5,t) + 1 - (iCluster-1) * nTransform;
		iRes = activations(4,t) + 1;
    	for j = 1:numElement
		    gaborX = floor(activations(2,t) + TransformedTemplate{iTransform,iCluster}.selectedRow(j));
		    gaborY = floor(activations(3,t) + TransformedTemplate{iTransform,iCluster}.selectedCol(j));
		    gaborO = TransformedTemplate{iTransform,iCluster}.selectedOri(j);
		    if gaborX > 0 && gaborX <= size(MAX1map{iRes,1},1) && gaborY > 0 && gaborY <= size(MAX1map{iRes,1},2)
		        trace = ARGMAX1map{iRes,gaborO+1}(gaborX,gaborY) + 1;
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
		    if gaborX > 0 && gaborX <= size(MAX1map{iRes,1},1) && gaborY > 0 && gaborY <= size(MAX1map{iRes,1},2)
		        val = SUM1map{iRes,gaborO+1}(gaborX,gaborY);
		    else
		        val = 0;
		    end
		    gaborMM = [gaborMM; max(0,sqrt(val)-.2)];
        end
		
		% render the sketches for this activated template
		tmpMatchedSym = displayMatchedTemplate([size(J{iRes},1) size(J{iRes},2)],gaborXX,...
			gaborYY,gaborOO,zeros(length(gaborXX),1,'single'),gaborMM,allSymbol,numOrient);
		
        
		scaling = double(size(J{end},1)) / double(size(J{iRes},1));
		tmpMatchedSym = imresize(tmpMatchedSym,scaling,'nearest');
		matchedSym = max( matchedSym, tmpMatchedSym );
		
        if showPartBoundingBox
		    margin = 2;
			largerPartSizeX = floor(partSizeX * scaling+.5);
			largerPartSizeY = floor(partSizeY * scaling+.5);
		    xx = repmat((1:largerPartSizeX),1,margin*2);
		    yy = [];
		    for y = [1:margin largerPartSizeY-margin+1:largerPartSizeY]
		        yy = [yy,ones(1,largerPartSizeX)*y];
		    end
		    yy = [yy,repmat((1:largerPartSizeY),1,margin*2)];
		    for x = [1:margin partSizeX-margin+1:partSizeX]
		        xx = [xx,ones(1,largerPartSizeY)*x];
		    end
		    inRow = single(xx-floor(largerPartSizeX/2)); inCol = single(yy-floor(largerPartSizeY/2));
		    tScale = 0; rScale = 1; cScale = 1; inO = zeros(numel(inRow),1,'single'); inS = zeros(numel(inRow),1,'single');
		    actualPartRotation = iTransform-1;
		    [outRow, outCol] = ...
		        mexc_TemplateAffineTransform(tScale,rScale,cScale,...
		            actualPartRotation,inRow,inCol,inO,inS,numOrient);

		    % directly overwrite the corresponding pixels
		    for p = 1:length(outRow)
		        x = floor(.5 + outRow(p) + activations(2,t)*scaling); y = floor(.5 + outCol(p) + activations(3,t)*scaling);
		        if x > 0 && x <= size(matchedBoundingBox,1) && y > 0 && y <= size(matchedBoundingBox,2)
		            matchedBoundingBox(x,y,:) = [.5 .9 .6];
		        end
		    end
		end
    end
    
	
    % overlay
    matchedSym = repmat(-single(matchedSym),[1 1 3]);
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
	
	tmp = single(repmat(J{1},[1 1 3])) / 255;
	alpha = .6;
    tmp = tmp * alpha + matchedSym * (1-alpha);
    imwrite( tmp, sprintf('%s/overlayed_image%d.png',destFolder,img) );
    imwrite( matchedSym, sprintf('%s/sketch_image%d.png',destFolder,img) );
end
fprintf(1,'\n');



disp(['===============> Getting parameters, mex-c code, and training images']);
mkdir('working'); % working directory to store .mat files. 
mkdir('output'); % directory that contains output from the code. 

%% parameters for the outer square bounding box
category = '1';
supressionModeInEStep = 'MatchingPursuit'; % matching pursuit in E step (other options: LocalSurroundSurpression)
templateSize = [60 60];
templateSize = single(templateSize);
partSize = floor(sqrt(templateSize(1)*templateSize(2))); % alias for the template size (radius)
numRandomStart = 3;
%% parameters for EM clustering
locationPerturbationFraction = .25; % the size of neighborhood for MAX2 pooling, as well as surround supression
rotationRange = 0:2:0; % allowed global rotation of the template
subsampleS2 = 1; % subsampling step size for computing SUM2 maps
maxNumClusterMember = 30; % maximum number of training examples in each cluster used in re-learning the template
S2Thres = 5; % cut-off value for detected instances of the template
S1softthres = 0.9; % soft thresholding cutting-off for S1 maps
numIter = 10; % number of EM iterations
%% parameters for active basis
numCluster = 10; % number of data clusters
numElement = 8; % number of Gabors in active basis at the first scale
epsilon = .1; % allowed correlation between selected Gabors 
subsample = 1; subsampleM1 = 1; % subsample in computing MAX1 maps
locationShiftLimit = 3; % shift in normal direction = locationShiftLimit*subsample pixels
orientShiftLimit = 1; % shift in orientation
%% parameters for Gabor filters
numScale = 1;
scales = 0.7; % scales of Gabor wavelets 
numOrient = 16;  % number of orientations
saturation = 6.; % saturation level for sigmoid transformation
doubleOrNot = -1;
%% parameters for exponential model
binSize = .2;  % binsize for computing histogram of q()
numStoredPoint = 50; % number of stored lambda values
spacing = .1; % spacing between two adjacent points
%% parameters for normalization
localOrNot = 1; % if we use local normalization or not. If not, set it to -1 
localHalfx = 20; localHalfy = 20; % the half range for local normalization, has to be quite large
windowNormalizeOrNot = -1; % whether normalize within the scanning window in detection 
if (localOrNot>0)
    windowNormalizeOrNot = -1;
end % if we use local normalization, we should not use window normalization in detection
thresholdFactor = .01;  % divide the response by max(average, maxAverage*thresholdFactor)

%% parameters for detection
inhibitFind = -1;  % whether to apply inhibition after detection for re-computing MAX2 score
resolutionGap = .1; % gap between consecutive resolutions in detection
numExtend = 0; % number of gaps extended both below and above zero
numResolution = numExtend*2 + 1;  % number of resolutions to search for in detection stage
originalResolution = numExtend + 1; % original resolution is the one at which the imresize factor = 1
allResolution = (-numExtend : numExtend)*resolutionGap + 1.;
%% read in positive images
sizeTemplatex = templateSize(1);
sizeTemplatey = templateSize(2);
halfTemplatex = floor(sizeTemplatex/2);
halfTemplatey = floor(sizeTemplatey/2);
imageFolder = 'positiveImage'; % folder of training images  
imageName = dir([imageFolder '/*.jpg']);
numImage = size(imageName, 1); % number of training images 
Ioriginal = cell(1, numImage);
for img = 1 : numImage
    tmpIm = imread([imageFolder '/' imageName(img).name]); 
    if size(tmpIm,3) == 3
        tmpIm = rgb2gray(tmpIm);
    end

    sx = size(tmpIm,1); sy = size(tmpIm,2);
    tmpIm = imresize( tmpIm, 600/sqrt(sx*sy), 'bilinear' );
    Ioriginal{img} = single(tmpIm);
    J0 = Ioriginal{img};
    J = cell(1, numResolution);
    for r=1:numResolution
       J{r} = imresize(J0, allResolution(r), 'nearest');  % scaled images
    end
    multipleResolutionImageName = ['working/multipleResolutionImage' num2str(img)]; 
    save(multipleResolutionImageName, 'J');
end

count = 0;
% generate the set of geometric transformations for each template
for templateScaleInd = 0:0 % large scale change
    for rotation = rotationRange 
        for rowScale = 2.^[0] % small scale change and aspect ratio change
            for colScale = 2.^[0]
                count = count + 1;
                templateTransform{count} = [templateScaleInd rowScale colScale rotation];
            end
        end
    end
end
nTransform = count;

% prepare morph-back mappings for geometric transforms
largestPartSizeX = templateSize(1); largestPartSizeY = templateSize(2);
denseX = -floor(largestPartSizeX/2) + (1:largestPartSizeX);
denseY = -floor(largestPartSizeY/2) + (1:largestPartSizeY);
count = 0;
inRow = zeros(length(denseX)*length(denseY),1,'single');
inCol = zeros(length(denseX)*length(denseY),1,'single');
inO = zeros(numel(inRow),1,'single'); 
inS = zeros(numel(inRow),1,'single');
for y = denseY
    for x = denseX
        count = count+1;
        inRow(count) = x;
        inCol(count) = y;
    end
end
outRow = cell(nTransform,1);
outCol = cell(nTransform,1);
for iT = 1:nTransform
	tScale = templateTransform{iT}(1); 
	rScale = templateTransform{iT}(2);
	cScale = templateTransform{iT}(3);
	[outRow{iT}, outCol{iT}] = ...
		mexc_TemplateAffineTransform(tScale,rScale,cScale,...
		    templateTransform{iT}(4),inRow,inCol,inO,inS,numOrient);
end

save('partLocConfig.mat','templateSize',...
    'category','numOrient','localOrNot','subsample','saturation','locationShiftLimit','orientShiftLimit',...
    'numElement','thresholdFactor','doubleOrNot', 'numCluster', 'numIter', 'subsampleS2', 'locationPerturbationFraction',...
    'partSize','rotationRange','nTransform','templateTransform','S1softthres','subsampleM1','localHalfx','localHalfy',...
    'inRow', 'inCol', 'outRow', 'outCol', 'largestPartSizeX', 'largestPartSizeY');
clear Ioriginal


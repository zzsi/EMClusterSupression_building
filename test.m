% test on new images

imname = dir('testImage/*.jpg');
names = cell(length(imname),1);
for i = 1:length(imname)
	names{i} = imname(i).name;
end
images = cell(length(names),1);
for i = 1:length(names)
	im = imread(sprintf('testImage/%s',names{i}));
	if size(im,3) == 3
		im = rgb2gray(im);
	end
	sx = size(im,1); sy = size(im,2);
    im = imresize( im, 600/sqrt(sx*sy), 'bilinear' );
	images{i} = single(im);
end

destfolder = 'output_test';
if ~exist(destfolder,'dir')
	mkdir(destfolder);
end
encode(images,names,destfolder);



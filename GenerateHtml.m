% GenerateHtml - Generates html documentation for the learning results.
%

clear
close all;
maxDisplayImg = 60;

% load the starting image number
load partLocConfig templateSize category locationShiftLimit orientShiftLimit numElement numCluster numIter

zipname = sprintf('EMScanSupression_%s.zip',date);
imFolder = 'EMScanSupression';

% delete the previous version
html_dir = 'document';
if ~exist(html_dir,'dir')
    mkdir(html_dir);
else
%     rmdir(html_dir,'s');  % TODO: find out why removing existing dir causes error (in GenerateHtml.m) 
%     mkdir(html_dir);
end

html_path = sprintf('%s/%s.html',html_dir,imFolder);
html = fopen(html_path,'w');

html_img_dir = sprintf('document/%s/',imFolder);
if ~exist(html_img_dir,'dir')
    mkdir(html_img_dir);
end

%% html header
% note: modified on Jan 13, 2010. Some url's are now absolute. Only 
%   copied files are linked with relative url.
tmp = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n';
tmp = [tmp '<html>\n'];
tmp = [tmp '<head>\n'];
tmp = [tmp '<title>Discovering a visual dictionary of active basis templates by EM clustering and sparsifiction</title>\n'];
tmp = [tmp '<link rel="stylesheet" href="http://www.stat.ucla.edu/~zzsi/plain_and_simple.css" type="text/css" media="screen" />\n'];
tmp = [tmp '<script type="text/javascript" src="http://www.stat.ucla.edu/~zzsi/1.js"></script>\n'];
tmp = [tmp '</head>\n'];
fprintf(html,tmp);
fprintf(html, '<body>\n');
fprintf( html, '<div id="header">\n');
fprintf( html, '<h1>Discovering a visual dictionary of active basis templates by EM clustering and sparsifiction</h1></div>\n' );
fprintf( html, '\n<div id="content">\n');

%% link to project page
fprintf( html, '\n<p><a href="http://www.stat.ucla.edu/~zzsi/HAB/exp2.html">Exp2 Home</a></p>\n' );

%% table of content
fprintf( html, '<div id="content">\n' );
fprintf( html, '<div id="TableOfContents">\n' );
fprintf( html, '<p>Contents</p>\n' );
fprintf( html, '<ul>\n' );
fprintf( html, '<li>\n' );
fprintf( html, '<A href="#download">Download</A>\n' );
fprintf( html, '</li>\n' );
fprintf( html, '<li>\n');
fprintf( html, '<a href="#traindata">Training examples</a>\n' );
fprintf( html, '</li>\n' );
fprintf( html, '<li>\n' );
fprintf( html, '<a href="#templates">Learned templates</a>\n' );
fprintf( html, '</li>\n' );
fprintf( html, '</ul>\n' );
fprintf( html, '</div>\n' );

%% for download
fprintf( html, '<div style="border-bottom:1 solid #dddddd; margin-top:0.3em;"></div>\n' );
fprintf( html, '<a name="download"></a> <table cellspacing="10" cellpadding="10" class="center" width="60%%">\n' );
fprintf( html, '\n<tr><td>\n' );
fprintf( html, '\n<b>Code and data: (<a href="%s">ZIP</a>).</b>\n', zipname );
fprintf( html, '\n</td>\n' );
fprintf( html, '\n<td>\n' );
fprintf( html, '\n<a href="http://www.stat.ucla.edu/~zzsi/HAB/hab_changelog.html">Change Log</a>\n' );
fprintf( html, '\n</td>\n' );
fprintf( html, '\n</tr>\n' );
fprintf( html, '\n<tr>\n' );
fprintf( html, '\n<td colspan=2 align=left>\n' );
fprintf( html, '\nRun StartFromHere.m in Matlab. You can monitor intermediate results in the folder: output/. \n' );
fprintf( html, '\n</td>\n' );
fprintf( html, '\n</tr>\n' );
fprintf( html, '\n</table>\n' );

%% training examples
fprintf(html, '<div style="border-top:1 solid #dddddd; margin-top:0.3em;"></div><h2> ');
fprintf(html, '<a name="traindata"></a>Training examples');
fprintf(html, ' </h2>');
fprintf( html, ['\n<p>A selection of the input images: </p>\n']);
% read the training examples
Iname = dir('positiveImage/*.jpg');
n = length(Iname);

if n > maxDisplayImg
	ind = randperm(n);
	selected_img = ind(1:maxDisplayImg);
	Iname = Iname(selected_img);
else
	selected_img = 1:n;
end

% render sketched images with activated partial templates
destFolder = sprintf('document/%s',imFolder);
displayActivations;

% move the images to img/ folder
for i = 1:length(Iname)
    new_img_name = Iname(i).name;
    src = sprintf('positiveImage/%s',new_img_name);
    dst = sprintf('document/%s/%s',imFolder, new_img_name);
    copyfile(src,dst);
end

% generate corresponding html
fprintf( html, '\n<p>' );
for i = 1:length(Iname)
    fprintf( html, '<img src="%s" alt="" height=80/>', sprintf('%s/%s',imFolder,Iname(i).name) );
end
fprintf( html, '\n</p>\n' );

%% show learned templates
fprintf( html, ['<div style="border-bottom:1 solid #dddddd; margin-top:0.3em;"></div>\n<a name="templates"></a><p>Learned templates for ' num2str(numCluster) ' clusters (after ' num2str(numIter) ' iterations)' ...
    ' for image patches randomly cropped/scanned from the input images. (some clusters may be empty): </p>\n']);
fprintf( html, '\n<p>' );
new_img_name = sprintf('template_sorted.png');
src = sprintf('./%s',new_img_name);
dst = sprintf('document/%s/%s', imFolder,new_img_name);
copyfile(src,dst);
fprintf( html, '\n<img src="%s" alt="" width="500"/> <br> <br>', sprintf('%s/%s',imFolder,new_img_name) );
fprintf( html, '\n</p>\n' );

% also show the intialized templates
fprintf( html, ['<div style="margin-top:0.3em;"></div>\n<a name="templates"></a><p>Initial templates learned from randomly initialized clusters:\n']);
fprintf( html, '\n<p>' );
new_img_name = sprintf('template_iter1.png');
src = sprintf('output/%s',new_img_name);
dst = sprintf('document/%s/%s', imFolder,new_img_name);
copyfile(src,dst);
fprintf( html, '\n<img src="%s" alt="" width="500"/> <br> <br>', sprintf('%s/%s',imFolder,new_img_name) );
fprintf( html, '\n</p>\n' );

%% show sketched images
fprintf( html, ['<div style="border-bottom:1 solid #dddddd; margin-top:0.3em;"></div>\n<a name="templates"></a><p>' ...
	'Sketching the observed images by overlaying the activated templates on them:' ...
	'</p>\n']);
fprintf( html, '\n<p>' );
for i = 1:length(Iname)
	fprintf( html, '<img src="%s" alt="" height=80/>', sprintf('%s/overlayed_image%d.png',imFolder,i) );
end
fprintf( html, '\n</p>\n' );

fprintf( html, '\n<p>Showing only the activated templates (for visual clarity): </p><p>' );
for i = 1:length(Iname)
	fprintf( html, '<img src="%s" alt="" height=80/>', sprintf('%s/sketch_image%d.png',imFolder,i) );
end
fprintf( html, '\n</p>\n' );


%% finishing off
fprintf( html, '\n\n\n</div> \n');
fprintf( html, '<div id="last" class="footer"></div>' );
fprintf(html, '</body> </html> \n');
fclose(html);
disp('finished generating Html.. check here:');
disp([pwd '/' html_path]);



% ZIPproject - create a zip file of the project

zipfilename = sprintf('EMScanSupression_%s.zip',date);

%% specify files to put in ZIP
filelist = {};
count = 1;

% .m files in current folder
names = dir('*.m');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% .c files in current folder
names = dir('*.c');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% .cpp files in current folder
names = dir('*.cpp');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% .h files in current folder
names = dir('*.h');
for i = 1:length(names)
    filelist{count} = sprintf('%s',names(i).name);
    disp([num2str(count) ': added ' filelist{count}]);
    count = count + 1;
end

% positive images
names = dir('positiveImage/*');
for i = 1:length(names)
    if names(i).name(1) ~= '.'
        filelist{count} = sprintf('positiveImage/%s',names(i).name);
        disp([num2str(count) ': added ' filelist{count}]);
        count = count + 1;
    end
end

% negative images
names = dir('negativeImage/*');
for i = 1:length(names)
    if names(i).name(1) ~= '.'
        filelist{count} = sprintf('negativeImage/%s',names(i).name);
        disp([num2str(count) ': added ' filelist{count}]);
        count = count + 1;
    end
end

names = dir('storedExponentialModel*.mat');
for i = 1:length(names)
    if names(i).name(1) ~= '.'
        filelist{count} = sprintf('%s',names(i).name);
        disp([num2str(count) ': added ' filelist{count}]);
        count = count + 1;
    end
end

names = dir('learning_result.mat');
for i = 1:length(names)
    if names(i).name(1) ~= '.'
        filelist{count} = sprintf('%s',names(i).name);
        disp([num2str(count) ': added ' filelist{count}]);
        count = count + 1;
    end
end

%% zip
zip(zipfilename,filelist,pwd);


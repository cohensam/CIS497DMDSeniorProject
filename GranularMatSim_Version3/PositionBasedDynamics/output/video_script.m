clear all
fileNames = dir(fullfile('*.png'));
N = length(fileNames);

tic

writerObj = VideoWriter('WaterSim-cow-1000-headon.mp4','MPEG-4');
writerObj.FrameRate = 30;
writerObj.Quality = 100; 

open(writerObj)

for k = 1:N
  filename = fileNames(k).name;
  images(:,:,:,k) = imread(filename);
end

writeVideo(writerObj,images)
close(writerObj)

toc
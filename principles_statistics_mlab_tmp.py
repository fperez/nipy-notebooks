# statstalk.m -*- matlab -*-
# Script for spm statistics tutorial
# If you want to run this script, I suggest you:
#  create a new directory somewhere, called (e.g.)  "statstalk";
#  copy this script to the new directory;
#  download the (11 megabyte) egscans.tar.gz archive into the same
#  directory.  The URL is:
#  http://imaging.mrc-cbu.cam.ac.uk/downloads/Tutscans/egscans.tar.gz
#  The archive contains 12 PET scans from one subject, which make up
#  an example dataset for the tutorial.
#  Then unpack the archive.  See the following webpage for instructions
#  on how to do this:
#  http://imaging.mrc-cbu.cam.ac.uk/imaging/UnpackCompressedFiles
#  Then start matlab, change directory to the statstalk directory
#  and type "statstalk" at the matlab prompt
# The script will run in matlab 4 or 5 with spm96 or higher, I believe
#
# Matthew Brett 21/8/99

# grand mean for data
GM = 50;


# transformation from voxel no to mm
M = spm_get_space(deblank(imgs(1, :)));
# mm to voxel no
iM = inv(M);
# coordinates of voxel of interest in mm (MNI space)
posmm = [-20 -42 34];
# coordinates in voxel space
posvox = iM * [posmm 1].T;
posvox = posvox(1:3);

# get data for this voxel
nimgs = size(imgs, 1);
vdata = zeros(nimgs, 1);
gdata = vdata;
disp('Getting image data and globals')
for i=1:nimgs
  iname = deblank(imgs(i, :));
  V = spm_vol(iname);
  vdata(i) = spm_sample_vol(V, posvox(1), posvox(2), posvox(3), 0);
  gdata(i) = spm_global(V);
end
disp('Done')

# demonstrate how global works on first image, without
# using the SPM global routine (which is in C)
# load data for first image using SPM utilities
disp(['Loading ' img1name ' for global calculation']);
V = spm_vol(img1name);
DIM = V.dim(1:3);
imgdata = spm_read_vols(V);
# Make one long vector, we're not going to use the array shape
imgdata = imgdata(:);
disp('Done')

# do global calculation
firstpassmean  = mean(imgdata);
tmp = find(imgdata > (firstpassmean/8));
globalimg1 = mean(imgdata(tmp));  # which = gdata(1)

# show how to find the data for the voxel we are interested in
# if we weren't using SPM routines
voxno = posvox(1)+(posvox(2)-1)*DIM(1)+(posvox(3)-1)*prod(DIM(1:2));
vdataimg1 = imgdata(voxno); # which = vdata(1)

# plot globals
figure
plot(gdata, vdata, 'x')
xlabel('TMVV for scan')
ylabel('Voxel value for -20 -42 34')
title('The relationship of overall signal to values for an individual voxel')

# proportionally scale data
pvdata = vdata ./ gdata

# scale to global mean of 50
Y = pvdata * GM

# our covariate
ourcov = np.array([5,4,4,2,3,1,6,3,1,6,5,2])

# plot data 
minx = 0 maxx = max(ourcov)+1 # min max x for plots
axs = [minx maxx min(Y)*0.98 max(Y)*1.02]
figure
plot(ourcov, Y, 'x')
axis(axs)
xlabel('Task difficulty')
ylabel('PS voxel value')
title('The relationship of task difficulty to voxel value')

# guess intercept, slope, estimated points
guessic   = 55
guesslope = 0.5
guessPts  = guessic + guesslope*ourcov

# plot data with guess intercept and slope
figure
plot(ourcov, Y, 'x')
axis(axs)
hold on 
plot([minx maxx], [minx maxx]*guesslope+guessic, '--') # guessed line
for i = 1:length(ourcov)     # guessed residuals
  plot(ourcov([i i]), [Y(i) guessPts(i)])
end
xlabel('Task difficulty')
ylabel('PS voxel value')
title('A first guess at a linear relationship, and its residuals')

# residuals from guess slope
guessRes = Y - guessPts

# design matrix, df
X  = [ourcov ones(nimgs,1)]
df = nimgs - rank(X)

# mean SoS for guess slope
guessRSS = sum(guessRes.^2) / df

# show design matrix
Xs = spm_DesMtx('sca', X)
figure
colormap('gray')
image((Xs+1)*32)
title('Design matrix for the first analysis')

# the analysis, giving slope and constant
B = pinv(X)*Y

# plot data with new least squares line
figure
plot(ourcov, Y, 'x')
hold on 
plot([minx maxx], [minx maxx]*B(1)+B(2), 'r')
axis(axs)
xlabel('Task difficulty')
ylabel('PS voxel value')
title('The least squares linear relationship')

# Contrast
C = [1 0]

# t statistic and significance test
RSS   = sum((Y - X*B).^2)
MRSS  = RSS / df
SE    = sqrt(MRSS*(C*pinv(X'*X)*C'))
t     = C*B./SE
ltp   = spm_Tcdf(t, df) # lower tail p
Z     = spm_invNcdf(ltp)
p     = 1-ltp           # upper tail p 

# print results to matlab window
fprintf('First TD analysis: t= #2.2f, p= #0.6f\n',...
	t, p)

# save data for analysis in e.g. SPSS
save voxdata.txt Y -ascii

# now analysis for added covariate
# note that this and all the other analyses here use
# exactly the same code as above

# design matrix, df
X = [ourcov (1:nimgs).T ones(nimgs,1)]
df = nimgs - rank(X)

# show design matrix
Xs = spm_DesMtx('sca', X)
figure
colormap('gray')
image((Xs+1)*32)
title('Design matrix for added covariate')

# the analysis, giving slopes and constant
B = pinv(X)*Y

# Contrast for TD, allowing for PR
C = [1 0 0]

# t statistic and significance test
RSS   = sum((Y - X*B).^2)
MRSS  = RSS / df
SE    = sqrt(MRSS*(C*pinv(X'*X)*C'))
t     = C*B./SE
p     = 1-spm_Tcdf(t, df) # upper tail p

# print results to matlab window
fprintf('Added covariate analysis: t= #2.2f, p= #0.6f\n',...
	t, p)

# Conditions analysis

# design matrix, df
X = [1 1 1 1 1 1 0 0 0 0 0 0...
     0 0 0 0 0 0 1 1 1 1 1 1]'
df = nimgs - rank(X)

# show design matrix
Xs = spm_DesMtx('sca', X)
figure
colormap('gray')
image((Xs+1)*32)
title('Design matrix for two conditions')

# the analysis, giving slopes (=means) 
B = pinv(X)*Y

# show that B contains the means
fprintf('\nBetas for condition analysis equal means by condition\n')
fprintf('Beta(R) #2.2f, mean(R) #2.2f\n'  , B(1), mean(Y(1:6)))
fprintf('Beta(A) #2.2f, mean(A) #2.2f\n', B(2), mean(Y(7:12)))

# Contrast for activation minus rest
C = [-1 1]

# t statistic and significance test
RSS   = sum((Y - X*B).^2)
MRSS  = RSS / df
SE    = sqrt(MRSS*(C*pinv(X'*X)*C'))
t     = C*B./SE
p     = 1-spm_Tcdf(t, df) # upper tail p

# print results to matlab window
fprintf('Analysis for two conditions: t= #2.2f, p= #0.6f\n',...
	t, p)

# Some extra analyses for good measure
fprintf('\nSome extra analyses to show how it all works...\n')

# do an ANCOVA for globals, instead of proportional scaling
# GM scale
Yanc = vdata * GM/mean(gdata)

# design matrix, df
X  = [ourcov gdata ones(nimgs,1)]
df = nimgs - rank(X)

# the analysis, giving slopes and constant
B = pinv(X)*Yanc

# Contrast for TD, allowing for ANCOVA for global
C = [1 0 0]

# t statistic and significance test
RSS   = sum((Yanc - X*B).^2)
MRSS  = RSS / df
SE    = sqrt(MRSS*(C*pinv(X'*X)*C'))
t     = C*B./SE
p     = 1-spm_Tcdf(t, df) # upper tail p

# print results to matlab window
fprintf('TD analysis with ANCOVA for global: t= #2.2f, p= #0.6f\n',...
	t, p)

# get data for another voxel

# coordinates of voxel of interest in mm (MNI space)
posmm  = [6 0 6]
# coordinates in voxel space
posvox = iM * [posmm 1]'
posvox = posvox(1:3)

# get data for this voxel
vdata2 = zeros(nimgs, 1)
for i=1:nimgs
  iname = deblank(imgs(i, :))
  V = spm_vol(iname)
  vdata2(i) = spm_sample_vol(V, posvox(1), posvox(2), posvox(3), 0)
end

# PS, GM scale
Y2 = vdata2./gdata * 50

# design matrix, df
X = [ourcov ones(nimgs,1)]
df = nimgs - rank(X)

# the analysis, giving slope and constant
B = pinv(X)*Y2

# Contrast
C = [1 0]

# t statistic and significance test
RSS   = sum((Y2 - X*B).^2)
MRSS  = RSS / df
SE    = sqrt(MRSS*(C*pinv(X'*X)*C'))
t     = C*B./SE
p     = 1-spm_Tcdf(t, df) # upper tail p

# print results to matlab window
fprintf('First TD analysis, second voxel: t= #2.2f, p= #0.6f\n',...
	t, p)

# two voxels at the same time, independently.
# Note that the code is just the same as for one voxel
Yboth = [Y Y2]
B     = pinv(X)*Yboth
RSS   = sum((Yboth - X*B).^2)
MRSS  = RSS / df
SE    = sqrt(MRSS*(C*pinv(X'*X)*C'))
t     = C*B./SE
p     = 1-spm_Tcdf(t, df) # upper tail p

# print results to matlab window
fprintf(['TD analysis for 1st and 2nd voxels '...
	'independently and simultaneously\n'])
fprintf('(which gives identical results to analysis one by one)\n')
fprintf('First  voxel: t= #2.2f, p= #0.6f\n',t(1), p(1))
fprintf('Second voxel: t= #2.2f, p= #0.6f\n',t(2), p(2))

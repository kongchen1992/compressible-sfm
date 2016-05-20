% Store directory of results folder
ROOT = '~/GitHub/Chen-CVPR-2016';
DATADIR = fullfile(ROOT, 'data');
RESULTDIR = fullfile(DATADIR, 'Results');

MOCAPDIR = fullfile(RESULTDIR, 'CMU-MoCap');
if ~exist(MOCAPDIR, 'dir')
    mkdir(MOCAPDIR);
end

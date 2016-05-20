function exp_dai_mocap(configI, configRank, Ks)

%[~, NumProc] = unix('nproc');
%NumCore = str2num(NumProc)/2;
NumCore = 3;
poolobj = parpool('local', NumCore);
parfor i = 1:numel(Ks)
    exp_dai_mocap_K(configI, configRank, Ks(i));
end
delete(gcp)

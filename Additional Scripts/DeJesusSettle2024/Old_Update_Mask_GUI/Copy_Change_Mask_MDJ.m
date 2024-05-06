%% Copy_Change_Mask_MDJ: Took apart Daan's old version of the code (pre-2021-07-08) to be compatible with new version.

warning('off','all');

for idx=1:length(MPStats)
    if idx==1
        MPStats_dummy = Update_Mask_GUI_MDJ(MPStats(idx));
        MPStats_dummy.FileName
    elseif idx>1
        MPStats_dummy(idx) = Update_Mask_GUI_MDJ(MPStats(idx));
        MPStats_dummy(idx).FileName
    end
end

clear MPStats;
MPStats = MPStats_dummy;
clear MPStats_dummy

disp('Done updating all masks for this MPStats file.');

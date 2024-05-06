%% Scatter plotting for Fig. 7-update Z-PCAs

load('IMG157_coeff-struct-mc.mat')
load('IMG157_coeff-struct-int.mat')

%% If run again, this is the celltype assignment

for idx=1:length(coeff_struct_int)
a = coeff_struct_int(idx).FileName;
if contains(a,'18k')
b = ' 18k';
else
b = ' 300';
end
coeff_struct_int(idx).Celltype = strcat(coeff_struct_int(idx).Celltype,b);
end
for idx=1:length(coeff_struct_mc)
a = coeff_struct_mc(idx).FileName;
if contains(a,'18k')
b = ' 18k';
else
b = ' 300';
end
coeff_struct_mc(idx).Celltype = strcat(coeff_struct_mc(idx).Celltype,b);
end
coeff_struct_int = sortStruct(coeff_struct_int,'Celltype');
coeff_struct_mc = sortStruct(coeff_struct_mc,'Celltype');

%% Plot Z-PCA for topography - run the calc line

celltype_list = {coeff_struct_mc.Celltype};
[u1,u2,~] = unique(celltype_list);

figure, hold on
% sgNT 18k
scatter(score([1:133],1),score([1:133],2),50,kulay(1,:),'x','LineWidth',2)
% sgNT 300
scatter(score([134:237],1),score([134:237],2),20,kulay(1,:),'filled','o')
% sgPTEN 18k
scatter(score([238:358],1),score([238:358],2),50,kulay(2,:),'x','LineWidth',2)
% sgPTEN 300
scatter(score([359:447],1),score([359:447],2),20,kulay(2,:),'filled','o')
pbaspect([1 1 1])

%% Plot Z-PCA for actin - run the calc line

celltype_list = {coeff_struct_int.Celltype};
[u1,u2,~] = unique(celltype_list);

figure, hold on
% sgNT 18k
scatter(score([1:133],1),score([1:133],2),50,kulay(1,:),'x','LineWidth',1.5)
% sgNT 300
scatter(score([134:237],1),score([134:237],2),20,kulay(1,:),'filled','o')
% sgPTEN 18k
scatter(score([238:358],1),score([238:358],2),50,kulay(2,:),'x','LineWidth',1.5)
% sgPTEN 300
scatter(score([359:447],1),score([359:447],2),20,kulay(2,:),'filled','o')
pbaspect([1 1 1])

%% gates_ZPCA overlay

for idx=1:length(gates_ZPCA)
 patch(gates_ZPCA(idx).Polygon_coords(:,1),gates_ZPCA(idx).Polygon_coords(:,2),[1 1 1],'FaceAlpha',0,'LineStyle','--');
end
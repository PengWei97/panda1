% 用于生成GBAnisotropy中的晶界能、晶界迁移率和晶界激活能
% 生成txt文件
%------------------参数处理-------------------------------

count = input('please input the number of grain:')
V=zeros(3*count,count);

if i = 1:3*count
    if 0<i<=count
        for j = 1:count
            V(i,j) = 0.708;
        end
    elseif count< i <= 2*count
        for j = 1:count
            V(i,j) = 0.708;
        end
    else
        for j = 1:count
            V(i,j) = 0.23;
        end


%-----------------------------输出---------------------------------------
% filename = strcat('grn_',num2str(count),'_rand_2D.tex')
filename = 'poly_anisotropy_mobility.txt'
f=fullfile('D:\Matlab2021A\bin\inputfile',filename);
[r,c]=size(V);
fid=fopen(filename,'w');
fprintf(fid,'>>>>>>>>> GB properties <<<<<<<<<\n');
% fprintf(fid,'\n');
fprintf(fid,'---GB energy (sigma_ij), mobility prefactor (mob0_ij), activation energy (Q_ij)---\n');
% fprintf(fid,'B ');
% fprintf(fid,'%d\n',r);

for i=1:r
    fprintf(fid,'   ');
    for j=1:c
        if j==c
            fprintf(fid,'%3.2f\n',V(i,j));%如果是最后一个，就换行
        else
            if V(i,j)>0 && V(i,j)<10
               fprintf(fid,' '); 
            end
            fprintf(fid,'%4.2f',V(i,j));%如果不是最后一个，就tab
            fprintf(fid,'   ');
        end
    end

end
fclose(fid);

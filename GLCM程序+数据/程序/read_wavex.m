%%  V2.01˵������Ӣ  2011-12-12���棩
%
%%  V2.0  ˵�������ͣ���Ӣ 2010-11-20���棩
%% ������ʵ�ֽ��ɼ������״����ݽ����˲���ֵ���õ�ѡ���ڵ�������Ϣ��
%   ���иó���֮ǰ��Ҫ������fun_rad_para.m�����ø�����ʼ�趨����
%   month_num  day_num ��Ҫ����ʵ�ʼ������ݵ����ڽ����޸�
%   path_data1Ϊ��Ҫ��������ݵĲ��ִ洢·��������ɸ�ԭʼ���ݾ���洢λ�����޸ģ�path_data2��ԭʼ�������ڵ��ļ���
%   �Ӻ���fun_rad_01_file_name_pt.m��path_name����ݵĵط�Ҫ���ݾ��������һ����������޸ģ�2009��2010��2011��
%   09���10���ԭʼ���ݴ洢��ʽ�иı䣬�Ӻ���fun_rad_02_read_data_pt.m
%   �еġ�data_all=bitand(data_all,16383);�������ڲ�ֵ10��11��Ĳ������ݵ�ʱ���õ�����ֵ09������ݵ�ʱ��Ͳ�����
%   fun_data_marine_get ��������ȡ������Ϣ��
%   ����ֵ��ɣ����õ�һ��������Ϣ�ļ�
%%
clear all
clc
close all
load para.mat

r0=600;               % the distence of the data chosen range.
phi=(75/180)*pi;      % the angle of the data chosen range.  2011-12-03

month_num=10;          % the number of the month.
day_num=22;           % the number of the day.

hour_beg=10;          % to judge wether the data exists, the hour begins.
min_beg=35;            % to judge wether the data exists, the minute begins.
min_int=1;            % the interval of the minutes.
time_i=0;             % the i in the for loop.
%%
month_name = fun_rad_011_name_gen(month_num,12);  % generate the month name.
day_name   = fun_rad_011_name_gen(day_num,31);    % generate the day name.
if month_num<10
    filename_write1=['0',num2str(month_num)]; %because the hour beg ins at 0. the mid night. so the main pro should consider this condition.
else
    filename_write1=[num2str(month_num)];
end
if day_num<10
    filename_write2=['0',num2str(day_num)]; %because the hour begins at 0. the mid night. so the main pro should consider this condition.
else
    filename_write2=[num2str(day_num)];
end

filename_write=['data_marine_',filename_write1,filename_write2,'.dat'];  %%%%  ������Ϣ
fd1_write=fopen(filename_write,'w');
fidt=fopen('GLCM22.txt','wt');

%%
[filename_ogr, path_name] = uigetfile('*.dat', 'X �����״������ļ�');

% path_data1='D:\2011����ʵ��\2011����ʵ������\';                       % the path of the data dir
% path_data2=['data_',month_name,'_',day_name,'\']; % the path of the certain date data.
% path_name=[path_data1,path_data2];                   % generate the whole data dir.

%%  begin the loop, to test the file exists.
for hour_i=1:24
    hour_num=hour_beg+(hour_i-1);
    for min_i=1:60                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        min_num=min_beg+min_int*(min_i-1);
        [file_name,file_flg]=fun_rad_01_file_name_pt(path_name,month_num,day_num,hour_num,min_num);%judge the file, if true, the file_flg=1,else file_flg=0
        if file_flg==1
            time_name =[path_name(end-5:end-4),path_name(end-2:end-1),',',...
                file_name(end-11:end-10),file_name(end-8:end-7),',']; % if the file exists, draw the file_name 4 the img output.
           [dirEffect,velReal_ave,vel2flow,shipDir]=fun_data_marine_get(file_name,time_name,fd1_write);
            if abs(dirEffect)<15
                tic;
                [out_test,filename,akx1,aky1,direction_end]=fun_wind_streaks_fft2(file_name);% the main function, to calculate the wave parameters.
               toc;
               t=toc;
                fprintf(fidt,'%s %.02f %.1f',file_name(end-22:end-7),direction_end,t);
                fprintf(fidt,'\r\n');
            end
        end
    end
end
fclose(fd1_write);
fclose(fidt);


%   [out_test,filename,akx1,aky1]=fun_data_read_915_32(O:\2010.10½����������\hrb\10_10_22_10\2010_10_22_10_35_32.dat)
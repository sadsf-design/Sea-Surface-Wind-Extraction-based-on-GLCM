function [out_test,filename,akx1,aky1,direction_end]=fun_wind_streaks_fft2(file_name)
ship_direction=93;%船艏向
R_num=600;
dr=7.5;
pic_num=1;
% file_name=['E:\实验仿真\test_MATLAB\扇形——测试\2010_10_22_10_35\2010_10_22_10_35_10.dat'];
r_beg=63;%半径起始位置
r_end=300;%半径终止位置
theta_beg_temp=13*100;
theta_end_temp=198*100;
theta_beg=13*100+ship_direction*100;%扇形起始角度
theta_end=198*100+ship_direction*100;%扇形终止角度
out_test_temp=[0];
for num=1:32
     imag_temp=[];
     if num<10
         filename=[file_name(1:end-6),'0',num2str(num),'.dat'];
     else
         filename=[file_name(1:end-6),num2str(num),'.dat'];
     end
    fid_0=fopen(filename,'r');
    status = fseek(fid_0,3114,'bof');
    temp_data=fread(fid_0,[1,3601],'ubit16');
    n_line=temp_data(1);%角度向总线数
    theta=temp_data(2:n_line+1)*0.01*pi/180;%全部角度的弧度值   
    theta0=temp_data(2:n_line+1)*0.01;%全部角度值
    thata1=temp_data(2:n_line+1);%全部角度值的100倍
    data_all=fread(fid_0,'ubit16');
    data_all=bitand(data_all,16383);%所有点的强度值
    
%     a = -0.2;
%     b = 0.03;
%     n_rayleigh = a + (-b .* log(1 - rand(M,N))).^0.5;
    
    fclose(fid_0);
    
%     data_all = data_all + data_all.*n_rayleigh;
    
    Gray=reshape(data_all,[R_num,n_line]);%原始图像全部半径和角度生成的灰度值矩阵
    Gray= medfilt2(Gray,[3,3]);  %中值滤波

%% ---------------------------生成1800个新角度并选取大小叠加平均---------------------------
   
   n_line_new=1800;
   thata1_temp=0;
   thata1_new=[];
   for i=1:1800
       thata1_new(i)=thata1_temp;%1800个新角度存放在thata1_new变量中
       thata1_temp=thata1_temp+20;    
   end
 %%插值，在3608个角度中寻找与1800个角度最相近的角度，从而形成新的600*1800的灰度矩阵
a=thata1;%全部角度值的100倍，原角度
b=thata1_new;%生成的新角度
aa=length(a);
bb=length(b);
Gray_new=zeros(600,1800);%生成0.2°为间隔的灰度矩阵
Gray_north=zeros(600,1800);%将新生成的矩阵转为正北向
c=1;
s=1;
for j=1:bb
    for i=c:aa      
        if a(i)>=b(j)
            Gray_new(:,s)=Gray(:,i);    % 赋值回波强度
%             Gray_new(:,s)=Gray_new(:,s)+Gray(:,i);    %%生成新矩阵
            s=s+1;
%             c=i+1; % 后加的部分
        break
%          c=i;
%         else
%             Gray_new(:,s)=Gray_new(:,s)+Gray(:,i);
        end    
    end
end
%%%%将新生成的图像转为正北向%%%%%%
gray_line=0;
line_north=ship_direction/0.2;%船艏向相对北的线数
for i=1:bb
    gray_line=i+line_north;
    if gray_line>1800
        gray_line=gray_line-1800;
    end
    Gray_north(:,gray_line)=Gray_new(:,i);   %转为正北向
end
dx=7.5;
dy=7.5;
dr=7.5;
%%在新的灰度值矩阵中截取扇形
max_thata1_new=max(thata1_new);
line_beg=ceil(theta_beg*n_line_new/max_thata1_new);% 新的角度起始
line_end=ceil(theta_end*n_line_new/max_thata1_new);% 新的角度终止   
angle=thata1_new(line_beg:line_end);
angle=angle*0.01*pi/180;%新的角度范围
[akx1,aky1]=fun_rad_022_Rtheta2XY(angle,r_beg,r_end,dr); %生成截取后的x,y坐标
image_north=Gray_north(r_beg:r_end,line_beg:line_end);
imag_temp=Gray_new(r_beg:r_end,line_beg:line_end);
image0=zeros(r_end-r_beg+1,line_end-line_beg+1,32);

% m=find(imag_temp>4500); %%%  hrb数据的时候改成 3500
% imag_temp(m)=4500;  %%%  hrb数据的时候改成 3500
% 
% n=find(imag_temp<500); %%%  hrb数据的时候改成 500
% imag_temp(n)=500;  %%%  hrb数据的时候改成 500
          
%--------------拟合去除距离干扰——————————————
% row_ave(:,num)=sum(image_north,2)./size(image_north,2);

image_north_32(:,:,num)=image_north;
% for k=1:size(image_north,2)
%     image_north_modify(:,k)=image_north(:,k)-row_ave(:,num);
% end
%----------------------------------绘制单幅图像----------------------------------

%   pcolor(akx1,aky1,image_north_modify(:,:,1))%转为正北方向
%                  hold on 
%                  axis([-2500 2500 -2500 2500])
%                pp=-2500:2500;
%                qq=-2500:2500;
%                pp_temp=zeros(1,length(pp));
%                qq_temp=zeros(1,length(qq));
%                plot(pp,qq_temp,'-.k')
%                hold on
%                plot(pp_temp,qq,'-.k')
% 
% 
%                  
%                   t=0:pi/20:2*pi;%加上从500开始到4500以1000为间隔的同心圆
%                   for i=600:300:2400
%                      xx=i.*cos(t);
%                      yy=i.*sin(t);
%                      axis([-2500 2500 -2500 2500]);
%                      plot(xx,yy,'-.k');
%                   end
%                   
%                  
%                  x=-2500:2500;%加上夹角辅助线
%                  y=-2500:2500;
%                  y=sqrt(3)/3*x;
%                  plot(x,y,'-.k');
%                  hold on 
%                  y=sqrt(3)*x;
%                  plot(x,y,'-.k');
%                  y=-sqrt(3)/3*x;
%                  plot(x,y,'-.k');
%                  y=-sqrt(3)*x;
%                  plot(x,y,'-.k');
%                  hold off
%                  
%                  axis([-2500 2500 -2500 2500])
%                  colorbar
%                  shading interp
%                  caxis auto
%                  axis equal
%                  axis([-2500 2500 -2500 2500])
%                  
%                  s1=filename(end-5:end-4);
%                  s2=filename(end-8:end-7);
%                  s3=filename(end-11:end-10);
%                  s4=filename(end-14:end-13);
%                  s5=filename(end-17:end-16);
%                  s6=filename(end-22:end-19);
%                  savename=strcat(s6,'.',s5,'.',s4,',',s3,':',s2,',Correct radar image');
%                  xlabel('distance[m]','FontName','Time New Roman');
%                  ylabel('distance[m]','FontName','Time New Roman');
%                  title(savename); 
%           
%                  x_p=[0.2 0.2];
%                  y_p=[0.7 0.9];
%                  annotation('arrow',x_p,y_p);
%                  x_pp=[0.18 0.22];
%                  y_pp=[0.8 0.8];
%                  annotation('line',x_pp,y_pp);
%                  text(-1950,2100,'N','FontName','Time New Roman');
%                  text(-1950,2100,'N','fontweight','bold');
%                  real_wind=(90-36)*pi/180;   %实际风向
%                  x_wind=[0.6+0.1*cos(real_wind) 0.6];
%                  y_wind=[0.6+0.1*sin(real_wind) 0.6];
%                  annotation('arrow',x_wind,y_wind);
%                  text(1100, 700,'参考风向','fontweight
% out_test_temp=out_test_temp+image_north_modify;%图像序列的叠加
end
% save bijiao.mat
% clear all
% clc
% load bijiao.mat
% %%---------------------3次多项式最小二乘拟合—————————————————-
% X=r_beg:r_end;
% for i=1:32
%    p=polyfit(X,row_ave(:,i)',3);
%    pp(i,:)=p;
% end
% ppp=sum(pp,1)./32;
% line_reduce=polyval(ppp,X);
% % hold on 
% % plot(X,row_ave(:,1),'.');
% % plot(X,line_reduce,'linewidth',2);
% % hold off
% % error=line_reduce'-row_ave(:,1);
% % corr=corrcoef(line_reduce,row_ave(:,1));%相关系数
% % bias=mean(error);%平均误差
% % standard=std(error);%标准偏差
% % y_ave=mean(row_ave(:,1));
% % r=sqrt(sum(line_reduce'-y_ave).^2/sum(row_ave(:,1)-y_ave).^2);%显著性检验
% % F=((238-2)*r^2)/1-r^2;
% 
% %%%%%%------拉格朗日多段拟合校正----------------
% for i=1:32
% X=r_beg:r_end;
% for h=1:size(X,2)
%     x(h,:)=[1 X(h) X(h)^2];%假设两段函数基函数都为1,x,x^2
% end
% x1=x(1:size(X,2)/2,:);
% x2=x(size(X,2)/2:size(X,2),:);    
% y=row_ave(:,i);
% y_size=size(y);
% y1=y(1:y_size/2);
% y2=y(y_size/2:y_size);
% z1=x(size(X,2)/2,:)';
% z2=x(size(X,2)/2,:)';
% a1_temp=inv(x1'*x1)*x1'*y1;
% a2_temp=inv(x2'*x2)*x2'*y2;
% lambda=(z1'*a1_temp-z2'*a2_temp)/(z1'*inv(x1'*x1)*z1+z2'*inv(x2'*x2)*z2);
% a1=a1_temp-lambda.*inv(x1'*x1)*z1;
% a2=a2_temp-lambda.*inv(x2'*x2)*z2;
% y1_count=x1*a1;
% y2_count=x2*a2;
% y_count(:,i)=[y1_count(1:end-1);y2_count];
% % error=y_count-y;
% % corr=corrcoef(y_count,y);%相关系数
% % bias=mean(error);%平均误差
% % standard=std(error);%标准偏差
% % hold on 
% % plot(X,y,'.');
% % plot(X,y_count,'linewidth',2);
% % hold off
% error_1=abs(y1_count-y1);%%计算绝对值误差及平均值
% error_2=abs(y2_count-y2);
% s_error=[size(error_1,1),size(error_2,1)];
% bias_1=mean(error_1);
% bias_2=mean(error_2);
% bias=[bias_1,bias_2];
% %%%%%———————————依据绝对值误差判断是否继续分段————————————
% M=max(s_error);
% error_1_temp=[error_1;zeros(M-s_error(1))];
% error_2_temp=[error_2;zeros(M-s_error(2))];
% error_temp=[error_1_temp,error_2_temp];
% aaa=1;
% K=1;%断点数
% y_count=[];
% while(aaa==1)
%     [error_temp,bias,aaa,y_count,s_error,K]=fun_revise(error_temp,bias,x,y,s_error,X,K,y_count);
% end
% y_count_end(:,i)=y_count;
% % figure
% %  plot(X,y_count_end,'r','LineWidth',2);
% %  hold on 
% %  plot(X,y,'.')
% %  xlabel('径向像元位置','fontweight');
% %  ylabel('NRCS','FontName','Time New Roman');
% %  legend('径向拟合 NRCS','雷达图像径向 NRCS');
% %  hold off
% end
% 
% %%%---------------------拟合结束-------------------
% for i=1:32
%     for j=1: size(image_north,2)
%     image_north_modify(:,j,i)=image_north_32(:,j,i)-y_count_end(:,i);%%-line_reduce';%%距离校正后的图像
%     end
% end
out_test_temp=sum(image_north_32,3);

    M = size(out_test_temp, 1);
    N = size(out_test_temp, 2);

    a = 0.01;    %高斯噪声均值
    b = 0.1; %高斯噪声标准差
    n_gaussian = a + b .* randn(M,N); 
    
    out_test_temp = out_test_temp + out_test_temp.*n_gaussian;
    
%     out_test_temp = medfilt2(out_test_temp,[3,3]);

%%------------------------------------------绘制校正后的图像--------------------------------
out_test=out_test_temp./32;
% out_test(2:5, 10:900) = 0;
% out_test(5:235, 900:903) = 0;
% out_test(235:238, 10:903) = 0;


 figure,pcolor(akx1,aky1,out_test)%转为正北方向
 colormap jet
                 hold on 
                 axis([-2500 2500 -2500 2500])
               pp=-2500:2500;
               qq=-2500:2500;
               pp_temp=zeros(1,length(pp));
               qq_temp=zeros(1,length(qq));
               plot(pp,qq_temp,'-.k')
               hold on
               plot(pp_temp,qq,'-.k')


                 
                  t=0:pi/20:2*pi;%加上从500开始到4500以1000为间隔的同心圆
                  for i=600:300:2400
                     xx=i.*cos(t);
                     yy=i.*sin(t);
                     axis([-2500 2500 -2500 2500]);
                     plot(xx,yy,'-.k');
                  end
                  
                 
                 x=-2500:2500;%加上夹角辅助线
                 y=-2500:2500;
                 y=sqrt(3)/3*x;
                 plot(x,y,'-.k');
                 hold on 
                 y=sqrt(3)*x;
                 plot(x,y,'-.k');
                 y=-sqrt(3)/3*x;
                 plot(x,y,'-.k');
                 y=-sqrt(3)*x;
                 plot(x,y,'-.k');
                 hold off
                 
                 axis([-2500 2500 -2500 2500])
                 colorbar
                 caxis([0 4000])
                 shading interp
%                  caxis auto
                 axis equal
                 axis([-2500 2500 -2500 2500])
      
                 xlabel('distance[m]','FontName','Time New Roman');
                 ylabel('distance[m]','FontName','Time New Roman');
                 
                 titl_name_temp = filename(end-22:end-6);
                 titl_name = [];
                 n_count = 1;
                 in22 = {'.'; '.'; '  ';  ':'; '  '};
                 for i = 1 : length(titl_name_temp)
                     if titl_name_temp(i) ~= '_'
                         titl_name = [titl_name, titl_name_temp(i)];
                     else
                         titl_name = [titl_name, in22{n_count}];
                         n_count = n_count + 1;
                     end
                 end
                 
                 title([titl_name, ' Overlay Wind Map'],'FontName','Time New Roman'); 
      
                 x_p=[0.2 0.2];
                 y_p=[0.7 0.9];
                 annotation('arrow',x_p,y_p);
                 x_pp=[0.18 0.22];
                 y_pp=[0.8 0.8];
                 annotation('line',x_pp,y_pp);
                 text(-1950,2100,'N','FontName','Time New Roman');
               real_wind=(90-39)*pi/180;   %实际风向
                  x_wind=[0.18+0.1*cos(real_wind) 0.15];
                 y_wind=[0.18+0.1*sin(real_wind) 0.15];
                 annotation('arrow',x_wind,y_wind);
                  text(-2200, -1900,'Wind','FontName','Time New Roman');
                axis([-2500 2500 -2500 2500]);
                 hold off
                 filename=[filename(end-17:end-4),'.jpg'];%保存图像
                 print( '-djpeg',filename);
                 close

%----------------------------------GLCM李诗语----------------------------------
%   [direction_end]=main3(filename);
[direction_end]=main3(filename);
  direction_end=0;
end

     



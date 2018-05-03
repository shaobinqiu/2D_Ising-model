%clear  
n=3;
eps=0.01;

%%%%%%%%%%%%Merge data and dot Minimum common multiple(1-25)
% apha=267711444;
%  A_20=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infs_20.txt');
% A_21=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infs_21.txt');
% A_22=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infs_22.txt');
% A_23=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infs_23.txt');
% A_24=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infs_24.txt');
% A_25=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infs_25.txt');
% A_1_25=[A_20;A_21;A_22;A_23;A_24;A_25];
% A=apha*A_1_25;
% A=round(100*A);
% save inf1_25.txt A -ascii -double
% step=1
%%%%%%%%%%%%%
% filename = '/home/qiusb/Documents/python_work/2D_Ising_model/1_25/all_inf1_25.txt';
% A_d = importdata(filename);
 %A=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/inf1_25.txt');
sits=[7 1 0 0; 
       7 1 2 0; 
       7 1 2 3;
       1 2 0 0;
       1 2 3 0];
   D=[];%situtations of considering sum 1st; sum 1st 2nd ;sum 1st 2nd 3rd; 1st 2nd ; 1st 2nd 3rd 
for c=5%1:size(sits,1)
%     new_folder1=['./',num2str(sits(c,1)*1000+sits(c,2)*100+sits(c,3)*10+sits(c,4)*1)];
%     mkdir(new_folder1); 
    inf=[];
    for ii=1:size(sits,2)
        if sits(c,ii)~=0
            inf=[inf A(:,sits(c,ii))] ;
        end
    end                                                                 %%%%%%%%%%set a information matrix
    if sits(c,1)==7
        inf=[inf;-inf(:,1) inf(:,2:size(inf,2))];
    end
%     plot(inf(:,1),inf(:,2),'*')
%     hold on 
    K = convhulln(inf,{'Qt','Pp'});
    dir=unique(K);
    
    
    x=[];y=[];
    for xx=1:size(dir,1)
       x=[x; inf(dir(xx,1),1)];
       y=[y;inf(dir(xx,1),2)];
    end
     k = convhull(x,y);
    plot(x(k)/26771144400,y(k)/26771144400,'r-')
    hold on 
    plot(x/26771144400,y/26771144400,'*')
xlabel('N\_1nd');
ylabel('N\_2nd');
    
    D=[D;size(dir)]%the number of convex p
%     for jj=1:size(dir)
%         d=inf(jj,:);
%         n=['C'];
%         for rr=1:size(d,2)
%             n=[n,'_',num2str(d(1,rr))];
%         end                                                             %convex point data(for folder name)
%         new_folder2=[new_folder1,'/',n];
%         mkdir(new_folder2); 
%         for kk=1:size(inf,1)                                     %%k :the dir of  structures which accord with convex point
%             diff=inf(kk,:)-d;
%             if sum(diff.^2)<eps
%                 file_dir= regexp(A_d{mod(kk,size(A,1)),1}, '5.54    ', 'split')
%                 copyfile(file_dir,new_folder2)
%             end
%         end
%     end
end
        

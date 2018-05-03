

eps=0.01;



%  filename = '/home/qiusb/Documents/python_work/2D_Ising_model/1_25/all_infev1_25.txt';
%  A_d = importdata(filename);
%A=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infev1_25.txt');


  
%B=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/infev1_25_unique.txt');
 %E=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/BC1_25.txt');

 %save infev1_25_unique.txt B -ascii -double
% F=load('/home/qiusb/Documents/python_work/2D_Ising_model/1_25/B_C1_25.txt');


 D=[A(:,[7,9,10,11,12])];


 tic
     K12 = convhulln(D(:,:),{'Qt','Qx','QbB'});                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
     dir12=unique(K12);
 toc
    
 W=[];
   
   for ii=1:size(A,1)
       for jj=1:size(dir12,1)
           if sum((A(ii,[7,9,10,11,12])-D(dir12(jj),1:5)).^2)<0.5 && sum(A(ii,[1,2,3,4,5,6,8]))==0
              W=[W;ii ];
           end
       end
   end
    
 dir=W;   
 R=E(dir,:);
%     for jj=1:size(dir)
%             
%         new_folder=['conv_stru_env12','/',n];
%         mkdir(new_folder); 
%                 file_dir= regexp(A_d{mod(dir(jj),size(A,1)),1}, '5.54    ', 'split')
%                 copyfile(file_dir{1,2},new_folder2)
%     end
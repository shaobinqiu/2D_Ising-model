clear

 filename = './1_25/all_inf1_25.txt';
 A_d = importdata(filename);
%A=load('./1_25/infev1_25.txt');
 1

dir=load('dir.txt dir' );    
        new_folder='conv_stru_env12';
        mkdir(new_folder);
    for jj=1:size(dir)
                jj 
                file_dir1= regexp(A_d{mod(dir(jj),size(A,1)),1}, '5.54    ', 'split');
                file_dir= file_dir1{1,2};
                copyfile(file_dir,new_folder)
                
    end

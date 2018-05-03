clear
A=load('d2C_B16_s_ij_s_i.txt');
inf=A(:,1:7);
X=[inf(:,1:2) inf(:,4) inf(:,5) inf(:,6) inf(:,7)]

plot(X(:,1),X(:,2),'*')
xlabel('\rho');
K = convhulln(inf)
patch('Vertice',X,'Faces',K,'FaceColor','y')
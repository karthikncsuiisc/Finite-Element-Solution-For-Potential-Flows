clear all
clc
fileID = fopen('1feflo.domn.cylinder.coarse_phi_u_v_p_h.dat','r');
formatSpec = '%d %lf %lf %lf %lf %lf';
sizeA = [6 Inf];

A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';
havg_mesh1=A(5,6);
unum=(A(:,3).^2+A(:,4).^2).^0.5;
fileID = fopen('1feflo.domn.cylinder.coarse_results.vtk','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
A=A';

coord=A(:,1:2);

[m n]=size(coord);
U=1;
a=0.5;

for i=1:1:m
    
    theta=tan(coord(i,2)/coord(i,1));
    r2=coord(i,1)^2+coord(i,2)^2;
    
    uax=U*(1-cos(2*theta)*a^2/r2);
    uay=-a^2/r2*U*sin(2*theta);
    uana(i)=sqrt(uax^2+uay^2);
end

E2=sum((unum'-uana).^2);

E2_mesh1=sqrt(E2/m)
havg_mesh1
%---------------------------------------------------------------------------------------------

% %mesh 2
fileID = fopen('2feflo.domn.cylinder.medium_phi_u_v_p_h.dat','r');
formatSpec = '%d %lf %lf %lf %lf %lf';
sizeA = [6 Inf];

A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';
havg_mesh2=A(5,6);
unum=(A(:,3).^2+A(:,4).^2).^0.5;

fileID = fopen('2feflo.domn.cylinder.medium_results.vtk','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
A=A';

coord=A(:,1:2);

[m n]=size(coord);
U=1;
a=0.5;

for i=1:1:m
    theta=tan(coord(i,2)/coord(i,1));
    r2=coord(i,1)^2+coord(i,2)^2;
    
    uax=U*(1-cos(2*theta)*a^2/r2);
    uay=-a^2/r2*U*sin(2*theta);
    uana(i)=sqrt(uax^2+uay^2);
end

E2=sum((unum'-uana).^2);

E2_mesh2=sqrt(E2/m)
havg_mesh2
%---------------------------------------------------------------------------------------------

%mesh 3
fileID = fopen('3feflo.domn.cylinder.fine_phi_u_v_p_h.dat','r');
formatSpec = '%d %lf %lf %lf %lf %lf';
sizeA = [6 Inf];

A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';
unum=(A(:,3).^2+A(:,4).^2).^0.5;
havg_mesh3=A(5,6);

fileID = fopen('3feflo.domn.cylinder.fine_results.vtk','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
A=A';

coord=A(:,1:2);

[m n]=size(coord);
U=1;
a=0.5;

for i=1:1:m
    theta=tan(coord(i,2)/coord(i,1));
    r2=coord(i,1)^2+coord(i,2)^2;
    
    uax=U*(1-cos(2*theta)*a^2/r2);
    uay=-a^2/r2*U*sin(2*theta);
    uana(i)=sqrt(uax^2+uay^2);
end

E2=sum((unum'-uana).^2);

E2_mesh3=sqrt(E2/m)
havg_mesh3
%--------------------------------------------------------------------------------
% %mesh 4
fileID = fopen('4feflo.domn.sphere.vfine_phi_u_v_p_h.dat','r');
formatSpec = '%d %lf %lf %lf %lf %lf';
sizeA = [6 Inf];

A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';
unum=(A(:,3).^2+A(:,4).^2).^0.5;
havg_mesh4=A(5,6)

fileID = fopen('4feflo.domn.sphere.vfine_results.vtk','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
A=A';

coord=A(:,1:2);

[m n]=size(coord);
U=1;
a=0.5;

for i=1:1:m
    theta=tan(coord(i,2)/coord(i,1));
    r2=coord(i,1)^2+coord(i,2)^2;
    
    uax=U*(1-cos(2*theta)*a^2/r2);
    uay=-a^2/r2*U*sin(2*theta);
    uana(i)=sqrt(uax^2+uay^2);
end

E2=sum((unum'-uana).^2);

E2_mesh4=sqrt(E2/m)
havg_mesh4

x= [havg_mesh4 havg_mesh3 havg_mesh2];
y= [E2_mesh4 E2_mesh3 E2_mesh2];



plot(x,y,'Linewidth',2)
hold on
plot(x,y,'o','Linewidth',2)
xlabel('Average mesh cell size');
ylabel('L2 norm of (u_{nummerical}-u_{accurate})');
title('L2 norm vs Average cell size');

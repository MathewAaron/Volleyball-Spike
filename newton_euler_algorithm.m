% Aaron Mathew 
% Code for Recursive Newton-Euler algorithm

clear all;
close all;
clc;

%% Declaring variables
% The Data provided is not made public
% The data-table provided is in format
% | Time(ms) , |Theta_shoulder| , |Theta dot_shoulder|, |Theta-double-dot_shoulder|, |Theta_elbow| , |Theta dot_elbow|, |Theta-double-dot_elbow|,\
% |Theta-double-dot_shoulder|, |Theta_wrist| , |Theta dot_wrist|, |Theta-double-dot_wrist|| 
data = xlsread('arm_spike.xlsx');
time = data(:,1);

theta_shoulder = data(:,2);
theta_dot_shoulder = data(:,3);
theta_d_dot_shoulder = data(:,4);

theta_elbow = data(:,5);
theta_dot_elbow = data(:,6);
theta_d_dot_elbow = data(:,7);

theta_wrist = data(:,8);
theta_dot_wrist = data(:,9);
theta_d_dot_wrist = data(:,10);

% Code for Recursive Newton Euler algoritm

%% Declaring variables
Lu = 0.47;
Ll = 0.45;
Lw = 0.15;
phi = [1;1;1];


%m = [2+Zm1,1+Zm2,1+Zm3] % Calculate COM

m = [4.964,2.369,1.2006];
g = 9.8; % gravity

s1 = [0.1239 ; 0 ; 0];
s2 = [0.1140 ; 0 ; 0];
s3 = [0.0106 ; 0 ; 0];

s = cat(2,s1,s2,s3);
syms t1 t2 t3

w = zeros(3,4);
w_dot = zeros(3,4);

v = zeros(3,4);
v(:,1) = [0,g,0];

v_c = zeros(3,3);
F_shoulder = zeros(3,length(theta_shoulder));
F_elbow = zeros(3,length(theta_shoulder));
F_wrist = zeros(3,length(theta_shoulder));

I1 = [0.0945 0 0; 0 0.0945 0; 0 0 0.0062];
I2 = [0.0407 0 0; 0 0.0407 0; 0 0 0.0015];
I3 = [0.0024 0 0; 0 0.0024 0; 0 0 0.0003];

I = cat(3,I1,I2,I3);
%% Getting Transform Matrices
% reference : https://www.mathworks.com/matlabcentral/fileexchange/103050-dh-table-solver
dh_shoulder = [t1 0 0 0];
dh_elbow = [t2 0 Lu 0];
dh_wrist = [t3 0 Ll 0];

% Computing 0T1, 1T2, 2T3 matrix
T01 = [cos(t1) -sin(t1) 0 0; sin(t1) cos(t1) 0 0; 0 0 1 0; 0 0 0 1];
T12 = [cos(t2) -sin(t2) 0 Lu; sin(t2) cos(t2) 0 0; 0 0 1 0; 0 0 0 1];
T23 = [cos(t3) -sin(t3) 0 Ll; sin(t3) cos(t3) 0 0; 0 0 1 0; 0 0 0 1];
T34 = [1 0 0 Lw;0 1 0 0; 0 0 1 0; 0 0 0 1];

% Computing 0R1, 1R2, 2R3 matrix
R01 = T01(1:3,1:3);
R12 = T12(1:3,1:3);
R23 = T23(1:3,1:3);
R34 = T34(1:3,1:3);
rot_matrix = cat(3,R01,R12,R23,R34);

% Computing inverse of rotation matrix
R10 = inv(R01);
R21 = inv(R12);
R32 = inv(R23);
R43 = inv(R34);
rot_matrix_inverse = cat(3,R10,R21,R32,R43);

%Computing position matrix
P01 = T01(1:3,4);
P12 = T12(1:3,4);
P23 = T23(1:3,4);
P34 = T34(1:3,4);
%position_matrix = cat(3,P01,P12,P23,P34);
pos = cat(3,P01,P12,P23,P34);

% outward

rot_matrix = cat(3,R01,R12,R23,R34);

n_var = zeros(3,4);
n_shoulder = zeros(3,length(time));
n_elbow = zeros(3,length(time));
n_wrist = zeros(3,length(time));

f = zeros(3,4);
n_final = zeros(1,length(time));
x = 1;
for j = 1:length(time)

    t1 = theta_shoulder(j);
    t2 = theta_elbow(j);
    t3 = theta_wrist(j);
    theta_dot = [theta_dot_shoulder(j);theta_dot_elbow(j);theta_dot_wrist(j)];
    theta_double_dot = [theta_d_dot_shoulder(j);theta_d_dot_elbow(j);theta_d_dot_wrist(j)];

    % Forward Propagation     
    for i=1:3 
       v_var = vpa(subs(rot_matrix_inverse(:,:,i)));
       
       w(:,i+1) = (v_var*w(:,i))+(phi(i)*[0; 0; theta_dot(i)]);
       w_dot(:,i+1) = v_var*w_dot(:,i) + (phi(i)*v_var*cross(w(:,i),[0; 0; theta_dot(i)]))+phi(i)*[0; 0; theta_double_dot(i)];

       v(:,i+1) = v_var*(v(:,i)+ cross(w_dot(:,i),pos(:,:,i)) + cross(w(:,i),cross(w(:,i),pos(:,:,i))));

       v_c(:,i) =  v(:,i+1)+cross(w_dot(:,i+1),s(:,i))+cross(w(:,i+1),cross(w(:,i+1),s(:,i)));
       F(:,i) = m(1,i)*v_c(:,i);
       N(:,i) = I(:,:,i)*w_dot(:,i+1) + cross(w(:,i+1),I(:,:,i)*w(:,i+1));

    end
    % Outward Propagation
    for i=3:-1:1
        
        v_var_inw = vpa(subs(rot_matrix(:,:,i+1)));
        f(:,i) = v_var_inw*f(:,i+1) + F(:,i);
        n_var(:,i) = v_var_inw*n_var(:,i+1) + N(:,i)+ cross(pos(:,:,i+1),v_var_inw*f(:,i+1))+cross(s(:,i),F(:,i));
    end
    

    for i=3:-1:1
        test = n_var(3,1:3);
        % saving the final torque values
       n_final(1,x:x+2) = test; 
    end
    x = x+3;
    
end


%% Extracting data

n_wrist = n_final(3:3:end);
n_elbow = n_final(2:3:end);
n_shoulder = n_final(1:3:end);

%% Export data

col1 = time(:,1);
col2 = transpose(n_shoulder(1,:));
col3 = transpose(n_elbow(1,:));
col4 = transpose(n_wrist(1,:));

T = [col1, col2,col3,col4];
writematrix(T,'torque_readings.txt', 'Delimiter', ' ');

% Rotational angles

T1 = [data(:,1),data(:,2),data(:,5),data(:,8)];
writematrix(T1,'jointangle.txt', 'Delimiter', ' ');

% Rotational velocity
T2 = [data(:,1),data(:,3),data(:,6),data(:,9)];
writematrix(T2,'jointanglevelocity.txt', 'Delimiter', ' ');

%%

%plots
figure(1)
plot(data(:,1),col2,data(:,1),col3,data(:,1),col4)
title('Torque Shoulder-Elbow-Wrist (without PID)')
legend('Shoulder','Wrist','Elbow');

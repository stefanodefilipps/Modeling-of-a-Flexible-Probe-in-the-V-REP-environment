% Define all the robots links given the DH table

l1 = Link ("d", 0 , "a", 0.01+0.005 , "alpha", 0);
links = [l1];
for(i=2:31)
    links(i) = Link ("d", 0 , "a", 0.01 , "alpha", 0);
end
links(32) = Link ("d", 0 , "a", 0.01 , "alpha", 0);
robot = SerialLink(links);
% robot.fkine(zeros(1,32))
% robot.jacob0(zeros(1,32))


% I suppose that i Know to which link I am applying the force

contacted_links = [11,32];

% I build a robot of n links up to the interested joint so to compute the
% partial jacobian
joints = parse_joint_file("/home/stefano/MedicalRoboticsProject/FINAL_Measured_Angles.txt");
torques = parse_torques_file("/home/stefano/MedicalRoboticsProject/FINAL_Measured_Torques.txt");
filtered_torques = Filter_torques(torques,5,5);
filtered_joints = Filter_joints(joints,5,5);
% plot(1:size(f_ext_i,1),filtered_torques(:,1));
% hold on
plot(1:size(f_ext_i,1),torques(:,1));
% plot(1:size(f_ext_i,1),filtered_joints(:,1));
% hold on
% plot(1:size(f_ext_i,1),joints(:,1));
measured_forces = parse_forces_file("/home/stefano/MedicalRoboticsProject/FINAL_Measured_Forces.txt");
lambda_fE=[20e-4 4e-6 0 0 0 0];
K = diag([8.5,ones(1,31)*1.5]);
[f_ext_i,S] = Estimate_Forces_ith(filtered_joints,filtered_torques,links,32,0.1,0.05,lambda_fE,K);
f_ext_i
x_filter_move=movmean(f_ext_i(:,1),[minindex maxindex]);
y_filter_move=movmean(f_ext_i(:,2),[minindex maxindex]);
%% FILTER
%Moving average (off-line)
minindex=5; %number of elements backward
maxindex=5; %number of elements forward
measured_x_filter_move=movmean(measured_forces(:,1),[minindex maxindex]);
measured_y_filter_move=movmean(measured_forces(:,2),[minindex maxindex]);
forces_filtered = [measured_x_filter_move,measured_y_filter_move];
size(forces_filtered)

%plot smoothed version
%First plot x
% plot(1:size(f_ext_i,1),x_filter_move);
% hold on
% plot(1:size(f_ext_i,1),measured_x_filter_move);
% legend("est fx","meas fx");
% Then i try to plot the y component
plot(1:size(f_ext_i,1),y_filter_move);
hold on
plot(1:size(f_ext_i,1),measured_y_filter_move);
hold on
legend("est fy","meas fy");
%  legend("est fy","meas fy","est fx","meas fx");

% smoothed = filter_signal(f_ext_i(:,1));
% plot(1:size(f_ext_i,1),f_ext_i(:,2));
% hold on
%plot(1:size(f_ext_i,1),measured_forces(:,1));
% hold on
% plot(1:size(f_ext_i,1),measured_forces(:,2));
% hold on
% plot(1:size(f_ext_i,1),smoothed);

function joints = parse_joint_file(file_name)
joints = zeros(1,32);
fid = fopen(file_name);
tline = fgetl(fid);
index = 1;
while ischar(tline)
    splitted = split(strtrim(tline),";");
    for(i=1:length(splitted)-1)
        number = splitted{i};
        if contains(number,"e-") == 1
            number = 0;
        else
            number = str2num(number);
        end
        joints(index,i) = number;
    end
    tline = fgetl(fid);
    index = index+1;
end
fclose(fid);
end

function measured_forces = parse_forces_file(file_name)
measured_forces = zeros(1,2);
fid = fopen(file_name);
tline = fgetl(fid);
index = 1;
while ischar(tline)
    splitted = split(strtrim(tline),";");
    for(i=1:length(splitted)-1)
        number = splitted{i};
        if contains(number,"e-") == 1
            number = 0;
        else
            number = str2num(number);
        end
        measured_forces(index,i) = number;
    end
    tline = fgetl(fid);
    index = index+1;
end
fclose(fid);
end

function torques = parse_torques_file(file_name)
torques = zeros(1,32);
fid = fopen(file_name);
tline = fgetl(fid);
index = 1;
while ischar(tline)
    splitted = split(strtrim(tline),";");
    for(i=1:length(splitted)-1)
        number = splitted{i};
        if contains(number,"e-") == 1
            number = 0;
        else
            number = str2num(number);
        end
        torques(index,i) = number;
    end
    tline = fgetl(fid);
    index = index+1;
end
fclose(fid);
end

%functon to estimate force at link i

function [f_ext_i,S] = Estimate_Forces_ith(joints,torques,links,idx,eps,lambda_max,lambda_fE,K)
f_ext_i = zeros(size(joints,1),6);
S = zeros(size(joints,1),6);

partial_robot = SerialLink(links(1:idx));

for(i=1:size(joints,1))
    current_joints = zeros(1,32);
    current_torques = zeros(32,1);
    for(j=1:size(joints,2))
        current_joints(j) = joints(i,j);
        current_torques(j) = torques(i,j);
    end
    J = partial_robot.jacob0(current_joints(1:idx));
    %jacobian_T = jacobian';
    size(torques);
    J_DLS=pinv(diag(lambda_fE)+J*J')*J;
    %current_torques = current_torques-K*current_joints';
    f_ext = J_DLS*current_torques(1:idx);
    jj = 1;
    f_ext_i(i,:) = f_ext;
end
end


% This function will smooth the x components that has a lot of variance
% apparently

function smoothed = filter_signal(signal)
std_ = std(signal);
mean_ = mean(signal);
smoothed = zeros(1,length(signal));
for(i=1:length(signal))
    value = signal(i);
    if(value > mean_)
        smoothed(i) = value - std_;
    else
        smoothed(i) = value + std_;
    end
end
end

% function error = compute_error(real,sig)
% error(size(real,1));
% for(i=1:size())
% end


function filtered_torques = Filter_torques(torques,min,max)
minindex=min; %number of elements backward
maxindex=max; %number of elements forward
filtered_torques = zeros(size(torques));
for(i=1:size(torques,2))
    filtered_torques(:,i) = movmean(torques(:,i),[minindex maxindex]);
end
end


function filtered_joints = Filter_joints(joints,min,max)
minindex=min; %number of elements backward
maxindex=max; %number of elements forward
filtered_joints = zeros(size(joints));
for(i=1:size(joints,2))
    filtered_joints(:,i) = movmean(joints(:,i),[minindex maxindex]);
end
end





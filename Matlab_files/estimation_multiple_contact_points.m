clear all
close all

% Define all the robots links given the DH table

l1 = Link ("d", 0 , "a", 0.01+0.005 , "alpha", 0);
links = [l1];
for(i=2:31)
    links(i) = Link ("d", 0 , "a", 0.01 , "alpha", 0);
end
links(32) = Link ("d", 0 , "a", 0.01 , "alpha", 0);
robot = SerialLink(links);


% I suppose that i Know to which link I am applying the force

contacted_links = [7,8,9,10,25,26];

% I build a robot of n links up to the interested joint so to compute the
% partial jacobian
joints = parse_joint_file("/home/stefano/MedicalRoboticsProject/force_experiment/exp6/FINAL_Measured_Angles.txt");
torques = parse_torques_file("/home/stefano/MedicalRoboticsProject/force_experiment/exp6/FINAL_Measured_Torques.txt");
filtered_torques = Filter_torques(torques,5,5);
filtered_joints = Filter_joints(joints,5,5);
measured_forces = parse_forces_file("/home/stefano/MedicalRoboticsProject/force_experiment/exp6/FINAL_Measured_Forces.txt");
%lambda_fE=[1.5e-1 2.5e-6 0 1.5e-1 2.5e-6 0 1.5e-1 2.5e-6 0];
%% FILTER
%Moving average (off-line)
minindex=5; %number of elements backward
maxindex=5; %number of elements forward
measured_x_filter_move=movmean(measured_forces(:,1),[minindex maxindex]);
measured_y_filter_move=movmean(measured_forces(:,2),[minindex maxindex]);
forces_filtered = [measured_x_filter_move,measured_y_filter_move];
size(forces_filtered)


lambda_fE=[8e-5 4e-6 0 8e-5 4e-6 0 8e-5 4e-6 0 8e-5 4e-6 0 8e-5 4e-6 0 8e-5 4e-6 0];
K = diag([8.5,ones(1,31)*1.5]);
f_estimated = Estimate_Forces(filtered_joints,filtered_torques,contacted_links,links,lambda_fE,K);
sum_forces = compute_sum_estimation(f_estimated)
sum_x_filter_move=movmean(sum_forces(:,1),[minindex maxindex]);
sum_y_filter_move=movmean(sum_forces(:,2),[minindex maxindex]);

t = tiledlayout(2,1);
nexttile
plot(1:size(sum_forces,1),sum_x_filter_move);
hold on
plot(1:size(sum_forces,1),measured_x_filter_move);
legend("est fx","meas fx");
nexttile
plot(1:size(sum_forces,1),sum_y_filter_move);
hold on
plot(1:size(sum_forces,1),measured_y_filter_move);
legend("est fy","meas fy");

% 
% %First i try to plot the y component
% plot(1:size(sum_forces,1),sum_forces(:,2));
% hold on
% plot(1:size(sum_forces,1),measured_forces(:,2));
% hold on
% %Now I also try to plot the force on the x direction
% %plot(1:size(sum_forces,1),sum_forces(:,1));
% plot(1:size(sum_forces,1),sum_x_filter_move);
% hold on
% plot(1:size(sum_forces,1),measured_forces(:,1));

%plot smoothed version
% 
% First i try to plot the y component
% plot(1:size(sum_forces,1),sum_y_filter_move);
% hold on
% plot(1:size(sum_forces,1),measured_y_filter_move);
% legend("est fy","meas fy");
% hold on
% Now I also try to plot the force on the x direction
% plot(1:size(sum_forces,1),sum_x_filter_move);
% hold on
% plot(1:size(sum_forces,1),measured_x_filter_move);
% legend("est fx","meas fx");

% legend("est fy","meas fy","est fx","meas fx");
% 
% axis([0 200 -0.1 0.1]);



% function to compute all the partial jacobians and put them in a cell
% structure

function jacobians = compute_jacobians(q,partial_robots,contacted_links)
jacobians = cell(1,length(partial_robots));
for(i=1:length(partial_robots))
    jacobians{i} = partial_robots(i).jacob0(q(1:contacted_links(i)));
end
end

% Now i need to extend the partial jacobians and take only the first 2 rows
% since I am interested in the force in the x,y plane

function Js = extend_jacobians(jacobians)
Js = cell(1,length(jacobians));
for(i=1:length(jacobians))
    jac = zeros(3,32);
    current_jacobian = jacobians{i};
    joint_number = size(current_jacobian,2);
    jac(1:3,1:joint_number) = current_jacobian(1:3,:);
    Js{i} = jac;
end
end

% I need a function that builds the correct J total matrix that needs to be
% pseudoinversed

function JT = force_mapping(Js)
JT = [];
for(i=1:length(Js))
    JT = [JT Js{i}'];
end
end

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

function f_estimated = Estimate_Forces(joints,torques,contacted_links,links,lambda_fE,K)
%f_estimated = [size(joints,1),length(contacted_links)*2];
f_estimated = cell(size(joints,1),length(contacted_links));
partial_robots = [SerialLink(links(1:contacted_links(1)))];

for(i=2:length(contacted_links))
    partial_robots(i) = SerialLink(links(1:contacted_links(i)));
end

for(i=1:size(joints,1))
    current_joints = zeros(1,32);
    current_torques = zeros(32,1);
    for(j=1:size(joints,2))
        current_joints(j) = joints(i,j);
        current_torques(j) = torques(i,j);
    end
    jacobians = compute_jacobians(current_joints,partial_robots,contacted_links);
    Js = extend_jacobians(jacobians);
    JT = force_mapping(Js);
    J = JT';
    J_DLS=pinv(diag(lambda_fE)+J*JT)*J;
    %current_torques = current_torques - K*current_joints';
    f_ext = J_DLS*current_torques;
    %f_ext = -J_DLS*K*current_joints';
    jj = 1;
    for(j=1:3:length(f_ext))
        f_estimated{i,jj} = f_ext(j:j+1);
        jj = jj + 1;
    end
end
end

% Function to compute the sum of the estimated forces in order to compaare
% with the resultant force measured from the force sensor at the base of
% the catheter

function sum_forces = compute_sum_estimation(estimated_forces)
sum_forces = zeros(size(estimated_forces,1),2);
for(i=1:size(estimated_forces,1))
    for(j=1:size(estimated_forces,2))
        sum_forces(i,1) = sum_forces(i,1) + estimated_forces{i,j}(1);
        sum_forces(i,2) = sum_forces(i,2) + estimated_forces{i,j}(2);
    end
end
end

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





function disp_traj(traj,slow,pts,NPro)

if nargin < 2
    slow = false;
    pts = false;
    NPro = size(traj,3);
end
if nargin < 3
    pts = false;
    NPro = size(traj,3);
end
if nargin < 4
    NPro = size(traj,3);
end

figure('Name','Trajectory Sanity Check')
view(-41.1,16.2);
hold on
for i = 1:NPro
    if ~pts
        plot3(traj(1,:,i),traj(2,:,i),traj(3,:,i))
    else
        plot3(traj(1,:,i),traj(2,:,i),traj(3,:,i),'.k');
    end
    if slow
        pause(0.1)
    end
end
hold off


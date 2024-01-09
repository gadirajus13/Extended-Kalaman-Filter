
test

function test
syms v w dt x y z roll pitch yaw real
s = [x y z roll pitch yaw]'
u = [v w]';

% Function 1
state_out = sys_evaluate_g(s,v,w,dt);
ccode(state_out)
% Function 2
G = sys_evaluate_G(s,v,w,dt);
ccode(G)
% Function 3
WMWt = sys_evaluate_WMWt(s,v,w,dt);
ccode(WMWt)
%% 

% Function 4: Find GPS FIT

% Format file of samples
fid = fopen('xyz_samples.txt', 'r');
% Read the header row
header = fgetl(fid);
nCols = numel(strfind(header, ',')) + 1;
format long
% Read the data from the file into a matrix
data = textscan(fid, '%f%f%f%f%f%f%f%f', 'delimiter', ',');
% Close the file
fclose(fid);
% Convert the cell array to a matrix
initialX = cell2mat(data);
nRows = size(initialX, 1) / nCols;

filename = 'latlong_samples.txt';
fid2 = fopen('latlong_samples.txt', 'r');
latlongsamples = importfile(filename, [1, Inf]);

X = [];
y = [];
count = 0;

for i=1:length(initialX)
    for j=1:length(latlongsamples)
        if(initialX(i,1) == latlongsamples(j,1))
            X = [X; initialX(i,2) initialX(i,3) initialX(i,4)];
            y = [y; latlongsamples(j,4) latlongsamples(j,5) latlongsamples(j,6)];
        end
    end
end
count
% X = [9.61683395889 4.10777798143 -1.86166330955;

Latreg = get_regression(X(:,1),y(:,1))
Longreg = get_regression(X(:,2),y(:,2))
plot(Longreg)
Altreg = get_regression(X(:,3),y(:,3))
end

function state_out = sys_evaluate_g(state,v,w,dt)
x = state(1,1);
y = state(2,1);
z = state(3,1);
roll = state(4,1);
pitch = state(5,1);
yaw = state(6,1);
Rx = ROTX(roll);
Ry = ROTY(pitch);
Rz = ROTZ(yaw);
R_xyz = Rx*Ry*Rz;
trans = [x; y; z];

E1 = [R_xyz trans; 0 0 0 1];
Rz_dt = ROTZ(w*dt);
E2 = [Rz_dt [v*dt; 0;0]; 0 0 0 1];
E = E1*E2;
E=simplify(E);

sx = E(1,4);
sy = E(2,4);
sz = E(3,4);
sroll = (atan2(E(2,3),E(3,3)));
spitch = (asin(E(1,3)));
syaw =atan2(E(1,2),E(1,1));

state_out = [sx sy sz sroll spitch syaw]';
state_out = (state_out);
end

function G = sys_evaluate_G(state,v,w,dt)

g = sys_evaluate_g(state,v,w,dt);
G = simplify(jacobian(g,state'));

end

function WMWt = sys_evaluate_WMWt(state,v,w,dt)
syms alpha1 alpha2 alpha3 alpha4
alphas = [alpha1 alpha2 alpha3 alpha4];
vals = [alphas(1)*v.^2+alphas(2)*w.^2 alphas(3)*v.^2+alphas(4)*w.^2];
M = diag(vals);
g = sys_evaluate_g(state,v,w,dt);
W = simplify(jacobian(g,[v w]));
WMWt = simplify(W*M*W');
end

function reg = get_regression(X,Y)

reg = fitlm(X,Y);
plot(reg);

end

function [nodes, element, BC, theta, Surf] = ReadMesh(Input)


fid = fopen(Input);

fgets(fid);
n = textscan(fid, '%f %f');

nodes = [n{1} n{2}];
fgets(fid);


e = textscan(fid, '%f %f %f');
element = [e{1} e{2} e{3}];

fgets(fid);

b = textscan(fid, '%f %f');

BC = [b{1} b{2}];

fgets(fid);

t = textscan(fid, '%f %f');

theta = [t{1} t{2}];

fgets(fid);

s = textscan(fid, '%f %f');

Surf = [s{1} s{2}];

fclose(fid);

% figure
% hold on
% axis equal
% plot(nodes(:,1),nodes(:,2),'.','Color','b','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2);
% hold off

% Solution plotter for IBVP output (surface plot)
clc, clear all

name = 'MATH4340_Prog06_Q1_Euler';
thex = textread('spacecoord.txt','%f');
thet = textread('timecoord.txt','%f');
[x,t] = meshgrid(thex,thet);
Uapprox = dlmread('approximatesol.txt');

figure(1)
mesh(x,t,Uapprox)
view(45,135);
title('Approximate solution')
saveas(gcf,strcat(name,'_ForEulApprox.png'))

if exist('truesol.txt','file')==2
    Utrue = dlmread('truesol.txt');
    figure(2)
    mesh(x,t,Utrue)
    view(45,135);
    title('True solution')
    saveas(gcf,strcat(name,'_true.png'))
end

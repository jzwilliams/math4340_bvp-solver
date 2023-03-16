% Solution plotter for MATH 4340 output
clc, clear all

name = 'MATH4340_Prog06_Q2-2-8_True-32.png';
[apprX,apprY] = textread('approximatesol.txt','%f %f');

if exist('truesol.txt','file') == 2
    [trueX,trueY] = textread('truesol.txt','%f %f');
    plot(apprX,apprY,trueX,trueY);
    legend('Approximate solution','Analytic Solution','Location',...
        'NorthEast')
else
    plot(apprX,apprY)
    title('Approximate solution')
end

saveas(gcf,name)
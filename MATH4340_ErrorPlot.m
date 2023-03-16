% Applies to program 05 as well.
name = 'MATH4340_Prog06_Q1TimeErr';

[h,err] = textread(strcat(name,'.txt'),'%f %f');

[sorth,index] = sort(h);
sorterr = err(index);

logh = log(h);
logerr = log(err);
[slope, yint] = polyfit(logh,logerr,1);
slope

loglog(sorth,sorterr)
legend('error','Location','South')
title('Error against h')

saveas(gcf,strcat(name,'.png'))
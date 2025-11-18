% A SOFT-KILL BESO CODE BY X. HUANG and Y.M. Xie
% Huang x, xie YM (2010) Topology optimization of continuum
%structures: methods and applications, Wiley, Chichester.
%doi:101002/9780470689486
function softkill(nelx,nely,volfrac,er,rmin)
% INITIALIZE
E=3e9;
nu=0.4;
xmin=1e-3;x(1:nely,1:nelx)=1.;
vol=1;loop=0;
dc = zeros(nely,nelx);
change=1;penal=1;
% START ITH ITERATION
while change>0.001&&loop<150
loop=loop+1;
if loop >1
    olddc=dc;
    vo1=max(vol*(1-er),volfrac); 
end
%FE-ANALYSIS
[U]=FEA(nelx,nely,x,penal,E,nu);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
[KE]=lk(E,nu);
c(loop)=0.;
for ely=1:nely
    for elx=1:nelx
        n1=(nely+1)*(elx-1)+ely;
        n2=(nely+1)*elx + ely;
        Ue=U([2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
        c(loop)=c(loop) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
        dc(ely,elx)=0.5 * penal*x(ely,elx)^penal*Ue'*KE*Ue;
    end
end
% FILTERING OF SENSITIVITIES
[dc]=check(nelx,nely,rmin,dc);
% STABLIZATION OF EVOLUTIONARY PROCESS
if loop>1
    dc=(dc+olddc)/2.;
                                                                                                end
% BESO DESIGN UPDATE
if loop>1
    [x]=ADDDEL(nelx,nely,vol,dc,x,xmin);
end
%DETERMINE CONVERGENCE FACTOR
if loop>10
    change=abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
end
% PRINT RESUTTS
fprintf('it.:%4i obj.:%10.4f Vol.:%6.3f ch.:%6.3f\n',loop,c(loop),sum(sum(x))/nelx/nely,change);
    colormap(gray);imagesc(1-x);caxis([0,1]);axis equal;axis off;
    drawnow;
end
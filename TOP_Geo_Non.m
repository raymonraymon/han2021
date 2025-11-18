function TOP_Geo_Non(nelx,nely,volfrac,er,rmin)
%initialize
L=nelx; H=nely;dx=L/nelx;dy=H/nely;a=dx/2;
E=3e9;nu=0.4;xmin=1e-3;x(1:nely,1:nelx)=1.;
penal = 1;vol=1;loop=0;change = 1;
%start ith iteration
while change >0.001 && loop<150
    loop = loop+1;
    if loop>1
        olddc = dc;
        vol = max(vol*(1-er),volfrac);
    end
    %FE-analysis and sensitivity analysis
    [U,dc,KE_matrix]= FEA9(nelx,nely,x,penal,E,nu,dx,dy,a);
    %objective function
    c(loop) = 0.;
    hv4 =1;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(nelx-1)+ely;
            n2 = (nely+1)* nelx   +ely;
            Ue=U([2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
            KE =KE_matrix(:,:,hv4);
            hv4 =hv4+1;
            c(loop) = c(loop)+0.5*x(ely,elx)^penal*Ue'*KE*Ue;
        end
    end
    %Filtering of sensitivities
    [dc] = check(nelx,nely,rmin,dc);
    %stablization of evolutionary process
    if loop >1 
        dc = (dc+olddc)/2;
    end
    %BESO design update
    if loop > 1
        x=ADDDEL(nelx,nely,vol,dc,x,xmin);
    end
    if loop>10
        change = abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
    end
    %print result
    fprintf('it.:%4i obj.:%10.4f Vol.:%6.3f ch.:%6.3f\n',loop,c(loop),sum(sum(x))/nelx/nely,change);
    colormap(gray);imagesc(1-x);caxis([0,1]);axis equal;axis off;
    drawnow;
end
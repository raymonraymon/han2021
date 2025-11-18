 function [U] = FEA(nelx,nely,x,penal,E,nu)
 % FE-ANALYSIS
 [KE] = lk(E,nu);
 K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
 U = zeros(2*(nely+1)*(nelx+1),1);
 for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx +ely;
        edof = [2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
 end
 K = (K+K')/2;
 % Loads
F=sparse(2*(nelx+1)*(nely+1)-nely,1,-100000,2*(nely+1)*(nelx+1),1);
 % Supports
 fixeddofs=1:2*(nely+1);
 alldofs =1:2*(nely+1)*(nelx+1);
 freedofs =setdiff(alldofs,fixeddofs);
 % Displacement field
 U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
 U(fixeddofs,:)= 0;
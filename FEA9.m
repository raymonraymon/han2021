function[U,dc,KE_matrix]=FEA9(nelx,nely,x,penal,E,nu,dx,dy,a)
%loads
F=sparse(2*(nelx+1)*(nely+1)-nely,1,-100000,2*(nely+1)*(nelx+1),1);
%supports
fixeddofs=1:2*(nely+1);
alldofs=1:2*(nelx+1)*(nely+1);
freedofs = setdiff(alldofs,fixeddofs);
%init global stiffness matrix,load and displacement vectors
Kt = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
U = sparse(2*(nelx+1)*(nely+1),1);
Fint = sparse(2*(nelx+1)*(nely+1),1);
BL=1/(4*a)*[-1 0 1 0 1 0 -1 0;
            0 1 0 -1 0 1 0 1 ;
            -1 -1 -1 1 1 1 1 -1];
BN0 = zeros(3,8);
B=BL+BN0;
DE = E/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
A0=dx*dy;detaU=1.;loop=0;
% Iteratively solve the geometric nonlinear fE balance equation
while detaU>0.001&&loop<50
% Convergence condition
loop=loop+1;
hv1=1;hv3=1;
for ely=1:nely
for elx=1:nelx
n1=(nely+1)*(elx-1)+ely;
n2=(nely+1)*elx+ely;
edof=[2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2];
Ue=U(edof);
estrain = x(ely,elx) * B * Ue;
estress = x(ely,elx) * DE * B * Ue;
Fe=B'*estress*A0;
Fint(edof)=Fint(edof)+Fe;
ex = estrain(1);ey =estrain(2);exy = estrain(3);
sx = estress(1);sy =estress(2);sxy = estress(3);
A=[ex 0.5*exy 0 0;0 0 0.5*exy ey;0.5*exy ey ex 0.5*exy];
G=1/(4*a)*...
[-1 0 1 0 1 0 -1 0;
0 -1 0 1 0 1 0 -1;
-1 0 -1 0 1 0 1 0;
0 -1 0 -1 0 1 0 1];
M=[sx 0 sxy 0;0 sx 0 sxy;sxy 0 sy 0;0 sxy 0 sy];
BN=A*G;
B=BL+BN;
KL=BL'*DE*BL*A0;
KS=G'*M*G*A0;
KN=(BL'*DE*BN+BN'*DE*BL+BN'*DE*BN)*A0;
KE=KL+KN+KS;
Kt(edof,edof)=Kt(edof,edof)+x(ely,elx)^penal*KE;
Bmatrix(:,:,hv1)=B;
KE_matrix(:,:,hv3)=KE;
hv1=hv1+1;hv3=hv3+1;
end
end
Kt=(Kt+Kt')/2;
% Displacement field
R=F-Fint;
Uk(freedofs,:) = Kt(freedofs,freedofs) \ R(freedofs,:);
Uk(fixeddofs,:)=0;
U=U+Uk;
detaU=norm(Uk,2)/norm(U,2);
end
%sensitivy analysis
hv2=1;
for ely=1:nely
for elx=1:nelx
n1=(nely+1)*(elx-1)+ely;
n2=(nely+1)*elx +ely;
edof=[2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2];
Ue=U(edof);
B=Bmatrix(:,:,hv2);
estress=x(ely,elx)*DE*B*Ue;
Fint=B'*estress*A0;
dc(ely,elx)=0.5*Ue'*Fint;
hv2=hv2+1;
end
end
%==_=_=_=
 %Bidirectional Evolutionary structural Optimization (BESO)
% stiffness Design with Geometric Nonlinearity
%2D Cantilever Beam
% Yongsheng Han.
%School of Mechanics, Civil Engineering and Architecture,
% Northwestern Polytechnical University, Xi'an, Shaanxi, China.
% Email:hanys0407@mail.nwpu.edu.cn
%September 25£¬2020
        
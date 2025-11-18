function [x]=ADDDEL(nelx,nely,volfra,dc,x,xmin)
%% BESO DESIGN UPDATE
 l1 = min(min(dc)); l2 = max(max(dc));
 while ((l2-l1)/l2 > 1.0e-5)
    th = (l1+l2)/2.0;
    x = max(xmin,sign(dc-th));
    if sum(sum(x))-volfra*(nelx*nely) > 0;
        l1 = th;
    else
        l2 = th;
    end
 end
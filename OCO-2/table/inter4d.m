function [data] = inter4d(x,t,p,absco,x_grid,t_grid1,t_grid2,p_grid)
%输出8组
    if x==x_grid(end)
        [k1,k2]=size(x_grid);
    else
        [k1,k2]=find((x_grid-x)>0);
    end
    nx2=k1(1);nx1=nx2-1;

    if p==p_grid(end)
        [k1,k2]=size(p_grid);
    else
        [k1,k2]=find((p_grid-p)>0);
    end
    np2=k1(1);np1=np2-1;

    if t==t_grid1(end)
        [k1,k2]=size(t_grid1);
    else
        [k1,k2]=find((t_grid1-t)>0);
    end
    nt12=k1(1);nt11=nt12-1;
    if t==t_grid2(end)
        [k1,k2]=size(t_grid2);
    else
        [k1,k2]=find((t_grid2-t)>0);
    end
    nt22=k1(1);nt21=nt22-1;
    profilenx1=[np1,nt11,nx1;np1,nt12,nx1;np2,nt21,nx1;np2,nt22,nx1];
    profilenx2=[np1,nt11,nx2;np1,nt12,nx2;np2,nt21,nx2;np2,nt22,nx2];
    profilen=[profilenx1;profilenx2];
    profilex1=[p_grid(np1),t_grid1(nt11),x_grid(nx1);p_grid(np1),t_grid1(nt12),x_grid(nx1);p_grid(np2),t_grid2(nt21),x_grid(nx1);p_grid(np2),t_grid2(nt22),x_grid(nx1)];
    profilex2=[p_grid(np1),t_grid1(nt11),x_grid(nx2);p_grid(np1),t_grid1(nt12),x_grid(nx2);p_grid(np2),t_grid2(nt21),x_grid(nx2);p_grid(np2),t_grid2(nt22),x_grid(nx2)];
    profile=[profilex1;profilex2];

    for i=1:8
        data1(:,i)=absco(:,profilen(i,3),profilen(i,2),profilen(i,1));%wave,x,t,p
    end
    data1=data1';
%插值T变4组%  3.5=interp1([1,2],[3,4],1.5) 
    for i=1:4
        data2(:,i)=interp1([profile(i*2-1,2),profile(i*2,2)],[data1(i*2-1,:);data1(i*2,:)],t);
    end
    data2=data2';
%插值p变2组%
    for i=1:2
        data3(:,i)=interp1([p_grid(np1),p_grid(np2)],[data2(i*2-1,:);data2(i*2,:)],p);
    end
    data3=data3';
%插值x结束%
    data=interp1([x_grid(nx1),x_grid(nx2)],[data3(1,:);data3(2,:)],x);
    data=data';
end%输出要是一行


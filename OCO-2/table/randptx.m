function [p,t,x,t_grid1,t_grid2] = randptx(p_grid,temperature,x_grid,pmodel)
    p=rand*(max(pmodel)-min(pmodel))+min(pmodel);x=rand*(max(x_grid)-min(x_grid))+min(x_grid);
    if p==p_grid(end)
        [k1,k2]=size(p_grid);
    else
        [k1,k2]=find((p_grid-p)>0);
    end
    k1=k1(1);
    t_grid1=temperature(:,k1-1);
    t_grid2=temperature(:,k1);
    t_grid_range=[t_grid1(1),t_grid1(end);t_grid2(1),t_grid2(end)];
    if t_grid_range(1,1)>t_grid_range(2,1)
        t_max=t_grid_range(2,2);
        t_min=t_grid_range(1,1);
    else
        t_max=t_grid_range(2,1);
        t_min=t_grid_range(1,2);
    end
    t=rand*(t_max-t_min)+t_min;
end


[az el]=view;
XYZlim=[get(gca,'Xlim');get(gca,'Ylim');get(gca,'Zlim')];
% cla;
% scatter3(x_clustered, y_clustered, z_clustered,36,clr(clustered,:),'Marker','.');
% scatter3(GS_seed(:,1),GS_seed(:,2),GS_seed(:,3),100,clr(1:nEle,:),'Marker','o','LineWidth',3);
% % scatter3(x_all, y_all, z_all,1,'Marker','o');
% view(az,el);axis vis3d, box on;rotate3d on;
% 
% xlabel('x'),ylabel('y'),zlabel('z');
set(gca,'Xlim',XYZlim(1,:));
set(gca,'Ylim',XYZlim(2,:));
set(gca,'Zlim',XYZlim(3,:));
view(az,el);
if strcmp(get(S.all,'Visible'),'off')
    set(S.all,'Visible','on');
else
    set(S.all,'Visible','off');
end

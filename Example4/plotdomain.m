x_coord = [0;400;800;1200;1600;1600;0];
y_coord = [150;100;150;100;150;-1000;-1000];
figure;
hold on;
fill(x_coord,y_coord,'w')
plot([400;1500],[100;-1000],'r','LineWidth',2)
plot([1200;1000],[100;-1000],'r','LineWidth',2)
axis equal; axis off; ax=axis;axis(ax*1.001);
text(0,150+25,'1','HorizontalAlignment','center','Color','k','FontSize',14)
text(400,100+25,'2','HorizontalAlignment','center','Color','k','FontSize',14)
text(800,150+25,'3','HorizontalAlignment','center','Color','k','FontSize',14)
text(1200,100+25,'4','HorizontalAlignment','center','Color','k','FontSize',14)
text(1600,150+25,'5','HorizontalAlignment','center','Color','k','FontSize',14)
text(1600,-1000-25,'6','HorizontalAlignment','center','Color','k','FontSize',14)
text(1500,-1000-25,'7','HorizontalAlignment','center','Color','k','FontSize',14)
text(1000,-1000-25,'8','HorizontalAlignment','center','Color','k','FontSize',14)
text(0,-1000-25,'9','HorizontalAlignment','center','Color','k','FontSize',14)

text(800,150-50,'H=height','HorizontalAlignment','center','Color','k','FontSize',14)
text(1600+90,-400,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',14)
text(800,-1000-50,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',14)
text(0-90,-400,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',14)


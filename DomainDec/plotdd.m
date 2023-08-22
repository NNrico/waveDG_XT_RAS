% plot decomposed domain

function [] = plotdd(DataDD,femregion)


X=DataDD.X;
T=DataDD.T;
m=DataDD.m;
n=DataDD.n;
Nx=DataDD.Nx;
NT=DataDD.NT;
nsub_t=DataDD.nsub_t;
nsub_x=DataDD.nsub_x;
list_sub=DataDD.list_sub;

figure
axis equal
%plot whole domain
rectangle_whole = @(x1,y1,x2,y2) rectangle('Position',[x1,y1,x2-x1,y2-y1],'LineWidth',3);
rectangle_whole(0,0,X,T)

%plot elements
rectangle_elements = @(x1,y1,x2,y2) rectangle('Position',[x1,y1,x2-x1,y2-y1],'LineWidth',1,'EdgeColor',"#C0C0C0");
for i =1:length(femregion.coords_element)
    x1=femregion.coords_element{i}(1,1);
    y1=femregion.coords_element{i}(1,2);
    x2=femregion.coords_element{i}(3,1);
    y2=femregion.coords_element{i}(3,2);
    rectangle_elements(x1,y1,x2,y2);
end

%plot subdomain
colors = ["r","b","y","g"];
rectangle_sub = @(x1,y1,x2,y2,i) rectangle('Position',[x1,y1,x2-x1,y2-y1],'LineWidth',2,'EdgeColor',colors(mod(i,4)+1));
for i =1:length(list_sub)
    elem1=list_sub(i);
    x1=femregion.coords_element{elem1}(1,1);
    y1=femregion.coords_element{elem1}(1,2);

    elem2=elem1+m-1+(n-1)*NT;
    x2=femregion.coords_element{elem2}(3,1);
    y2=femregion.coords_element{elem2}(3,2);
    rectangle_sub(x1,y1,x2,y2,i);
end
% xlim([0,1])
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);

%plot caratteristica c=1
% line([0,T],[0,T])



end
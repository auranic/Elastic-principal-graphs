function drawSimple2DIllustration(data,np,ed,varargin)

color = 'r';
linewidth = 3;

for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'Color')
            color = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'LineWidth')
            linewidth = varargin{i + 1};
        end
end


plot(data(:,1),data(:,2),'ko','MarkerEdgeColor',[0.5;0.5;0.5]); hold on;
drawGraph2D(np,ed,'LineWidth',linewidth,'ShowClusterNumbers',0,'LineColor',color);

axis off; axis equal;

end


function SubPlot
% -----------------------------------------------------------
% General usage: This function inserts various MATLAB figure (.fig) files
%                into one figure with multiple subplots
%
% Notes:
%      All desired .fig files should contain only one figure inside it (no
%      subplot inside the .fig files).
%      All desired .fig files should be defined in the 2D space.
%
% Author: Farhad Sedaghati
% Contact: <farhad_seda@yahoo.com>
% Updated: 10/08/2019 to include the multiselect option, and get the legend,
% x and y limit, scale, tick, and ticklabels
% -----------------------------------------------------------
% disp('---------------------------------------------------------------------');
% disp('note that all desired .fig files should contain only one figure inside it.');
% disp('note that all desired .fig files should be defined in the 2D space.');
% disp('Please press enter to continue');
% pause
% Number of the run
N=0; % Initialize 'N'
% answer='y';
% Nf = input('How many figures are you trying to merge?');
% K  = input('Number of figures in a row?');
% while strcmpi(answer,'y')    
%     % Get the path and filename of the desired fig file
% for i = 1:Nf
    filename=0;
    run=0;
    while isequal(filename,0)
        if run==2
            error('please choose your fig file');
        end
        disp('---------------------------------------------------');
        disp('Select the desired fig file(files) which you want to insert in the subplot: ');
        [filename,pathname]=uigetfile('MultiSelect', 'on',{'*.fig';'*.FIG'},'Select the .fig file(files) you want to insert in the subplot');
        run=run+1;
    end
    clear FN % initialize 'FN' container
    if ~iscell(filename) % If single-select is used
        N=N+1;
        FN=[pathname filename];
        % open figure
        h(N) = openfig(FN,'new');
        % get handle to axes of figure
        ax(N)=gca;
    else % If multi-select is used
        for figureIndex = 1:size(filename,2)
            N=N+1;
            FN=[pathname filename{1,figureIndex}];
            % open figure
            h(N) = openfig(FN,'new');
            % resize the figure
            pbaspect([2 1 1]);
            % get handle to axes of figure
            ax(N)=gca;
            disp('pause for 1 seconds to copy figure objects');
            pause(1); % Critical line
        end
    end
%     answer=input('Do you have more .fig files to read? \n','s');
% end
% K=input('How many figures in a row do you want to have? \n');
% if isempty(K)
%     K = 2;
% end
K = 2; % I assumed that the number of columns is always 2.
subplotName = ["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"];
% Creat a figure
figure;
set(gcf,'units','centimeters','position',[10 10 24 12]);
% set(gcf,'units','centimeters','position',[10 10 12 18]);
% set(gcf,'units','centimeters','position',[10 10 12 6]);
for i=1:N
    % create and get handle to the subplot axes
    s(i) = subplot(ceil(N/K),K,i); 
    % get handle to all the children in the figure
    aux=get(ax(i),'children');
    for j=1:size(aux)
        fig(i) = aux(j);
        copyobj(fig(i),s(i)); 
        hold on
    end
    % copy children to new parent axes i.e. the subplot axes
    xlab = get(get(ax(i),'xlabel'),'string');
    ylab = get(get(ax(i),'ylabel'),'string');
%     tit = get(get(ax(i),'title'),'string');
%     Legend = get(get(ax(i),'legend'),'string');
    xLimits = get(ax(i),'XLim');
    yLimits = get(ax(i),'YLim');
    ScaleX = get(ax(i),'Xscale');
    ScaleY = get(ax(i),'Yscale');
    TickX = get(ax(i),'Xtick');
    TickY = get(ax(i),'Ytick');
    TicklabelX = get(ax(i),'Xticklabel');
    TicklabelY = get(ax(i),'Yticklabel');
    set(gca, 'XScale', ScaleX, 'YScale', ScaleY, 'Xtick', TickX, 'Ytick',...
        TickY, 'Xticklabel', TicklabelX, 'Yticklabel', TicklabelY);
    xlabel(xlab,'interpreter','latex');
    ylabel(ylab,'interpreter','latex');
    xlim(xLimits);
    ylim(yLimits);
%     legend(Legend)
    legend('off')
    box on;
    % Additional control on the figure
    figureTitle = subplotName(i);
    title(figureTitle,'Units', 'normalized', 'Position', [-0.22, 1.0]);
    FontSize = 12;
    figureLineWidth = 1;
    boxLineWidth = 1;
    set(findall(gca, 'Type', 'Line'),'LineWidth',figureLineWidth);
    set(gca,'fontsize',FontSize,'linew',1,'fontname','times');
    set(gca,'linew',boxLineWidth);
    pbaspect([2 1 1])
end
return
%% Save file
% /Users/hanulhwang/Documents/RT_instability_connor/Connor_figures/2Dfigures/edited
figurePath = '~/Documents/RT_instability_connor/Connor_figures/merged_figures/';
figureName = 'new_2d_A_kZ';
saveas(gcf,fullfile(figurePath,figureName),'fig');
% saveas(gcf,fullfile(figurePath,figureName),'epsc');
disp('plot saved!');
function plot_network(areas,line,n_areas)
    

    % angles = linspace(0, 2*pi, n_areas); 
    % circle = randi([-5,5],n_areas,2);
    % circle(:,1) = 1.5*cos(angles);
    % circle(:,2) = 1.5*sin(angles);

    NodeColor = zeros(size(areas,1),3);
    % XData = zeros(n_areas,1);
    % YData = zeros(n_areas,1);

    %according to a figure i found online
    coordinates = [57,24
                    121,24
                    77,72
                    89,106
                    87,177
                    158,177
                    200,177
                    104,312
                    104,371
                    104,417
                    158,105
                    236,82
                    286,167
                    335,129
                    349,166
                    227,225
                    314,264
                    394,264
                    441,226
                    409,351
                    421,391
                    435,423
                    439,466
                    566,433
                    391,514
                    334,479
                    170,456
                    148,408
                    148,359
                    322,312
                    210,361
                    250,389
                    584,95
                    629,164
                    551,231
                    650 226
                    674,139
                    689,371
                    660,79
                    682,42
                    740,43
                    776,43
                    734,168
                    751,134
                    760,216
                    746,247
                    777,258
                    808,202
                    861,265
                    902,150
                    923,134
                    818,90
                    825,43
                    867,43
                    976,43
                    916,43
                    902,106
                    922,103
                    1018,60
                    1018,72
                    1018,219
                    1010,397
                    992,96
                    984,211
                    933,454
                    920,397
                    933,368
                    838,427
                    803,426
                    656 434
                    629,483
                    604,462
                    606,502
                    636,551
                    660,573
                    765,537
                    797,571
                    815,515
                    828,492
                    864,537
                    854,514
                    754,617
                    655,643
                    640,661
                    639,706
                    638,733
                    657,749
                    682,706
                    735,705
                    732,750
                    814,750
                    859,705
                    850,680
                    872,642
                    836,642
                    847,617
                    843,580
                    897,654
                    930,538
                    923,643
                    906,732
                    851,733
                    970,732
                    973,643
                    1013,643
                    1000,589
                    1060,643
                    1021,684
                    1020,714
                    1012,760
                    974,783
                    1053,760
                    258,294
                    228,438
                    289,427
                    865,433
                    328,82
                    709,537];
    coordinates(:,2) = -coordinates(:,2);
    coordinates = [coordinates ; coordinates(line(end-13:end,1),:)];

    XData(:,1) = coordinates(:,1);
    XData(end-13:end,1) = XData(end-13:end,1) + ones(14,1)*0;
    YData(:,1) = coordinates(:,2);
    YData(end-13:end,1) = YData(end-13:end,1) + ones(14,1)*40;


     matlab_colors =        [     0         0    1.0000
    1.0000         0         0
         0    1.0000         0
         0         0    0.1724
    1.0000    0.1034    0.7241
    1.0000    0.8276         0
         0    0.3448         0
    0.5172    0.5172    1.0000
    0.6207    0.3103    0.2759
         0    1.0000    0.7586
         0    0.5172    0.5862
         0         0    0.4828
    0.5862    0.8276    0.3103
    0.9655    0.6207    0.8621
    0.8276    0.0690    1.0000
    0.4828    0.1034    0.4138
    0.9655    0.0690    0.3793
    1.0000    0.7586    0.5172
    0.1379    0.1379    0.0345
    0.5517    0.6552    0.4828
    75/255    220/255   59/255  %asdasda
    0.5172    0.4483         0
    0.4483    0.9655    1.0000
    0.6207    0.7586    1.0000
    0.4483    0.3793    0.4828
    0.6207         0         0
         0    0.3103    1.0000
         0    0.2759    0.5862
    0.8276    1.0000         0
    0.7241    0.3103    0.8276
    0.2414         0    0.1034
    0.4660 0.6740 0.1880
    1.0000    0.4828    0.3793
    0.2759    1.0000    0.4828
    0.0690    0.6552    0.3793
    0.8276    0.6552    0.6552
    0.8276    0.3103    0.5172
    0.4138         0    0.7586
    0.1724    0.3793    0.2759
         0    0.5862    0.9655
    0.0345    0.2414    0.3103
    0.6552    0.3448    0.0345
    0.4483    0.3793    0.2414
    0.0345    0.5862         0
    0.6207    0.4138    0.7241
    1.0000    1.0000    0.4483
    0.6552    0.9655    0.7931
    0.5862    0.6897    0.7241
    0.6897    0.6897    0.0345
    0.1724         0    0.3103];



    matlab_colors = [
        255, 0, 0;       % Red
        0, 255, 0;       % Green
        0, 0, 255;       % Blue
        255, 255, 0;     % Yellow
        255, 0, 255;     % Magenta
        0 0.4470*255 0.7410*255     % Cyan
        128, 0, 0;       % Maroon
        255, 215, 0;     % Gold (Replaced Dark Green)
        0, 0, 128;       % Navy
        128, 128, 0;     % Olive
        0.4940*255 0.1840*255 0.5560*255
        0, 128, 128;     % Teal asdasdasd
        192, 192, 192;   % Silver
        128, 128, 128;   % Gray
        255, 165, 0;     % Orange
        255, 192, 203;   % Pink
        139, 69, 19;     % Saddle Brown
        210, 105, 30;    % Chocolate
        0, 255, 127;     % Spring Green
        255, 215, 0;     % Gold
        75, 220, 59;      % Red-Orange
        255, 20, 147;    % Deep Pink
        70, 130, 180;    % Steel Blue
        0, 255, 255;     % Aqua
        255, 140, 0;     % Dark Orange
        0, 0, 0;         % Black
        128 128 128;   % White
        70, 130, 180;    % Sky Blue (Replaced Ivory)
        128, 0, 0;       % Dark Red
        0, 128, 0;       % Dark Green
        0, 0, 128;       % Dark Blue
        255, 99, 71;     % Tomato
        0, 128, 128;     % Dark Cyan sadas
        128, 0, 128;     % Dark Magenta
        255, 140, 105    % Salmon
    ];


    matlab_colors = matlab_colors./255;

    %matlab_colors = lines(7);

    % K=3;
    % matlab_colors=reshape(matlab_colors(:,perms(1:K)),[],K);

    g = graph(line(:,1),line(:,2));
    %g.Nodes.NodeColors = degree(g);
    

    figure('Position',4*[0 0 412 266]);
    hold on


    R = 0.1;
    %node_names = {};
    for i = 1:n_areas
        
        mask = areas == i;
        n = sum(mask);
        
        % cells = cell(1, n);
        % area_name = sprintf('Area %d',i);
        % cells(:) =  {area_name};
        % node_names = {node_names,cells};
        
        NodeColor(mask,:) = repmat(matlab_colors(i,:),n,1);
        
    end
    xlim([min(XData)-10 max(XData)+200])
    ylim([min(YData)-20 max(YData)+20])
    set(gca, 'XTick', [], 'XTickLabel', []);
    set(gca, 'YTick', [], 'YTickLabel', []);
    axis(gca,'equal')
    set(get(gca, 'XAxis'), 'Visible', 'off');
    set(get(gca, 'YAxis'), 'Visible', 'off');


    gscatter(XData,YData,areas,matlab_colors)

    h = plot(g);
    lgd = legend({'Area 1','Area 2','Area 3','Area 4','Area 5','Area 6','Area 7','Area 8','Area 9','Area 10','Area 11','Area 12','Area 13','Area 14','Area 15','Area 16','Area 17','Area 18','Area 19','Area 20','Area 21','Area 22','Area 23','Area 24','Area 25','Area 26','Area 27','Area 28','Area 29','Area 30','Area 31','Area 32','Area 33','Area 34','Area 35'},'Location','east');
    lgd.FontSize = 12;
    hold off;
    


    h.LineStyle = '--';
    h.LineWidth = 1;
    h.NodeColor = NodeColor;
    h.XData = XData;
    h.YData = YData;
    h.MarkerSize = 5;
    set(gcf,'renderer','Painters');
    title='./fig/graph_areas.eps';
    saveas(gca,title,'epsc');




    
end




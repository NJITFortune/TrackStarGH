function thetaplayer(in, im)

figure(1); 

framestoskip = 8; % Higher numbers for faster, less precise plaback.

for j=1:framestoskip:length(im)-1 
    clf;
    imshow(im(j).gray, 'InitialMagnification', 300); hold on;
    
% Fish 1
    if in(1).theta(j) < 90
        plot(in(1).x(j), in(1).y(j), 'g*', 'MarkerSize', 10);
    end
    if in(1).theta(j) >= 90 && in(1).theta(j) < 180
        plot(in(1).x(j), in(1).y(j), 'c*', 'MarkerSize', 10);
    end
    if in(1).theta(j) >= 270
        plot(in(1).x(j), in(1).y(j), 'r*', 'MarkerSize', 10);
    end
    if in(1).theta(j) >= 180 && in(1).theta(j) < 270
        plot(in(1).x(j), in(1).y(j), 'm*', 'MarkerSize', 10);
    end
    plot([in(1).x(j) in(1).x(j)+(20*sin((in(1).theta(j)/360)*2*pi))], ...
        [in(1).y(j) in(1).y(j)-(20*cos((in(1).theta(j)/360)*2*pi))], '-g');
% Fish 2
if length(in) == 2
    if in(2).theta(j) < 90
        plot(in(2).x(j), in(2).y(j), 'go', 'MarkerSize', 10);
    end
    if in(2).theta(j) >= 90 && in(2).theta(j) < 180
        plot(in(2).x(j), in(2).y(j), 'co', 'MarkerSize', 10);
    end
    if in(2).theta(j) >= 270
        plot(in(2).x(j), in(2).y(j), 'ro', 'MarkerSize', 10);
    end
    if in(2).theta(j) >= 180 && in(2).theta(j) < 270
        plot(in(2).x(j), in(2).y(j), 'mo', 'MarkerSize', 10);
    end
    plot([in(2).x(j) in(2).x(j)+(20*sin((in(2).theta(j)/360)*2*pi))], ...
        [in(2).y(j) in(2).y(j)-(20*cos((in(2).theta(j)/360)*2*pi))], '-m');
end

    text(30,30, num2str(in(1).theta(j)), 'Color', 'w');
    text(100,30, num2str(j), 'Color', 'g');
    
    text(30,40, '*', 'Color', 'r', 'FontSize', 18);
    text(40,40, '*', 'Color', 'g', 'FontSize', 18);
    text(30,50, '*', 'Color', 'm', 'FontSize', 18);
    text(40,50, '*', 'Color', 'c', 'FontSize', 18);

    pause(0.05); 
end

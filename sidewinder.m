function out = sidewinder(in, im)

out = in;

%% Initial plot

% Cutoff for what is a potential jump in theta.
maxcut = 330; % Usually 330 
mincut = 75; % Usually 25

    figure(1); clf;
    plot(abs(diff(in.theta)), '-*'); 
    hold on;
    plot([1 length(in.theta)], [maxcut maxcut], 'm');
    plot([1 length(in.theta)], [mincut mincut], 'c');
    ylim([0 360]);

%% Find potential reversals
    idx = find(abs(diff(in.theta)) > mincut & abs(diff(in.theta)) < maxcut);
    % idx = idx(2:end); % Omit first because why?

    length(idx)
    
    suebee = -15:5:15;
    
    j=1;

while j < length(idx)+1
    
    figure(2); clf;
    imshow(im(idx(j)).gray, 'InitialMagnification', 300); 
    hold on;
    
            % Plot trajectory of the fish
        for pp = 1:length(suebee)
            ff(pp) = suebee(pp) + idx(j);
            if ff(pp) <= 0; ff(pp) = 1; end
            if ff(pp) > length(out.x); ff(pp) = length(out.x); end
        end    
        for yy = 2:length(ff)
            plot([out.x(ff(yy-1)), out.x(ff(yy))], [out.y(ff(yy-1)), out.y(ff(yy))], 'b-', 'LineWidth', yy/2);
        end
        
    for kk = 1:length(suebee)
        if suebee(kk) > 0; colr = '-c'; end
        if suebee(kk) < 0; colr = '-m'; end
        
        cdx = idx(j)+suebee(kk);
        if cdx <= 0; cdx = 1; end
        if cdx > length(out.x); cdx = length(out.x); end
        
    plot([out.x(cdx) out.x(cdx)+(20*sin((out.theta(cdx)/360)*2*pi))],... 
        [out.y(cdx) out.y(cdx)-(20*cos((out.theta(cdx)/360)*2*pi))], colr); % Line showing direction    
    plot(out.x(cdx), out.y(cdx), 'g*'); % Dot where the fish is located

    end
            
    % Replot current direction in white
    plot([out.x(idx(j)) out.x(idx(j))+(40*sin((out.theta(idx(j))/360)*2*pi))],... 
        [out.y(idx(j)) out.y(idx(j))-(40*cos((out.theta(idx(j))/360)*2*pi))], '*-w');
    plot(out.x(idx(j)), out.y(idx(j)), 'g*');
    
    text(300, 10, num2str(idx(j)), 'Color', 'w');
    % text(250, 10, num2str(in.theta(idx(j))), 'Color', 'y');
    asdf = input('1=do nothing, 2=rotate, 3=go back one: ');
    
    if asdf == 2 % Rotate!
        out.theta(idx(j)+1:end) = out.theta(idx(j)+1:end) + 180;
        out.theta(out.theta > 359) = out.theta(out.theta > 359)-360;
    end
    
    if asdf == 3 % Made a mistake, go back 1
       j = j-2; 
    end
    
    j=j+1;
    
end

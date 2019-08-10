function scrollerdoller(s, im, rago)

lingth = 10;

    f = figure;
    % axf = axes('Units','pixels');
    frms = rago(2) - rago(1);

sld = uicontrol('Style', 'slider','Min',lingth+rago(1),'Max',rago(2),'Value',lingth+rago(1),'SliderStep', [2/frms, 2/frms],...
    'Position', [100 20 300 20], 'Callback', @pozz); 

% 'Position', [400 20 120 20],

    imshow(im(lingth+rago(1)).gray); hold on;
    plot(s(1).x(rago(1):lingth+rago(1)), s(1).y(rago(1):lingth+rago(1)), 'g-*');
    plot(s(2).x(rago(1):lingth+rago(1)), s(2).y(rago(1):lingth+rago(1)), 'm-*');    
    
    f.Visible = 'on'; 

function pozz(source,~)
        % For R2014a and earlier: get(source,'Value');
        cla;
        curr = round(source.Value);
        imshow(im(curr).gray);
        text(10,10, num2str(curr), 'Color', 'w');
        plot(s(1).x(curr-lingth:curr), s(1).y(curr-lingth:curr), 'g-*');
        plot(s(2).x(curr-lingth:curr), s(2).y(curr-lingth:curr), 'm-*');        
end

end


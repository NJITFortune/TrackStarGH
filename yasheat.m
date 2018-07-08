function yasheat(in)

stepsize = 5;

HotA = zeros(64,64);

    for xx = 0:stepsize:320-stepsize
        for yy = 0:stepsize:320-stepsize
            
HotA(1+(xx/stepsize), 1+(yy/stepsize)) = HotA(1+(xx/stepsize), 1+(yy/stepsize)) + length(find(in.x > xx & in.x < xx+stepsize & in.y > yy & in.y < yy+stepsize));        
        
        end
    end
figure; surf(1:64, 1:64, HotA); view(0,90);




%heat map heatmap

function vidplay(im, spadata, rango)

figure(1);

spadata(1).x = medfilt1(spadata(1).x, 5);
spadata(1).y = medfilt1(spadata(1).y, 5);
    if length(spadata) ==2
spadata(2).x = medfilt1(spadata(2).x, 5);
spadata(2).y = medfilt1(spadata(2).y, 5);
    end
for j=rango(1)+10:1:rango(2)
   
    clf;    
    imshow(im(j).gray, 'InitialMagnification', 300);
    hold on;
    plot(spadata(1).x(j-9:j), spadata(1).y(j-9:j), 'g-*');
        text(spadata(1).x(j)+5, spadata(1).y(j)+5, num2str(j), 'Color', 'g');
    if length(spadata) ==2
        plot(spadata(2).x(j-9:j), spadata(2).y(j-9:j), 'm-*');
        text(spadata(2).x(j)+5, spadata(2).y(j)+5, num2str(j), 'Color', 'm');
    end       
    text(10,10, num2str(j), 'Color', 'w');
    pause(0.01);
end

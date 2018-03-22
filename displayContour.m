function [] = displayContour(phi,I,mag,color)
%DISPLAYCONTOUR Show the contour boundary on the image

imshow(I,'InitialMagnification',mag);
hold on; [c,h] = contour(phi,[0 0],color,'Linewidth',3); hold off;
	delete(h);
    test = isequal(size(c,2),0);
	while (test==false)
        s = c(2,1);
        if ( s == (size(c,2)-1) )
            t = c;
            hold on; plot(t(1,2:end)',t(2,2:end)',color,'Linewidth',3);
            test = true;
        else
            t = c(:,2:s+1);
            hold on; plot(t(1,1:end)',t(2,1:end)',color,'Linewidth',3);
            c = c(:,s+2:end);
        end
	end    

end


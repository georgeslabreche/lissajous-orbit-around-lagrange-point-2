function animate_2d_plot(fig, X, Y, plot_title, x_label, y_label, animation_start_index, x_lim, y_lim, filename)
    %ANIMATE_2D_PLOT Animate the 2D orbit plot
    axis tight manual; % this ensures that getframe() returns a consistent size
    
    for i = animation_start_index : length(X)
        
        % Plot.
        plot(X(animation_start_index:i), Y(animation_start_index:i));
        
        % Set axis limits.
        xlim(x_lim);
        ylim(y_lim);
        
        % Set axis labels.
        xlabel(x_label); 
        ylabel(y_label);
        
        % Set title.
        title(plot_title);
        
        % Draw.
        drawnow;
        
        % Capture the plot as an image 
        frame = getframe(fig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        
        % Write to the GIF File 
        if i == animation_start_index
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.0005); 
        else 
          imwrite(imind,cm,filename,'gif', 'WriteMode','append', 'DelayTime',0.0005); 
        end 
    end
end


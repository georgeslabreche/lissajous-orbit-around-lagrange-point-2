function animate_3d_plot(fig, X, Y, Z, plot_title, x_label, y_label, z_label, animation_start_index, x_lim, y_lim, z_lim, view_az, view_el, include_plane_projections, filename)
    %ANIMATE_3D_PLOT Animate the 3D orbit plot.
    axis tight manual; % this ensures that getframe() returns a consistent size
    for i = animation_start_index : length(X)
        % Plot.
        plot3(X(animation_start_index:i),Y(animation_start_index:i),Z(animation_start_index:i), 'LineWidth', 1);
        
        % Set perspective.
        view(view_az,view_el);

        % Set axis limits.
        xlim(x_lim);
        ylim(y_lim);
        zlim(z_lim);

        % Set axis labels.
        xlabel(x_label);
        ylabel(y_label);
        zlabel(z_label);
        
        % Set title.
        title(plot_title)
        
        % Show grid.
        grid on
        
        % Plot projections on the planes.
        if include_plane_projections
            hold on

            % Project on yx-plane.
            plot3(X(animation_start_index:i),...
                Y(animation_start_index:i),...
                zeros(length(Z(animation_start_index:i)),1)+z_lim(1),...
                'r', 'LineWidth', 1, 'LineStyle', ':');

            % Project on yz-plane.
            plot3(zeros(length(X(animation_start_index:i)),1)+x_lim(1),...
                Y(animation_start_index:i),...
                Z(animation_start_index:i),...
                'r', 'LineWidth', 1, 'LineStyle', ':');

            % Project on zx-plane.
            plot3(X(animation_start_index:i),...
                zeros(length(Y(animation_start_index:i)),1)+y_lim(1),...
                Z(animation_start_index:i),...
                'r', 'LineWidth', 1, 'LineStyle', ':');
            hold off
        
        end

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


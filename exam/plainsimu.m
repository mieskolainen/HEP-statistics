% Plain boarding simulation
% ------------------------------------------------------------------------
% In the animation:  blue ~ searching for the seat or sitting
%                     red ~ busy with bag
%                   green ~ busy with seat interference
%
% INPUT:             N  =   Number of simulation runs
%            simu_type  =   'outside-in' or 'random'
%               visual  =   Animation (true / false)
%          createvideo  =   Create video (true / false)
%
% OUTPUT:           bt  =   Boarding times (N x 1)
%
% Mikael Mieskolainen, 2014

function bt = plainsimu(N, simu_type, visual, createvideo)

if (visual == true)
    close all;
end

% Boarding times
bt = zeros(N,1);

if (strcmp(simu_type, 'outside-in'))
    fprintf('Simulation type: outside-in boarding \n');
elseif (strcmp(simu_type, 'random'))
    fprintf('Simulation type: random boarding \n');
else
    fprintf('Error: Give a valid simulation type string! \n');
    return; 
end

% Cell array of seat positions, seat positions are a matrix 25 x 4

% Seat index list, we have 100 seats
seat_ind = {};
for i = 1:25
    for j = 1:4
        seat_ind{end+1} = [i,j];
    end
end

% The simulation is discrete time,
% people move forward on each k-th time step amount of d_step meters
% So our space-time discretization is the following:
d_step = 0.5;   % space discretization in meters
                % => fixes also the minimum distance between persons
k_to_s = 1;     % time conversion from steps to seconds

% e.g. d_step = 0.5 and k_to_s = 1   gives average walking speed 0.5 m/s
%  or  d_step = 1   and k_to_s = 0.5 gives average walking speed 2.0 m/s
%  or  d_step = 1   and k_to_s = 1   gives average walking speed 1.0 m/s
%  or  d_step = 1   and k_to_s = 2   gives average walking speed 0.5 m/s

% Gate to plane in discrete steps (15 in meters)
gate2plane = 15 * d_step;

if (createvideo == true)
    N = 1; % Force only one run
end

% Simulation runs
for n = 1:N

    fprintf('Simulation run %d/%d \n', n, N);
    
    if (visual == true)
        fig = figure;
        if (createvideo)
            filename = sprintf('boarding%d_%s.avi', n, simu_type);
            writerObj = VideoWriter(filename);
            open(writerObj);
        end
    end

    % Finite state machine datatype for each people, omega is our set of
    % people, 100 in this case

    % .position denotes an index in the corridor vector
    %   (minus indices denote position in the queue from gate)
    % .busy is a timer how many time steps the people will be stuck at his
    %   position (because he's putting his luggage or seat interference)
    % .seated means if the person is already sitting

    omega = cell(100,1);
    for i = 1:length(omega)
        omega{i} = struct('position', -1, 'busy', 0, 'if_bit', false, ...
                          'seated', false, 'seat', seat_ind{i}, ...
                          'xy', [0 0]);
    end

    % Outside inside boarding
    if (strcmp(simu_type, 'outside-in'))
    
        % First draw outside (columns 1 or 4) people queue numbers
        q_order_first = randperm(50);
        q_order_second = randperm(50) + 50;
        j = 1;
        for i = 1:length(omega)
            if (omega{i}.seat(2) == 1 || omega{i}.seat(2) == 4)
                omega{i}.position = - q_order_first(j) * d_step - gate2plane;
                j = j + 1;
            end
        end
        % Then draw inside (columns 2 or 3) people queue numbers
        j = 1;
        for i = 1:length(omega)
            if (omega{i}.seat(2) == 2 || omega{i}.seat(2) == 3)
                omega{i}.position = - q_order_second(j) * d_step - gate2plane;
                j = j + 1;
            end
        end

    % Random boarding
    else 
        q_order = randperm(100);
        for i = 1:length(omega)
            omega{i}.position = - q_order(i) * d_step - gate2plane;
        end
    end
    
    % Update their initial xy positions for visualization
    for i = 1:length(omega)
       omega{i}.xy = [omega{i}.position 2.5]; 
    end
    
    % Start simulation
    for k = 1:1e6 % Just some big number, doesn't matter really
        
        % First find the current order of people
        positions = zeros(length(omega),1);
        for i = 1:length(omega)
            positions(i) = omega{i}.position;
        end
        % Now sort people
        [~,sort_hash] = sort(positions, 'descend');
        
        % Loop over people in descending order (first in queue treated
        % first), this is a must because otherwise we cannot move
        % the queue right. Loop only non-seated people
        N_seated = 0;
        for i = 1:length(omega)
            N_seated = N_seated + omega{i}.seated;
        end
        if (N_seated == length(omega))
            fprintf('Boarding took k = %d steps ( t = %0.0f sec ) \n', ...
                    k-1, (k-1)*k_to_s);
            bt(n) = (k-1)*k_to_s; % Save boarding time
            break;
        end
        
        for p = 1:length(omega) - N_seated

            % Choose the highest in the queue
            i = sort_hash(p);
            
            % # A. Still in the queue
            if (omega{i}.seated == false && omega{i}.busy == 0)
                
                % Is somebody standing in front of him
                next_spot_is_free = true;
                for j = 1:length(omega)
                    if (omega{j}.position == omega{i}.position + d_step)
                        next_spot_is_free = false;
                        break;
                    end
                end
                % Can move on
                if (next_spot_is_free)
                    omega{i}.position = omega{i}.position + d_step;
                    
                    % If found the right seat row
                    if (omega{i}.position == omega{i}.seat(1))
                    
                        % Make him busy for 20 seconds on average
                        omega{i}.busy = randi([round(18/k_to_s),...
                                               round(22/k_to_s)]);
                    end
                end
                
            % # B. He is busy at his row
            elseif (omega{i}.busy > 0)
                omega{i}.busy  = omega{i}.busy - 1;
                
                % He's now ready to sit down
                if (omega{i}.busy == 0)                    
                    
                    seat_interference = false;

                    % Haven't had interference before
                    if (omega{i}.if_bit == false)
                        
                        % It's a window spot
                        if (omega{i}.seat(2) == 1 || omega{i}.seat(2) == 4)
                        
                            % Find out if somebody is sitting on that row
                            for j = 1:length(omega)
                            
                                % Is it interfering, i.e. right next to him
                                if (omega{j}.seated == true && ...
                                omega{j}.seat(1) == omega{i}.seat(1) && ...
                                abs(omega{j}.seat(2)-omega{i}.seat(2)) == 1)
                                
                                % Make him busy for average 10 seconds
                                omega{i}.busy = randi([round(8/k_to_s), ...
                                                       round(12/k_to_s)]);
                                omega{i}.if_bit = true;
                                seat_interference = true;
                                break;
                                end
                            end
                        end
                    end
                    
                    % Finally make him sit
                    if (seat_interference == false)
                        omega{i}.seated = true;
                        omega{i}.position = -inf;
                        omega{i}.xy = [omega{i}.seat(1) omega{i}.seat(2)];
                    end
                end
            
            % # C. He has seated
            elseif (omega{i}.seated == true)
                % Do nothing
            end
        end
        
        % Update people positions for non-seated (for visualization)
        for i = 1:length(omega)
            if (omega{i}.seated == false)
                omega{i}.xy = [omega{i}.position 2.5];
            end
        end
        
        % Draw animation
        if (visual == true)

            XY = zeros(length(omega),2);
            blue = [];
            red = [];
            green = [];

            for i = 1:length(omega)
                XY(i,:) = omega{i}.xy;
                if (omega{i}.busy == 0)
                    blue(end+1) = i;
                elseif (omega{i}.busy > 0 && omega{i}.if_bit == false)
                    red(end+1) = i;
                elseif (omega{i}.busy > 0 && omega{i}.if_bit == true)
                    green(end+1) = i;
                end
            end
            clf; % Clear figure, a must!
            plot(XY(blue,1), XY(blue,2), 'b.', 'MarkerSize', 12); hold on;
            plot(XY(red,1), XY(red,2), 'r.', 'MarkerSize', 12);
            plot(XY(green,1), XY(green,2), 'g.', 'MarkerSize', 12);
            axis([-1 26 0 5]);
            set(gca,'XTick', 1:2:25);
            set(gca,'YTick', 1:4);
            set(gca,'ZTick', []);
            set(gca,'YTickLabel',{'A', 'C', 'D', 'F'});
            xlabel('Seat row # (also meters)');
            %view([-20,68]);
            title(sprintf('Discrete time step k = %d ( t = %0.0f sec )', ...
                  k, k*k_to_s));
            drawnow;

            % Capture video
            if (createvideo)

                % Opengl setup
                set(gca,'nextplot','replacechildren');
                set(gcf,'Renderer','zbuffer');

                % Capture videoframe
                frame = getframe;
                writeVideo(writerObj,frame);
            end
        end
    end
    
    if (visual == true)
        close all;
        if (createvideo)
            close(writerObj);
        end
    end
end


end
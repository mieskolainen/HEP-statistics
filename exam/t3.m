%% Task 3
% Statistical Methods
% Home Exam, Spring 2014
%
% Mikael Mieskolainen

clear; close all;

% Let's do this task with Fisher's linear discriminant which is optimal in
% variance min/max sense. It minimizes the intraclass (inside one class)
% variance and maximises the interclass (between two classes) variance.

% Read in data for class A
fid = fopen('./data/exam_sample_cata.txt');
XA = fscanf(fid, '%f %f', [2 inf]);
fclose(fid);

% Read in data for class B
fid = fopen('./data/exam_sample_catb.txt');
XB = fscanf(fid, '%f %f', [2 inf]);
fclose(fid);

% Class targets, 0 == class A, 1 == class B
T = [zeros(size(XA,2),1); ones(size(XB,2),1)];

% Train the LDA, get the weight vector w
X = [XA XB];
[w,mu,C] = LDAtrain(X, T);

% The weight vector
w

%% i.)

% Plot data points in 2D
figure;
plot(XA(1,:), XA(2,:), 'ro'); hold on;
plot(XB(1,:), XB(2,:), 'xb');
xlabel('x_1'); ylabel('x_2');
legend('Class A', 'Class B');

figure;
set(gcf,'visible','off');

% Calculate the projection (test statistic of LDA)
% for each event, i.e. "the new variable" x3
x3 = X'*w;

% Plot distributions
x3_bins = -17:0.075:-2;
[nA, ~] = hist(x3(T == 0), x3_bins);
pdfA = nA / sum(nA); % Normalized to a valid pdf
bar(x3_bins, pdfA, 'r'); hold on;

[nB, ~] = hist(x3(T == 1), x3_bins);
pdfB = nB / sum(nB); % Normalized to a valid pdf
bar(x3_bins, pdfB, 'b');

xlabel('$x_3 = \langle \mathbf{x},\mathbf{w} \rangle$ (Fisher variable)', ...
       'interpreter', 'latex');
ylabel('pdf($x_3$)', 'interpreter', 'latex');
axis tight;


%% ii.)

% Calculate cumulative distributions of class A and class B in x_3 space
CDF_A = cumsum(pdfA);
CDF_B = cumsum(pdfB);

% Calculate the x3_cut necessary to get 95% rejection from class A
% (find the bin which is the closest)
[~, ind] = min(abs(CDF_A - 0.05));
x3_cut = x3_bins(ind);
fprintf('Task 3. ii.) Requirement x_3 < %0.1f, gives 95%% rejection of class A \n', ...
        x3_cut);

% Calculate the efficiency of class B with this t_cut
fprintf('Task 3. ii.) Efficiency of class B with x_{3,cut} = %0.1f is %0.2f \n', ...
        x3_cut, CDF_B(ind));

% Plot x3_cut in the previous figure (not asked for, but good to plot)
set(gcf,'visible','on');
plot([x3_cut x3_cut], [0 0.06], 'k');

legend('Class A', 'Class B', 'x_{3,cut}'); axis tight;
    
% Lets plot the events with this selection rule x3 < x3_cut
%{
figure;
X_cutA = X(:, x3 < x3_cut & T == 0);
X_cutB = X(:, x3 < x3_cut & T == 1);

plot(X_cutA(1,:), X_cutA(2,:), 'ro'); hold on;
plot(X_cutB(1,:), X_cutB(2,:), 'xb');
xlabel('x_1'); ylabel('x_2');
legend('Class A', 'Class B');
%}

% Plot data points in 2D
figure;
plot(XA(1,:), XA(2,:), 'ro'); hold on;
plot(XB(1,:), XB(2,:), 'xb');
xlabel('x_1'); ylabel('x_2');

hold on;

% Affine decision boundary plot
x = linspace(min(X(1,:)), max(X(1,:)), 2);
% <x,w> = x3_cut <=> x2 = (x3_cut - x1w1)w2
d = [x; (x3_cut - x*w(1))/w(2)];
plot(d(1,:), d(2,:), '-k'); hold on;

% Plot orthogonal projection plane (line) (not asked for, but good to plot)
theta = pi/2; % rotate 90 deg
R_xy = [cos(theta) -sin(theta); 
     sin(theta) cos(theta)];
d_ = R_xy*d;
xy_off = [1, -0.5]; % Arbitrary for visual purposes
plot(d_(1,:) + xy_off(1), d_(2,:) + xy_off(2), '--k');

legend('Class A', 'Class B', ...
       'Decision boundary: <x,w> = x_{3,cut}', ...
       'Fisher projection plane (x_3-space)', 'Location', 'SouthEast');
axis equal;


%% iii.)

% Observed event x = [x1; x2]
x = [0.50; 0.25];

% Now train the so-called Naive Gaussian-Bayes classifier with independent,
% uncorrelated component x1,x2 distributions for both classes.

% Maximum likelihood estimates of means and variances

muA = mean(XA,2);        % Mean vector of class A
muB = mean(XB,2);        % Mean vector of class B

covA = cov(XA');         % Covariance matrix of class A
covB = cov(XB');         % Covariance matrix of class B
covA = diag(diag(covA)); % Make diagonal, i.e. only variances left
covB = diag(diag(covB)); % Make components, i.e. only variances left

% Now calculate the probabilities using multivariate Gaussian distributions,
% Because we have _independent_ component marginal pdfs, we
% could also use 1-dimensional Gaussian distributions and just multiply
% those, i.e. pdf(\vec{x}|class Y) = pdf(x_1|class Y) x pdf(x_2|class Y),
% however, using multidimensional Gaussians is a more compact and general
% notation.

% Then, using Bayes formula with equal class priors,
% i.e. P(A) = P(B) = 0.5, gives us the (posteriori) probabilities:

denom = mvnpdf(x,muA,covA) + mvnpdf(x,muB,covB);
P_A_x = mvnpdf(x,muA,covA) / denom;
P_B_x = mvnpdf(x,muB,covB) / denom;

fprintf('Task 3. iii.) P(class A | \\vec{x}) = %0.3f \n', P_A_x);
fprintf('Task 3. iii.) P(class B | \\vec{x}) = %0.3f \n', P_B_x);

% So most probably the event vector x came from the class B.


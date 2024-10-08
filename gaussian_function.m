function gaussian_integration()

    mu = 0;  % Mean of the distribution
    sigma = 1;  % Standard deviation of the distribution
    a = -5;  % Lower bound of the integration
    b = 5;  % Upper bound of the integration
    n = 1000;  % Number of sub-intervals

    % Calculate areas using different methods
    riemann_area = riemann_sum(@gaussian, a, b, mu, sigma, n);
    trapezoidal_area = trapezoidal_approx(@gaussian, a, b, mu, sigma, n);
    exact_area = definite_integral(@gaussian, a, b, mu, sigma);

    % Display results
    fprintf('Estimated area using Riemann sum: %.5f\n', riemann_area);
    fprintf('Estimated area using trapezoidal approximation: %.5f\n', trapezoidal_area);
    fprintf('Exact area using integration: %.5f\n', exact_area);

    % Calculate margin of error
    riemann_sum_margin_of_error = (abs(riemann_area - exact_area) / exact_area) * 100;
    trapezoidal_approx_margin_of_error = (abs(trapezoidal_area - exact_area) / exact_area) * 100;

    fprintf('Margin of error for Riemann sum: %.5f%%\n', riemann_sum_margin_of_error);
    fprintf('Margin of error for trapezoidal approximation: %.5f%%\n', trapezoidal_approx_margin_of_error);

    % Plot Gaussian function
    plot_gaussian(@gaussian, mu, sigma, a, b);
end

function y = gaussian(x, mu, sigma)
    % Gaussian function
    y = 1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((x - mu) / sigma).^2);
end

function area = riemann_sum(gaussian_func, a, b, mu, sigma, n)
    % Riemann sum approximation
    dx = (b - a) / n;  % Width of each sub-interval
    area = 0;
    for i = 0:n-1
        x = a + i * dx;  % Left endpoint of each sub-interval
        area = area + gaussian_func(x, mu, sigma) * dx;
    end
end

function area = trapezoidal_approx(gaussian_func, a, b, mu, sigma, n)
    % Trapezoidal approximation
    x = linspace(a, b, n);
    y = gaussian_func(x, mu, sigma);
    area = trapz(x, y);
end

function area = definite_integral(gaussian_func, a, b, mu, sigma)
    % Exact integral using integral function
    area = integral(@(x) gaussian_func(x, mu, sigma), a, b);
end

function plot_gaussian(gaussian_func, mu, sigma, a, b)
    % Plot the Gaussian function and shaded area
    x = linspace(-5, 5, 1000);
    y = gaussian_func(x, mu, sigma);

    figure('Position', [100, 100, 800, 400]);
    plot(x, y, 'b', 'LineWidth', 1.5);
    hold on;

    % Fill the area under the curve from a to b
    fill_x = linspace(a, b, 100);
    fill_y = gaussian_func(fill_x, mu, sigma);
    fill([fill_x, fliplr(fill_x)], [fill_y, zeros(size(fill_y))], 'skyblue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    title('Gaussian Function and Area Estimations');
    xlabel('x');
    ylabel('y');
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid on;
    legend('Gaussian Function', sprintf('Area from %.2f to %.2f', a, b));
    hold off;
end
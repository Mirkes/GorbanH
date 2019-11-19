function DrawTriangle(eq, c0)
% DrawTriangle draws triange and draw lines of partial equilibrium for all
% reactions.
%
% Inputs:
%   eq is row vector of equilibrium.
%   c0 is row vector of point to show projections

    figure('Color', 'w');
    P = [0, 0; 1, 0; 0.5, sind(60)];
    T = [1, 2, 3];
    TR = triangulation(T,P);
    triplot(TR, 'k', 'Linewidth', 1.5);
    hold on;

    % Define shift in points
    ps = 5;
    % define horizontal size of arrow
    has = 0.1;

    % Decoration
    fontSize = 20;
    axis off;
    axis equal;
    text(P(1, 1), P(1, 2), '\it A_1', "FontSize", fontSize, 'HorizontalAlignment', 'right');
    text(P(2, 1), P(2, 2), '\it A_2', "FontSize", fontSize);
    text(P(3, 1), P(3, 2) + 0.02, '\it A_3', "FontSize", fontSize,...
        'HorizontalAlignment', 'center');

    % Define balance and renormalise equilibrium for unit sum
    eq = eq(:)';
    b = sum(eq);
    eq = eq / b;

    % Define line for A1<->A2
    c3 = linspace(0, 1)';
    c = [eq(1)/(eq(1)+eq(2))*(1-c3), eq(2)/(eq(1)+eq(2))*(1-c3), c3];
    conv = barycentricToCartesian(TR, ones(length(c3),1), c);
    plot(conv(:, 1), conv(:, 2), '-m', 'Linewidth', 1.5);
    quiver(conv(ps, 1) + has, conv(ps, 2) - has, -has, has, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
    text(conv(ps, 1) + has, conv(ps, 2) - has, '\it A_1\leftrightarrow A_2', "FontSize", fontSize);

    % Define line for A2<->A3
    c1 = linspace(0, 1)';
    c = [c1, eq(2)/(eq(3)+eq(2))*(1-c1), eq(3)/(eq(3)+eq(2))*(1-c1)];
    conv = barycentricToCartesian(TR, ones(length(c1),1), c);
    plot(conv(:, 1), conv(:, 2), '-m', 'Linewidth', 1.5);
    quiver(conv(ps, 1) + has - 0.02, conv(ps, 2) - 0.02, -has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
    text(conv(ps, 1) + has - 0.02, conv(ps, 2) - 0.02, '\it A_2\leftrightarrow A_3', "FontSize", fontSize);

    % Define line for A3<->A1
    c2 = linspace(0, 1)';
    c = [eq(1)/(eq(3)+eq(1))*(1-c2), c2, eq(3)/(eq(3)+eq(1))*(1-c2)];
    conv = barycentricToCartesian(TR, ones(length(c2),1), c);
    plot(conv(:, 1), conv(:, 2), '-m', 'Linewidth', 1.5);
    if eq(2) < 0.1
        ps = ps - 3;
    end

    quiver(conv(ps, 1) - has + 0.02, conv(ps, 2) - 0.02, has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
    text(conv(ps, 1) - has + 0.01, conv(ps, 2) - 0.02, '\it A_3\leftrightarrow A_1',...
        "FontSize", fontSize, 'HorizontalAlignment', 'right');
    
    % Draw equilibrium
    conv = barycentricToCartesian(TR, 1, eq);
    plot(conv(1, 1), conv(1, 2), '*r', 'Linewidth', 1.5);
    if eq(2) < 0.1
        quiver(conv(1, 1) + has, conv(1, 2), - has, 0, 'r', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
        text(conv(1, 1) + has, conv(1, 2) + 0.01, '\it c^{eq}', "FontSize", fontSize);
        % Draw point c0
        convC = barycentricToCartesian(TR, 1, c0);
        plot(convC(1, 1), convC(1, 2), '*b', 'Linewidth', 1.5);
        text(convC(1, 1) + 0.01, convC(1, 2), '\it c', "FontSize", fontSize);
    else
        quiver(conv(1, 1) - has, conv(1, 2), has, 0, 'r', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
        text(conv(1, 1) - has, conv(1, 2) + 0.01, '\it c^{eq}', "FontSize", fontSize, 'HorizontalAlignment', 'right');
        % Draw point c0
        convC = barycentricToCartesian(TR, 1, c0);
        plot(convC(1, 1), convC(1, 2), '*b', 'Linewidth', 1.5);
        text(convC(1, 1) - 0.01, convC(1, 2), '\it c', "FontSize", fontSize, 'HorizontalAlignment', 'right');
    end

    % Draw arrow to the first partial quilibrium
    % Define line for A1<->A2
    c3 = c0(3);
    c = [eq(1)/(eq(1)+eq(2))*(1-c3), eq(2)/(eq(1)+eq(2))*(1-c3), c3];
    conv = barycentricToCartesian(TR, ones(length(c3),1), c);
    shiftX = conv(1, 1) - convC(1, 1);
    shiftY = conv(1, 2) - convC(1, 2);
    len = sqrt(shiftX ^ 2 + shiftY ^ 2);
    q = quiver(convC(1, 1), convC(1, 2), shiftX, shiftY, ':b', ...
        'LineWidth', 2, 'AutoScale', 'off');
    q.Head.LineStyle = 'solid';
    q.MaxHeadSize = 0.9 * 0.09 / len;

    % Define line for A2<->A3
    c1 = c0(1);
    c = [c1, eq(2)/(eq(3)+eq(2))*(1-c1), eq(3)/(eq(3)+eq(2))*(1-c1)];
    conv = barycentricToCartesian(TR, ones(length(c1),1), c);
    shiftX = conv(1, 1) - convC(1, 1);
    shiftY = conv(1, 2) - convC(1, 2);
    len = sqrt(shiftX ^ 2 + shiftY ^ 2);
    q = quiver(convC(1, 1), convC(1, 2), shiftX, shiftY, ':b', ...
        'LineWidth', 2, 'AutoScale', 'off');
    q.Head.LineStyle = 'solid';
    q.MaxHeadSize = 0.9 * 0.09 / len;

    % Define line for A3<->A1
    c2 = c0(2);
    c = [eq(1)/(eq(3)+eq(1))*(1-c2), c2, eq(3)/(eq(3)+eq(1))*(1-c2)];
    conv = barycentricToCartesian(TR, ones(length(c2),1), c);
    shiftX = conv(1, 1) - convC(1, 1);
    shiftY = conv(1, 2) - convC(1, 2);
    len = sqrt(shiftX ^ 2 + shiftY ^ 2);
    q = quiver(convC(1, 1), convC(1, 2), shiftX, shiftY, ':b', ...
        'LineWidth', 2, 'AutoScale', 'off');
    q.Head.LineStyle = 'solid';
    q.MaxHeadSize = 0.9 * 0.09 / len;
end
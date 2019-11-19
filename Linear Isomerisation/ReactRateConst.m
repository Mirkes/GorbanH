% This script depicts Reaction Rate Polyhedron.
vertices = [0, 0, 0.333333333;...
            0, 0, 0.5;...
            0, 0.5, 0;...
            0.5, 0, 0;...
            0.333333333, 0.333333333, 0];
labels = {' 1/3 (0, 0, 1/3)',...
          ' -1 (0, 0, 1/2)',...
          ' -1 (0, 1/2, 0)',...
          ' -1 (1/2, 0, 0)',...
          ' 1/3 (1/3, 1/3, 0)'};
edges = [1, 2;
         1, 3;
         1, 4;
         2, 3;
         2, 4;
         2, 5;
         3, 4;
         3, 5;
         4, 5];

figure;
% Draw vertices
hold on;
for k = 1:size(vertices, 1)
    plot3(vertices(k, 1), vertices(k, 2), vertices(k, 3), 'bo', 'MarkerFaceColor', 'b');
    text(vertices(k, 1), vertices(k, 2), vertices(k, 3), labels{k});
end
% Draw edges
for k = 1:size(edges, 1)
    plot3([vertices(edges(k, 1), 1), vertices(edges(k, 2), 1)],...
          [vertices(edges(k, 1), 2), vertices(edges(k, 2), 2)],...
          [vertices(edges(k, 1), 3), vertices(edges(k, 2), 3)], 'b-');
end
view(77, 32);

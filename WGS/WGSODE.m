function dy = WGSODE(~, y, kp, km, b)
    % Transform to complete list of substances
    % y contains H2O and CO only
    c = [y(1), b(1) - y(1), y(2), b(2) - y(2), 0,  0];
    %    H2O    H2          CO    CO2          red Ox
    c(6) = b(3) - c(1) - c(3) - 2 * c(4); % red
    c(5) = b(4) - c(6);                   % Ox
    dy = [-kp(1) * c(1) * c(5) + km(1) * c(2) * c(6);...
        -kp(2) * c(3) * c(6) + km(2) * c(4) * c(5)];
end
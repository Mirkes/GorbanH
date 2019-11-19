function dy = modelODE(~, y, kp, km)
    % add one additional coordinate
    c = [y(1), 1 - sum(y), y(2)];
    dy = [-kp(1) * c(1) + km(1) * c(2) - 2 * kp(3) * c(1) ^ 2 + 2 * km(3) * c(2) * c(3);...
        kp(2) * c(2) - km(2) * c(3) + kp(3) * c(1) ^ 2 - km(3) * c(2) * c(3)];
end
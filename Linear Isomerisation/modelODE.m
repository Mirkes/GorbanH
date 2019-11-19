function dy = modelODE(~, y, kp, km)
    % add one additional coordinate
    c = [y(1), y(2), 1 - sum(y)];
    dy = [-kp(1) * c(1) + km(1) * c(2) + kp(3) * c(3) - km(3) * c(1);...
           kp(1) * c(1) - km(1) * c(2) - kp(2) * c(2) + km(2) * c(3)];
end
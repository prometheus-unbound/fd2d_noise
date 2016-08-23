
function [GY] = gaussblur2d(Y, x1, x2, sigma)


    [n1, n2] = size(Y);
    Mtmp = zeros(size(Y));
    GY = zeros(size(Y));
    w1 = zeros(1, n1);
    w2 = zeros(1, n2);

    sqrt14 = sqrt(14.);

    
    for i = 1:n1
        Mtmp(i,:) = 0;
        for j = 1:n1
            dist = abs(x1(i) - x1(j));
            if (dist < sqrt14 * sigma(1))
                w = exp(- (dist / sigma(1)) ^ 2 / 2.0);
                w1(j) = w1(j) + w;
                Mtmp(i,:) = Mtmp(i,:) + w * Y(j,:);
            end
        end
    end
    
    
    for i = 1:n2
        for j = 1:n2
            dist = abs(x2(i) - x2(j));
            if (dist < sqrt14 * sigma(2))
                w = exp(- (dist / sigma(2)) ^ 2 / 2.0);
                w2(j) = w2(j) + w;
                GY(:, i) = GY(:, i) + w * Mtmp(:, j);
            end
        end
    end
    
    
    Mtmp = w1' * w2;
    GY = GY ./ Mtmp;


end


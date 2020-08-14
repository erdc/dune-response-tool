function mouseMove (~, ~, handle)
C = handle.CurrentPoint;
title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
drawnow
end
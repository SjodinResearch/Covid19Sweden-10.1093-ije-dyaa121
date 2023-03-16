function expo = exports_of_X(T, X, alpha)

    expo = alpha * [sum(T' .* X(1, :)',2)'; sum(T' .* X(2, :)',2)'; sum(T' .* X(3, :)',2)'];

end

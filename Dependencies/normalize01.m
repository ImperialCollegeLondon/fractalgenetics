function g = normalize01(f)

% Normalize to the range of [0, 1]
fmin = min(f(:));
fmax = max(f(:));

TOLERANCE = 0.1;

if (fmax > TOLERANCE)
  g = (f - fmin)/(fmax - fmin);  
else
  g(:) = 0.0;
end

end
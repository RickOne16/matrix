function [X, alpha, condest] = sqrtm2(A)
##Copyright (C) 2015-2016 Marco Caliari
##Copyright (C) 2015-2016 Mudit Sharma
##Copyright (C) 2011-2015 Awad H. Al-Mohy
##Copyright (C) 2011-2015 Nicholas J. Higham

##This file is part of Octave.

##Octave is free software; you can redistribute it and/or modify it
##under the terms of the GNU General Public License as published by the
##Free Software Foundation; either version 3 of the License, or (at your
##option) any later version.

##Octave is distributed in the hope that it will be useful, but WITHOUT
##ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
##for more details.

##You should have received a copy of the GNU General Public License
##along with Octave; see the file COPYING.  If not, see
##<http://www.gnu.org/licenses/>.

##This file incorporates work covered by the following copyright and
##permission notice:

##Redistribution and use in source and binary forms, with or without
##modification, are permitted provided that the following conditions are met:
##* Redistributions of source code must retain the above copyright notice, this
##  list of conditions and the following disclaimer.
##* Redistributions in binary form must reproduce the above copyright notice,
##  this list of conditions and the following disclaimer in the documentation
##  and/or other materials provided with the distribution.

##SQRTM2    Matrix square root.
##         X = SQRTM2(A) is a square root of the matrix A (A = X*X).
##         X is the unique square root for which every eigenvalue has
##         nonnegative real part (the principal square root).
##         If A is real with a negative real eigenvalue then a complex
##         result will be produced.
##         If A is singular then A may not have a square root.
##         [X, ALPHA, CONDEST] = SQRTM2(A) returns a stability factor ALPHA and
##         an estimate CONDEST for the matrix square root condition number of X.
##         The residual NORM(A-X^2,'fro')/NORM(A,'fro') is bounded by
##         (n+1)*ALPHA*EPS and the Frobenius norm relative error in X is
##         bounded approximately by N*ALPHA*CONDEST*EPS, where N = MAX(SIZE(A)).

##         References:
##         N. J. Higham, Computing real square roots of a real
##            matrix, Linear Algebra and Appl., 88/89 (1987), pp. 405-430.
##         A. Bjorck and S. Hammarling, A Schur method for the square root of a
##            matrix, Linear Algebra and Appl., 52/53 (1983), pp. 127-140.

   n = max (size(A));
  [Q, T] = schur (A);        % T is real/complex according to A.
  [Q, T] = rsf2csf (Q, T);   % T is now complex Schur form.

## Compute upper triangular square root R of T, a column at a time.

  nzeig = length(find(diag(T)==0));
  if (nzeig)
    warning('Matrix is singular and may not have a square root.')
  endif

  R = zeros(n);
  for j=1:n
    R(j,j) = sqrt(T(j,j));
    for i=j-1:-1:1
      s = 0;
    for k=i+1:j-1
      s = s + R(i,k)*R(k,j);
    endfor
    if (R(i,i) + R(j,j) ~= 0)
      R(i,j) = (T(i,j) - s)/(R(i,i) + R(j,j));
    endif
    endfor
  endfor

  if (nargout > 1)
    alpha = norm(R,'fro')^2 / norm(T,'fro'); 
  endif

  X = Q*R*Q';

  if (nargout > 2)
    if (nzeig)
      condest = inf;
    else

      ##Power method to get condition number estimate.
      warns = warning;
      warning('off');

      tol = 1e-2;
      x = ones(n^2,1);    ##Starting vector.
      cnt = 1;
      e = 1;
      e0 = 0;
      while (abs(e-e0) > tol*e & cnt <= 6)
         x = x/norm(x);
         x0 = x;
         e0 = e;
         Sx = solve(R, x);
         x = solve(R, Sx, 'T');
         e = sqrt(real(x0'*x));  ## sqrt of Rayleigh quotient.
         fprintf('cnt = %2.0f, e = %9.4e\n', cnt, e)
         cnt = cnt+1;
      endwhile

      condest = e*norm(A,'fro')/norm(X,'fro');
      warning(warns);

    endif

  endif
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = solve(R, b, tran)
##SOLVE       Solves Kronecker system.
##           x = SOLVE(R, b, TRAN) solves
##                 A*x = b  if TRAN = '',
##                A'*x = b  if TRAN = 'T',
##           where A = KRON(EYE,R) + KRON(TRANSPOSE(R),EYE).
##           Default: TRAN = ''.

  if (nargin < 3) 
    tran = ''; 
  endif

  n = max(size(R));
  x = zeros(n^2,1);

  I = eye(n);

  if (isempty(tran))

   # Forward substitution.
   for i = 1:n
       temp = b(n*(i-1)+1:n*i);
       for j = 1:i-1
           temp = temp - R(j,i)*x(n*(j-1)+1:n*j);
       endfor
       x(n*(i-1)+1:n*i) = (R + R(i,i)*I) \ temp;
   endfor

  elseif strcmp(tran,'T')

   ## Back substitution.
   for i = n:-1:1
       temp = b(n*(i-1)+1:n*i);
       for j = i+1:n
           temp = temp - conj(R(i,j))*x(n*(j-1)+1:n*j);
       endfor
       x(n*(i-1)+1:n*i) = (R' + conj(R(i,i))*I) \ temp;
   endfor

  endif

return
endfunction
#####################################################

%!assert( sqrtm2([3 4;5 8]) ,[1.2910 1.0328;1.2910 2.5820], 6e-5);

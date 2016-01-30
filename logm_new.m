function [X, s, m] = logm_new (A, maxsqrt)

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

##Copyright (c) 2011, Awad H. Al-Mohy and Nicholas J. Higham
##All rights reserved.

##Redistribution and use in source and binary forms, with or without
##modification, are permitted provided that the following conditions are met:
##* Redistributions of source code must retain the above copyright notice, this
##  list of conditions and the following disclaimer.
##* Redistributions in binary form must reproduce the above copyright notice,
##  this list of conditions and the following disclaimer in the documentation
##  and/or other materials provided with the distribution.

##THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
##AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
##IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
##DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
##FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
##DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
##SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
##CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
##OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
##OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
##LOGM_NEW  Matrix logarithm by Schur-based inverse scaling and squaring.
##  X = LOGM_NEW(A,MAXSQRT) computes the logarithm of A, for a matrix
##  with no nonpositive real eigenvalues, using the inverse scaling and
##  squaring method with Pade approximation and a Schur decomposition.
##  [X, S, M] = LOGM_NEW(A) returns the number S of square roots
##  computed and the degree M of the Pade approximant.
##  At most MAXSQRT matrix square roots are  computed.
    
##  This code is intended for double precision.

##  Reference: A. H. Al-Mohy and N. J. Higham, Improved Inverse Scaling 
##  and Squaring Algorithms for the Matrix Logarithm, MIMS EPrint 2011.83,
##  The University of Manchester, October 2011.  
##  Name of corresponding algorithm in that paper: 
##  Algorithm 4.1/iss_schur_new.

##  Awad H. Al-Mohy and Nicholas J. Higham, October 19, 2011.

  if (nargin < 2 || isempty (maxsqrt))
    maxsqrt = 100; 
  endif

  xvals = [1.586970738772063e-005
           2.313807884242979e-003
           1.938179313533253e-002
           6.209171588994762e-002
           1.276404810806775e-001
           2.060962623452836e-001
           2.879093714241194e-001];

  n = length (A);
  mmax = 7; 
  foundm = false;
## First form complex Schur form (if A not already upper triangular).
  if (isequal (A, triu (A)))
    T = A; 
    Q = eye (n);
  else
    [Q,T] = schur (A, 'complex');
  endif
  T0 = T;
  if (any ( imag (diag (T)) == 0 & real (diag (T)) <= 0 ))
    warning ('A must not have nonpositive real eigenvalues!')
  endif
  p = 0;
  s0 = opt_cost (diag (T)); 
  s = s0;
  for (k = 1:min (s, maxsqrt))
    T = sqrtm_tri (T);
  endfor

  d2 = normAm (T-eye (n), 2)^(1/2); 
  d3 = normAm (T-eye (n), 3)^(1/3);
  alpha2 = max (d2, d3);
  if (alpha2 <= xvals (2))
    m = find (alpha2 <= xvals (1:2), 1); 
    foundm = true;
  endif

  while (~foundm)
    more = 0; ## more square roots
    if (s > s0)
      d3 = normAm (T-eye (n), 3)^(1/3); 
    endif
    d4 = normAm (T-eye (n), 4)^(1/4);
    alpha3 = max (d3, d4);
    if (alpha3 <= xvals (mmax))
      j = find ( alpha3 <= xvals (3:mmax), 1) + 2;
      if (j <= 6)
        m = j; 
        break
      else
        if (alpha3/2 <= xvals (5) && p < 2)
          more = 1; 
          p = p+1;
        endif
      endif
    endif
    if (~more)
      d5 = normAm (T-eye (n), 5)^(1/5);
      alpha4 = max (d4, d5);
      eta = min (alpha3, alpha4);
      if (eta <= xvals (mmax))
        m = find (eta <= xvals (6:mmax), 1) + 5;
        break
      endif
    endif
    if (s == maxsqrt)
      m = mmax; 
      break
    endif
    T = sqrtm_tri(T); s = s + 1;
  endwhile

## Compute accurate superdiagonal of T^(1/2^s).
  for (k = 1:n-1)
    ## Tkk = T0(k:k+1,k:k+1);
    ## Tkk = powerm2by2(Tkk,1/2^s);
    ## T(k:k+1,k:k+1) = Tkk;
    T (k:k+1, k:k+1) = powerm2by2 (T0 (k:k+1, k:k+1), 1/2^s);
  endfor

## Compute accurate diagonal of T^(1/2^s) - I.
  d = sqrt_power_1 (diag (T0), s);
## T = triu(T,1) + diag(d);
  T(1:n+1:end) = d;
  Y = logm_pf (T, m);
  X = 2^s*Y;

## Compute accurate diagonal and superdiagonal of log(T).
  for (j = 1:n-1)
    X(j:j+1, j:j+1) = logm2by2 (T0 (j:j+1, j:j+1));
  endfor

  X = Q*X*Q';
  
## Nested functions

#######################
function s = opt_cost (d)
    ## for i = 0:double(intmax)
    ##     if max( abs(d-1)) <= xvals(mmax), s = i; return, end
    ##     d = sqrt(d);
    ## end
  s = 0;
  while (norm(d-1,inf) > xvals(mmax))
    d = sqrt(d); 
    s = s+1;
  endwhile
endfunction

##################################
function S = logm_pf (A, m)
##LOGM_PF   Pade approximation to matrix log by partial fraction expansion.
##   LOGM_PF(A,m) is an [m/m] Pade approximant to LOG(EYE(SIZE(A))+A).

  [nodes,wts] = gauss_legendre (m);
## Convert from [-1,1] to [0,1].
  nodes = (nodes + 1)/2;
  wts = wts/2;

  S = zeros (n);
  for (j=1:m)
    S = S + wts (j)*(A/(eye (n) + nodes (j)*A));
  endfor
endfunction

endfunction

## Subfunctions

####################################
function R = sqrtm_tri (T)
## Compute upper triangular square root R of T, a column at a time.
  n = length (T);
  R = zeros (n);
  for (j=1:n)
    R (j, j) = sqrt (T (j, j));
    for (i=j-1:-1:1)
      R (i, j) = (T (i, j) - R (i, i+1:j-1)*R (i+1:j-1, j))/(R (i, i) + R (j, j));
    endfor
  endfor
endfunction
######################################
function X = powerm2by2 (A, p)
##POWERM2BY2    Power of 2-by-2 upper triangular matrix.
##   POWERM2BY2(A,p) is the pth power of the 2-by-2 upper
##   triangular matrix A, where p is an arbitrary real number.

  a1 = A (1, 1);
  a2 = A (2, 2);
  a1p = a1^p;
  a2p = a2^p;

  loga1 = log (a1);
  loga2 = log (a2);

  X = diag ([a1p a2p]);

  if (a1 == a2)
    X (1, 2) = p*A (1, 2)*a1^(p-1);
  elseif (abs(a1) < 0.5*abs(a2) || abs(a2) < 0.5*abs(a1))
    X (1, 2) =  A (1, 2) * (a2p - a1p) / (a2 - a1);
  else ## Close eigenvalues.
   w = atanh ((a2-a1)/(a2+a1)) + 1i*pi*unwinding (loga2-loga1);
   dd = 2 * exp (p*(loga1+loga2)/2) * sinh (p*w) / (a2-a1);
   X (1, 2) = A (1, 2)*dd;

  endif
endfunction

####################################
function r = sqrt_power_1 (a, n)
##SQRT_POWER_1    Accurate computation of a^(2^n)-1.
##  SQRT_POWER_1(A,N) computes a^(2^n)-1 accurately.

## A. H. Al-Mohy.  A more accurate Briggs method for the logarithm.
## Numer. Algorithms, DOI: 10.1007/s11075-011-9496-z.

  if (n == 0)
    r = a-1; 
    return
  endif
  n0 = n;
  if (angle(a) >= pi/2)
    a = sqrt(a); 
    n0 = n-1;
  endif
  z0 = a - 1;
  a = sqrt (a);
  r = 1 + a;
  for (i=1:n0-1)
    a = sqrt (a);
    r = r.*(1+a);
  endfor
  r = z0./r;
endfunction

############################
function X = logm2by2(A)
%LOGM2BY2    Logarithm of 2-by-2 upper triangular matrix.
%   LOGM2BY2(A) is the logarithm of the 2-by-2 upper triangular matrix A.

  a1 = A (1, 1);
  a2 = A (2, 2);

  loga1 = log (a1);
  loga2 = log (a2);
  X = diag ([loga1 loga2]);

  if (a1 == a2)
    X (1, 2) = A (1, 2)/a1;

  elseif (abs (a1) < 0.5*abs (a2) || abs (a2) < 0.5*abs (a1))
    X (1, 2) =  A (1, 2) * (loga2 - loga1) / (a2 - a1);

  else ## Close eigenvalues.
   dd = (2*atanh ((a2-a1)/(a2+a1)) + 2*pi*1i*unwinding (loga2-loga1)) / (a2-a1);
   X (1, 2) = A (1, 2)*dd;

  endif
endfunction

###################################
function [x,w] = gauss_legendre (n)
##GAUSS_LEGENDRE  Nodes and weights for Gauss-Legendre quadrature.
##   [X,W] = GAUSS_LEGENDRE(N) computes the nodes X and weights W
##   for N-point Gauss-Legendre quadrature.

## Reference:
## G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
## rules, Math. Comp., 23(106):221-230, 1969.

  i = 1:n-1;
  v = i./sqrt ((2*i).^2-1);
  [V,D] = eig (diag (v, -1)+diag (v, 1));
  x = diag (D);
  w = 2*(V(1, :)'.^2);
endfunction

##################################
function u = unwinding (z, k)
##UNWINDING    Unwinding number.
##   UNWINDING(Z,K) is the K'th derivative of the
##   unwinding number of the complex number Z.
##   Default: k = 0.

  if (nargin == 1 || k == 0)
    u = ceil ((imag (z) - pi)/(2*pi));
  else
    u = 0;
  endif
endfunction

function [c,mv] = normAm (A, m)
##NORMAM   Estimate of 1-norm of power of matrix.
##NORMAM(A,m) estimates norm(A^m,1).
##If A has nonnegative elements the estimate is exact.
##[C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
##matrix-vector products computed involving A or A^*.

##Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
##Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
##970-989, 2009.

##Awad H. Al-Mohy and Nicholas J. Higham, April 19, 2010.

  n = length (A);
  if (isequal (A, abs (A)))
    e = ones (n, 1);
    for j=1:m         ## for positive matrices only
      e = A'*e;
    endfor
    c = norm (e, inf);
    mv = m;
  else
    afun_power = @(flag, X) afun_power_A (flag, X, A, m);
    [c, v, w, it] = normest1 (afun_power);
    mv = it (2)*2*m; ## Since t = 2.
  endif

function [est, v, w, k] = normest1 (A, t, x0)
##est = mynormest1(A)
##est = mynormest1(Afun)
##based on Algorithm 2.4 of 
##A BLOCK ALGORITHM FOR MATRIX 1-NORM ESTIMATION,
##WITH AN APPLICATION TO 1-NORM PSEUDOSPECTRA
##NICHOLAS J. HIGHAM AND FRANCOISE TISSEUR
  if (nargin <= 1)
    t = 2;
  endif
  if (isnumeric (A))
## A is a matrix
    n = length (A);
    Afun = @(x) A*x;
    A1fun = @(x) A'*x;
    realm = isreal (A);
  else
    n = A ('dim', []);
    realm = A ('real', []);
    Afun = @(x) A ('notransp', x);
    A1fun = @(x) A ('transp', x);
  endif
  
  if (nargin <= 2)
    X = 2*rand (n, t)-1;
    X = X/diag (sum (abs (X)));
  endif
  
  ind_hist = [];
  estold = 0;
  ind = zeros (n, 1);
  S = zeros (n, t);
  k = 1;
  k2 = 0;
  conv = false;
  
  while (~conv)
    Y = Afun (X);
    k2 = k2+1;
    [est, ind_best] = max (sum (abs (Y)));
    if ((est > estold) || (k == 2))
      w = Y (:, ind_best);
    endif
    if ((est <= estold) && (k >= 2))
      est = estold;
      break;
    endif
    estold = est;
    Sold = S;
    S = sign (Y);
    S (S==0) = 1;
    possiblebreak = false;
    if (realm)
    ## test parallel
      p = false (1, t);
      p (1) = true;
      i = 1;
      possiblebreak = true;
      while (i <= t)
        for j = i+1:t
	        if (abs (S (:, i)'*S (:, j)) == n) 
	        ## parallel column
            p (j) = true;
	        endif
        endfor
        if (possiblebreak && all (p))
 	      ## check if all the parallel columns of S are parallel to Sold
	        possiblebreak = false;
	        
          for i = 1:t
            if (abs (S (:,1)'*Sold (:, i)) == n) %'      
              possiblebreak = true;
              conv = true;
              break ## for loop
            endif      
	        endfor
        
        else
	        possiblebreak = false;
        endif
        if (possiblebreak)
	        break ## while loop
        endif
        if (any (p (2:t)))
	      ## some parallel in S, unless t=1
	        S (:, p) = 2*round (rand (n, sum (p)))-1;
	        p = false (1, t);
        else
	        i = i+1;
        endif
      endwhile
    endif
    if (~possiblebreak)
      Z = A1fun (S);
      k2 = k2+1;
      h = max (abs (Z), [], 2);
      ind = 1:n;
      if ((max (h) == h (ind_best)) && (k >= 2))
        break;
      endif
      [h, ind] = sort (h, 'descend');
      if (t > 1)
        if (all (ismember (ind (1:t), ind_hist)))
          break
        endif
        ind = ind (~ismember (ind (1:n), ind_hist));
        ind = ind (1:t);
      endif
      for j = 1:t
        X (:,j) = [zeros(ind (j)-1, 1);1;zeros(n-ind (j), 1)];
      endfor
      k = k+1;
    
    endif
  
  endwhile
  
  k (2) = k2;
  k = k (:);
  v = zeros (n, 1);    
  v (ind_best) = 1;

endfunction

endfunction

function Z = afun_power_A (flag, X, A, m)
  ##AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.
	   
  if (isequal (flag, 'dim'))
    Z = length (A);
  elseif isequal (flag, 'real')
    Z = isreal (A);
  else
    [p,q] = size (X);
    if (p ~= length (A))
      error ('Dimension mismatch');
    endif
    if (isequal (flag, 'notransp'))
      for i = 1:m
      X = A*X; 
      endfor
    elseif isequal (flag, 'transp')
      for i = 1:m
        X = A'*X; 
      endfor
    endif
    
    Z = X;    
  endif
endfunction
%!assert (logm_new ([1 -1 -1;0 1 -1; 0 0 1]), [0 -1 -1.5; 0 0 -1; 0 0 0], 1e-5)
%!assert (logm_new (10), log (10),2*eps)
%!assert (full (logm_new (eye (3))), logm_new (full (eye (3))))
%!assert (full (logm_new (10*eye (3))), logm_new (full (10*eye (3))), 8*eps)
%!assert (logm_new (expm ([0 1i; -1i 0])), [0 1i; -1i 0], 10 * eps)

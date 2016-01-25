function [F, n_swaps, n_calls, terms, ind, T] = funm (A, fun, delta, tol, prnt, m)

##Copyright (C) 2015 Marco Caliari
##Copyright (C) 2015 Mudit Sharma
##Copyright (C) 2010 Nicholas J. Higham
##Copyright (C) 2010 Awad H. Al-Mohy
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

##Part of this code was originally distributed in the accompanying files
##of "2009.9: Awad H. Al-Mohy and Nicholas J. Higham (2009) A New Scaling
##and Squaring Algorithm for the Matrix Exponential. SIAM Journal On Matrix
##Analysis and Applications., 31 (3). pp. 970-989. ISSN 1095-7162"

##Copyright (c) Awad H. Al-Mohy and Nicholas J. Higham, April 20, 2010.

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
##OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##FUNM  Evaluate general matrix function.
##       F = FUNM(A,FUN), for a square matrix argument A, evaluates the
##       function FUN at the square matrix A.
##       FUN(X,K) must return the K'th derivative of
##       the function represented by FUN evaluated at the vector X.
##       The MATLAB functions COS, SIN, EXP, LOG can be passed as FUN,
##       i.e. FUNM(A,@COS), FUNM(A,@SIN), FUNM(@EXP), FUNM(A,@LOG).
##       For matrix square roots use SQRTM(A) instead.
##       For matrix exponentials, either of EXPM(A) and FUNM(A,@EXPM)
##       may be the faster or the more accurate, depending on A.
##
##       F = FUNM(A,FUN,DELTA,TOL,PRNT,M) specifies a tolerance
##       DELTA used in determining the blocking (default 0.1),
##       and a tolerance TOL used in a convergence test for evaluating the
##       Taylor series (default EPS).
##       If PRNT is nonzero then information describing the
##       behaviour of the algorithm is printed.
##       M, if supplied, defines a blocking.
##
##       [F,N_SWAPS,N_CALLS,TERMS,IND,T] = FUNM(A,FUN,...) also returns
##       N_SWAPS:  the total number of swaps in the Schur re-ordering.
##       N_CALLS:  the total number of calls to ZTREXC for the re-ordering.
##       TERMS(I): the number of Taylor series terms used when evaluating
##                 the I'th atomic triangular block.
##       IND:      a cell array specifying the blocking: the (I,J) block of
##                 the re-ordered Schur factor T is T(IND{I},IND{j}),
##       T:        the re-ordered Schur form.

  if (isequal (fun, @cos) || isequal (fun, 'cos'))
    fun = @fun_cos; 
  endif
  if (isequal (fun, @sin) || isequal (fun, 'sin'))
    fun = @fun_sin; 
  endif
  if (isequal (fun, @exp) || isequal (fun, 'exp'))
    fun = @fun_exp; 
  endif

  if (nargin < 3 || isempty (delta))
    delta = 0.1; 
  endif
  if (nargin < 4 || isempty (tol))
    tol = eps; 
  endif
  if (nargin < 5 || isempty (prnt))
    prnt = 0;  
  endif
  if (nargin < 6)
    m = []; 
  endif

  n = length (A);

## First form complex Schur form (if A not already upper triangular).
  if (isequal (A,triu (A)))
    T = A; 
    U = eye (n);
  else
    [U, T] = schur (A, 'complex');
  endif

  if (isequal (T, tril (T))) ## Handle special case of diagonal T.
    F = U*diag (feval (fun, diag (T)))*U';
    n_swaps = 0; 
    n_calls = 0; 
    terms = 0; 
    ind = {1:n};
    return
  endif

## Determine reordering of Schur form into block form.
  if (isempty (m))
    m = blocking (T, delta, abs (prnt)>=3); 
  endif

  if (prnt)
    fprintf('delta (blocking) = %9.2e, tol (TS) = %9.2e\n', delta, tol)
  endif

  [M, ind, n_swaps] = swapping (m);
  n_calls = size (M,1);
  if (n_calls > 0)            ## If there are swaps to do...
    [U, T] = swap (U, T, M);  ## MEX file
  endif

  m = length (ind);

## Calculate F(T)
  F = zeros (n);

  for col=1:m
    j = ind {col};
    [F(j,j), n_terms] = funm_atom (T (j, j), fun, tol, abs (prnt)*(prnt ~= 1));
    terms (col) = n_terms;

    for row=col-1:-1:1
      i = ind{row};
      if (length(i) == 1 && length(j) == 1)
         ## Scalar case.
        k = i+1:j-1;
        temp = T (i, j)*(F (i, i) - F (j, j)) + F (i, k)*T (k, j) - T (i, k)*F (k, j);
        F (i, j) = temp/(T (i, i)-T (j, j));
      else
        k = cat (2, ind{row+1:col-1});
        rhs = F (i, i)*T (i, j) - T (i, j)*F (j, j) + F (i, k)*T (k, j) - T (i, k)*F (k, j);
        F (i, j) = sylv_tri (T (i, i),-T (j, j), rhs);
      endif
    endfor
  endfor

F = U*F*U';

## As in FUNM:
  if (isreal (A) && norm (imag (F), 1) <= 10*n*eps*norm (F, 1))
    F = real (F);
  endif
endfunction
###############################
###############################
function f = fun_cos (x, k)
##FUN_COS
  if (nargin < 2 | k == 0)
    f = cos(x);
  else
    g = mod (ceil(k/2),2);
    h = mod (k,2);
  if (h == 1)
    f = sin(x)*(-1)^g; 
  else
    f = cos(x)*(-1)^g; 
  endif
  endif 
endfunction

function f = fun_exp(x,k)
##FUN_EXP
f = exp(x);
endfunction

function f = fun_sin(x,k)
##FUN_SIN
  if (nargin < 2 | k == 0)
    f = sin(x);
  else
    k = k - 1;
    g = mod (ceil (k/2),2);
    h = mod(k,2);
  if (h == 1)
    f = sin(x)*(-1)^g; 
  else
    f = cos(x)*(-1)^g; 
  endif
  endif
endfunction
function f = fun_log (x,k)
##FUN_LOG
##Only to be called for plain log evaluation.
  f = log(x);
endfunction

function m = blocking(A,delta,showplot)
##BLOCKING  Produce blocking pattern for block Parlett recurrence.
##          M = BLOCKING(A, DELTA, SHOWPLOT) accepts an upper triangular matrix
##          A and produces a blocking pattern, specified by the vector M,
##         for the block Parlett recurrence.
##         M(i) is the index of the block into which A(i,i) should be placed.
##         DELTA is a gap parameter (default 0.1) used to determine the
##         blocking.
##         Setting SHOWPLOT nonzero produces a plot of the eigenvalues
##         that indicates the blocking:
##           - Black circles show a set of 1 eigenvalue.
##           - Blue circles show a set of >1 eigenvalues.
##             The lines connect eigenvalues in the same set.
##             Red squares show the mean of each set.

##         For A coming from a real matrix it should be posible to take
##         advantage of the symmetry about the real axis.  This code does not.

  a = diag(A); 
  n = length(a);
  m = zeros(1,n); 
  maxM = 0;

  if (nargin < 2 | isempty(delta))
    delta = 0.1; 
  endif
  if (nargin < 3) 
    showplot = 0;  
  endif

  if (showplot)
   clf
   hold on
  endif

  for i = 1:n

    if (m(i) == 0)
    
      m(i) = maxM + 1; ## If a(i) hasn`t been assigned to a set
        maxM = maxM + 1; ## then make a new set and assign a(i) to it.
    endif

    for j = i+1:n
      if (m(i) ~= m(j))    ## If a(i) and a(j) are not in same set.
        if (abs (a (i)-a (j)) <= delta)
          if (showplot)
            plot(real([a(i) a(j)]),imag([a(i) a(j)]),'o-')
          endif

          if (m(j) == 0)
            m(j) = m(i); ## If a(j) hasn`t been assigned to a
        ## set, assign it to the same set as a(i).
          else
          p = max(m(i),m(j)); 
          q = min(m(i),m(j));
          m(m==p) = q; ## If a(j) has been assigned to a set
                       ## place all the elements in the set
                       ##containing a(j) into the set
                       ## containing a(i) (or vice versa).
          m(m>p) = m(m>p) -1;
          maxM = maxM - 1;
                                 ## Tidying up. As we have deleted set
                                 ## p we reduce the index of the sets
                                 ## > p by 1.
           endif
         endif
       endif
    endfor
  endfor

  if (showplot)
    for i = 1:max(m)
      a_ind = a(m==i);
        if (length(a_ind) == 1)
          plot(real(a_ind),imag(a_ind),'ok' )
        else
##plot(real(mean(a_ind)),imag(mean(a_ind)),'sr' )
        endif
    endfor
    grid
    hold off
    box on
  endif
endfunction

function [M, ind, n_swaps] = swapping (m)
##SWAPPING  Confluent permutation by swapping adjacent elements.
##         [ISWAP,IND,N_SWAPS] = SWAPPING(M) takes a vector M containing
##         and constructs a swapping scheme that produces
##         a confluent permutation, with elements ordered by ascending
##         average position. The confluent permutation is obtained by using
##         the LAPACK routine ZTREX to move m(ISWAP(i,2)) to m(ISWAP(i,1))
##         by swapping adjacent elements, for i = 1:SIZE(M,1).
##         The cell array vector IND defines the resulting block form:
##         IND{i} contains the indices of the i'th block in the permuted form.
##         N_SWAPS is the total number of swaps required.

  mmax = max(m); 
  M = []; 
  ind = {}; 
  h = zeros(1, mmax);
  g = zeros(1, mmax);

  for i = 1:mmax
    p = find (m==i);
    h(i) = length (p);
    g(i) = sum (p)/h (i);
  endfor

  [x, y] = sort (g);
  mdone = 1;

  for i = y
    if (any (m (mdone:mdone+h (i)-1) ~= i))
      f = find(m==i); 
      g = mdone:mdone+h(i)-1;
      ff = f (f~=g); 
      gg = g (f~=g);

      ## Create vector v = mdone:f(end) with all elements of f deleted.
      v = mdone-1 + find (m (mdone:f (end)) ~= i);

      ## v = zeros(1,f(end)-g(1)+1);
      ## v(f-g(1)+1) = 1; v = g(1)-1 + find(v==0);

      M (end+1:end+length (gg),:) = [gg' ff'];

      m (g (end)+1:f (end)) = m (v);
      m (g) = i*ones (1,h (i));
      ind = cat (2,ind,{mdone:mdone+h (i)-1}); 
      mdone = mdone + h (i);
    else
      ind = cat (2,ind,{mdone:mdone+h(i)-1}); 
      mdone = mdone + h (i);
    endif
  endfor

  n_swaps = sum (abs (diff (M')));
endfunction

function [F, n_terms] = funm_atom (T, fun, tol, prnt)
##FUNM_ATOM  Function of triangular matrix with nearly constant diagonal.
##          [F, N_TERMS] = FUNM_ATOM(T, FUN, TOL, PRNT) evaluates function
##          FUN at the upper triangular matrix T, where T has nearly constant
##          diagonal.  A Taylor series is used.
##          FUN(X,K) must return the K'th derivative of
##          the function represented by FUN evaluated at the vector X.
##          TOL is a convergence tolerance for the Taylor series,
##          defaulting to EPS.
##          If PRNT ~= 0 trace information is printed.
##          N_TERMS is the number of terms taken in the Taylor series.
##          N_TERMS  = -1 signals lack of convergence.

  if (nargin < 3 | isempty (tol))
    tol = eps; 
  endif
  if (nargin < 4)
    prnt = 0; 
  endif

  if (isequal(fun,@fun_log))   ## LOG is special case.
    [F, n_terms]  = logm_isst (T, prnt);
  return
  endif

  itmax = 500;

  n = length (T);
  if (n == 1)
    F = feval (fun, T, 0); 
    n_terms = 1; 
    return
  endif

  lambda = sum (diag (T))/n;
  F = eye (n)*feval (fun, lambda, 0);
  f_deriv_max = zeros (itmax+n-1, 1);
  N = T - lambda*eye (n);
  mu = norm ( (eye(n)-abs(triu(T,1)))\ones(n,1),inf );

  P = N;
  max_d = 1;

  for k = 1:itmax
    f = feval(fun,lambda,k);
    F_old = F;
    F = F + P*f;
    rel_diff = norm(F - F_old,inf)/(tol+norm(F_old,inf));
    if (prnt)
      fprintf('%3.0f: coef = %5.0e', k, abs(f)/factorial(k));
      fprintf('  N^k/k! = %7.1e', norm(P,inf));
      fprintf('  rel_d = %5.0e',rel_diff);
      fprintf('  abs_d = %5.0e',norm(F - F_old,inf));
    endif
    P = P*N/(k+1);

    if (rel_diff <= tol)

      ## Approximate the maximum of derivatives in convex set containing
      ## eigenvalues by maximum of derivatives at eigenvalues.
      for j = max_d:k+n-1
          f_deriv_max(j) = norm(feval(fun,diag(T),j),inf);
      endfor
      max_d = k+n;
      omega = 0;
      for j = 0:n-1
        omega = max(omega,f_deriv_max(k+j)/factorial(j));
      endfor

      trunc = norm(P,inf)*mu*omega;  ## norm(F) moved to RHS to avoid / 0.
      if (prnt)
        fprintf('  [trunc,test] = [%5.0e %5.0e]', ...
        trunc, tol*norm(F,inf))
      endif
      if (prnt == 5)
        trunc = 0; 
      endif ## Force simple stopping test.
      if (trunc <= tol*norm(F,inf))
        n_terms = k;
      if (prnt)
        fprintf('\n')
      endif 
        return
      endif
    endif

    if (prnt)
      fprintf('\n')
    endif

  endfor
  n_terms = -1;
endfunction



function [X, iter] = logm_isst(T, prnt)
##LOGM_ISST   Log of triangular matrix by Schur-Pade method with scaling.
##        X = LOGM_ISST(A) computes the logarithm of an upper triangular
##        matrix A, for a matrix with no nonpositive real eigenvalues,
##        using the inverse scaling and squaring method with Pade
##        approximation.  TOL is an error tolerance.
##        [X, ITER] = LOGM_ISST(A, PRNT) returns the number ITER of square
##        roots computed and prints this information if PRNT is nonzero.

## References:
##S. H. Cheng, N. J. Higham, C. S. Kenney, and A. J. Laub, Approximating the
##   logarithm of a matrix to specified accuracy, SIAM J. Matrix Anal. Appl.,
##   22(4):1112-1125, 2001.
##N. J. Higham, Evaluating Pade approximants of the matrix logarithm,
##   SIAM J. Matrix Anal. Appl., 22(4):1126-1135, 2001.

  if (nargin < 2)
    prnt = 0; 
  endif
  n = length(T);

  if (any( imag(diag(T)) == 0 & real(diag(T)) <= 0 ))
    error('A must not have nonpositive real eigenvalues!')
  endif

  if (n == 1)
    X = log(T); 
    iter = 0; 
    return
  endif

  R = T;
  maxlogiter = 50;

  for iter = 0:maxlogiter

    phi  = norm(T-eye(n),'fro');

    if (phi <= 0.25)
      if (prnt)
      fprintf('LOGM_ISST computed %g square roots. \n', iter)
      endif
      break
    endif
    if (iter == maxlogiter)
      error('Too many square roots in LOGM_ISST.\n')
    endif

    ## Compute upper triangular square root R of T, a column at a time.
    for j=1:n
      R(j,j) = sqrt(T(j,j));
      for i=j-1:-1:1
        R(i,j) = (T(i,j) - R(i,i+1:j-1)*R(i+1:j-1,j))/(R(i,i) + R(j,j));
      endfor
    endfor
    T = R;
  endfor

  X = 2^(iter)*logm_pf(T-eye(n),8);
endfunction
#########################################
function S = logm_pf(A,m)
##LOGM_PF   Pade approximation to matrix log by partial fraction expansion.
##         Y = LOGM_PF(A,m) approximates LOG(I+A).

  [nodes,wts] = gauss_legendre(m);
##Convert from [-1,1] to [0,1].
  nodes = (nodes + 1)/2;
  wts = wts/2;

  n = length(A);
  S = zeros(n);

  for j=1:m
    S = S + wts(j)*(A/(eye(n) + nodes(j)*A));
  endfor
endfunction
###############################
function [x,w] = gauss_legendre(n)
##GAUSS_LEGENDRE  Nodes and weights for Gauss-Legendre quadrature.

##Reference:
##G. H. Golub and J. H. Welsch, Calculation of Gauss quadrature
##rules, Math. Comp., 23(106):221-230, 1969.

  i = 1:n-1;
  v = i./sqrt((2*i).^2-1);
  [V,D] = eig( diag(v,-1)+diag(v,1) );
  x = diag(D);
  w = 2*(V(1,:)'.^2);
endfunction

function X = sylv_tri(T,U,B)
##SYLV_TRI    Solves triangular Sylvester equation.
##           x = SYLV_TRI(T,U,B) solves the Sylvester equation
##           T*X + X*U = B, where T and U are square upper triangular matrices.

  m = length(T);
  n = length(U);
  X = zeros(m,n);

##Forward substitution.
  for i = 1:n
    X(:,i) = (T + U(i,i)*eye(m)) \ (B(:,i) - X(:,1:i-1)*U(1:i-1,i));
  endfor
endfunction
%!assert (funm(10,@log), log (10),eps)
%!assert (funm (([1 -1 -1;0 1 -1; 0 0 1]),@log), [0 -1 -1.5; 0 0 -1; 0 0 0], 1e-5)
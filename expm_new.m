function [F, s] = expm_new (A, schur_fact)

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

##EXPM_NEW  Matrix exponential.
##   EXPM_NEW(A) is the matrix exponential of A computed using
##   an improved scaling and squaring algorithm with a Pade approximation.
##   It exploits triangularity (if any) of A.
##   EXPM_NEW(A,1) uses an initial transformation to complex Schur form.
##   EXPM_NEW(A,2) uses an initial transformation to real Schur form if A
##   is real.

##   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
##   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
##   970-989, 2009.

##   Awad H. Al-Mohy and Nicholas J. Higham, April 20, 2010.

  upper_triang = isequal (A, triu (A));
  lower_triang = isequal (A, tril (A));

  if (nargin < 2 || isempty (schur_fact) || upper_triang || lower_triang)
    schur_fact = false;
  endif
 
  if (schur_fact == 1)
    [Q,A] = schur (A, 'complex');
  elseif (schur_fact == 2)
    [Q,A] = schur (A);  ## Complex Schur form will result if A is complex.
  endif

## Coefficients of leading terms in the backward error functions h_{2m+1}.
  Coeff = [1/100800, 1/10059033600, 1/4487938430976000,...
           1/5914384781877411840000, 1/113250775606021113483283660800000000];
  u = eps/2;
  n = length (A);
  have_A4 = 0;   ## prior evaluation of A4
  have_A6 = 0;   ## prior evaluation of A6

  m_vals = [3 5 7 9 13];
  ## theta_m for m=1:13.
  theta = [##3.650024139523051e-008
           ##5.317232856892575e-004
            1.495585217958292e-002   ##m_vals = 3
           ##8.536352760102745e-002
            2.539398330063230e-001   ##m_vals = 5
           ##5.414660951208968e-001
            9.504178996162932e-001   ##m_vals = 7
           ##1.473163964234804e+000
            2.097847961257068e+000   ##m_vals = 9
           ##2.811644121620263e+000
           ##3.602330066265032e+000
           ##4.458935413036850e+000
            4.250000000000000e+000]; ## m_vals = 13

  s = 0;

  A2 = A*A;
  eta1 = max (normAm (A2, 2)^(1/4), normAm (A2, 3)^(1/6));
  t = eval_alpha (A, 1);

  if (eta1 <= theta(1) && t == 0)
    F = PadeApproximantOfDegree (m_vals(1));
    schur_backtrans
    return
  endif

  A4 = A2*A2; 
  have_A4 = 1;
  eta2 = max (norm (A4, 1)^(1/4), normAm (A2, 3)^(1/6));
  t = eval_alpha (A, 2);

  if (eta2 <= theta(2) && t == 0)
    F = PadeApproximantOfDegree (m_vals(2));
    schur_backtrans
    return
  endif

  A6 = A2*A4 ; 
  have_A6 = 1;
  eta3 = max (norm (A6, 1)^(1/6), normAm (A4, 2)^(1/8));
  h = zeros (4, 1); 

  h (3) = eval_alpha (A, 3);
  h (4) = eval_alpha (A, 4);

  for i = 3:4
    if (eta3 <= theta(i) && h (i) == 0)
      F = PadeApproximantOfDegree (m_vals(i));
      schur_backtrans
      return
    endif
  endfor

  eta4 = max (normAm (A4, 2)^(1/8), normAm (A2, 5)^(1/10));
  eta5 = min (eta3, eta4);
  s = max (ceil (log2 (eta5/theta(end))), 0); ## Zero must be here
  t = eval_alpha (A/2^s, 5);
  s = s + t;
  A = A/2^s;  A2 = A2/2^(2*s); A4 = A4/2^(4*s); A6 = A6/2^(6*s); ## Scaling
  F = PadeApproximantOfDegree (m_vals(end));
  
  if (lower_triang)
    A = A'; 
    F = F'; 
  endif

  if (upper_triang || lower_triang || schur_fact)
    F = expm_sqtri (A, F, s);
      if (lower_triang)
        F = F'; 
      endif
  else
    for k= 1:s
      F = F*F;
    endfor
  endif
  schur_backtrans

##Nested Functions

function t = eval_alpha (A, k)
  alpha = Coeff(k)*normAm (abs (A), 2*m_vals(k)+1)/norm (A, 1);
  t = max (ceil (log2 (alpha/u)/(2*m_vals(k))), 0);
endfunction

function F = PadeApproximantOfDegree (m)
##PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
##F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
##Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
##Series are evaluated in decreasing order of powers, which is
##in approx. increasing order of maximum norms of the terms.

  c = getPadeCoefficients;
  ## Evaluate Pade approximant.
    switch (m)
      case {3, 5, 7, 9}
        Apowers = cell (ceil ((m+1)/2), 1);
        Apowers {1} = eye (n);
        Apowers {2} = A2; 
        start = 3;
      
      if (have_A4)
        Apowers {3} = A4;  
        start = 4;
      endif
      
      if (have_A6)
        Apowers {4} = A6;
        start = 5;
      endif
      
      for j = start:ceil ((m+1)/2)
        Apowers {j} = Apowers {j-1}*Apowers {2};
      endfor

      U = zeros (n); V = zeros (n);

      for j = m+1:-2:2
        U = U + c (j)*Apowers {j/2};
      endfor
      U = A*U;
      
      for j = m:-2:1
        V = V + c (j)*Apowers {(j+1)/2};
      endfor
      case 13
      ## For optimal evaluation need different formula for m >= 12.
      U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) ...
          + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye (n) );

      V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) ...
          + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye (n);

    endswitch
    F = (-U+V)\(U+V);
  function c = getPadeCoefficients
##GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
##C = GETPADECOEFFICIENTS returns coefficients of numerator
##of [M/M] Pade approximant, where M = 3,5,7,9,13.
    switch m
      case (3)
        c = [120, 60, 12, 1];
      case (5)
        c = [30240, 15120, 3360, 420, 30, 1];
      case (7)
        c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
      case (9)
        c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
             2162160, 110880, 3960, 90, 1];
      case (13)
        c = [64764752532480000, 32382376266240000, 7771770303897600, ...
             1187353796428800,  129060195264000,   10559470521600, ...
             670442572800,      33522128640,       1323241920,...
             40840800,          960960,            16380,  182,  1];
    endswitch
  
  endfunction

endfunction

function schur_backtrans
  if (schur_fact)
    F = Q*F*Q'; 
  endif
endfunction

##Nested Functions
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

function X = expm_sqtri (T, F, s)
## EXPM_SQTRI   Squaring phase of scaling and squaring method.
##   X = EXPM_SQTRI(T/2^s,F,s) carries out the squaring phase
##   of the scaling and squaring method for an upper quasitriangular T,
##   given T/2^s and a Pade approximant F to e^{T/2^s}.
##   It corrects the diagonal blocks blocks at each step.

##   This M-file exploits Code Fragment 2.1 and Code Fragment 2.2 of the
##   reference below.

##   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
##   Algorithm for the Matrix Exponential,SIAM J. Matrix Anal. Appl. 31(3):
##   970-989, 2009.

##   Awad H. Al-Mohy and Nicholas J. Higham, April 19, 2010.

  n = length (T);
  k = 1;
## To turn off exact superdiagonal computation force "istriangular = 0".
  istriangular = isequal (T, triu (T));

  if (n > 1)
    c = abs (diag (T, -1)) > 0;    ## sum(c) = number of 2by2 full blocks
   ## NUMBLK blocks with i'th block in rows/cols INDX{i}.
    numblk = n - sum (c);         ## The number of blocks
    indx = cell (numblk, 1);
    if (c (end) == 0)
      indx {end} = n; c = [c ; 0];
    endif
    for j = 1:numblk
      if c (k)
        indx {j} = k:k+1; k = k+2;
      else
        indx {j} = k; k = k+1;
      endif
    endfor
  endif

  for i = 0:s
    if (i > 0)
    F = F*F;
    endif
    if (istriangular)
    ## Compute diagonal and first superdiagonal.
      for j = 1:2:n
        if (j < n)
          F (j:j+1, j:j+1) = expmT2by2 ( 2^i * T (j:j+1, j:j+1) );
        else
          F (n, n) = exp (2^i * T (n, n));
        endif
      endfor
    else
       ## Quasitriangular case: compute (block) diagonal only.
      for j = 1:numblk
        F (indx {j}, indx {j}) = expm2_by_2 ( 2^i * T (indx{j},indx {j}) );
      endfor
    endif
  endfor

  X = F;

#####################
function X = expm2_by_2 (A)
## EXPM2_BY_2  Exponential for a general 2-by-2 matrix A.

  if (length (A) == 1)
    X = exp (A);
  else
    X = expm2by2full (A);
  endif
endfunction

#####################
function X = expm2by2full (A)

## EXPM2BY2FULL   Exponential of 2-by-2 full matrix.

  a = A (1, 1);
  b = A (1, 2);
  c = A (2, 1);
  d = A (2, 2);

  delta = sqrt ((a-d)^2 + 4*b*c);

  X = exp ((a+d)/2)  * ...
      [ cosh (delta/2) + (a-d)/2*sinch (delta/2),  b*sinch (delta/2)
        c*sinch (delta/2),  cosh (delta/2) + (d-a)/2*sinch (delta/2) ];
endfunction

######################
function y = sinch (x)
  if (x == 0)
    y = 1;
  else
    y = sinh (x)/x;
  endif
endfunction

#####################
function X = expmT2by2 (A)
##EXPMT2BY2    Exponential of 2-by-2 upper triangular matrix.
##   EXPMT2BY2(A) is the exponential of the 2-by-2 upper triangular matrix A.

## Modified from FUNM (EXPM2by2).

  a1 = A (1, 1);
  a2 = A (2, 2);

  ave = (a1+a2)/2; df  = abs (a1-a2)/2;

  if (max (ave, df) < log (realmax))
   ## Formula fine unless it overflows.
    x12 = A (1, 2)*exp ((a1+a2)/2 ) * sinch ((a2-a1)/2 );
  else
   ## Formula can suffer cancellation.
    x12 = A (1, 2)*(exp (a2)-exp (a1))/(a2-a1);
  endif

  X = [exp(a1)  x12
         0      exp(a2)];

endfunction

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
#############################################
#############################################
## 'involuntary matrix A'
%!test 
%! A = [-4 0.5 1/3 0.25;-120 20 15 12;240 -45 -36 -30;-140 28 23+1/3 20];
%! EA = [-3.15772413976001, 0.587600596821908, 0.391733731214599, 0.293800298410948;...
%!       -141.024143237258, 25.0471045076916, 17.628017904657, 14.1024143237255;...
%!        282.048286474513, -52.884053713971, -40.7641623363609, -35.2560358093133;...
%!       -164.528167110131, 32.9056334220263, 27.4213611850215, 25.0471045076907];
%! assert (expm_new (A), EA, 1e4*eps);

## 'cauchy matrix A'
%!test 
%! A = [0.5 1/3;1/3 0.25];
%! EA = [1.73390927994721 0.49530669334232;0.49530669334232 1.36242925994047]; 
%! assert (expm_new (A), EA, 20*eps);

## 'chebyshev spectral differentiation matrix A'
%!test 
%! A = [1.5 -2 0.5;0.5 -3.06151588455594e-017 -0.5;-0.5 2 -1.5];
%! EA = [3 -3 1;1 -6.24500451351651e-016 8.32667268468867e-017;
%!       -1.11022302462516e-016 0.999999999999999 5.55111512312578e-017];
%! assert (expm_new (A), EA, 5*eps);

%!test 
%! A=[0 1e-8 0;-(2e10+4e8/6) -3 2e10;200/3 0 -200/3];
%! EA = [0.446849468283145, 1.54044157383964e-09, 0.46281145355881;...
%!       -5743067.77947979, -0.0152830038686872, -4526542.71278561;...
%!         0.447722977849464, 1.54270484519604e-09, 0.463480648837687];
%! assert (expm_new (A), EA, 5e8*eps)

%!test
%! A = [ -0.04999999999999998 -1.110223024625157e-17 0 0 0 0 0 0 -0.05000000000000002 0;...
%!      -1.110223024625157e-17 -0.1 0 0 0 0 0 0 0.1 0;... 
%!       0 0 -0.02 0.01 0 0.5555555555555556 0 0 0 0.05;...
%!       0 0 -0.01 -0.025 0 1.388888888888889 0 0 0 0.125;...
%!       0 0 0 0 0 -1.333333333333333 0 0 0 0;...
%!       0 0 0 0 1 -2 0 0 0 -0.36;... 
%!       0 0 0 0 0 0 0 -0.4800000000000001 0 0; ...
%!       0.3 3.33066907387547e-17 -0.12 0 0 0 1 -1.2 -0.3 0;...
%!       0 0 0 0 0 0 0 0 0 0;...
%!       0 0 0 0 0 0 0 0 0 0];
%! EA = [ 0.951229424500714 -1.03010947471459e-017 0 0 0 0 0 0 -0.048770575499286 0;...                  
%!       -1.03010947471459e-017 0.904837418035959 0 0 0  0 0 0 0.0951625819640404 0;...
%!        0 0 0.98014974536191 0.00977735959901868 0.143591830497017  0.193881201026749 0 0 0 -0.00157486331694516;...
%!        0 0 -0.00977735959901867 0.9752610655624 0.354213062326703  0.472736022670418 0 0 0 -0.00431458931636983;...
%!       -0 -0 -0 -0 0.656030167547536  -0.463706176832214 -0 -0 -0 0.123829139682887;...
%!        0 0 0 0 0.34777963262416  -0.0395290977007842 0 0 0 -0.125200667744698;...
%!       -0.0474340852766379 -4.98063943411835e-018 0.0191836344938302 7.05264225051023e-005 0.000702022263725654  0.00235697456047202 0.838951554331466 -0.258192515000521 0.0491949821244824 0.000104596353349409;...
%!        0.156429271325676 1.5780937851063e-017 -0.0637473413457518 -0.000395985800782656 -0.00491036366765005 -0.0106365609947363 0.537901072917753 0.193470266830163 -0.166311372424975 -0.000248930595280668;...
%!        0 0 0 0 0 0 0 0 1 0;...
%!        0 0 0 0 0 0 0 0 0 1];
%! assert(expm_new(A), EA, 3*eps)


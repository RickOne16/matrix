function [] = gmres ( A, b , x0 , m )

## GMRES algorithm 6.9
## made by Mudit Sharma

  n = length(A);
  if( nargin < 3 )
     m = n; 
     x0 = [0 ; 0];
  endif
  
  if ( nargin < 4 )
     m = n;
  endif

 [size_m, size_n]  =  size (A);
 if ( size_m ~= size_n ) error ( ' enter a square matrix A ');
 endif
  
 V = zeros( n, m+1 );
 H = zeros( m+1, m );

 r0 = b - (A * x0);
 be = norm( r0 );
 v = r0 / be;

 V(:,1) = v;

   for j = 1:m
      
      w = A * V(:,j);
     
      for i = 1:j
        
        H(i,j) = V(:,i)' * w;
        w = w - V(:,i) * H(i,j);

      endfor
    
      H(j+1,j) = norm(w);
      
      if(!H(j+1,j)) 
        m = j;
        break; 
      endif

      V(:,j+1) = w / H(j+1,j);
   
  endfor
  
  ## Calculate argmin(y) for the function ||b - Ax0 - AVy|| (equation 6.28)
  ## which algorithm to be used 5.3 or 5.4 ?
  
endfunction 

 

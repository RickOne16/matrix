function [H,V] = arnoldi_gram_schmidt(A,v,m)

## Arnoldi Gram - Schmidt , algorithm 6.2
## made by : Mudit Sharma
  
  n = length(A);
  if( nargin < 3 ) m = n; 
  endif

 [size_m, size_n]  =  size (A);
 if ( size_m != size_n ) error ( ' enter a square matrix A ');
 endif

 if ( norm(v) != 1 ) error (' enter unity norm vector v ' );
 endif
  
 V = zeros( n, m+1 );
 H = zeros( m+1, m );

  V(:,1) = v;

   for j = 1:m
      
      w = A * V(:,j);
     
      for i = 1:j
        
        H(i,j) = V(:,i)' * w;
        w = w - V(:,i) * H(i,j);

      endfor
    
      H(j+1,j) = norm(w);
      if(!H(j+1,j)) error (' H is not an upper hessenberg matrix! ');
      endif

      V(:,j+1) = w / H(j+1,j);

   endfor
   V
endfunction



%--------------------------------------------------------------------------
%---- Code to generate a low discrepency sequence for scanning angles -----
%--------------------------------------------------------------------------

% based on the code vdc.m developped by John Burkardt (http://people.sc.fsu.edu/~jburkardt/m_src/van_der_corput/van_der_corput.html)


% parameters:
%    Np: number of projections
%    Ns: number of segments in the circle partition
%    p : a prime integer used as the base
%    Nz: projection of the angle 0 each Nz projectios


function Seq = Low_discrepency_squence(Np,Ns,p,Nz)

Seq = [];
i = 0;
while length(Seq) < Np
   if (length(Seq)>0 && mod (length(Seq),Nz)==0)
       Seq = [Seq,0];
   end
    
   b = i;
   i = i+1;
   r = 0;
   q = 1/p;
   
   while (b~=0)
        a = mod(b,p);
        r = r + a*q;
        q = q/p;
        
        b = floor(b/p);
   end
   
   r = 360*r/Ns;
   
   k = 0;
   while (length(Seq) < Np && k < Ns)
       Seq =[Seq,r+k*360/Ns];
       k = k+1;
   end   
end
end

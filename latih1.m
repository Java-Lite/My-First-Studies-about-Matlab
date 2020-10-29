G= [53.8516   22.3607   31.6228;
    22.3607   58.1378   35.7771;
    31.6228   35.7771   67.2309];
B=[ -68.1986  116.5651  108.4349;
     116.5651  -63.4349  116.5651;
     108.4349  116.5651  -67.2490];
 V=[1.05; 1.00; 1.04];
del=0;
% carip2(1,1,1)
  n = 1;
  while n <= 3 %ambil baris 2 kolom 2 kebawah
      for i = 1:3 %baris berapa dan kolom ke kanan
         nilaiG = G(n,i)
         nilaiB = B(n,i)%tergantung kebutuhan
         nilaiV = V;
         nilaiD=del;
         carip2(nilaiG,nilaiB,nilaiV,nilaiD);
      end 
      n = n+1;
  end
  
  
 function P = carip2(G,B,V,D)
    P(i)=P(i)+V(i)*V(j)*G(i,j)*cos(B(i,j)-D(i)+D(j))
 end

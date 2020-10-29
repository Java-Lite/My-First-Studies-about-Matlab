%% General Program For Newton Raphson Load flow 
% Enter the busdata, and Loaddata in mention form
% Bus data (Bus	 Bus   Vol	   Vol	    Generating	  Load	   Reactive Power limit
%         	no   type  Mag(pu) angle     Pg   QG    Pl     Ql    Qmax   Qmin   
 
% Load Data(From   To       R         X          B            Tap
%            Bus    Bus     (pu)      (pu)       (pu)         Ratio)

clc
clear all
%% Bus data Bus type 1(slack Bus),type2(PV_Bus Bus), type3(PQ_Bus Bus)
tic
		  % Bus	 Bus   Vol	 Vol	  Generating	    Load	   Reactive Power limit
          % no   type  Mag   angle     PG   QG        Pl     Ql        Qmax   Qmin   
busdata=  [ 1     3     0      0       0     0       0.20   0.097       0       0;
            2     3     0      0       0     0       0.30   0.145       0       0;
            3     3     0      0       0     0       0.20   0.097       0       0;
            4     3     0      0       0     0       0.30   0.145       0       0;
            5     3     0      0       0     0       0.20   0.097       0       0;
            6     2    1.00    0      0.30  0.145      0      0         0       0;
            7     2    1.00    0      0.15  0.0726     0      0         0       0;
            8     2    1.00    0      0.20  0.097      0      0         0       0;
            9     2    1.00    0      0.20  0.097      0      0         0       0;
            10    1    1.05    0      0.20  0.097      0      0         0       0];
%            From   To       R         X          B           Tap
%            Bus    Bus     (pu)      (pu)       (pu)         Ratio
linedata=[   1      2       0.02      0.06        0            0
             1      6       0.06      0.24        0            0
             1      9       0.04      0.16        0            0
             2      3       0.06      0.24        0            0
             2      6       0.06      0.24        0            0
             3      7       0.06      0.24        0            0 
             4      7       0.04      0.16        0            0 
             4      8       0.06      0.24        0            0
             5      6       0.04      0.16        0            0
             5     10       0.06      0.24        0            0
             6      9       0.01      0.04        0            0
             8     10       0.04      0.16        0            0
             9     10       0.08      0.32        0            0];

%% Data arranged for Linedata in the different vector 
fb=linedata(:,1);tb=linedata(:,2);
r=linedata(:,3);x=linedata(:,4);
b=linedata(:,5);a=linedata(:,6);
z=r+1i*x; 					% Impedance of branch/eg : 0.0200+0.0600i 
y=1./z;                     % admittance of branch/eg : 5.0000-15.0000i
b=1i*b;  				% admittance of branch 
nl=length(fb);					% No of branch 
No_of_Bus=max(max(fb),max(tb));		% No of Bus 

%% Formation of YBus matrix 

Y=zeros(No_of_Bus,No_of_Bus);				% Initialize of YBus 
for k=1:nl
    Y(fb(k),tb(k))=Y(fb(k),tb(k))-y(k)/a(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end
for m=1:No_of_Bus
    for n=1:nl
        if fb(n)==m
            Y(m,m)=Y(m,m)+y(n)/a(n)^2+b(n);
        elseif tb(n)==m
            Y(m,m)=Y(m,m)+y(n)+b(n);
        end
    end
end
G=real(Y);B=imag(Y);			% Separation of YBus 
%% Data arranged for Linedata in the different vector
BMva=100;
busNo=busdata(:,1);type=busdata(:,2);V=busdata(:,3);del=busdata(:,4);
Pg=busdata(:,5)/BMva;Qg=busdata(:,6)/BMva;Pl=busdata(:,7)/BMva;
Ql=busdata(:,8)/BMva;Qmin=busdata(:,9)/BMva;Qmax=busdata(:,10)/BMva;
PV_Bus=find(type==2|type==1);PQ_Bus=find(type==3);  	% type1(Slack),type2(PV_Bus Bus),type3(PQ_Bus Bus )
No_of_PQ_Bus=length(PQ_Bus);No_of_PV_Bus=length(PV_Bus);
Active_Power_specified=Pg-Pl;Reactive_Power_specified=Qg-Ql; % Net Power flow through different node 
Iter=1;Tol=1; % Iterantion And tolerance 
%% Newton Raphson Load Flow 
while Tol>1e-5
    P=zeros(No_of_Bus,1);
    Q=zeros(No_of_Bus,1);
    for i=1:No_of_Bus
        for j=1:No_of_Bus
 P(i)=P(i)+V(i)*V(j)*(G(i,j)*cos(del(i)-del(j))+B(i,j)*sin(del(i)-del(j)));
 Q(i)=Q(i)+V(i)*V(j)*(G(i,j)*sin(del(i)-del(j))-B(i,j)*cos(del(i)-del(j)));
        end
    end
%Matrica koja sadr�i visine faznih vodi�a i gromobranskog u�eta u odnosu na zemlju%
y = [20 20 20 20 20 20 25]';

%Matrica koja sadr�i razmak faznih vodi�a i gromobranskog u�eta od stuba%
x = [2 4 6 -2 -4 -6 0]';

%Matrica koja sadr�i polupre�nike faznih vodi�a i gromobranskog u�eta%
r = [7.45 7.45 7.45 7.45 7.45 7.45 3.5]';
r = r/1000;

%Matrica koja sadr�i provjes faznih vodi�a i gromobranskog u�eta%
f = [3 3 3 3 3 3 2]';

%Matrica otpornosti faznih vodi�a i gromobranskog u�eta%
Rdc = [0.2304 0.2304 0.2304 0.2304 0.2304 0.2304 0.7444]';

%Definiranje konstante%
mi0 = 4e-7*pi;
ro = 100;
fr = 50;
omega = 2*pi*fr;




%Uzimanje provjesa u obzir%
ym = y-(3/2)*f;

%Modificiranje polupre�nika%
rp = 0.7788*r;

%Udaljenost izme�u i-tog i j-tog vodi�a%
for i=1:7
    for j=1:7
        broj = (x(i)-x(j))^2+(y(i)-y(j))^2;
        d(i,j) = sqrt(broj);
    end
end

%Udaljenosti uvedene radi Carson-ovog korekcionog faktora%
for i=1:7
    for j=1:7
        broj = (x(i)-x(j))^2+(y(i)+y(j))^2;
        H(i,j) = sqrt(broj);
    end
end

%Ra�unanje matrice induktiviteta%
for i=1:7
    for j=1:7
        if i==j
           L(i,i) = (mi0/(2*pi))*(log((H(i,i)/rp(i)))/log(exp(1)))*10^6;
        else
           L(i,j) = (mi0/(2*pi))*(log((H(i,j)/d(i,j)))/log(exp(1)))*10^6;
        end
    end
end

%Ra�unanje matrice reaktansi%
X = omega*L;
X = X/1000;

%Ra�unanje faktora snage%
for i=1:7
    for j=1:7
        if i==j
           cosfi(i,i) = 1;
        else
           cosfi(i,j) = (H(i,i)+H(j,j))/(2*H(i,j));
        end
    end
end

%Ra�unanje Carsonovih korekcionih faktora%
for i=1:7
    for j=1:7
        s = H(i,j);
        a = 0.00281*s*sqrt(fr/ro);
        
        deltaX(i,j) = 4e-4*omega*(0.5*(0.6159315-(log(a)/log(exp(1))))+(sqrt(2)/6)*a*cosfi(i,j));
        deltaR(i,j) = 4e-4*omega*((pi/8)-(sqrt(2)/6)*a*cosfi(i,j));
    end
end

%Formiranje matrice aktivnih otpornosti%
for i=1:7
    for j=1:7
        if i==j
           R(i,i) = Rdc(i);
        else
           R(i,j) = 0;
        end
    end
end

%Kona�ne matrice aktivnih otpornosti i reaktansi%
Ru = R + deltaR;
Xu = X + deltaX;

%Formiranje matrice serijskih impedansi%
Z = complex(Ru,Xu);

%Formiranje matrice admitansi%
Y = Z^-1;

%Eliminacija gromobranskog u�eta%
for i=1:6
    for j=1:6
        Z11(i,j) = Z(i,j);
    end
end

for i=1:6
    Z12(i) = Z(i,7);
end

for j=1:6
    Z21(j) = Z(7,j);
end

Z22 = Z(7,7);

%Matrica impedanci sa izba�enim gromobranskim u�etom%
Zn = Z11 - Z12'*Z22^-1*Z21;

%Matrica admitansi sa izba�enim gromobranskim u�etom%
Yn=Zn^-1;




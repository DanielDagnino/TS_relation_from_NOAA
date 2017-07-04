%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
pos1 = [520, 100, 1800, 1000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
lines = load('lines.txt');
nl = length(lines);

% Read.
fidz = fopen( 'z2.txt', 'r' );
fidT = fopen( 'T2.txt', 'r' );
fidS = fopen( 'S2.txt', 'r' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
zI = (0:5:2000)';

Tz_TS = zeros(length(zI),nl);
Sz_TS = zeros(length(zI),nl);

Tmin = 100*ones(length(zI),1);
Tmax = zeros(length(zI),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = figure(3); clf(h); set(h,'OuterPosition',pos1);

cont = 0;
for k=1:nl

  %---------------------------------------------%
  % 
  z = zeros(lines(k)+2,1);
  T = zeros(lines(k)+2,1);
  S = zeros(lines(k)+2,1);
  
  % 
  z(2:(lines(k)+1)) = fscanf( fidz, '%f', lines(k) );
  T(2:(lines(k)+1)) = fscanf( fidT, '%f', lines(k) );
  S(2:(lines(k)+1)) = fscanf( fidS, '%f', lines(k) );
  
  % 
  z(1) = -0.01;
  T(1) = T(2);
  S(1) = S(2);
  
  z(lines(k)+2) = max( 2000, z(lines(k)+1) ) +0.00001;
  T(lines(k)+2) = T(lines(k)+1);
  S(lines(k)+2) = S(lines(k)+1);
  
  %---------------------------------------------%
%   correct = true;
%   iz=2;
%   while iz < lines(k)+2
%     if ( z(iz-1) >= z(iz) )
%       [k iz]
% %       [lines(k) length(z)]
% %       [z(iz-1) z(iz) z(iz+1) z(iz+2)]
% %       [T(iz-1) T(iz) T(iz+1) T(iz+2)]
%       correct = false;
%       
%       z = z(z~=iz);
%       T = T(T~=iz);
%       S = S(S~=iz);
%       lines(k) = lines(k)-1;
%       iz = iz-1;
%       
% %       break;
% %       return;
%     end
%     iz = iz+1;
%   end

  %---------------------------------------------%
  correct = true;
  for iz=2:(lines(k)+2)
    if ( abs(T(iz-1)-T(iz)) >= 2 )
      correct = false;
      break;
    end
  end

  %---------------------------------------------%
  % Eliminamos XBT con cambios exagerados (sin sentido).
  for iz=1:(lines(k)+2)
    if ( z(iz) > 2001 )
      z = z(1:(iz-1));
      T = T(1:(iz-1));
      S = S(1:(iz-1));
      lines(k) = iz-1-2;
      break;
    end
  end
  
  %---------------------------------------------%
  if ( correct )
    % Cont.
    cont = cont+1;
    
    % Interpolamos.
    TI = interp1( z, T, zI, 'linear' );
    SI = interp1( z, S, zI, 'linear' );
    
    % Save.
    Tz_TS(:,cont) = TI;
    Sz_TS(:,cont) = SI;
  
    % Corregimos la primera medida que sale una mierda muchas veces.
    TI(1) = TI(2);
    SI(1) = SI(2);
    
    Tmin = min( TI, Tmin );
    Tmax = max( TI, Tmax );

%     % Ploteamos.
%     plot( zI, SI, '-' ); hold on;
% %     plot( zI, TI ); hold on;
  end
  
end

fclose(fidz);
fclose(fidT);
fclose(fidS);

% Remove empty entries.
Tz_TS = Tz_TS(:,1:cont);
Sz_TS = Sz_TS(:,1:cont);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h = figure(4); clf(h); set(h,'OuterPosition',pos1);
% % plot( zI, Tmin, zI, Tmax );
% 
% h = figure(5); clf(h); set(h,'OuterPosition',pos1);
% % for k=1:cont
% % for k=1:1
% %   plot( Tz_TS(:,k), Sz_TS(:,k), 'o' ); hold on;
% % end

iz=390;
p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
plot( Sz_TS(iz,:), Tz_TS(iz,:), 'o' ); hold on;
Tr = p(2) + p(1)*Sz_TS(iz,:);
plot( Sz_TS(iz,:), Tr, '-r' ); hold on;

% iz=10;
% p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
% plot( Sz_TS(iz,:), Tz_TS(iz,:), 'o' ); hold on;
% Tr = p(2) + p(1)*Sz_TS(iz,:);
% plot( Sz_TS(iz,:), Tr, '-r' ); hold on;


Tr = p(2) + p(1)*Sz_TS(iz,:);
Tr = Tz_TS(iz,:)-Tr;
error = sqrt(sum(Tr.^2))/length(Tr)
h = figure(2); clf(h);
plot( Sz_TS(iz,:), Tr, 'o' ); hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iz1 = 1 + floor( 150/5);
iz2 = 1 + floor( 600/5);
iz3 = 1 + floor(1100/5);
iz4 = 1 + floor(1600/5);

e1 = 0;
for iz = iz1:iz2
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Tr = p(2) + p(1)*Sz_TS(iz,:);
  Tr = Tz_TS(iz,:)-Tr;
  e1 = e1 + sqrt(sum(Tr.^2))/length(Tr);
end
e1/(iz2-iz1+1)

e2 = 0;
for iz = iz2:iz3
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Tr = p(2) + p(1)*Sz_TS(iz,:);
  Tr = Tz_TS(iz,:)-Tr;
  e2 = e2 + sqrt(sum(Tr.^2))/length(Tr);
end
e2/(iz3-iz2+1)

e3 = 0;
for iz = iz3:iz4
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Tr = p(2) + p(1)*Sz_TS(iz,:);
  Tr = Tz_TS(iz,:)-Tr;
  e3 = e3 + sqrt(sum(Tr.^2))/length(Tr);
end
e3/(iz4-iz3+1)

ee = (e1+e2+e3)/(iz4-iz1+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e1 = 0;
for iz = iz1:iz2
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Sr = (Tz_TS(iz,:)-p(2))/p(1);
  Sr = Sz_TS(iz,:)-Sr;
  e1 = e1 + sqrt(sum(Sr.^2))/length(Sr);
end
e1/(iz2-iz1+1)

e2 = 0;
for iz = iz2:iz3
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Sr = (Tz_TS(iz,:)-p(2))/p(1);
  Sr = Sz_TS(iz,:)-Sr;
  e2 = e2 + sqrt(sum(Sr.^2))/length(Sr);
end
e2/(iz3-iz2+1)

e3 = 0;
for iz = iz3:iz4
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Sr = (Tz_TS(iz,:)-p(2))/p(1);
  Sr = Sz_TS(iz,:)-Sr;
  e3 = e3 + sqrt(sum(Sr.^2))/length(Sr);
end
e3/(iz4-iz3+1)

ee = (e1+e2+e3)/(iz4-iz1+1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
h = figure(6); clf(h); set(h,'OuterPosition',pos1);

% File.
fid = fopen( './TSz_regres_NOAA.txt', 'w' );

for iz=2:length(zI)
  %  
  p = polyfit( Sz_TS(iz,:), Tz_TS(iz,:), 1 );
  Tr = p(2) + p(1)*Sz_TS(iz,:);
  
  %  
  plot( Sz_TS(iz,:), Tr, '-r' ); hold on;
  
  %  
  fprintf( fid, '%16.9e %16.9e  %16.9e\n', zI(iz), p(2), p(1) );

end

fclose(fid);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   C�lculo do Fluxo de Pot�ncia Linearezado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Desenvolvido por: Claudio Siervi M. Jr. 24/05/15

%Para a disciplina de An�lise e Opera��o de Sistemas El�tricos de Pot�ncia,
%do Programa de P�s-Gradua��o em Eng. El�trica da Universidade Federal do Paran�.

%Ministrada pela Prof� Dr�. Elizete Maria Louren�o.

function LPF

clc;
clear all;

% Os arquivos de entrada devem estar na mesma pasta que a fun��o

% ArquivoBarras = 'DadosBarras.txt';
% ArquivoLinhas = 'DadosLinhas.txt';

ArquivoBarras = 'DadosBarras_14b.txt';
ArquivoLinhas = 'DadosLinhas_14b.txt';

% ArquivoLinhas = 'DadosLinhas_30b.txt';
% ArquivoBarras = 'DadosBarras_30b.txt';

DadosBarras = importdata(ArquivoBarras);
DadosLinhas = importdata(ArquivoLinhas);

%08/19/93 UW ARCHIVE           100.0  1962 W IEEE 14 Bus Test Case
Sbase = 100;    % Base para transforma��o

%% Leitura dos dados de Barra
[numb, numItensB] = size(DadosBarras);

Num  = DadosBarras(:, 1); Tipo  = DadosBarras(:, 2); V  = DadosBarras(:, 3); 
Teta = DadosBarras(:, 4); Pd  = DadosBarras(:, 5); Qd  = DadosBarras(:, 6);
Pg  = DadosBarras(:, 7); Qg  = DadosBarras(:, 8); %bshb  = j*Barras(:, 9);

% Transforma pot�ncia em pu e os �ngulos em radianos

Teta   = (pi/180)*Teta;    % �ngulo em radianos
Pd = (1/Sbase)*Pd; Qd = (1/Sbase)*Qd; % Pot�ncias em p.u.
Pg = (1/Sbase)*Pg; Qg = (1/Sbase)*Qg;

%% Leitura dos dados de Linha
[numl, NumItensL] = size(DadosLinhas);

na = DadosLinhas(:, 1); nb = DadosLinhas(:, 2); r = DadosLinhas(:, 3); 
x = DadosLinhas(:, 4); bshl = DadosLinhas(:, 5)*j; 
    %a = DadosLinhas(:, 6);
    %phi = DadosLinhas(:, 7);
a = ones(numl,1); % Rela��o de transforma��o
phi = zeros(numl,1); % Angulo de defasagem
phi   = (pi/180)*phi;    % �ngulo do TD em radianos

%% Imped�ncia Linha
z = r + j*x;

%% Adimit�ncia Linha
y = 1./z;

%% a) Determina��o da matriz de admit�ncias Y

Y = zeros(numb, numb);
for l = 1:numl
     % Admit�ncia Barra
     Y(na(l),na(l)) = Y(na(l),na(l)) + (a(l).^2) * y(l) - bshl(l);
     Y(nb(l),nb(l)) = Y(nb(l),nb(l)) + (a(1).^2) * y(l) - bshl(l);
     
     % Admit�ncia Linha
     Y(na(l),nb(l)) = - a(l) * y(l) * exp(-j*phi(l));
     Y(nb(l),na(l)) = - a(l) * y(l) * exp(-j*phi(l));
end

%% Vetor de Condut�ncia
G = real(Y);

%% Vetor de Suscept�ncia
B = imag(Y);

%% 2-b) Matriz de Suscept�ncias
  
Bl = zeros(numb, numb); % barras k , m

for l = 1:numl
    Bl(na(l),na(l)) = Bl(na(l),na(l)) + (1/x(l));
    Bl(nb(l),nb(l)) = Bl(nb(l),nb(l)) + (1/x(l));
    
    Bl(na(l),nb(l)) = -(1/x(l));
    Bl(nb(l),na(l)) = -(1/x(l));
end

%% 2)-   Fluxo Linearizado
%
perdas = sum(Pg) - sum(Pd);

pos_BR = find(Tipo==3); % Barra de Refer�ncia: IEEE -> 3, Prof� -> 2
if size(pos_BR,1) == 0
    disp('N�o existe barra de refer�ncia no sistema.');
    exit 
end
Pg(pos_BR) = Pg(pos_BR) - perdas;

P = Pg - Pd;    % Pot�ncia injetada

Bl(pos_BR,pos_BR) = 1000000000; % Zera coluna e linha pivo quando invertida   

P_TD = P;   % Pot�ncia injetada na presen�a de TD

pos_TD = find(phi ~= 0);    % Posi��o dos TD's

cont_null = nnz(pos_TD);    % Verifica se existe TD
if (cont_null ~= 0)
    TD = (phi(pos_TD)/x(pos_TD)); % Valor TD
    
    P_TD(na(pos_TD)) = P_TD(na(pos_TD)) - TD;
    P_TD(nb(pos_TD)) = P_TD(na(pos_TD)) + TD;
end

Teta_FPL_rad = Bl\P_TD; % �ngulos das tens�es

Teta_FPL_gr = Teta_FPL_rad *(180/pi); % radiano -> grau

% Teta = Teta*(180/pi);
% 
% %   Verifica��o
%  disp([(1:numb)' Teta_FPL_gr Teta]);

% Distribui��o dos FPL
for l = 1:numl
    Dist_FPL(l) = Teta_FPL_gr(na(l)) - Teta_FPL_gr(nb(l));
end

% Verifica��o
Dist_FPL = -Dist_FPL';
disp(['    Distridui��o dos Fluxos de Pot�ncia']);
disp(['      nb    ' '    na    ' '  Dist_FPL']);
disp([nb na Dist_FPL]);

Escreve_Resultados(DadosBarras, DadosLinhas, Teta_FPL_gr)

end

function Escreve_Resultados(Barras, Linhas, Teta_Final)

 [linBarras,colBarras]=size(Barras);

NomeArquivo = strcat('ieee' , int2str(linBarras) ,'b.txt');
fid = fopen(NomeArquivo, 'wt');

    % Imprime cabe�alho 
    titulo.Data = datestr(today,'dd/mm/yy'); 
    titulo.OriginatorName = 'UW ARCHIVE         ';
    titulo.MVABase = sprintf('%3.1f',100);
    titulo.Year = datestr(today,'yyyy') ;
    titulo.Season = 'S';
    titulo.CaseIdentification =  strcat('IEEE ' , int2str(linBarras) , ...
        ' Bus Test Case');
    titulo.NumItems = sprintf('%-5i ITEMS',linBarras);

    fprintf(fid,'%8s  ', titulo.Data);
    fprintf(fid,'%s  ', titulo.OriginatorName);
    fprintf(fid,'%s  ', titulo.MVABase);
    fprintf(fid,'%s ', titulo.Year);
    fprintf(fid,'%s ', titulo.Season);
    fprintf(fid,'%s\n', titulo.CaseIdentification);

    fprintf(fid,'BUS DATA FOLLOWS                            ');
    fprintf(fid,'%s\n\n', titulo.NumItems);

    % Imprime tabela BUS DATA FOLLOWS
    for i=1:linBarras
        % Estruturas no formato IEEE
        BusDataFollows(i).BusNumber = Barras(i,1);
        BusDataFollows(i).Name = sprintf('%-4s%-4i%-2s','Bus', i);
        %BusDataFollows(i).Name = sprintf('%-4s%-4i%-2s','Bus', i, 'HV');
        BusDataFollows(i).FlowAreaNum = 1;
        BusDataFollows(i).LossZoneNumber = 1;
        BusDataFollows(i).Type = Barras(i,2);
        BusDataFollows(i).FinalVoltage = Barras(i,3);
        BusDataFollows(i).FinalAngle = sprintf('%-2.2f',Teta_Final(i,1)); 
        BusDataFollows(i).LoadMw = sprintf('%-2.1f',Barras(i,5)); 
        BusDataFollows(i).LoadMVAR = sprintf('%-2.1f',Barras(i,6)); 
        BusDataFollows(i).GenerationMW = sprintf('%7.1f',Barras(i,7)); 
        BusDataFollows(i).GenerationMVAR = sprintf('%7.1f',Barras(i,8)); 
        BusDataFollows(i).BaseKV = sprintf('%6.1f',0);
        BusDataFollows(i).DesiredVolts = sprintf('%5.1f',0);
        BusDataFollows(i).MaximumMVAR = sprintf('%7.1f',0);
        BusDataFollows(i).MinimumMVAR = sprintf('%7.1f',0);
        BusDataFollows(i).ShuntConductance = sprintf('%7.1f',0);
        %dados(i).ShuntSusceptamce = sprintf('%7.2f',Barras(i,9));
        BusDataFollows(i).ShuntSusceptamce = sprintf('%7.2f', 0);
        BusDataFollows(i).RemoteContrBusNo = sprintf('%3s','  0');
        
        %   Impress�o de dados
        fprintf(fid,'%3i', BusDataFollows(i).BusNumber);        
        fprintf(fid, '%s','   ');               
        fprintf(fid, '%-10s',BusDataFollows(i).Name);           
        fprintf(fid, '%s','  ');
        fprintf(fid, '%-1i',BusDataFollows(i).FlowAreaNum);     
        fprintf(fid, '%s',' ');
        fprintf(fid, '%2i',BusDataFollows(i).LossZoneNumber);   
        fprintf(fid, '%s','  ');
        fprintf(fid, '%1i',BusDataFollows(i).Type);             
        fprintf(fid, '%s','  ');
        fprintf(fid, '%4.3f',BusDataFollows(i).FinalVoltage);   
        fprintf(fid, '%s',' ');
        fprintf(fid, '%6s',BusDataFollows(i).FinalAngle);       
        fprintf(fid, '%s',' ');
        fprintf(fid, '%8s',BusDataFollows(i).LoadMw);           
        fprintf(fid, '%s',' ');
        fprintf(fid, '%9s',BusDataFollows(i).LoadMVAR);         
        fprintf(fid, '%s',' ');
        fprintf(fid, '%-7s',BusDataFollows(i).GenerationMW);    
        fprintf(fid, '%s',' ');
        fprintf(fid, '%-7s',BusDataFollows(i).GenerationMVAR);  
        fprintf(fid, '%s','  ');
        fprintf(fid, '%-6s',BusDataFollows(i).BaseKV);          
        fprintf(fid, '%s','  ');
        fprintf(fid, '%-5s',BusDataFollows(i).DesiredVolts);    
        fprintf(fid, '%s',' ');
        fprintf(fid, '%-7s',BusDataFollows(i).MaximumMVAR);     
        fprintf(fid, '%s',' ');
        fprintf(fid, '%-7s',BusDataFollows(i).MinimumMVAR);     
        fprintf(fid, '%s',' ');
        fprintf(fid, '%-7s',BusDataFollows(i).ShuntConductance);    
        fprintf(fid, '%s',' ');
        fprintf(fid, '%-7s',BusDataFollows(i).ShuntSusceptamce);    
        fprintf(fid, '%s','  ');
        fprintf(fid, '%-3s',BusDataFollows(i).RemoteContrBusNo);    
        fprintf(fid,'\n');
    end

    fprintf(fid,'-999\n');
    
    [linLinhas,colLinhas]=size(Linhas); 
    
    % Imprime tabela BRANCH DATA FOLLOWS
    titulo.NumItemsBranch = sprintf('%-5i ITEMS',linLinhas);

    fprintf(fid,'BRANCH DATA FOLLOWS                            ');
    fprintf(fid,'%s\n\n', titulo.NumItemsBranch);

    for i=1:linLinhas       
        BranchDataFollows(i).BusNumber = Linhas(i,1);
        BranchDataFollows(i).ZBusNumber = Linhas(i,2);
        BranchDataFollows(i).LoadFlowArea = 1;
        BranchDataFollows(i).LossZone = 1;
        BranchDataFollows(i).Circuit = 1;
        BranchDataFollows(i).Type = 0;
        BranchDataFollows(i).BranchResistanceR = Linhas(i,3);
        BranchDataFollows(i).BranchReactanceX = Linhas(i, 4);
        BranchDataFollows(i).LineChargingB = Linhas(i, 5);
        BranchDataFollows(i).LineMVARatingNo1 = 0;
        BranchDataFollows(i).LineMVARatingNo2 = 0;
        BranchDataFollows(i).LineMVARatingNo3 = 0;
        BranchDataFollows(i).ControlBusNumber = 0;
        BranchDataFollows(i).Side = 0;
        BranchDataFollows(i).TransformerRatio = 0;
        BranchDataFollows(i).TransformerAngle = 0;
        BranchDataFollows(i).MinimumTap = 0;
        BranchDataFollows(i).MaximumTap = 0;
        BranchDataFollows(i).StepSize = 0;
        BranchDataFollows(i).MinimumVoltage = 0;
        BranchDataFollows(i).MaximumVoltage = 0;

        fprintf(fid,'%3i  ', BranchDataFollows(i).BusNumber);
        fprintf(fid,'%3i  ', BranchDataFollows(i).ZBusNumber);
        fprintf(fid,'%i ', BranchDataFollows(i).LoadFlowArea);
        fprintf(fid,'%i  ', BranchDataFollows(i).LossZone);
        fprintf(fid,'%i ', BranchDataFollows(i).Circuit);
        fprintf(fid,'%i ', BranchDataFollows(i).Type);
        fprintf(fid,'%-9.5f ', BranchDataFollows(i).BranchResistanceR);
        fprintf(fid,'%-9.5f  ', BranchDataFollows(i).BranchReactanceX);
        fprintf(fid,'%-9.5f ', BranchDataFollows(i).LineChargingB);
        fprintf(fid,'%-4i  ', BranchDataFollows(i).LineMVARatingNo1);
        fprintf(fid,'%-4i  ', BranchDataFollows(i).LineMVARatingNo2);
        fprintf(fid,'%-4i  ', BranchDataFollows(i).LineMVARatingNo3);
        fprintf(fid,'%-3i ', BranchDataFollows(i).ControlBusNumber);
        fprintf(fid,'%i   ', BranchDataFollows(i).Side);
        fprintf(fid,'%5.3f  ', BranchDataFollows(i).TransformerRatio);
        fprintf(fid,'%6.1f ', BranchDataFollows(i).TransformerAngle);
        fprintf(fid,'%-6.1f ', BranchDataFollows(i).MinimumTap);
        fprintf(fid,'%-6.1f  ', BranchDataFollows(i).MaximumTap);
        fprintf(fid,'%-5.1f  ', BranchDataFollows(i).StepSize);
        fprintf(fid,'%-5.1f  ', BranchDataFollows(i).MinimumVoltage);
        fprintf(fid,'%-5.1f\n', BranchDataFollows(i).MaximumVoltage);
    end
    fprintf(fid,'-999\n');

    % Imprime tabela LOSS ZONES FOLLOWS
    titulo.NumItemsZones = sprintf('%-5i ITEMS',1);

    fprintf(fid,'LOSS ZONES FOLLOWS                     ');
    fprintf(fid,'%s\n', titulo.NumItemsZones);

    for i=1:1       
        LossZones(i).Number = i;
        LossZones(i).Name = strcat('IEEE ' , int2str(linBarras) ,' BUS');
        
        fprintf(fid,'%3i ', LossZones(i).Number);      
        fprintf(fid,'%s\n', LossZones(i).Name);
    end
    fprintf(fid,'-99\n');

    % Imprime tabela INTERCHANGE DATA FOLLOWS
    titulo.NumItemsInterchange = sprintf('%-5i ITEMS',1);

    fprintf(fid,'INTERCHANGE DATA FOLLOWS               ');
    fprintf(fid,'%s\n', titulo.NumItemsInterchange);

    for i=1:1
        InterchangeData(i).Number = i;
        InterchangeData(i).SlackBus = 2;
        InterchangeData(i).AlternateSwing = 'Bus 2     HV';
        InterchangeData(i).AreaExport = 0;
        InterchangeData(i).AreaTolerance = 999.99;
        InterchangeData(i).AreaCode =  strcat('IEEE',int2str(linBarras)); 
        InterchangeData(i).AreaName = strcat('IEEE',int2str(linBarras),' Bus Test Case');

        fprintf(fid,'%2i ', InterchangeData(i).Number);
        fprintf(fid,'%4i ', InterchangeData(i).SlackBus);
        fprintf(fid,'%s', InterchangeData(i).AlternateSwing);
        fprintf(fid,'%7.1f  ', InterchangeData(i).AreaExport);
        fprintf(fid,'%6.2f  ', InterchangeData(i).AreaTolerance);
        fprintf(fid,'%s  ', InterchangeData(i).AreaCode);
        fprintf(fid,'%s\n', InterchangeData(i).AreaName);

    end
    fprintf(fid,'-9\n');

    fprintf(fid, 'TIE LINES FOLLOWS                     0 ITEMS\n');
    fprintf(fid, '-999\n');
    fprintf(fid, 'END OF DATA\n');

fclose(fid);

disp(['Os resultados est�o no arquivo:' NomeArquivo]);
end

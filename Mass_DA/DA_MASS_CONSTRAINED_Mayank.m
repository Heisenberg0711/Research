function [C,Y,D_Act,P] = DA_MASS_CONSTRAINED_Mayank
close all
%     data = load('Data/D31.txt');
%     X = data(:,1:2);

    % Data for comparison with the simultaneous appraoch
%     Centroids = [2, 2; 3, 5; 5, 5; 4.25, 2];
%     rng('default');
%     rng(1);
%     X = Centroids;
%     for i = 1 : length(Centroids)
%         for j = 1 : 30
%             R = 0.01 + 0.5*rand(1);
%             theta = 2*pi*rand(1);
%             randx = Centroids(i,1) + R*cos(theta);
%             randy = Centroids(i,2) + R*sin(theta);
%             X = [X; randx, randy];
%         end
%     end
% 
%     C = [6.5, 4];
%     [M, N] = size(X);
%     %Normalize data 
%     Centroids = [2.5, 3; 3, 5; 4.25, 2];
%     rng('default');
%     rng(1);
%     X = Centroids;
%     for i = 1 : length(Centroids)
%         for j = 1 : 30
%             R = 0.01 + 0.5*rand(1);
%             theta = 2*pi*rand(1);
%             randx = Centroids(i,1) + R*cos(theta);
%             randy = Centroids(i,2) + R*sin(theta);
%             X = [X; randx, randy];
%         end
%     end
% 
%     C = [6.5, 4];
%     [M, N] = size(X);

Centroids = [1.5 2.5; 3 1; 5 1; 8 1; 10, 2.5; 9, 5];
rng('default');
rng(1);
X = Centroids;
for i = 1 : length(Centroids)
    for j = 1 : 20+randi([0,40],1,1)
        R = 0.1 + 0.5*rand(1);
        theta = 2*pi*rand(1);
        randx = Centroids(i,1) + R*cos(theta);
        randy = Centroids(i,2) + R*sin(theta);
        X = [X; randx, randy]; %Adding datapoints to center X
    end
end

[M,N] = size(X);
C = [4, 7];
Xmean = mean(X);
Xstd = std(X);

%Normalizing the data so that it has zero mean and unit std.
X = (X-ones(M,1)*Xmean)./(ones(M,1)*Xstd);
    Kmax = 5; Tmin = 0.0005; alpha = 0.9; PERTURB = 0.01; STOP = 1e-7;
    
% Set initial data
% T is the inital temperature
% Px is the probablity of each data point
% Y represents cluster locations

    T = 120; K = 1; Px = (1/M)*ones(M,1); Y = Px'*X; Py = 1; 
    
    SPLIT = 0.4; inds = []; flag = 0;
    while(T > Tmin)
        if(flag == 0)
%             initialize and perturb
            if(K < Kmax) %phase 1, duplicate codevectors
				YY = [Y;Y] + PERTURB*randn(2*K,N); %Introducing new cluster
                Pyy = 0.5*[Py;Py];
                KK = 2*K;
            else %phase 2, do nothing
                YY = Y;
                Pyy = Py;
                KK = K;
            end
            L_old = inf;
            
%             iterate until convergence
            while 1
                [D,D_Act] = distortion(X,YY,M,N,KK);
                num = repmat(Pyy',[M 1]).*exp(-D/T);
                den = repmat(sum(num,2),[1 KK]);
                P = num./den;
                Pyy = P'*Px;
                YY = P'*(X.*repmat(Px,[1 N]))./repmat(Pyy,[1 N]);
                L = -T*Px'*log(sum(exp(-D/T),2));
                if(norm(L-L_old) < STOP)
                    break;
                end
                L_old = L;
            end
            
%determine the distince codevectors
            Y = [];
            Py = [];
            dist=2*SPLIT;
            
            for i = 1:KK
                for j = 1:size(Y,1)
                    dist = norm(YY(i,:) - Y(j,:));
                    if(dist < SPLIT)
                        break;
                    end
                end
                if(dist > SPLIT)
                    Y = [Y;YY(i,:)]; 
                    inds = [inds;i];
                    Py = [Py;Pyy(i)];
                else
                    Py(j) = Py(j) + Pyy(i);
                end

                K = size(Y,1);
                if(K > Kmax)
                    [sortPy, sortPy_ind] = sort(Py);
                    Y = Y(sortPy_ind(1:Kmax),:);
                    Py = sortPy(1:Kmax);
                    flag = 1;
                    K = Kmax;
                end
            end
        elseif(flag == 1)
            while(1)
                [D,D_Act] = distortion(X,Y,M,N,K);
                num = repmat(Py',[M 1]).*exp(-D/T);
                den = repmat(sum(num,2),[1 K]);
                P = num./den;
                Py = P'*Px;
                Y = P'*(X.*repmat(Px,[1 N]))./repmat(Py,[1 N]);
                L = -T*Px'*log(sum(exp(-D/T),2));
                if(norm(L-L_old) < STOP)
                    break;
                end
                L_old = L;
            end
        end
        
        T = alpha*T;
    end
    
    X = X.*(ones(M,1)*Xstd) + ones(M,1)*Xmean;
    Y = Y.*(ones(K,1)*Xstd) + ones(K,1)*Xmean;
    
%     figure();
%     %h1_1 =plot(X(:,1),X(:,2),'+b','MarkerSize',8);
%     h1 =plot(X(:,1),X(:,2),'+b','MarkerSize',16);hold on;
%     h2 = plot(Y(:,1),Y(:,2),'*r','MarkerSize',20);
%     %h2_1 = plot(Y(:,1),Y(:,2),'*r','MarkerSize',8);
%     h6 = plot([Y(1,1),C(1,1)],[Y(1,2),C(1,2)],':r','LineWidth',2);
%     h7 = plot([Y(2,1),C(1,1)],[Y(2,2),C(1,2)],':r','LineWidth',2);
%     h8 = plot([Y(3,1),C(1,1)],[Y(3,2),C(1,2)],':r','LineWidth',2);
%     text(Y(1,1)+0.25,Y(1,2)-0.25,'Y1','VerticalAlignment','top','FontSize',25);
%     text(Y(2,1)+0.25,Y(2,2)-0.25,'Y2','VerticalAlignment','top','FontSize',25);
%     text(Y(3,1)+0.30,Y(3,2)-0.30,'Y3','VerticalAlignment','top','FontSize',25);
%     %plot([Y(end,1), C(1,1)],[Y(end,2), C(1,2)],':r');
%     h3 = plot(C(1,1),C(1,2),'square','MarkerSize',22,'MarkerFaceColor','k');
%     text(C(1,1)-0.05,C(1,2)-0.1,'r','VerticalAlignment','top','FontSize',25);
%     %h3_1 = plot(C(1,1),C(1,2),'dk','MarkerSize',8,'MarkerFaceColor','k');
%     h_legend = legend([h1, h2, h3],{'Sensors', 'CRUs','Common Destination'},'Location', 'NorthEast');
%     set(h_legend,'FontSize',14,'Interpreter', 'latex');
%     a=get(h_legend,'children');
%     disp(a);
%     set(a,'markersize',2); 
%     %xlim([2 6.7]);
%     %ylim([1.3 5.5]);
%     hold off;
%     %disp(Y);   
%figure();
figure();
h1 = plot(X(:,1),X(:,2),'^b','MarkerSize',6,'MarkerFaceColor','b'); hold on;
h2 = plot(Y(:,1),Y(:,2),'or','MarkerSize',10,'MarkerFaceColor','r');
h3 = plot(C(1,1),C(1,2),'diamond','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','k');
text(Y(1,1),Y(1,2),'$r_1$','VerticalAlignment','top','HorizontalAlignment','right','FontSize',20,'Interpreter','Latex');
text(Y(2,1),Y(2,2),'$r_2$','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',20,'Interpreter','Latex');
text(Y(3,1),Y(3,2),'$r_3$','VerticalAlignment','top','HorizontalAlignment','right','FontSize',20,'Interpreter','Latex');
text(Y(4,1),Y(4,2),'$r_4$','VerticalAlignment','top','HorizontalAlignment','right','FontSize',20,'Interpreter','Latex');
text(Y(5,1),Y(5,2),'$r_5$','VerticalAlignment','top','HorizontalAlignment','right','FontSize',20,'Interpreter','Latex');
text(C(1,1),C(1,2),' $\delta$ ','VerticalAlignment','middle','HorizontalAlignment','right','FontSize',25,'Interpreter','Latex');
h_legend = legend([h2, h1, h3],{ 'Resources','Sensors','Destination'}, 'Location', 'NorthEast');
xlim([0 11]);
ylim([0 7.5]);
set(h_legend,'FontSize',14,'Interpreter', 'latex');

pos = zeros(M,1);
for i = 1 : M
    [c, pos(i)] = max(P(i,:));
end
R1_sensors = X(pos == 1,:);
R2_sensors = X(pos == 2,:);
R3_sensors = X(pos == 3,:);
R4_sensors = X(pos == 4,:);
R5_sensors = X(pos == 5,:);
plot(R1_sensors(:,1),R1_sensors(:,2),'^r','MarkerSize',18,'MarkerFaceColor','r');
plot(Y(1,1),Y(1,2),'or','MarkerSize',25,'MarkerFaceColor','r');
plot(R2_sensors(:,1),R2_sensors(:,2),'^g','MarkerSize',18,'MarkerFaceColor','g');
plot(Y(2,1),Y(2,2),'og','MarkerSize',25,'MarkerFaceColor','g');
plot(R3_sensors(:,1),R3_sensors(:,2),'^r','MarkerSize',18,'MarkerFaceColor','r');
plot(Y(3,1),Y(3,2),'or','MarkerSize',25,'MarkerFaceColor','r');
plot(R4_sensors(:,1),R4_sensors(:,2),'^b','MarkerSize',18,'MarkerFaceColor','b');
plot(Y(4,1),Y(4,2),'ob','MarkerSize',25,'MarkerFaceColor','b');
plot(R5_sensors(:,1),R5_sensors(:,2),'^r','MarkerSize',18,'MarkerFaceColor','r');
plot(Y(5,1),Y(5,2),'or','MarkerSize',25,'MarkerFaceColor','r');
%h_legend = legend([h2, h4, h1, h3],{ 'Resources','Sensors','Destination'}, 'Location', 'NorthEast');
h3 = plot(C(1,1),C(1,2),'diamond','MarkerSize',25,'MarkerFaceColor','k','MarkerEdgeColor','k');
xlim([0 11]);
ylim([0 7.5]);
%set(h_legend,'FontSize',14,'Interpreter', 'latex');
hold off;
end

%compute distance between every datapoint and every resource
function [D,D_Act] = distortion(X,Y,M,N,K)
    X = X'; X = X(:); X = X';
    Xar = repmat(X,[K 1]);
    Yar = repmat(Y,[1 M]);
    D2 = (Xar - Yar).^2;
    
    D = zeros(K,M);
    for i = 1:N
        D = D + D2(:,i:N:end);
    end
    D = D';
    D_Act = D;
    Dm = min(D')';
    D = D - repmat(Dm,[1 K]);
    
end

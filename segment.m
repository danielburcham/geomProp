function [N] = segment(A,B,R1,R2,input1,imageType)
%{
DESCRIPTION
Segments tomogram into solid and decayed regions and returns the associated 
boundary coordinates
NOTES
This version computes estimates for trees with one more more columns of 
decay and one or more open cavities
---------------------------------------------------------------------------
INPUTS
A: Geometry image file showing only the blue trunk boundary line.
B: Sonic tomogram image file showing visualized decay pattern.
R1: Coordinate systen reference object for trunk geometry file.
R2: Coordinate systen reference object for tomogram file.
input1: Choose liberal (0) or conservative (1) estimate.
OUTPUTS
N: 1x2 cell array with Cartesian coordinates (m) for solid and decayed
regions contained in the first and second cell, respectively. Coordinate
pairs are recorded as row-wise entries in a mX2 matrix.
---------------------------------------------------------------------------
%}

switch imageType
    case 'pit'
        % Choose histogram thresholds for solid region
        C = zeros(6,1);
        C(1,1) = 0.014;
        C(2,1) = 0.776;
        C(3,1) = 0.160;
        C(4,1) = 1.000;
        C(5,1) = 0.000;
        C(6,1) = 1.000;
        
        % Convert RGB image to HSV color space
        D = rgb2hsv(A);
        
        % Create mask based on histogram thresholds for solid region
        E = (D(:,:,1) >= C(1,1) ) & (D(:,:,1) <= C(2,1)) & ...
            (D(:,:,2) >= C(3,1) ) & (D(:,:,2) <= C(4,1)) & ...
            (D(:,:,3) >= C(5,1) ) & (D(:,:,3) <= C(6,1));
        
        % Find connected regions, select maximum, and trace boundary
        F = bwconncomp(E,4);
        numPixels = cellfun(@numel,F.PixelIdxList);
        [~,idx] = max(numPixels);
        G = zeros(size(A,1),size(A,2));
        G(F.PixelIdxList{idx}) = 1;
        [H,~,~,I] = bwboundaries(G,'holes',8);
        
        switch input1
            case 0
                %Convert RGB image to HSV color space
                Z = rgb2hsv(B);
                
                % Choose histogram thresholds for decayed (blue) region
                Y = zeros(6,1);
                Y(1,1) = 0.497;
                Y(2,1) = 0.638;
                Y(3,1) = 0.087;
                Y(4,1) = 1.000;
                Y(5,1) = 0.592;
                Y(6,1) = 1.000;
                X = (Z(:,:,1) >= Y(1,1) ) & (Z(:,:,1) <= Y(2,1)) & ...
                    (Z(:,:,2) >= Y(3,1) ) & (Z(:,:,2) <= Y(4,1)) & ...
                    (Z(:,:,3) >= Y(5,1) ) & (Z(:,:,3) <= Y(6,1));
                W = bwareaopen(X,30);
                
            case {1, 2}
                %Convert RGB image to LAB color space
                Z = rgb2lab(B);
                
                % Choose histogram thresholds for decayed (blue and purple) region
                Y = zeros(6,1);
                Y(1,1) = 20.922;
                Y(2,1) = 100.000;
                Y(3,1) = -26.902;
                Y(4,1) = 81.975;
                Y(5,1) = -63.537;
                Y(6,1) = -3.780;
                X = (Z(:,:,1) >= Y(1,1) ) & (Z(:,:,1) <= Y(2,1)) & ...
                    (Z(:,:,2) >= Y(3,1) ) & (Z(:,:,2) <= Y(4,1)) & ...
                    (Z(:,:,3) >= Y(5,1) ) & (Z(:,:,3) <= Y(6,1));
                W = bwareaopen(X,30);
                
                if input1==2
                    % Choose histogram thresholds for decayed (green) region
                    V=zeros(6,1);
                    V(1,1)=0.000;
                    V(2,1)=100.000;
                    V(3,1)=-68.602;
                    V(4,1)=-10.461;
                    V(5,1)=-99.915;
                    V(6,1)=62.944;
                    U = (Z(:,:,1) >= V(1,1) ) & (Z(:,:,1) <= V(2,1)) & ...
                        (Z(:,:,2) >= V(3,1) ) & (Z(:,:,2) <= V(4,1)) & ...
                        (Z(:,:,3) >= V(5,1) ) & (Z(:,:,3) <= V(6,1));
                    T = bwareaopen(U,30);
                    S = logical(W+T);
                    W = S;
                end
        end
        R = bwconncomp(W,6);
        numPixels = cellfun(@numel,R.PixelIdxList);
        Q = find(numPixels>100); % Select decayed areas with > 100 pixels
        P = zeros(size(B,1),size(B,2));
        for i = 1:size(Q,2)
            P(R.PixelIdxList{Q(i)}) = 1;
        end
        P(1:20,:) = 0; % Exclude tomogram colorbar from decayed areas
        P=logical(P);
        O = bwboundaries(P,'noholes',8);
        
        % Convert to registered Cartesian coordinate system
        ind1 = find(I(:,1));
        [~,ind2] = max(cellfun('length',H(ind1)));
        ind3 = ind1(ind2);
        [x1, y1] = intrinsicToWorld(R1,H{ind3,1}(:,2),H{ind3,1}(:,1)); %Solid region
        ymax = max(y1);
        y1 = ymax - y1;
        N = cell(2,1);
        N{1,1} = [x1, y1];
        N{2,1}{1,1} = [x1, y1];
        if size(O,1)==1
            [x2, y2] = intrinsicToWorld(R2,O{1,1}(:,2),O{1,1}(:,1));
            y2 = ymax - y2;
            N{2,1}{2,1} = [x2, y2];
        elseif size(O,1)>1
            for i=1:size(O,1)
                [x2, y2] = intrinsicToWorld(R2,O{i,1}(:,2),O{i,1}(:,1));
                y2 = ymax - y2;
                N{2,1}{i+1,1} = [x2, y2];
            end
        end
    case 'bin'
        C = imbinarize(rgb2gray(A));
        D = imbinarize(rgb2gray(B));
        
        % Select solid region
        % Select areas with >n pixels
        % Find connected regions, select maximum, and trace boundary
        F = bwconncomp(C,4);
        numPixels = cellfun(@numel,F.PixelIdxList);
        [~,idx] = max(numPixels);
        G = zeros(size(A,1),size(A,2));
        G(F.PixelIdxList{idx}) = 1;
        H = logical(G);
        I = bwboundaries(H,'holes',8);
        I = I(cellfun(@numel,I)>100,:);
        
        % Select decayed areas
        % Select areas with >n pixels
        J = bwconncomp(D,6);
        numPixels = cellfun(@numel,J.PixelIdxList);
        K = find(numPixels>100); % Select decayed areas with > 100 pixels
        L = zeros(size(B,1),size(B,2));
        for i = 1:size(K,2)
            L(J.PixelIdxList{K(i)}) = 1;
        end
        M = bwboundaries(L,'noholes',8);
        
        % Convert to registered Cartesian coordinate system
        N = cell(2,1);
        if size(I,1)==1
            [x1, y1] = intrinsicToWorld(R1,I{1,1}(:,2),I{1,1}(:,1)); %Solid region
            ymax = max(y1);
            y1 = ymax - y1;
            N{1,1} = [x1, y1];
            N{2,1}{1,1} = [x1, y1];
        elseif size(I,1)>1
            [~,O] = sort(cellfun(@length,I),'descend');
            I = I(O);
            for i=1:size(I,1)
                [x1,y1] = intrinsicToWorld(R1,I{i,1}(:,2),I{i,1}(:,1));
                if i==1
                    ymax = max(y1);
                end
                y1 = ymax - y1;
                N{1,1}{i,1} = [x1, y1];
                N{2,1}{i,1} = [x1, y1];
            end
        end
        
        if size(M,1)==1
            [x2, y2] = intrinsicToWorld(R2,M{1,1}(:,2),M{1,1}(:,1));
            y2 = ymax - y2;
            N{2,1}{size(I,1)+1,1} = [x2, y2];
        elseif size(M,1)>1
            for i=1:size(M,1)
                [x2, y2] = intrinsicToWorld(R2,M{i,1}(:,2),M{i,1}(:,1));
                y2 = ymax - y2;
                N{2,1}{i+size(I,1),1} = [x2, y2];
            end
        end

end
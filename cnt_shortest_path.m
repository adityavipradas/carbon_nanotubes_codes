function [mindist] = cnt_shortest_path(CNT_num, film_length, CNT_length)
% 1. Variable initialization
mindist = 0;
% -------------------------------------------------------------------------
% Get number of CNTs and thin film dimensions
% -------------------------------------------------------------------------
CNT_lengthset = CNT_length + randn(CNT_num,1); % length + variant term
% disp('Randomly creating metallic and semiconducting CNTs...');
%randomly generate CNT origins
CNT_originX = film_length * rand(CNT_num, 1);
CNT_originY = film_length * rand(CNT_num, 1);

% Cells to hold each CNT's matrix of points 
X = cell(1,CNT_num);
Y = cell(1,CNT_num);

%end point of line 1
CNT_angle1 = -90 + 180 * rand(CNT_num, 1);
CNT_endX1 = CNT_originX + (CNT_lengthset./2) .* cosd(CNT_angle1);
CNT_endY1 = CNT_originY + (CNT_lengthset./2) .* sind(CNT_angle1);

%declare arrays
CNT_endX2 = [];
CNT_endY2 = [];

%end point of line 2
CNT_angle2 = -30 + 60 * rand(CNT_num, 1);
for i=1:CNT_num
    if (CNT_angle1(i) > 0 && CNT_angle2(i) > 0) || (CNT_angle1(i) < 0 && CNT_angle2(i) < 0)
        CNT_endX2(i,1) = CNT_endX1(i) + (CNT_lengthset(i)/2) * cosd(CNT_angle1(i) + CNT_angle2(i));
        CNT_endY2(i,1) = CNT_endY1(i) + (CNT_lengthset(i)/2) * sind(CNT_angle1(i) + CNT_angle2(i));
    elseif (CNT_angle1(i) > 0 && CNT_angle2(i) < 0) || (CNT_angle1(i) < 0 && CNT_angle2(i) > 0) 
        CNT_endX2(i,1) = CNT_endX1(i) + (CNT_lengthset(i)/2) * cosd(CNT_angle1(i) - CNT_angle2(i));
        CNT_endY2(i,1) = CNT_endY1(i) + (CNT_lengthset(i)/2) * sind(CNT_angle1(i) - CNT_angle2(i));
    elseif CNT_angle2(i) == 0
        CNT_endX2(i,1) = CNT_endX1(i) + CNT_lengthset(i)/2;
        CNT_endY2(i,1) = CNT_endY1(i);
    end
end

%CNT_typeset = randint(CNT_num,1,[1,3]); % Matlab R2008a
 CNT_typeset = randi(3,CNT_num,1); % Matlab R2010
% -------------------------------------------------------------------------
% 1,2 = Semiconducting, 3 = Metallic
% -------------------------------------------------------------------------
% CHECK - If the types are equally spread out, equally probable 
% hist(CNT_typeset);
% Checked, the types are equally distributed
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Metallic CNT matrices - preprocessing
% a) Setting the mCNT counter - works only for 2/3 ratio
% b) Creating matrix containing references to metallic elements in main
%    CNT matrix
% ---------------------- ---------------------------------------------------
mCNT_count = tabulate(CNT_typeset);
mCNT_count = mCNT_count(3,2);
mCNT_ref = zeros(mCNT_count,1);
counter = 0;
for i = 1:CNT_num 
   if (CNT_typeset(i,1) == 3)
       counter = counter + 1;
       mCNT_ref(counter) = i;
   end
end
clear counter i;

if mCNT_count == 0
    return 
end

% vid = 'cnt';
% 
% %plot all CNTs
% h = figure();
% grid on;
% hold on;
% title('All CNTs(blue-metallic, red-semiconducting)');
% for i =1 :CNT_num
%     X{i} = [CNT_originX(i), CNT_endX1(i), CNT_endX2(i)];
%     Y{i} = [CNT_originY(i), CNT_endY1(i), CNT_endY2(i)];
%     if CNT_typeset(i) == 3
%         plot(X{i}, Y{i}, 'blue');
%     else
%         plot(X{i}, Y{i}, 'red');
%     end
%     hold on;
% end
% u = num2str(1);
% vname = [vid u 'jpg'];
% saveas(h, vname, 'jpg')
% set(clf,'visible','off')
% 
%slopes and intercepts of both the lines
CNT_slope1 = (CNT_endY1 - CNT_originY)./(CNT_endX1 - CNT_originX);
CNT_slope2 = (CNT_endY2 - CNT_endY1)./(CNT_endX2 - CNT_endX1);
CNT_intercept1 = CNT_endY1 - CNT_slope1 .* CNT_endX1;
CNT_intercept2 = CNT_endY2 - CNT_slope2 .* CNT_endX2;

%intersection with bottom and left edges
startX = [0];
startY = [0];
counterX = 0;
counterY = 0;

% disp('Checking intersections with domain boundaries...');
%intersections with left and bottom edges
for i = 1:CNT_num
    %left edge
    if CNT_endX1(i)<0
        CNT_endX1(i) = 0;
        CNT_endY1(i) = CNT_intercept1(i);
        CNT_endX2(i) = CNT_endX1(i);
        CNT_endY2(i) = CNT_endY1(i);
        if(CNT_typeset(i) == 3)
            counterX = counterX + 1;
            startX(counterX,1) = i;
        end
    elseif CNT_endX2(i)<0
        CNT_endX2(i) = 0;
        CNT_endY2(i) = CNT_intercept2(i);
        if(CNT_typeset(i) == 3)
            counterX = counterX + 1;
            startX(counterX,1) = i;
        end
    end
    %bottom edge
    if CNT_endY1(i)<0
        CNT_endY1(i) = 0;
        CNT_endX1(i) = -1 * CNT_intercept1(i)/CNT_slope1(i);
        CNT_endY2(i) = CNT_endY1(i);
        CNT_endX2(i) = CNT_endX1(i);
        if(CNT_typeset(i) == 3)
            counterY = counterY + 1;
            startY(counterY,1) = i;
        end
    elseif CNT_endY2(i)<0
        CNT_endY2(i) = 0;
        CNT_endX2(i) = -1 * CNT_intercept2(i)/CNT_slope2(i);
        if(CNT_typeset(i) == 3)
            counterY = counterY + 1;
            startY(counterY,1) = i;
        end
    end
end

if length(startX) == 0 || length(startY) == 0
    return
end

%intersection with top and right edges
endX = [0];
endY = [0];
counterX2 = 0;
counterY2 = 0;

for i = 1:CNT_num
    %right edge
    if CNT_endX1(i) > film_length
        CNT_endX1(i) = film_length;
        CNT_endY1(i) = CNT_slope1(i)*CNT_endX1(i) + CNT_intercept1(i);
        CNT_endY2(i) = CNT_endY1(i);
        CNT_endX2(i) = CNT_endX1(i);
        if(CNT_typeset(i) == 3)
            counterX2 = counterX2 + 1;
            endX(counterX2, 1) = i;
        end
    elseif CNT_endX2(i) > film_length
        CNT_endX2(i) = film_length;
        CNT_endY2(i) = CNT_slope2(i)*CNT_endX2(i) + CNT_intercept2(i);
        if(CNT_typeset(i) == 3)
            counterX2 = counterX2 + 1;
            endX(counterX2, 1) = i;
        end
    end
    if CNT_endY1(i) > film_length
        CNT_endY1(i) = film_length;
        CNT_endX1(i) = (CNT_endY1(i) - CNT_intercept1(i))/CNT_slope1(i);
        CNT_endY2(i) = CNT_endY1(i);
        CNT_endX2(i) = CNT_endX1(i);
        if CNT_typeset(i) == 3
            counterY2 = counterY2 + 1;
            endY(counterY2, 1) = i;
        end
    elseif CNT_endY2(i) > film_length
        CNT_endY2(i) = film_length;
        CNT_endX2(i) = (CNT_endY2(i) - CNT_intercept2(i))/CNT_slope2(i);
        if CNT_typeset(i) == 3
            counterY2 = counterY2 + 1;
            endY(counterY2, 1) = i;
        end
    end
end

if length(endX) == 0 || length(endY) == 0
    return
end

%intersection free plot
Xn = cell(1, CNT_num);
Yn = cell(1, CNT_num);

% h1 = figure();
% hold on;
% grid on;
% title('All CNTs within domain(blue-metallic, red-semiconducting)');
 for i =1 :CNT_num
     Xn{i} = [CNT_originX(i), CNT_endX1(i), CNT_endX2(i)];
     Yn{i} = [CNT_originY(i), CNT_endY1(i), CNT_endY2(i)];
%     if CNT_typeset(i) == 3
%         plot(Xn{i}, Yn{i}, 'blue');
%     else
%         plot(Xn{i}, Yn{i}, 'red');
%     end
%     hold on;
 end
% u = num2str(2);
% vname = [vid u 'jpg'];
% saveas(h1, vname, 'jpg')
% set(clf,'visible','off')

%metallic intersection free plot
% h2 = figure();
% hold on;
% grid on;
% title('All metallic CNTs(with edge numbers)'); 
% for i = 1: mCNT_count
%     plot(Xn{mCNT_ref(i)}, Yn{mCNT_ref(i)}, 'blue');
%     text(CNT_originX(mCNT_ref(i)), CNT_originY(mCNT_ref(i)), num2str(mCNT_ref(i)), 'FontSize', 8);
%     hold on;
% end
% u = num2str(3);
% vname = [vid u 'jpg'];
% saveas(h2, vname, 'jpg')
% set(clf,'visible','off')
% 
% disp('Searching all the intersecting points...'); 
%find all CNT intersections
CNT_map = zeros([mCNT_count, mCNT_count]);
tempx = cell(mCNT_count, mCNT_count);
tempy = cell(mCNT_count, mCNT_count);
for i = 1:mCNT_count
    i1Xo = min(Xn{mCNT_ref(i)}(1,1), Xn{mCNT_ref(i)}(1,2));
    i1Xe = max(Xn{mCNT_ref(i)}(1,1), Xn{mCNT_ref(i)}(1,2));
    i2Xo = min(Xn{mCNT_ref(i)}(1,2), Xn{mCNT_ref(i)}(1,3));
    i2Xe = max(Xn{mCNT_ref(i)}(1,2), Xn{mCNT_ref(i)}(1,3));
    for j = 1:mCNT_count
        j1Xo = min(Xn{mCNT_ref(j)}(1,1), Xn{mCNT_ref(j)}(1,2));
        j1Xe = max(Xn{mCNT_ref(j)}(1,1), Xn{mCNT_ref(j)}(1,2));
        j2Xo = min(Xn{mCNT_ref(j)}(1,2), Xn{mCNT_ref(j)}(1,3));
        j2Xe = max(Xn{mCNT_ref(j)}(1,2), Xn{mCNT_ref(j)}(1,3));
        if (i~=j && CNT_slope1(mCNT_ref(i))~=CNT_slope1(mCNT_ref(j)))
            intersect = (CNT_intercept1(mCNT_ref(i)) - CNT_intercept1(mCNT_ref(j)))/(CNT_slope1(mCNT_ref(j)) - CNT_slope1(mCNT_ref(i)));
            if intersect >= i1Xo && intersect >= j1Xo && intersect <= i1Xe && intersect <= j1Xe
                tempx{i,j}(1,1) = intersect;
                tempy{i,j}(1,1) = CNT_slope1(mCNT_ref(i))*intersect + CNT_intercept1(mCNT_ref(i));
                CNT_map(i,j) = mCNT_ref(i);
            else
                tempx{i,j}(1,1) = NaN;
                tempy{i,j}(1,1) = NaN;
            end
        else
                tempx{i,j}(1,1) = NaN;
                tempy{i,j}(1,1) = NaN;
        end
        if (i~=j && CNT_endX1(mCNT_ref(j))~= CNT_endX2(mCNT_ref(j)) && CNT_endY1(mCNT_ref(j))~= CNT_endY2(mCNT_ref(j)) && CNT_slope1(mCNT_ref(i))~=CNT_slope2(mCNT_ref(j)))
            intersect = (CNT_intercept1(mCNT_ref(i)) - CNT_intercept2(mCNT_ref(j)))/(CNT_slope2(mCNT_ref(j)) - CNT_slope1(mCNT_ref(i)));
            if intersect >= i1Xo && intersect >= j2Xo && intersect <= i1Xe && intersect <= j2Xe
                tempx{i,j}(1,2) = intersect;
                tempy{i,j}(1,2) = CNT_slope1(mCNT_ref(i))*intersect + CNT_intercept1(mCNT_ref(i));
                CNT_map(i,j) = mCNT_ref(i);
            else
                tempx{i,j}(1,2) = NaN;
                tempy{i,j}(1,2) = NaN;
            end
        else
            tempx{i,j}(1,2) = NaN;
            tempy{i,j}(1,2) = NaN;
        end
        if (i~=j && CNT_endX1(mCNT_ref(i))~= CNT_endX2(mCNT_ref(i)) && CNT_endY1(mCNT_ref(i))~= CNT_endY2(mCNT_ref(i)) && CNT_slope2(mCNT_ref(i))~=CNT_slope1(mCNT_ref(j)))
            intersect = (CNT_intercept2(mCNT_ref(i)) - CNT_intercept1(mCNT_ref(j)))/(CNT_slope1(mCNT_ref(j)) - CNT_slope2(mCNT_ref(i)));
            if intersect >= i2Xo && intersect >= j1Xo && intersect <= i2Xe && intersect <= j1Xe
                tempx{i,j}(1,3) = intersect;
                tempy{i,j}(1,3) = CNT_slope2(mCNT_ref(i))*intersect + CNT_intercept2(mCNT_ref(i));
                CNT_map(i,j) = mCNT_ref(i);
            else
                tempx{i,j}(1,3) = NaN;
                tempy{i,j}(1,3) = NaN;
            end
        else
            tempx{i,j}(1,3) = NaN;
            tempy{i,j}(1,3) = NaN;
        end
        if (i~=j && CNT_endX1(mCNT_ref(j))~= CNT_endX2(mCNT_ref(j)) && CNT_endY1(mCNT_ref(j))~= CNT_endY2(mCNT_ref(j)) && CNT_endX1(mCNT_ref(i))~= CNT_endX2(mCNT_ref(i)) && CNT_endY1(mCNT_ref(i))~= CNT_endY2(mCNT_ref(i)) && CNT_slope2(mCNT_ref(i))~=CNT_slope2(mCNT_ref(j)))
            intersect = (CNT_intercept2(mCNT_ref(i)) - CNT_intercept2(mCNT_ref(j)))/(CNT_slope2(mCNT_ref(j)) - CNT_slope2(mCNT_ref(i)));
            if intersect >= i2Xo && intersect >= j2Xo && intersect <= i2Xe && intersect <= j2Xe
                tempx{i,j}(1,4) = intersect;
                tempy{i,j}(1,4) = CNT_slope2(mCNT_ref(i))*intersect + CNT_intercept2(mCNT_ref(i));
                CNT_map(i,j) = mCNT_ref(i);
            else
                tempx{i,j}(1,4) = NaN;
                tempy{i,j}(1,4) = NaN;
            end
        else
            tempx{i,j}(1,4) = NaN;
            tempy{i,j}(1,4) = NaN;
        end
    end
end

if length(CNT_map) == 0
    return
end

% disp(tempx);
% disp(tempy);
% if isempty(CNT_map) == 1
%     error('No connection with top or bottom boundaries');
% end
[A,B] = find(CNT_map);
CNT_map = [B A];
% disp(CNT_map);

if length(CNT_map) == 0
    return
end

% disp('Removing duplicate entries...');
%remove duplicates
pair_count = length(CNT_map);
for i = 1:pair_count
    if CNT_map(i,1) ~= 0 || CNT_map(i,2) ~= 0
        pair_1 = CNT_map(i,1);
        pair_2 = CNT_map(i,2);
        l1 = find(CNT_map(:,1) == pair_2)';
        l2 = find(CNT_map(:,2) == pair_1)';
        for k = 1:length(l1)
            for l = 1:length(l2)
                if(l1(k) == l2(l))
                    duplicate = l1(k);
                end
            end
        end
        CNT_map(duplicate,1) = 0; %Zero this duplicate
        CNT_map(duplicate,2) = 0;
    end
end
% disp(CNT_map);

if length(CNT_map) == 0
    return
end

%disp('Removing all dead-ends...');
%dead end removal
dead_check = 0;
while dead_check == 0
    dead_check = 1;
    for i = 1:mCNT_count
        counter = [find(CNT_map(:,1) == i); find(CNT_map(:,2) == i)];
        if length(counter) == 1
            if isempty(find(endY == mCNT_ref(i))) == 1 && isempty(find(startY == mCNT_ref(i)))== 1 
                coord = counter(1);
                CNT_map(coord,1) = 0;
                CNT_map(coord,2) = 0;
                dead_check = 0;
            end
        end
    end
%     disp(CNT_map);
end


if length(CNT_map) == 0
    return
end

% disp(CNT_map);

%plot paths
% h3 = figure();
% hold on;
% grid on;
% title('CNTs with dead ends removed(with edge numbers)');
% for i = 1:length(CNT_map)
%     if CNT_map(i,1) ~= 0 || CNT_map(i,2) ~= 0
%         plot(Xn{mCNT_ref(CNT_map(i,1))}, Yn{mCNT_ref(CNT_map(i,1))});
%         hold on;
%         plot(Xn{mCNT_ref(CNT_map(i,2))}, Yn{mCNT_ref(CNT_map(i,2))});
%         hold on;
%         text(CNT_originX(mCNT_ref(CNT_map(i,1))), CNT_originY(mCNT_ref(CNT_map(i,1))), num2str(mCNT_ref(CNT_map(i,1))), 'FontSize', 8);
%         text(CNT_originX(mCNT_ref(CNT_map(i,2))), CNT_originY(mCNT_ref(CNT_map(i,2))), num2str(mCNT_ref(CNT_map(i,2))), 'FontSize', 8);
%     end
% end
% u = num2str(4);
% vname = [vid u 'jpg'];
% saveas(h3, vname, 'jpg')
% set(clf,'visible','off')

%disp('Creating master lists...');
%master list of all non-dead end intersections
node = 0;
counter = 0;
MasterList = [];
CNT_map_clean = [];
for i = 1:length(CNT_map)
    if CNT_map(i,1) ~= 0
        node = node + 1;
        counter = counter + 1;
        MasterList(counter, 1) = mCNT_ref(CNT_map(i,1));
        MasterList(counter, 2) = mCNT_ref(CNT_map(i,2));
        for j = 1:4
            if isnan(tempx{CNT_map(i,1), CNT_map(i,2)}(1,j)) == 0
                MasterList(counter, 3) = tempx{CNT_map(i,1), CNT_map(i,2)}(1,j);
                MasterList(counter, 4) = tempy{CNT_map(i,1), CNT_map(i,2)}(1,j);
            end
        end
        MasterList(counter, 5) = node;
        CNT_map_clean(counter, 1) = CNT_map(i,1);
        CNT_map_clean(counter, 2) = CNT_map(i,2);
    end
end
CNT_map = CNT_map_clean;
pair_count = length(CNT_map);
% disp(CNT_map);
% disp(MasterList);

if length(MasterList) == 0
    return 
end

if length(CNT_map) == 0
    return
end

if length(endY) == 0
    return
end

%top and bottom boundaries representations
Master_Top = [];
counter = 0;
for i = 1:length(endY)
    if endY(i) == 0
        return
    end
    if CNT_endY2(endY(i))==film_length && (any_equal(MasterList(:,1), endY(i)) == 1 || any_equal(MasterList(:,2), endY(i)) == 1)
        node = node + 1;
        counter = counter + 1;
        Master_Top(counter,1) = endY(i);
        Master_Top(counter,2) = -2;
        Master_Top(counter,3) = CNT_endX2(endY(i));
        Master_Top(counter,4) = CNT_endY2(endY(i));
        Master_Top(counter,5) = node;
    end
end

if length(Master_Top) == 0
    return 
end

if length(startY) == 0
    return
end

Master_Bottom = [];
counter = 0;
for i = 1:length(startY)
    if startY(i) == 0
        return
    end
    if CNT_endY2(startY(i))==0 && (any_equal(MasterList(:,1), startY(i)) == 1 || any_equal(MasterList(:,2), startY(i)) == 1)
        node = node + 1;
        counter = counter + 1;
        Master_Bottom(counter,1) = startY(i);
        Master_Bottom(counter,2) = -1;
        Master_Bottom(counter,3) = CNT_endX2(startY(i));
        Master_Bottom(counter,4) = CNT_endY2(startY(i));
        Master_Bottom(counter,5) = node;
    end
end

% if isempty(Master_Top) == 1 || isempty(Master_Bottom) == 1 || isempty(MasterList) == 1
%     error('No connection with top or bottom boundaries');
% end

%plot with node numbers
% h4 = figure();
% hold on;
% grid on;
% title('All nodes(with node numbers)');
% for i = 1:length(MasterList)
%     plot(Xn{MasterList(i,1)}, Yn{MasterList(i,1)}, 'blue');
%     plot(Xn{MasterList(i,2)}, Yn{MasterList(i,2)}, 'blue');
%     text(MasterList(i,3), MasterList(i,4), num2str(MasterList(i,5)), 'FontSize', 8);
%     if any_equal(Master_Bottom(:,1), MasterList(i,1)) == 1
%         [seek, col] = find(Master_Bottom(:,1) == MasterList(i,1));
%         text(Master_Bottom(seek, 3), Master_Bottom(seek, 4), num2str(Master_Bottom(seek, 5)), 'FontSize', 8);
%     end
%     if any_equal(Master_Bottom(:,1), MasterList(i,2)) == 1
%         [seek, col] = find(Master_Bottom(:,1) == MasterList(i,2));
%         text(Master_Bottom(seek, 3), Master_Bottom(seek, 4), num2str(Master_Bottom(seek, 5)), 'FontSize', 8);
%     end
%     if any_equal(Master_Top(:,1), MasterList(i,1)) == 1
%         [seek, col] = find(Master_Top(:,1) == MasterList(i,1));
%         text(Master_Top(seek, 3), Master_Top(seek, 4), num2str(Master_Top(seek, 5)), 'FontSize', 8);
%     end
%     if any_equal(Master_Top(:,1), MasterList(i,2)) == 1
%         [seek, col] = find(Master_Top(:,1) == MasterList(i,2));
%         text(Master_Top(seek, 3), Master_Top(seek, 4), num2str(Master_Top(seek, 5)), 'FontSize', 8);
%     end
%     hold on;
% end
% % disp(Master_Top);
% % disp(Master_Bottom);
% u = num2str(5);
% vname = [vid u 'jpg'];
% saveas(h4, vname, 'jpg')
% set(clf,'visible','off')
% 
% disp('Making sparse matrix for Dijkstra algorithm implementation...');
%form sparse matrix
vect1 = [];
vect2 = [];
cnt = 0;
W = [];

if length(Master_Bottom) == 0
    return
end

numRow = size(MasterList);
%disp(MasterList);
%disp(numRow(1));
%disp('---------------------');
    
for i = 1:numRow(1)
    if any_equal(Master_Bottom(:,1), MasterList(i,1)) == 1
        cnt = cnt + 1;
        a = MasterList(i,3);
        b = MasterList(i,4);
        [search col] = find(Master_Bottom(:,1) == MasterList(i,1));
        c = Master_Bottom(search, 3);
        d = Master_Bottom(search, 4);
        if abs(a - CNT_originX(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_originX(MasterList(i,1))) && abs(b - CNT_originY(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_originY(MasterList(i,1)))
            weight = sqrt((a - CNT_endX1(MasterList(i,1)))^2 + (b - CNT_endY1(MasterList(i,1)))^2) + sqrt((c - CNT_endX1(MasterList(i,1)))^2 + (d - CNT_endY1(MasterList(i,1)))^2);
        elseif abs(a - CNT_endX1(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_endX2(MasterList(i,1))) && abs(b - CNT_endY1(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_endY2(MasterList(i,1)))
            weight = sqrt((a-c)^2 + (b-d)^2);
        end
        vect1(1,cnt) = MasterList(i,5);
        vect2(1,cnt) = Master_Bottom(search,5);
        W(1,cnt) = weight;
        
        cnt = cnt + 1;
        vect2(1,cnt) = MasterList(i,5);
        vect1(1,cnt) = Master_Bottom(search,5);
        W(1,cnt) = weight;
    end
    
    if any_equal(Master_Bottom(:,1), MasterList(i,2)) == 1
        cnt = cnt + 1;
        a = MasterList(i,3);
        b = MasterList(i,4);
        [search col] = find(Master_Bottom(:,1) == MasterList(i,2));
        c = Master_Bottom(search, 3);
        d = Master_Bottom(search, 4);
        if abs(a - CNT_originX(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_originX(MasterList(i,2))) && abs(b - CNT_originY(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_originY(MasterList(i,2)))
            weight = sqrt((a - CNT_endX1(MasterList(i,2)))^2 + (b - CNT_endY1(MasterList(i,2)))^2) + sqrt((c - CNT_endX1(MasterList(i,2)))^2 + (d - CNT_endY1(MasterList(i,2)))^2);
        elseif abs(a - CNT_endX1(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_endX2(MasterList(i,2))) && abs(b - CNT_endY1(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_endY2(MasterList(i,2)))
            weight = sqrt((a-c)^2 + (b-d)^2);
        end
        vect1(1,cnt) = MasterList(i,5);
        vect2(1,cnt) = Master_Bottom(search,5);
        W(1,cnt) = weight;
        
        cnt = cnt + 1;
        vect2(1,cnt) = MasterList(i,5);
        vect1(1,cnt) = Master_Bottom(search,5);
        W(1,cnt) = weight;
    end
    
    if any_equal(Master_Top(:,1), MasterList(i,1)) == 1
        cnt = cnt + 1;
        a = MasterList(i,3);
        b = MasterList(i,4);
        [search col] = find(Master_Top(:,1) == MasterList(i,1));
        c = Master_Top(search, 3);
        d = Master_Top(search, 4);
        if abs(a - CNT_originX(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_originX(MasterList(i,1))) && abs(b - CNT_originY(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_originY(MasterList(i,1)))
            weight = sqrt((a - CNT_endX1(MasterList(i,1)))^2 + (b - CNT_endY1(MasterList(i,1)))^2) + sqrt((c - CNT_endX1(MasterList(i,1)))^2 + (d - CNT_endY1(MasterList(i,1)))^2);
        elseif abs(a - CNT_endX1(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_endX2(MasterList(i,1))) && abs(b - CNT_endY1(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_endY2(MasterList(i,1)))
            weight = sqrt((a-c)^2 + (b-d)^2);
        end
        vect1(1,cnt) = MasterList(i,5);
        vect2(1,cnt) = Master_Top(search,5);
        W(1,cnt) = weight;
        
        cnt = cnt + 1;
        vect2(1,cnt) = MasterList(i,5);
        vect1(1,cnt) = Master_Top(search,5);
        W(1,cnt) = weight;
    end
    
    if any_equal(Master_Top(:,1), MasterList(i,2)) == 1
        cnt = cnt + 1;
        a = MasterList(i,3);
        b = MasterList(i,4);
        search = find(Master_Top(:,1) == MasterList(i,2));
        c = Master_Top(search, 3);
        d = Master_Top(search, 4);
        if abs(a - CNT_originX(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_originX(MasterList(i,2))) && abs(b - CNT_originY(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_originY(MasterList(i,2)))
            weight = sqrt((a - CNT_endX1(MasterList(i,2)))^2 + (b - CNT_endY1(MasterList(i,2)))^2) + sqrt((c - CNT_endX1(MasterList(i,2)))^2 + (d - CNT_endY1(MasterList(i,2)))^2);
        elseif abs(a - CNT_endX1(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_endX2(MasterList(i,2))) && abs(b - CNT_endY1(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_endY2(MasterList(i,2)))
            weight = sqrt((a-c)^2 + (b-d)^2);
        end
        vect1(1,cnt) = MasterList(i,5);
        vect2(1,cnt) = Master_Top(search,5);
        W(1,cnt) = weight;
        
        cnt = cnt + 1;
        vect2(1,cnt) = MasterList(i,5);
        vect1(1,cnt) = Master_Top(search,5);
        W(1,cnt) = weight;
    end
    
    for j = 1:numRow(1)
        if i~=j && any_equal([MasterList(j,1) MasterList(j,2)], MasterList(i,1))==1
            cnt = cnt + 1;
            a = MasterList(i,3);
            b = MasterList(i,4);
            c = MasterList(j,3);
            d = MasterList(j,4);                    
            if abs(a - CNT_originX(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_originX(MasterList(i,1))) && abs(c - CNT_originX(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_originX(MasterList(i,1))) && abs(b - CNT_originY(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_originY(MasterList(i,1))) && abs(d - CNT_originY(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_originY(MasterList(i,1)))
                weight = sqrt((a-c)^2 + (b-d)^2);
            elseif abs(a - CNT_endX1(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_endX2(MasterList(i,1))) && abs(c - CNT_endX1(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_endX2(MasterList(i,1))) && abs(b - CNT_endY1(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_endY2(MasterList(i,1))) && abs(d - CNT_endY1(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_endY2(MasterList(i,1)))
                weight = sqrt((a-c)^2 + (b-d)^2);
            elseif abs(a - CNT_endX1(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_endX2(MasterList(i,1))) && abs(c - CNT_originX(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_originX(MasterList(i,1))) && abs(b - CNT_endY1(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_endY2(MasterList(i,1))) && abs(d - CNT_originY(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_originY(MasterList(i,1)))
                weight = sqrt((a - CNT_endX1(MasterList(i,1)))^2 + (b - CNT_endY1(MasterList(i,1)))^2) + sqrt((c - CNT_endX1(MasterList(i,1)))^2 + (d - CNT_endY1(MasterList(i,1)))^2);
            elseif abs(c - CNT_endX1(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_endX2(MasterList(i,1))) && abs(a - CNT_originX(MasterList(i,1))) < abs(CNT_endX1(MasterList(i,1)) - CNT_originX(MasterList(i,1))) && abs(d - CNT_endY1(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_endY2(MasterList(i,1))) && abs(b - CNT_originY(MasterList(i,1))) < abs(CNT_endY1(MasterList(i,1)) - CNT_originY(MasterList(i,1)))
                weight = sqrt((a - CNT_endX1(MasterList(i,1)))^2 + (b - CNT_endY1(MasterList(i,1)))^2) + sqrt((c - CNT_endX1(MasterList(i,1)))^2 + (d - CNT_endY1(MasterList(i,1)))^2);
            end
            vect1(1,cnt) = MasterList(i,5);
            vect2(1,cnt) = MasterList(j,5);
            W(1,cnt) = weight;
        elseif i~=j && any_equal([MasterList(j,1) MasterList(j,2)], MasterList(i,2))==1
            cnt = cnt + 1;
            a = MasterList(i,3);
            b = MasterList(i,4);
            c = MasterList(j,3);
            d = MasterList(j,4);                          
            if abs(a - CNT_originX(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_originX(MasterList(i,2))) && abs(c - CNT_originX(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_originX(MasterList(i,2))) && abs(b - CNT_originY(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_originY(MasterList(i,2))) && abs(d - CNT_originY(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_originY(MasterList(i,2)))
                weight = sqrt((a-c)^2 + (b-d)^2);
            elseif abs(a - CNT_endX1(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_endX2(MasterList(i,2))) && abs(c - CNT_endX1(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_endX2(MasterList(i,2))) && abs(b - CNT_endY1(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_endY2(MasterList(i,2))) && abs(d - CNT_endY1(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_endY2(MasterList(i,2)))
                weight = sqrt((a-c)^2 + (b-d)^2);
            elseif abs(a - CNT_endX1(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_endX2(MasterList(i,2))) && abs(c - CNT_originX(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_originX(MasterList(i,2))) && abs(b - CNT_endY1(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_endY2(MasterList(i,2))) && abs(d - CNT_originY(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_originY(MasterList(i,2)))
                weight = sqrt((a - CNT_endX1(MasterList(i,2)))^2 + (b - CNT_endY1(MasterList(i,2)))^2) + sqrt((c - CNT_endX1(MasterList(i,2)))^2 + (d - CNT_endY1(MasterList(i,2)))^2);
            elseif abs(c - CNT_endX1(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_endX2(MasterList(i,2))) && abs(a - CNT_originX(MasterList(i,2))) < abs(CNT_endX1(MasterList(i,2)) - CNT_originX(MasterList(i,2))) && abs(d - CNT_endY1(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_endY2(MasterList(i,2))) && abs(b - CNT_originY(MasterList(i,2))) < abs(CNT_endY1(MasterList(i,2)) - CNT_originY(MasterList(i,2)))
                weight = sqrt((a - CNT_endX1(MasterList(i,2)))^2 + (b - CNT_endY1(MasterList(i,2)))^2) + sqrt((c - CNT_endX1(MasterList(i,2)))^2 + (d - CNT_endY1(MasterList(i,2)))^2);
            end
            vect1(1,cnt) = MasterList(i,5);
            vect2(1,cnt) = MasterList(j,5);
            W(1,cnt) = weight;
        end
    end    
end

mat = sparse(vect1, vect2, W);
%bio = view(biograph(mat,[],'ShowWeights','on'));

% disp('Plotting all the paths...');
% mkdir('paths');
% cd('paths');
% vid1 = 'path';
% plotnum = 0;

allpaths = [];
pathtrace = cell([],1);
count = 0;
for i = 1:length(Master_Bottom(:,1))
    for j = 1:length(Master_Top(:,1))
        [dist, path, pred] = graphshortestpath(mat, Master_Bottom(i,5), Master_Top(j,5), 'Directed', 'true', 'Method', 'Dijkstra');
        if dist ~= Inf
            count = count + 1;
%             plotnum = plotnum + 1;
            allpaths(count, 1) = Master_Bottom(i,5);
            allpaths(count, 2) = Master_Top(j,5);
            allpaths(count, 3) = dist;
            pathtrace{count,1} = path;
%             h5 = figure();
%             title(sprintf('Path %d',count));
%             xlim([0, film_length]);
%             ylim([0, film_length]);
%             grid on;
%             hold on;
%             for k = 1:length(path)
%                 if k == 1
%                     [seek col] = find(Master_Bottom(:,5) == path(k));
%                     plot(Xn{Master_Bottom(seek,1)}, Yn{Master_Bottom(seek,1)});
%                     text(Master_Bottom(seek, 3), Master_Bottom(seek, 4), num2str(Master_Bottom(seek, 5)), 'FontSize', 8);
%                     hold on;
%                 end
%                 if k == length(path)
%                     [seek col] = find(Master_Top(:,5) == path(k));
%                     plot(Xn{Master_Top(seek,1)}, Yn{Master_Top(seek,1)});
%                     text(Master_Top(seek, 3), Master_Top(seek, 4), num2str(Master_Top(seek, 5)), 'FontSize', 8);
%                     hold on;
%                 end
%                 if k > 1 && k < length(path)-1
%                     m = [];
%                     n = [];
%                     [seek1 col] = find(MasterList(:,5) == path(k));
%                     m = [MasterList(seek1, 1), MasterList(seek1, 2)];
%                     [seek2 col] = find(MasterList(:,5) == path(k+1));
%                     n = [MasterList(seek2, 1), MasterList(seek2, 2)];
%                     for p = 1:2
%                         for q = 1:2
%                             if m(p) == n(q)
%                                 store = m(p);
%                             end
%                         end
%                     end
%                     plot(Xn{store}, Yn{store});
%                     text(MasterList(seek1, 3), MasterList(seek1, 4), num2str(MasterList(seek1, 5)), 'FontSize', 8);
%                     text(MasterList(seek2, 3), MasterList(seek2, 4), num2str(MasterList(seek2, 5)), 'FontSize', 8);                    
%                 end
%             end
%             u = num2str(plotnum);
%             vname = [vid1 u 'jpg'];
%             saveas(h5, vname, 'jpg');
%             set(clf, 'visible', 'off');
        end
    end
end

if length(allpaths) == 0
    return
end

% mkdir('shortest_path');
mindist = min(allpaths(:,3));
% minrow = find(allpaths(:,3) == min(allpaths(:,3)));
% filename = strcat('path', num2str(minrow), 'jpg.jpg');
% copyfile(filename, 'shortest_path\');

%fprintf('Distance of shortest path is %f units \n', mindist);
%disp('Done');



function tatti()
    a = textread('matrix1.txt','%s');
    r = [];
    c = [];
    v = [];
    format long
    for i=1:3:length(a)
        r = [r , str2num(char(a(i)))];       
         c = [c , str2num(char(a(i+1)))];      
         v = [v , str2double(char(a(i+2)))];
    end
    lhs = sparse(r,c,v);
    save('matrix1.mat', 'lhs')
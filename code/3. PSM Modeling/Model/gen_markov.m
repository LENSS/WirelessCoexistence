clear all;
s = 2; %packet size (slots)

w = 3; %contention window size (slots)
a = 5; %ATIM window size (slots)
n = 3; %number of nodes besides myself

fprintf('backing off states.\n');

for t=0:1
    for j=a:-1:1
        for i=n:-1:1
            if j<10
                if j-s>=0
                    fprintf('-1,%d ,%d   \t', j-s,i);
                else
                    fprintf('            \t', j-s,i);
                end
            else
                if j-s>=0
                    fprintf('-1,%d,%d   \t', j-s,i);
                else
                    fprintf('            \t', j-s,i);
                end
            end
            for k=0:w-1
                if k~=w-1
                    if j<10
                        fprintf('%d,%d,%d ,%d  \t', t, k,j,i);
                    else
                        fprintf('%d,%d,%d,%d  \t', t, k,j,i);
                    end
                else
                    if j<10
                        fprintf('%d,%d,%d ,%d', t, k,j,i);
                    else
                        fprintf('%d,%d,%d,%d', t, k,j,i);
                    end
                end
    %             fprintf('0,%d,%d   \t\t', k,j);
            end
            fprintf('\n');
        end
    end
    fprintf('\n');
    w=w*2-1;
end

fprintf('waiting for window ends states.\n');

for j=a:-1:1
    for i=n:-1:1
        
        for k=1:6
            if j<10
                fprintf('-1,%d ,%d   \t', j,n-i+1);
            else
                fprintf('-1,%d,%d   \t', j,n-i+1);
            end
        end

        
%             fprintf('0,%d,%d   \t\t', k,j);
        fprintf('\n');
    end
end


fprintf('tx states during the waiting for queue becoming not empty.\n');

for j=a:-1:0
    for i=n:-1:0
        
        for k=1:6
            if j<10
                fprintf('-3,%d ,%d   \t', j,i);
            else
                fprintf('-3,%d,%d   \t', j,i);
            end
        end

        
%             fprintf('0,%d,%d   \t\t', k,j);
        fprintf('\n');
    end
end


% syms a b c p1 p2 p3 p4 p5 T22 T21 T20 T11 T10 x1 x2 x3 x4 x5 x6 Pw;
% p1=Pw*Pw;
% p2=Pw*(1-Pw);
% p3=Pw^3;
% p4=2*Pw^2*(1-Pw);
% p5=Pw*(1-Pw)^2;
% 
% [x1, x2, x3, x4, x5, x6] = solve(...
%                                 x1-x4-x5*p2-x6*p4==a, ...
%                                 x2-x5*p1-x6*p4==b,...
%                                 x3-x6*p3==c,...
%                                 x4-x1*T22==0,...
%                                 x5-x1*T21-x2*T11==0,...
%                                 x6-x1*T20-x2*T10-x3==0,...
%                                 x1, x2, x3, x4, x5, x6);
% 
% disp(simplify(x1));
% disp(simplify(x2));
% disp(simplify(x3));
% disp(simplify(x4));
% disp(simplify(x5));
% disp(simplify(x6));







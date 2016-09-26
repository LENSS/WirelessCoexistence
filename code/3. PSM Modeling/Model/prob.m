
clear all;
syms p;
syms tau;

s = 4; %packet size (slots)
w0 = 5; %contention window size (slots)

k = 10; %ATIM window size (slots)
nn = 1; %meanningless here
no = 4; %nodes number, including myself

% k = w+s+1; %for computing purpose 
%(i,j,k) represents one state in Markov Chain


%stage 0


i = 0; %counter stage 
j = w0*(i+1); %contention window size (slots)

A=zeros(j,k*nn); %column, row
A = sym(A);

B=zeros(j,(s)*nn); %column, row
B = sym(B);



% A(:,:,2)=0;
% A(:,:,3)=0;

% for n=1:j*k
%     
%     A(n) = ps+n;
%     
% end
% 
% 
% disp(A(5,14));

init = 1; %init state stationary prob.


% for q=j:-1:1
%     At=zeros(j,k*nn); %column, row
%     At = sym(At);
%     At(q,1) = init/w0;
%     
%     Bt=zeros(j,(s)*nn); %column, row
%     Bt = sym(Bt);
    
for m=j:-1:1
    A(m,1) = init/j;
    for n=1:k*nn
%             disp(m);
%             disp(n);
        if A(m,n)~=0
            if m>1 
%                     disp('m>1');
                if n+1*nn<=k*nn
                    A(m-1,n+1*nn) = A(m,n)*(1-p)+A(m-1,n+1*nn);
                end          
                if n+s*nn<=k*nn
                    A(m-1,n+s*nn) = A(m,n)*p+A(m-1,n+s*nn);
                end

                if n+1*nn>k*nn
                    B(m-1,n+1*nn-k*nn) = A(m,n)*(1-p)+B(m-1,n+1*nn-k*nn);
%                         disp(Bt(m-1,n+1-k));
                end

                if n+s*nn>k*nn
                    B(m-1,n+s*nn-k*nn) = A(m,n)*p+B(m-1,n+s*nn-k*nn);
%                         disp(Bt(m-1,n+s-k));
                end
             end
        end
    end
end

%     At.'
%     Bt.'
%     A=A+At;
%     B=B+Bt;
% end
A = simplify(A);
B = simplify(B);


pr=zeros(1,j); %column, row
pr = sym(pr);

for q=1:j
        
    for n=1:(s)*nn
        pr(1,q) = B(q,n) + pr(1,q);
    end
%     pr(1,q) = simplify(pr(1,q));
end
pr = simplify(pr);

% Asum=zeros(j,1); %column, row
% Asum = sym(Asum);
% for q=1:j
%         
%     for n=1:(k)*nn
%         Asum(q,1) = A(q,n) + Asum(q,1);
%     end           
%     Asum(q,1) = simplify(Asum(q,1));
% end
% Asum = Asum.';%simple transpose, not conjugate transpose
% disp('Asum:');
% disp(Asum);




% result=zeros(j,1); %column, row
% result = sym(result);



result = [pr(j) pr(1:(j-1))];
for q=1:j
%     Asum=zeros(j,1); %column, row
%     Asum = sym(Asum);
% 
%     for n=1:(k)*nn
%         Asum(q,1) = A(q,n) + Asum(q,1);
%     end           
%     Asum(q,1) = simplify(Asum(q,1));
    result(q) = result(q)/((j+1-q)/j-pr(q));
end

disp('A:');
disp(A.');%simple transpose, not conjugate transpose


disp('B:');
disp(B.');

disp('pr:');
disp(pr);

disp('result:');
disp(result);

tau = 0.2;
p = 1-(1-tau)^(no-1);
disp(eval(result));

syms tau;
syms p;

At=A(1,:);
% disp(At.');
% return;

%stage 1 to m

for i = 1:1 %counter stage 
    j = w0*(i+1); %contention window size (slots)

    A=zeros(j,k*nn); %column, row
    A = sym(A);

    B=zeros(j,(s)*nn); %column, row
    B = sym(B);

    for m=j:-1:1
        for n=1:k*nn
            if At(1,n)~=0
                if n+s*nn<=k*nn
                    A(m,n+s*nn) = At(1,n)*p/j+A(m,n+s*nn);
                end
            end
        end
    end
    % A.'



    % return;


    % A(:,:,2)=0;
    % A(:,:,3)=0;

    % for n=1:j*k
    %     
    %     A(n) = ps+n;
    %     
    % end
    % 
    % 
    % disp(A(5,14));

    % init = 1; %init state stationary prob.


    % for q=j:-1:1
    %     At=zeros(j,k*nn); %column, row
    %     At = sym(At);
    %     At(q,1) = init/w0;
    %     
    %     Bt=zeros(j,(s)*nn); %column, row
    %     Bt = sym(Bt);

    for m=j:-1:1
    %     A(m,1) = 1/j;
        for n=1:k*nn
    %             disp(m);
    %             disp(n);
            if A(m,n)~=0
                if m>1 
    %                     disp('m>1');
                    if n+1*nn<=k*nn
                        A(m-1,n+1*nn) = A(m,n)*(1-p)+A(m-1,n+1*nn);
                    end          
                    if n+s*nn<=k*nn
                        A(m-1,n+s*nn) = A(m,n)*(p)+A(m-1,n+s*nn);
                    end

                    if n+1*nn>k*nn
                        B(m-1,n+1*nn-k*nn) = A(m,n)*(1-p)+B(m-1,n+1*nn-k*nn);
    %                         disp(Bt(m-1,n+1-k));
                    end

                    if n+s*nn>k*nn
                        B(m-1,n+s*nn-k*nn) = A(m,n)*(p)+B(m-1,n+s*nn-k*nn);
    %                         disp(Bt(m-1,n+s-k));
                    end
                 end
            end
        end
    end

    %     At.'
    %     Bt.'
    %     A=A+At;
    %     B=B+Bt;
    % end
    A = simplify(A);
    B = simplify(B);



    pr=zeros(1,j); %column, row
    pr = sym(pr);

    for q=1:j

        for n=1:(s)*nn
            pr(1,q) = B(q,n) + pr(1,q);
        end
    %     pr(1,q) = simplify(pr(1,q));
    end
    pr = simplify(pr);


    % Asum=zeros(j,1); %column, row
    % Asum = sym(Asum);
    % for q=1:j
    %         
    %     for n=1:(k)*nn
    %         Asum(q,1) = A(q,n) + Asum(q,1);
    %     end           
    %     Asum(q,1) = simplify(Asum(q,1));
    % end
    % Asum = Asum.';%simple transpose, not conjugate transpose
    % disp('Asum:');
    % disp(Asum);




    % result=zeros(j,1); %column, row
    % result = sym(result);



    result = [pr(j) pr(1:(j-1))];
    for q=1:j
    %     Asum=zeros(j,1); %column, row
    %     Asum = sym(Asum);
    % 
    %     for n=1:(k)*nn
    %         Asum(q,1) = A(q,n) + Asum(q,1);
    %     end           
    %     Asum(q,1) = simplify(Asum(q,1));
        result(q) = result(q)/((j+1-q)/j-pr(q));
    end

    disp('A:');
    disp(A.');%simple transpose, not conjugate transpose


    disp('B:');
    disp(B.');

    disp('pr:');
    disp(pr);

    disp('result:');
    disp(result);

    tau = 0.2;
    p = 1 - (1-tau)^(no-1);
    disp(eval(result));

    syms tau;
    syms p;


    At=A(1,:);


end


% for q=1:j
%     disp(eval(pr(q)));
% end


% 
% if A(1)=='2'
%     disp('ffffff');
% else
%     disp('dddddd');
% end


% syms a;
% 
% b = a*ps;
% c = b/a;


% 
% disp(c);
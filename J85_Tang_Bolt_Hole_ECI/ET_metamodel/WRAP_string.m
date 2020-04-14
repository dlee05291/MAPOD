function Y=WRAP_string(string1,n_columns)

% WRAP_string wrap a long syms expression. 
%   WRAP_string(string1) rewrites a long syms equation (string1) and wraps 
%   it so that it becomes more readable. It automatically divides the 
%   expression into individual lines with three dots "..." at the end of 
%   each line. Moreover, it detects the existence of operators +,-,*,/ so 
%   that it doesn't cut a long name of a variable (such as Lamda_w).
%
% WRAP_string(string1,n_columns) allows you to specify the maximum width 
%   of the line (default value is 72 columns + the 3 dots "..." = 75 which 
%   is the default format for autowrap commenting
%
% For example: If you have the following syms expression:
% syms Lamda_w1 Lamda_u1 Lamda_w2 Lamda_u2 L1 L2 E A1 I1 A2 I2 mu1 mu2 wn
% string1='1/4*(mu1*wn^2*I1*sin(Lamda_u1*L1)*I2*(mu2*wn^2/E/I2)^(3/4)*cos(Lamda_u2*L2)*cos(Lamda_w2*L2)*sinh(Lamda_w2*L2)*cos(alpha)^2
% This code will output the same string in a more readable and compact way 

% Author: Mohannad Hakeem 5-28-2009

if (nargin==1) 
n_columns=72 % default width of each line is 72
end
Y=[];



N_lines=floor(length(string1)/n_columns);

FF=find(string1=='+' | string1=='-' | string1=='*' | string1=='/');
aa=[1 FF];
    % aa finds the indices of the locations of operators   

% N_lines=4;

diff1=zeros(1,N_lines);
b=zeros(1,N_lines);
b(1)=1;
i=2;
j=1;

while j==1  
proposed_line=n_columns*(i-1)-diff1(i-1);

if (proposed_line>length(string1))
    j=2;
last=strcat(string1(aa(b(i-1))+1:length(string1)));
Y=strvcat(Y,last);

else
    b(i)=max(find(aa<=proposed_line));

%     if (b(i)<70*(i-1))
        diff1(i)=n_columns*(i-1)-aa(b(i));
        
    if(i==2)
        k=0;
    else
        k=1;
    end
    
      X=strcat(string1(aa(b(i-1))+k:aa(b(i))),'...');
      Y=strvcat(Y,X);  
        %     end
i=i+1;
end   

end      
%a class containing static methods used in the fun_eval functions of the
%systems and in the system_standalone.m to create the fun_eval
classdef commonFunctions
    methods (Static)
        % Output:
        % x - N+1 Chebyshev nodes on [a,b] (x_0=b, x_N=a),
        % w - weights of the quadrature formula in [a,b],
        % D - differentiation matrix
        % q - row vector of the barycentric weights
        % see Trefethen
        function [w,x,D,q]=cheb(N,a,b)
        if N==0
            x=1;
            D=0;
            return
        end
        p=pi*(0:N)'/N;
        x=((b-a)*cos(p)+b+a)/2;

        c=[2;ones(N-1,1);2].*(-1).^(0:N)';
        X=x(:,ones(1,N+1)); %X=repmat(x,1,N+1);
        dX=X-X';
        D=(c*(1./c)')./(dX+(eye(N+1)));
        D=D-diag(sum(D,2)); %D=D-diag(sum(D'));

        % Quadrature weights
        w=zeros(1,N+1);
        ii=2:N;
        v=ones(N-1,1);
        if mod(N,2)==0
            w(1)=1/(N^2-1);
            w(N+1)=w(1);
            for k=1:N/2-1
                v=v-2*cos(2*k*p(ii))/(4*k^2-1);
            end
            v=v-cos(N*p(ii))/(N^2-1);
        else
            w(1)=1/N^2;
            w(N+1)=w(1);
            for k=1:(N-1)/2
                v=v-2*cos(2*k*p(ii))/(4*k^2-1);
            end
        end
        w(ii)=2*v/N;
        w=w*abs(b-a)/2;

        % Barycentric weights
        q=1./prod(dX'+eye(N+1)); %q=1./prod(dX'+eye(N+1)); % row vector of the barycentric weights
        end
        
        % Computes the value of the interpolating polynomial (Nodes,Values) in theta,
        % using barycentric interpolation with Weights.
        % (see Berrut, Trefethen, 2004)
        %%%%%%% ATTENZIONE ALLA DIMENSIONE d1 o d2
        function PP = interpoly(Theta,Nodes,Values,Weights)
            
            
            TOLL=eps;
            %getting the index and closest val to Theta
            [closest_val,ii]=commonFunctions.closest_value(Nodes,Theta);
            %if is sufficiently close approximate the value
            closest_val=abs(Theta-closest_val);
            if(closest_val<TOLL)
                PP=Values(ii);
                return;
            end
            
            n=size(Theta);
            numer=zeros(n);
            denom=zeros(n);
            exact=zeros(n);
            for j=1:length(Nodes)
                xdiff=Theta-Nodes(j);
                temp=Weights(j)./xdiff;
                numer=numer+temp*Values(j);
                denom=denom+temp;
                exact(xdiff==0)=j;
            end
            PP=numer./denom;
            jj=find(exact);
            PP(jj)=Values(exact(jj));
        end
        
        
        
        %da eliminare, solo per test
        function eqIn = substituteIntegralVars(eqIn)
        expression = "\\int_{[^}]+}\^{[^}]+}{[^}]+}{[^}]+}";

        %getting the arrays of the beginning and ending positions of each match found with the reg exp 
        [inizio,fine]=regexp(eqIn,expression);

        %getting the value of how many matches we have found with the reg exp 
        [~,matches]=size(inizio); %strings begin from 1...

        %foreach match (i.e. for each integral) (from the last)
        for l=matches:-1:1
            %setting the parameters and the regular expression to extract the
            %contents between {}
            integral=extractBetween(eqIn,inizio(l)+5,fine(l));
            expression="{[^}]+}";

            %inizio1 and fine1 mark the beginning and the end of each part of
            %the integral, their len is 4 (a,b,function to integrate, delta)
            [inizio1,fine1]=regexp(integral,expression);


            %getting the various parts
            funzione=extractBetween(integral,inizio1(3)+1,fine1(3)-1);
            diff=extractBetween(integral,inizio1(4)+1,fine1(4)-1);

            %getting the whole integral (i.e. \int...{dx})
            integral=(extractBetween(eqIn,inizio(l),fine(l)));


            %find in the function to integrate all the instances of the
            %variable 
            [inizio2,fine2]=regexp(funzione,"\W"+diff+"\W");
            substitute=integral_var_+diff;
            [~,matches2]=size(inizio2); %strings begin from 1...

            %foreach match (i.e. for each integral)
            for ll=matches2:-1:1
                funzione=strrep(funzione,extractBetween(funzione,inizio2(ll),fine2(ll)),substitute);
                %susbstituing the intergral with dot(f^,weights)*(b-a)
            end

        end
        end
        
        
        
        %CITE:
        %Benjamin Bernard (2022). Binary search for closest value in an array (https://www.mathworks.com/matlabcentral/fileexchange/37915-binary-search-for-closest-value-in-an-array), MATLAB Central File Exchange. Retrieved April 11, 2022. 
        function [v, inf] = closest_value(arr, val)
        % Returns value and index of arr that is closest to val. If several entries
        % are equally close, return the first. Works fine up to machine error (e.g.
        % [v, i] = closest_value([4.8, 5], 4.9) will return [5, 2], since in float
        % representation 4.9 is strictly closer to 5 than 4.8).
        % ===============
        % Parameter list:
        % ===============
        % arr : increasingly ordered array
        % val : scalar in R
        len = length(arr);
        inf = 1;
        sup = len;
        % Binary search for index
        while sup - inf > 1
            med = floor((sup + inf)/2);

            % Replace >= here with > to obtain the last index instead of the first.
            %modified with <= (array ordered from bigger to smallest)
            if arr(med) <= val 
                sup = med;
            else
                inf = med;
            end
        end
        % Replace < here with <= to obtain the last index instead of the first.
        if sup - inf == 1 && abs(arr(sup) - val) < abs(arr(inf) - val)
            inf = sup;
        end  
        v = arr(inf);
        end
        
        %get the different parameters from an equation containing an
        %integral, format: \int_{a}^{b}{expression}{integration variable}
        function out = getIntegral(eqIn)
            out="";
            expression = "\\int_{[^}]+}\^{[^}]+}{[^}]+}{[^}]+}";

            %getting the arrays of the beginning and ending positions of each match found with the reg exp 
            [inizio,fine]=regexp(eqIn,expression);

            %getting the value of how many matches we have found with the reg exp 
            [~,matches]=size(inizio); %strings begin from 1...

            %foreach match
            for l=1:matches
                integral=extractBetween(eqIn,inizio(l)+5,fine(l));
                expression="{[^}]+}";
                [inizio1,fine1]=regexp(integral,expression);
                [~,matches1]=size(inizio1); %strings begin from 1...
                for kk=1:matches1
                    disp(extractBetween(integral,inizio1(kk)+1,fine1(kk)-1));
                end
                disp("Nuovo");
                disp(extractBetween(eqIn,inizio(l),fine(l)));
            end
        end
    end
end
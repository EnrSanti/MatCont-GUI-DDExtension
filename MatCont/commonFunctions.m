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
        function delayFunctionsns = getDelayFunctions(eqsIn,coords,tempi,dim)

        delayFunctionsns="";
        %getting the coordinates (x,y...)
        [coords,~]=split(coords,",");
        %getting the times
        [tempi,~]=split(tempi,",");

        %getting how many time variables we have
        [no_times,~]=size(tempi);

        %for each equation in the system
        for eqIndex=1:dim
            %getting the rhs of the current equation considered
            [eq,~]=split(eqsIn(eqIndex,:),"=");
            eq=eq(2);

            %getting the string itself
            eq=cell2mat(eq);

            %for each coordinate substitute
            for i=1:dim
                %for each time var
                for j=1:no_times

                    %if we have multiple time variables, extract one at the time
                    times(j)=string(tempi(j));

                    %the regular expression of cor_xyz[t(i)...]
                    expression = coords(i)+"\["+times(j)+"\W(\["+times(j)+"\W[^\]]*\]|[^\]])*\]";


                    %getting the arrays of the beginning and ending positions of each match found with the reg exp 
                    [inizio,fine]=regexp(eq,expression);

                    %getting the value of how many matches we have found with the reg exp 
                    [~,matches]=size(inizio); %strings begin from 1...

                    %foreach match substitute the expression:
                    for l=1:matches
                        %we extract the delay from the string we have found
                        %(e.g. [t-g(x)] -> g(x))
                        replace=extractBetween(eq,inizio(l)+1+strlength(coords(i))+strlength(times(j)),fine(l)-1);
                        replace=replace{1};
                        delayFunctionsns(end+1)=replace;

                    end
                end

            end
        end
        delayFunctionsns=delayFunctionsns(2:end);
        end
 
    end
end
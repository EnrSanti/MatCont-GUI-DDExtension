classdef commonFunctions
    methods (Static)
        function [w,x,D,q]=cheb(N,a,b)
        % Output:
        % x - N+1 Chebyshev nodes on [a,b] (x_0=b, x_N=a),
        % w - weights of the quadrature formula in [a,b],
        % D - differentiation matrix
        % q - row vector of the barycentric weights
        % see Trefethen

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
    end
end
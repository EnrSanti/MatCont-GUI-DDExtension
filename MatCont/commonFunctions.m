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
        q=[1/2,ones(1,N-1),1/2].*(-1).^(0:N); %q=1./prod(dX'+eye(N+1)); % row vector of the barycentric weights
        end
        
        % Computes the value of the interpolating polynomial (Nodes,Values) in theta,
        % using barycentric interpolation with Weights.
        % (see Berrut, Trefethen, 2004)
        %%%%%%% ATTENZIONE ALLA DIMENSIONE d1 o d2
        function PP = interpoly(Theta,Nodes,Values,Weights)
            
            %{
            TOLL=eps;
            %getting the index and closest val to Theta
            [closest_val,ii]=commonFunctions.closest_value(Nodes,Theta);
            %if is sufficiently close approximate the value
            closest_val=abs(Theta-closest_val);
            if(closest_val<TOLL)
                PP=Values(ii);
                return;
            end
            %} 
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
        

        %function that given a row vector, returns the string that identifies it
        %with a ";" at the end
        %e.g. [1,2,3] -> "[1,2,3];"
        %vett: a row vector
        function strRow = RowVett2Str(vett)
            strRow="["+RowVett2StrAux(vett)+"];";   
        end
        %function that given a column vector returns a string identifying the
        %vector with a ";" at the end (e.g. [1;2;3] -> "[1;2;3]";)
        %vett: a column vector
        function strCol = ColVett2Str(vett)
            %getting the size of the vector
            [rows, ~]=size(vett);
            %starting the string
            strCol="["+vett(1,:);
            %for each element, concatenate it with ";"
            for i=2:(rows)
                strCol=strCol+";"+vett(i,:);
            end
            %add at the end ];"
            strCol=strCol+"];";
        end
        %function that given a row vector, returns the string of all its elements
        %divided by ","
        %e.g. [1,2,3] -> "1,2,3"
        %vett: a row vector
        function strRow = RowVett2StrAux(vett)
            %getting the 1st element
            strRow=vett(1);

            vett=vett(2:end);
            %for every element (not the first, concatenate it with a "," to the
            %string"
            for el=vett
                strRow=strRow+","+el;
            end
        end
        %function that given a diff equation, the coordinates, times and the
        %dimension of the systems (i.e. how many total coordinates we have) returns
        %the string identifying the rhs of the original equation ready to be
        %insered in the systems file (according to Francesca's notation)
        %eqIn: a string containing the differential equation according to the
        %matcont format (typed in by the user)
        %coords: a string containing the coordinates of the system, divided by ","
        %(e.g "x,y,z")
        %tempi: a string containing the time variables, divided by "," (e.g
        %"t1,t2")
        %dim: integer, the number of coordinates in the system
        function eq = parseDDE(eqIn,coords,tempi,dim,REcoords,DDEcoords)

            REno=length(REcoords);
            DDEno=length(DDEcoords);
            %getting the coordinates (x,y...)
            [coords,~]=split(coords,",");
            %getting the times
            [tempi,~]=split(tempi,",");

            %getting how many time variables we have
            [no_times,~]=size(tempi);

            %getting the rhs of the current equation considered
            [eq,~]=split(eqIn,"=");
            eq=eq(2);

            %getting the string itself
            eq=cell2mat(eq);

            %for each coordinate substitute
            for i=1:dim
                %aux=dim;
                %aux=aux-i; ultima modifica
                %for each time var
                for j=1:no_times

                    %if we have multiple time variables, extract one at the time
                    times(j)=string(tempi(j));


                    %the regular expression of
                    %cor_xyz[t(i)..cor_xyz[...]...] can recognise a nested
                    %delay
                    expression = coords(i)+"\["+times(j)+"\W(\["+times(j)+"\W[^\]]*\]|[^\]])*\]";


                    %getting the arrays of the beginning and ending positions of each match found with the reg exp 
                    [inizio,fine]=regexp(eq,expression);

                    %getting the value of how many matches we have found with the reg exp 
                    [~,matches]=size(inizio); %strings begin from 1...

                    %foreach match substitute the expression:
                    for l=matches:-1:1
                        %we extract the delay from the string we have found
                        %(e.g. [t-g(x)] -> g(x))
                        replace=extractBetween(eq,inizio(l)+1+strlength(coords(i))+strlength(times(j)),fine(l)-1);

                        %it's not needed to eval(replace) to actually evaluate the
                        %expression, anyway, creating the string to substitute to
                        %x[t-...] which corresponds to the call to the
                        %interpolation function

                        %da portare fuori dal ciclo
                        approx="";
                        if(~(i>DDEno))%%dde coord
                            approx="commonFunctions.interpoly("+replace+",ScaledNodes,[yM("+(i)+");VM("+(i)+":d2:end)],BaryWeights)";
                        else%re coord
                            %modifica qui
                            approx="commonFunctions.interpoly("+replace+",ScaledNodes,"+"[0;derState("+(i-DDEno)+":d1:end)],BaryWeights)";
                        end
                        %in the rhs we substitute the coordinate with a delay with
                        %the function that will compute its value (e.g y[t-2*TAU]
                        %-> interpoly(-2*TAU,...))
                        delayfound=extractBetween(eq,inizio(l),fine(l));
                        eq = convertStringsToChars(strcat(eq(1:inizio(l) - 1),strcat( approx, eq(inizio(l)+(fine(l)+1-inizio(l)):end))));
                        %eq=strrep(eq,delayfound,approx); 
                    end
                end

            end
            %now each coordinate with a delay has been substitued by the
            %approximating value (that has to be computed)

            %for each coordinate (without delay) substitute yM(i)
            for i=1:dim
                eq=strrep(eq,coords(i),strcat("yM(",strcat(sprintf("%d",(i)),")")));
            end
            %var con yM(i) sostituite
        end


        % function that given the equations (the whole system), the coordinates,
        % the time variables and the system dimension, returns  array containing
        % all the delay functions used in the system (with duplicates) TODO:
        % optimize
        %p eqsIn: an array of strings describing the whole system (typed by the
        % user) in the matcont format
        %tempi: a string containing the time variables, divided by "," (e.g
        %"t1,t2")
        %dim: integer, the number of coordinates in the system
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
                    for i=1:length(coords)
                        %for each time var
                        for j=1:no_times

                            %if we have multiple time variables, extract one at the time
                            times(j)=string(tempi(j));
                            for kk=1:length(coords)
                                %the regular expression of cor_xyz[t(i)...] not
                                %between {...}
                                exp2=coords(kk)+"\["+times(j)+"\W[^\]]*\]";
                                expression = coords(i)+"\["+times(j)+"\W("+exp2+"|[^\]\[]*)\]";


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
                end
                delayFunctionsns=checkIntegralVars(unique(delayFunctionsns(2:end)),eqsIn,dim,tempi);
        end



        % function that given the delay functions of the system, the eqautions in the system 
        % and the number of the equation in the system removes from the delays
        % (vector of strings) the strings that contain a variable present in an
        % integral (i.e. from all the delays in the system we keep only the ones non
        % depending on differentiation variables) and adds to the vector the delays
        % found in the first extreme of integration
        % delays
        % eqsIn: an array of strings describing the whole system (typed by the
        % user) in the matcont format
        % delays: a string vector containing the delays in the system (e.g.
        % x'=...x[t-g(t)]...y[t-h(t)] ->[g(t),h(t)]
        % dim: an integer denotin the number of equations in the system
        % tempi: a string containing the time variables, divided by "," (e.g
        %"t1,t2")
        function delays= checkIntegralVars(delays,eqsIn,dim,tempi)

            %getting the differetiation variables in the system
            [integralVars,delaysToAdd]=getIntegralVars(eqsIn,dim,tempi);

            %the number of deleted delays
            deleted=0;
            %for each delay
            for i=1:(length(delays))
                %if the delay contains one of the differentiation variables
                if(contains(delays(i-deleted),integralVars))
                    %delete element and increment the number of deleted elements
                    delays(i-deleted)=[]; %elimina elemento
                    deleted=deleted+1;
                end
            end 
            %concat to create a unique vector of delays
            delays=horzcat(delays,delaysToAdd);
        end



        % function that given the equation in the system, and the dimension of the
        % system returns the array of the differentiation variables used and the delays arising from the first extreme of integraion .
        % eqsIn: an array of strings describing the whole system (typed by the
        % user) in the matcont format
        % dim: an integer denotin the number of equations in the system
        %tempi: a string containing the time variables, divided by "," (e.g
        %"t1,t2")
        function [integralVars,delaysToAdd] = getIntegralVars(eqsIn,dim,tempi)

            delaysToAdd=[""];
            %initializing empty string array
            integralVars=[""];
            %expression recognising the structure of an integral
            integralRegExp = "\\int_{[^}]+}\^{[^}]+}{[^}]+}{[^}]+}"; 
            integralPartsRegExp="{[^}]+}";

            %for each equation
            for eqIndex=1:dim

                %getting the rhs of the current equation considered
                [eqIn,~]=split(eqsIn(eqIndex,:),"=");
                eqIn=eqIn(2);

                %getting the string itself
                eqIn=cell2mat(eqIn);
                %getting the arrays of the beginning and ending positions of each match found with the reg exp 
                [inizio,fine]=regexp(eqIn,integralRegExp); 

                %getting the value of how many matches we have found with the reg exp 
                [~,matches]=size(inizio); %strings begin from 1...

                %foreach match (i.e. for each integral)
                for l=1:matches
                    %setting the parameters and the regular expression to extract the
                    %contents between {}
                    integral=extractBetween(eqIn,inizio(l)+5,fine(l));

                    %inizio1 and fine1 mark the beginning and the end of each part of
                    %the integral, their len is 4 (a,b,function to integrate, delta)
                    [inizio1,fine1]=regexp(integral,integralPartsRegExp);

                    inizio1=cell2mat(inizio1);
                    fine1=cell2mat(fine1);

                    %getting the various parts (the first and second extreme of intregration
                    %and the integration variable)
                    a=extractBetween(integral,inizio1(1)+1,fine1(1)-1);
                    b=extractBetween(integral,inizio1(2)+1,fine1(2)-1);
                    diff=extractBetween(integral,inizio1(4)+1,fine1(4)-1);

                    %adding the integration variable to the list to return
                    integralVars(end+1)=string(diff);

                    %extracting the string from the cells
                    a=string(a);
                    b=string(b);



                    %{
                    for timeNo=1:length(tempi)
                        if(contains(a,tempi(timeNo)))
                            a=extractAfter(a,length(tempi(timeNo)));
                        end
                    end
                    %}

                    %adding the extremese to the delay list to return (since
                    %depending on a,b >0 or a,b<0 \int_{a}^{b}{... x[t +/- s]}{s}
                    %in the first case (a,b>0 and x[t-s]) the maximum value of s
                    %will be in b, and in the second case (a,b<0 and x[t+s]) it
                    %will be in a
                    delaysToAdd(end+1)=a;
                    delaysToAdd(end+1)=b;
                end
            end
            integralVars=integralVars(2:end); %removing the starting ""
            delaysToAdd=delaysToAdd(2:end); %removing the starting ""
        end



        %function that given a diff equation, gets the different parameters from an equation containing an
        %integral, format: \int_{a}^{b}{expression}{integration variable} and
        %returns the actual equation to insert in the .m file
        function eqIn = parseIntegral(eqIn) 
            global gds;
            expression = "\\int_{[^}]+}\^{[^}]+}{[^}]+}{[^}]+}";

            %getting the arrays of the beginning and ending positions of each match found with the reg exp 
            [inizio,fine]=regexp(eqIn,expression);

            %getting the value of how many matches we have found with the reg exp 
            [~,matches]=size(inizio); %strings begin from 1...

            %the expression to get the different integral components    
            expression="{[^}]+}";

            %foreach match (i.e. for each integral) (loop "backwards")
            for l=matches:-1:1
                %setting the parameters to extract the
                %contents between {}
                integral=extractBetween(eqIn,inizio(l)+5,fine(l));


                %inizio1 and fine1 mark the beginning and the end of each part of
                %the integral, their len is 4 (a,b,function to integrate, dx)
                [inizio1,fine1]=regexp(integral,expression);


                %getting the various parts
                a=extractBetween(integral,inizio1(1)+1,fine1(1)-1);      
                b=extractBetween(integral,inizio1(  2)+1,fine1(2)-1);
                funzione=extractBetween(integral,inizio1(3)+1,fine1(3)-1);
                funzione=funzione{1};
                diff=extractBetween(integral,inizio1(4)+1,fine1(4)-1);

                %getting the whole integral (i.e. \int...{dx})
                integral=(extractBetween(eqIn,inizio(l),fine(l)));

                %creating the string *(a-(b)) to not recompute it every time
                b_a="*("+a+"-("+b+"))"; %le doppie pararentesi, a può avere un segno


                %checking if funzione begins with the intergration varaible
                %if(funzione inizia con var, sostituisci var con new_
                if(length(cell2mat(regexp(extractBetween(funzione,1,length(diff{1})),diff)))>0)
                    eqIn=strrep(eqIn,integral,"\int_{"+a+"}^{"+b+"}"+"{var_int_NEW"+diff+extractBetween(funzione,length(diff{1})+1,length(funzione))+"}{"+diff+"}");

                    integral=extractBetween(eqIn,inizio(l)+5,fine(l)+11);
                    expression="{[^}]+}";

                    %inizio1 and fine1 mark the beginning and the end of each part of
                    %the integral, their len is 4 (a,b,function to integrate, delta)
                    [inizio1,fine1]=regexp(integral,expression);


                    %getting the various parts again
                    a=extractBetween(integral,inizio1(1)+1,fine1(1)-1);      
                    b=extractBetween(integral,inizio1(2)+1,fine1(2)-1);
                    funzione=extractBetween(integral,inizio1(3)+1,fine1(3)-1);
                    funzione=funzione{1};
                    diff=extractBetween(integral,inizio1(4)+1,fine1(4)-1);

                    %getting the whole integral (i.e. \int...{dx})
                    integral=(extractBetween(eqIn,inizio(l),fine(l)+11));
                    fine(l)=fine(l)+11;

                end


                %checking if funzione ends with diff
                lenDiff=length(diff{1});
                if(length(cell2mat(regexp(extractBetween(funzione,length(funzione)-lenDiff,length(funzione)),diff)))>0)
                    eqIn=strrep(eqIn,integral,"\int_{"+a+"}^{"+b+"}"+"{"+extractBetween(funzione,1,length(funzione)-lenDiff)+"var_int_NEW"+diff+"}{"+diff+"}");

                    integral=extractBetween(eqIn,inizio(l)+5,fine(l)+11);
                    expression="{[^}]+}";

                    %inizio1 and fine1 mark the beginning and the end of each part of
                    %the integral, their len is 4 (a,b,function to integrate, delta)
                    [inizio1,fine1]=regexp(integral,expression);


                    %getting the various parts
                    a=extractBetween(integral,inizio1(1)+1,fine1(1)-1);      
                    b=extractBetween(integral,inizio1(2)+1,fine1(2)-1);
                    funzione=extractBetween(integral,inizio1(3)+1,fine1(3)-1);
                    funzione=funzione{1};
                    diff=extractBetween(integral,inizio1(4)+1,fine1(4)-1);

                    %getting the whole integral (i.e. \int...{dx})
                    integral=(extractBetween(eqIn,inizio(l),fine(l)+11));
                end


                [inizio2,fine2]=regexp(funzione,"\W"+diff+"(?=\W)"); %lookahead

                [~,matches2]=size(inizio2); %strings begin from 1...


                %foreach match (i.e. for each var in the integral)
                before="";
                after="";
                diff="var_int_NEW"+diff;
                while(length(inizio2)>0)
                    ll=length(inizio2);
                    if(length(regexp(funzione(inizio2(ll)),"\W"))>0)
                        before=funzione(inizio2(ll));

                    end
                    if(length(regexp(funzione(fine2(ll)+1),"\W"))>0)
                        after=funzione(fine2(ll)+1);

                    end
                    substituteString=extractBetween(funzione,inizio2(ll),fine2(ll)+1);
                    findsubstitute=substituteString;
                    if(before=="^")
                        findsubstitute="\"+findsubstitute; %escaping ^
                        findsubstitute=cellstr(findsubstitute);
                    end
                    if(after=="^")
                        findsubstitute=extractBefore(substituteString,strlength(substituteString))+"\^"; %escaping ^
                        findsubstitute=cellstr(findsubstitute);
                    end

                    [deleteInizio,deleteFine]=regexp(funzione,findsubstitute);

                    %deleteing the initial positions of the same string (there will
                    %all be susbstitued)
                    [X,Y] = ismember(deleteInizio{1},inizio2);
                    inizio2(Y(X)) = [];

                    %deleteing the final positions of the same string (there will
                    %all be susbstitued)
                    [X,Y] = ismember(deleteFine{1},fine2);
                    fine2(Y(X)) = [];


                    funzione=strrep(funzione,substituteString,before+diff+after);
                    funzione=funzione{1};

                end

                %riscala
                fCap=strrep(funzione,diff,"thetaCap"+b_a+"+"+b);

                %susbstituing the intergral with dot(f^,weights)*(b-a)
                %must substitute only the last
                %eqIn=strrep(eqIn,integral,"dot("+fCap+",wCap)"+b_a);

                posizioni = strfind(eqIn, integral);
                last = posizioni(end);

                actual_b_a="*("+b+"-("+a+"))"; %le doppie pararentesi, a può avere un segno
                eqIn = strcat(strcat(extractBetween(eqIn,1,last - 1), "dot("+fCap+",wCap)"+actual_b_a), extractBetween(eqIn,last+strlength(integral),strlength(eqIn)));
            end
        end

        %function that given an equation returns true iff the LHS doesn't contain
        %the derivative of a coordinate (i.e. if the equation is RE)
        %eq: equation of the system with the substituted coordinates
        function answ=isRE(eq)
            answ=true;
            [eqComponents,~]=strsplit(eq,"=");
            if(contains(eqComponents{1},"'"))
                answ=false;
            end
        end
        %function that given the equations in the system returns the number of the
        %renewal equations in the system
        %equations: the system of equations with the substituted coordinates
        function no=getREno(equations)
            global gds;
            no=0
            %for each equation in the system check if is an RE
            for i=1:gds.dim
               %the equation considered 
               eq=equations(i,:) ;
               if(isRE(eq))
                   no=no+1; %no++
               end
            end
        end
        %function that given the file identificator of the file (opened) to write,
        %the format, the current content of the file and what to write, wites
        %content on the file and updates fileContent
        %fid_write: the file identificator of the openend file
        %format: a string representing the format of how to write the string on the
        %file
        %fileContent: a string containing the current content of the file, updated
        %after the file has been written
        %content: a string containing the new content to write on the file
        function fileContent = write_M_and_File_Content(fid_write,format,fileContent,content) 
            fprintf(fid_write,format,content);   
            fileContent = [fileContent,  sprintf(format,content)];
        end
         %a function that given an equation (LHS=RHS) substitutes all the products,
        %exponents and divisions with the respetive component wise operation
        function eqIn = parseREDot(eqIn) 
            eqIn=strrep(eqIn,"*",".*");
            eqIn=regexprep(eqIn,"\^(?!{)",".^"); %non strrep, the integral has {}^{}...
            eqIn=strrep(eqIn,"/","./");
        end
        function string=replace_sys_inputDDE(string)
            global gds;
            if isempty(string)
                string = '';
                return
            end
            string_sys = cellstr(string);
            string=''; temp = '';eq = '';
            for j = 1:(gds.dim)
                sys{j} = strcat(gds.coordinates{j,1},char(39));
            end
            [temp,eq] = parse_input(string_sys,sys);
            num_temp = size(temp,1);
            for i = 1:num_temp
                string{i,1} = strcat(temp{i,1},';');
            end
        end
        %-_-_-_-_-_-_%



    end
end

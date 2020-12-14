function [str, fun, span, FittedVarsCell] = C_schot_fit_A(x,y,boxin,A,es,N,Vb,n,varargin)
f = 10^ceil(abs(log10(1/y(end-1))));
if ~isempty(varargin)
    if isa(varargin{1}, 'numeric')
        fit_cond_val = cell2mat(varargin(1:2:end)');
        fit_cond_names = varargin(2:2:end);
    else
        fit_cond_val = cell2mat(varargin(2:2:end)');
        fit_cond_names = varargin(1:2:end);
    end
else
    fit_cond_val = [];
    fit_cond_names = {};
end
if length(boxin)==2
    box = [boxin -inf inf];
elseif isempty(boxin)
    box = [-inf inf -inf inf];
end

syms c V
independent = {'V'};
T = 295;
k = 8.617e-5;
q = 1.6e-19;
e0 = 8.854e-14;
sqef = sqrt(q*e0)*f;
cschot = c/sqrt(2*n*(n*(Vb-0.0257)-V));
coefficients = {};
sym_factors = {};
sym_titles = {};
sym_title_find = {'N'; 'Doping [cm^{-3}]'};
coefficients_lim = [];
prod_values = [];


if isa(Vb,'sym')
    coefficients{end+1} = 'Vb';
    cmp = strcmp('Vb',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0 10 0.5];
    end
end
if isa(n,'sym')
    coefficients{end+1} = 'n';
    cmp = strcmp('n',fit_cond_names);
    if any(cmp)
        coefficients_lim(end+1,:) = fit_cond_val(cmp,:);
        fit_cond_val(cmp,:) = [];
        fit_cond_names(cmp) = [];
    else
        coefficients_lim(end+1,:) = [0.01 1000 1];
    end
end



sqrt_param = {'es','N'};
char_param = {'A','es','N'};

% only linear param 
param = {A,es,N};
flag = 1;
for i=1:length(param)
    p = param{i};
    if isa(p,'sym')
        if flag==1
            coefficients{end+1} = 'c';
            coefficients_lim(end+1,:) = [sqef sqef sqef];
            flag = 0;
        end
        sym_factors{end+1} = char(p);
        sym_titles{end+1} = sym_title_find(2,strcmpi(sym_title_find(1,:),char(p)));
        if ~isempty(fit_cond_val(strcmp(char(p),fit_cond_names),:))
            if any(strcmp(char(p),sqrt_param))
                coefficients_lim(end,:) = coefficients_lim(end,:).*sqrt(abs(fit_cond_val(strcmp(char(p),fit_cond_names),:)));
            else
                coefficients_lim(end,:) = coefficients_lim(end,:).*fit_cond_val(strcmp(char(p),fit_cond_names),:);
            end
        else
        	coefficients_lim(end,:) = [-inf inf 0];
        end
    else
        if any(strcmp(char_param{i},sqrt_param))
            if flag==1
                coefficients{end+1} = 'c';
                coefficients_lim(end+1,:) = [sqef sqef sqef]*sqrt(p);
                flag = 0;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:)*sqrt(p);
            end
            prod_values(end+1) = sqrt(p);
        else
            if flag==1
                coefficients{end+1} = 'c';
                coefficients_lim(end+1,:) = [sqef sqef sqef]*p;
                flag = 0;
            else
                coefficients_lim(end,:) = coefficients_lim(end,:)*p;
            end
            prod_values(end+1) = p;
        end
    end
end
% only linear param ^

fo = fitoptions('Method','NonlinearLeastSquares','Lower',coefficients_lim(:,1),'Upper',coefficients_lim(:,2),'StartPoint',coefficients_lim(:,3),'Robust','Bisquare');
str=['Fit Range= '  num2str(boxin)  newline];
if coefficients{end} ~= 'c'
    c = A*sqrt(N)*sqef;
end
ft = fittype(char(cschot), 'independent',independent, 'coefficients',coefficients, 'options',fo);
dvec = ~excludedata(x,y,'box',box);
try
    [fitt, gof] = fit(x(dvec),y(dvec)*f,ft);
catch ME
    if contains(ME.identifier, 'complex','IgnoreCase',true)
        warning('Complex value computed by model function, fitting cannot continue. Assisgning NaN values to outputs.');
        fitt = [];
        gof = [];
    end
end
        
if ~isempty(fitt)
    FittedVarsCell = {};
    if coefficients{end} == 'c' && ~isempty(sym_factors)
        if any(strcmp(sym_factors{1},sqrt_param))
            str = [str sym_factors{:} '=' num2str((fitt.c/sqef/prod(prod_values))^2,2) ' '];
            eval([sym_factors{1} '=' num2str((fitt.c/sqef/prod(prod_values))^2) ';']);
            FittedVarsCell(:, end+1) = {sym_factors{1}; (fitt.c/sqef/prod(prod_values))^2; sym_titles{1}};
        else
            str = [str sym_factors{:} '=' num2str(fitt.c/sqef/prod(prod_values),2) ' '];
            eval([sym_factors{1} '=' num2str(fitt.c/sqef/prod(prod_values)) ';']);
            FittedVarsCell(:, end+1) = {sym_factors{1}; fitt.c/sqef/prod(prod_values); sym_titles{1}};
        end
    end
    if any(strcmp(coefficients, 'Vb'))
        str = [str 'Vb=' num2str(fitt.Vb,2)];
        eval(['Vb=' num2str(fitt.Vb) ';']);
        FittedVarsCell(:, end+1) = {'Vb'; Vb; 'V Built in [V]'};
    end
    if any(strcmp(coefficients, 'n'))
        str = [str ' n=' num2str(fitt.n,2)];
        eval(['n=' num2str(fitt.n) ';']);
        FittedVarsCell(:, end+1) = {'n'; n; 'Ideality Factor'};
    end
    str = [str ' R^{2}=' num2str(gof.rsquare,3) ];
    FittedVarsCell(:, end+1) = {'R2'; gof.rsquare; 'R^{2}'};
    fun = @(V) A.*sqrt(q.*es.*e0.*N./(2.*n.*(n.*(Vb-k*T)-V)));
    span = [min(x), min(max(x), Vb-2*k*T)];
else
    if coefficients{end} == 'c' && ~isempty(sym_factors)
        coefficients = [sym_factors(1) coefficients];
    end
    str = ['Complex Values Computed.' newline 'Fitting aborted'];
    fun = [];
    span = [];
    FittedVarsCell = [coefficients {'R2'}; num2cell(nan(1,length(coefficients)+1)); coefficients {'R2'}]; 
end



function [pftcon] = read_xml_params(xmlfile)


xmlroot = xml2struct(xmlfile);
pnames  = fieldnames(xmlroot.all); % Names of parameters, last is pft specialized
n_pfts    = numel(xmlroot.all.pfts.pft);

% Initialize

% Get the name of the pft's from the attribute field
for ip=1:n_pfts
    pftcon.tag{ip} = xmlroot.all.pfts.pft{ip}.Attributes.tag;
end

% These are the numerical floating point parameters the model needs
expt_par_dp = {'c2b','eclim','llspan','bl_min','h_max','h_min','slatop', ...
    'd_adult','d_sap','l2f_ratio','agb_fraction','latosa_int', ...
    'latosa_slp','d2h1_ad','d2h2_ad','d2h3_ad','d2h1_sap', ...
    'd2h2_sap','d2bl1_ad','d2bl2_ad','d2bl3_ad','d2bl1_sap', ...
    'd2bl2_sap','d2bag1','d2bag2','wood_density'};

% These are the string parameters the model needs
expt_par_str = {'hallom_mode','lallom_mode','fallom_mode','aallom_mode', ...
    'callom_mode','sallom_mode'};

% Load up the parameter structure with the floating point values
for iv = 1:numel(expt_par_dp)
    % Check to see if this is a cross-pft parameter designation
    id = find(strcmp(pnames,expt_par_dp{iv}));
    if(isempty(id))
        % Turns out this is probably a pft specific parameter
        for ip=1:n_pfts
            id = find(strcmp(expt_par_dp{iv},fieldnames(xmlroot.all.pfts.pft{ip})),1);
            if(isempty(id))
                display('Expected parameter cannot be found in XML');
                display(sprintf('%s, pft %d',expt_par_dp{iv},ip));
                return;
            else
                pftcon.(expt_par_dp{iv})(ip) = str2double(xmlroot.all.pfts.pft{ip}.(expt_par_dp{iv}).Text);
            end
        end
    else
        for ip=1:n_pfts
            pftcon.(expt_par_dp{iv})(ip) = str2double(xmlroot.all.(pnames{id}).Text);
        end
    end
end

% Load up the parameter structure with the string values
for iv = 1:numel(expt_par_str);
    % Check to see if this is a cross-pft parameter designation
    id = find(strcmp(pnames,expt_par_str{iv}));
    if(isempty(id))
        % Turns out this is probably a pft specific parameter
        for ip=1:n_pfts
            id = strcmp(expt_par_str{iv},fieldnames(xmlroot.all.pfts.pft{ip}));
            if(isempty(id))
                display('Expected parameter cannot be found in XML');
                display(sprintf('%s, pft %d',expt_par_str{iv},ip));
                return;
            else
                pftcon.(expt_par_str{iv}){ip} = strtrim(xmlroot.all.pfts.pft{ip}.(expt_par_str{iv}).Text);
            end
        end
    else
        for ip=1:n_pfts
            pftcon.(expt_par_str{iv}){ip} = strtrim(xmlroot.all.(pnames{id}).Text);
        end
    end
end

end
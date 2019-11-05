%% User name variables
hysys_name = "\spreadsheet linked new.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";

%% Connect to HYSYS
hysys_path = pwd + hysys_name;
hysys = HYSYSFile(hysys_path, spreadsheet_input, spreadsheet_output);

% To get a variable:
inputs = hysys.feasible_point;
inputs = cell2mat(struct2cell(inputs))';
[outputs, outputs_struct] = hysys.get_output(inputs);
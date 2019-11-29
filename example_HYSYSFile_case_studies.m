%% User name variables
hysys_name = "\case study.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";

%% Connect to HYSYS
tic
hysys_path = pwd + hysys_name;
hysys = HYSYSFile_case_studies(hysys_path, spreadsheet_input, spreadsheet_output);
toc

%% Initialising Points for Training
input1_steps = 21;
input2_steps = 21;
lb = [40, 80000];                          % Lower bounds
ub = [300, 200000];                        % Upper bounds
sample_points = zeros(input1_steps*input2_steps,2);
input1 = (lb(1):((ub(1)-lb(1))/(input1_steps-1)):ub(1));
input2 = (lb(2):((ub(2)-lb(2))/(input2_steps-1)):ub(2));
data_grid = zeros(input1_steps*input2_steps,27);

for i = 1:input1_steps
    sample_points(1+(i-1)*input2_steps:i*input2_steps,1) = input1(i);
    sample_points(1+(i-1)*input2_steps:i*input2_steps,2) = input2';
end

toc

for point = 1:size(sample_points, 1)
    [outputs] = hysys.get_output(sample_points(point,:));
    data_grid(point,:) = outputs;
end

data_grid = [sample_points, data_grid]; 
toc

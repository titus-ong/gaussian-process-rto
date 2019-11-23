%% User name variables
hysys_name = "\case study - 3D.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";

%% Connect to HYSYS
tic
hysys_path = pwd + hysys_name;
hysys = HYSYSFile_case_studies_3D(hysys_path, spreadsheet_input, spreadsheet_output);
toc

%% Initialising Points for Training
input1_steps = 8;
input2_steps = 9;
input3_steps = 9;
lb = [60000, 200, 40];                          % Lower bounds
ub = [130000, 600, 120];                        % Upper bounds
sample_points = zeros(input1_steps*input2_steps*input3_steps,3);
data_grid = zeros(input1_steps*input2_steps*input3_steps,27);

input1 = (lb(1):((ub(1)-lb(1))/(input1_steps-1)):ub(1));
input2 = (lb(2):((ub(2)-lb(2))/(input2_steps-1)):ub(2));
input3 = (lb(3):((ub(3)-lb(3))/(input3_steps-1)):ub(3));
[x, y, z] = meshgrid(input1, input2, input3);


for i = 1:input1_steps
    for j = 1:input2_steps
        for k = 1:input3_steps
            sample_points((i-1)*input2_steps*input3_steps+(j-1)*input3_steps+(k-1)+1,:) = [input1(i), input2(j), input3(k)];
        end
    end
end

toc

for point = 1:size(sample_points, 1)
    [outputs] = hysys.get_output(sample_points(point,:));
    data_grid(point,:) = outputs;
end

data_grid = [sample_points, data_grid]; 
toc

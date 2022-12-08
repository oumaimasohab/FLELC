function problems = get_list_problems(files)

n = length(files);
problems = [];

for i=1:n
    files(i)
    fid = fopen(files(i), 'rt');
    datacell = textscan(fid, '%[^,]%*[^\n]', 'HeaderLines', 1);
    fclose(fid);
    list_problems = datacell{1};
    problems = [problems; list_problems];
end
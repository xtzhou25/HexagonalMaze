% splice kilosort recordings
% add two-minute spacing between two recordings
clear

chunkSize = 128 * 30000 * 60 * 5;
% 5-minute chunk


spacing = zeros(30000*60*2*128, 1);
spacing = int16(spacing);

Nrecordings = 2;
% how many recordings you want to splice

fidout = fopen('OFC2-2_20230619S1S2.dat', 'w');
rootpath = pwd;

fidin = zeros(Nrecordings,1);
for i = 1:Nrecordings
    [file, path] = uigetfile('*.dat', 'Select a dat file', [rootpath, '\example.dat']);
    fidin(i) = fopen([path, file], "r");
end



for i = 1:Nrecordings

    while ~feof(fidin(i))
        % Read a chunk of data from input_file1.dat
        data = fread(fidin(i), chunkSize, 'int16');

        % Write the chunk of data to the output file
        fwrite(fidout, data, 'int16');
    end
    fclose(fidin(i));
    if i < Nrecordings %避免在最后加上spacing
        fwrite(fidout, spacing, 'int16');
    end

end
fclose(fidout);

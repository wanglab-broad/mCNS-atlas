function Npz2Tif( input_path, file_name, output_name )
%

    labels_file = unzip(fullfile(input_path, file_name), input_path);

    labels = readNPY(labels_file{1});

    imwrite(uint16(labels), fullfile(input_path, output_name))



end


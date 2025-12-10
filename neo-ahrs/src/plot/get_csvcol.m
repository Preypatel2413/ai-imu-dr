function [col_index] = get_csvcol(filename, header)
  fid = fopen(filename, 'r');
  headers_line = fgetl(fid);

  if ~ischar(headers_line)
    error("Unable to read the header line.");
  end

    headers = strsplit(headers_line, ',');
    col_index = find(strcmp(headers, header));

  if isempty(col_index)
    error("Invalid header name.");
  end

  fclose(fid);
end

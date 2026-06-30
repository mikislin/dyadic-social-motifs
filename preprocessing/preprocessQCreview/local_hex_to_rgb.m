function rgb = local_hex_to_rgb(hex)
hex = char(strtrim(string(hex)));
rgb = sscanf(hex(2:end), '%2x%2x%2x', [1 3]) ./ 255;
end
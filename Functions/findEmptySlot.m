function outmap = findEmptySlot(map,insertsize)

x_dim = size(map,1);
y_dim = size(map,2);

spacing_kernel = ones(insertsize,insertsize);
spacing_map = conv2(map,spacing_kernel,'same');

outmap = spacing_map == 0 ;


end
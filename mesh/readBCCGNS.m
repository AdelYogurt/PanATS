function [point_list,marker_index_list]=readMarkerCGNS(mesh_filestr,marker_name_list)
% read BC data and point data from CGNS file
%
idx_base=1;
idx_zone=1;

% open file
[idx_file, ierr] = cg_open(mesh_filestr, CG_MODE_READ);

% read point
[~, ~, iphysdim, ierr] = cg_base_read(idx_file, idx_base);
[~, size, ierr] = cg_zone_read(idx_file, idx_base, idx_zone);

% Define the range of vertices
rmin = 1; % lower range idx of vertices
rmax = size(1); % upper range idx of vertices

% Read mesh_data coordinates (must use SIDS-standard names here)
point_list = zeros(rmax(1), iphysdim);
[point_list(:, 1), ierr] = cg_coord_read(idx_file, idx_base, idx_zone, ...
    'CoordinateX', RealDouble, rmin, rmax, point_list(:, 1));
[point_list(:, 2), ierr] = cg_coord_read(idx_file, idx_base, idx_zone, ...
    'CoordinateY', RealDouble, rmin, rmax, point_list(:, 2));

if (iphysdim == 3)
    [point_list(:, 3), ierr] = cg_coord_read(idx_file, idx_base, idx_zone, ...
        'CoordinateZ', RealDouble, rmin, rmax, point_list(:, 3));
end

% read BC point index
[out_nbocos, ierr] = cg_nbocos(idx_file, idx_base, idx_zone);
marker_index_list=struct();
for idx_BC=1:out_nbocos

    io_NormalList=0;
    [out_boconame, ~, ~, out_npnts, ~, ~, ~, ~, ierr] =...
        cg_boco_info(idx_file, idx_base, idx_zone, idx_BC, io_NormalList);

    [boolean,idx]=judgeBCName(out_boconame,marker_name_list);
    if boolean
        io_pnts=zeros(out_npnts,1,'int64');
        ierr = cgnslib_mex(int32(141), idx_file, idx_base, idx_zone, idx_BC, io_pnts, io_NormalList);
        marker_index_list.(marker_name_list{idx})=io_pnts;
    end
end

ierr = cg_close(idx_file);

    function [boolean,idx]=judgeBCName(name,BC_name_list)
        % judge if the same BC
        %
        boolean=false(1);
        for idx=1:length(BC_name_list)
            if strcmpi(['BC1_on_',BC_name_list{idx}],name) % ignore capital and small letter
                boolean=true(1);
                return;
            end
        end
    end

end

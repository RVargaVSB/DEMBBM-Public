function elements = ApplyNecking(elements, cords, r, neckRatio)

    x_middle = cords(1)/2;
    y_middle = cords(2)/2;
    z_middle = cords(3)/2;

    y_half = (cords(2) * neckRatio) / 2;
    z_half = (cords(3) * neckRatio) / 2;

    neck_width = 2 * r;

    nodes = elements.node;

    in_x_range = abs(nodes(1,:) - x_middle) < neck_width;

    outside_neck_y = abs(nodes(2,:) - y_middle) > y_half;
    outside_neck_z = abs(nodes(3,:) - z_middle) > z_half;

    removeMask = in_x_range & (outside_neck_y | outside_neck_z);
    keepMask = ~removeMask;

    elements.node     = nodes(:, keepMask);
    elements.radius   = elements.radius(keepMask);
    elements.prop     = elements.prop(keepMask);
    elements.propMat  = elements.propMat(keepMask);
end
function [delta] = get_best_delta(nclusters)
    angle_between_clusters = 2*pi/nclusters;
    radius = 1;
    half_distance = sin(angle_between_clusters/2)*radius;
    delta = 0.2 * half_distance;
end